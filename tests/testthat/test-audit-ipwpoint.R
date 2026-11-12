# Audit tests for ipwpoint(): reference agreement, invariances, input checks.
# Tests encode desired behaviour; failures are linked to the findings table.

make_point <- function(n = 400, seed = 1) {
  set.seed(seed)
  d <- data.frame(l = stats::rnorm(n), f = factor(sample(c("a", "b", "c"), n, TRUE)))
  d$a <- stats::rbinom(n, 1, stats::plogis(-0.3 + d$l + 0.5 * (d$f == "b")))
  d$m <- cut(d$l + stats::rnorm(n), 3, labels = c("1", "2", "3"))
  d$x <- 0.5 * d$l + stats::rnorm(n)
  d
}
bin_w <- function(d, ...) {
  ipwpoint(exposure = a, family = "binomial", link = "logit",
           numerator = ~ 1, denominator = ~ l + f, data = d, ...)$ipw.weights
}

test_that("binomial weights match closed form for ATE, ATT, ATC and ATO", {
  d <- make_point()
  p <- stats::fitted(stats::glm(a ~ l + f, binomial, d))
  expect_equal(ipwpoint(a, "binomial", "logit", denominator = ~ l + f, data = d)$ipw.weights,
               unname(ifelse(d$a == 1, 1 / p, 1 / (1 - p))))
  expect_equal(ipwpoint(a, "binomial", "logit", denominator = ~ l + f, data = d,
                        estimand = "ATT")$ipw.weights,
               unname(ifelse(d$a == 1, 1, p / (1 - p))))
  expect_equal(ipwpoint(a, "binomial", "logit", denominator = ~ l + f, data = d,
                        estimand = "ATC")$ipw.weights,
               unname(ifelse(d$a == 1, (1 - p) / p, 1)))
  expect_equal(ipwpoint(a, "binomial", "logit", denominator = ~ l + f, data = d,
                        estimand = "ATO")$ipw.weights,
               unname(ifelse(d$a == 1, 1 - p, p)))
})

test_that("binomial weights agree with WeightIt", {
  skip_if_not_installed("WeightIt")
  d <- make_point(1000, 2)
  for (est in c("ATE", "ATT", "ATC", "ATO")) {
    expect_equal(
      ipwpoint(a, "binomial", "logit", denominator = ~ l + f, data = d, estimand = est)$ipw.weights,
      unname(WeightIt::weightit(a ~ l + f, data = d, method = "glm", estimand = est)$weights),
      tolerance = 1e-8, label = est)
  }
  expect_equal(bin_w(d),
               unname(WeightIt::weightit(a ~ l + f, data = d, method = "glm", stabilize = TRUE)$weights),
               tolerance = 1e-8)
})

test_that("weights do not depend on the data frame class", {
  d <- make_point()
  ref <- bin_w(d)
  skip_if_not_installed("tibble")
  expect_equal(bin_w(tibble::as_tibble(d)), ref)
  skip_if_not_installed("data.table")
  expect_equal(bin_w(data.table::as.data.table(d)), ref)
})

test_that("weights are invariant to coding, column order, reference level and scale", {
  d <- make_point()
  ref <- bin_w(d)
  d2 <- d; d2$a <- as.integer(d2$a)
  expect_equal(bin_w(d2), ref)
  expect_equal(bin_w(d[, rev(names(d))]), ref)
  d3 <- d; d3$f <- stats::relevel(d3$f, "c")
  expect_equal(bin_w(d3), ref)
  d4 <- d; d4$l <- d4$l * 1000
  expect_equal(bin_w(d4), ref, tolerance = 1e-6)
  expect_equal(ipwpoint("a", "binomial", "logit", ~ 1, ~ l + f, d)$ipw.weights, ref)
})

test_that("a factor 0/1 exposure gives the same weights whatever its level order", {
  d <- make_point()
  ref <- bin_w(d)
  d1 <- d; d1$a <- factor(d1$a, levels = c(0, 1))
  expect_equal(bin_w(d1), ref)
  d2 <- d; d2$a <- factor(d2$a, levels = c(1, 0))
  expect_equal(bin_w(d2), ref)
})

test_that("ordinal exposure without a defined order is rejected", {
  d <- make_point()
  d$m <- c("low", "mid", "high")[as.integer(d$m)]
  expect_error(ipwpoint(m, "ordinal", "logit", ~ 1, ~ l, d), regexp = "order")
})

test_that("ordinal: ordered factor and its integer codes give the same weights", {
  d <- make_point()
  o <- d; o$m <- factor(c("low", "mid", "high")[as.integer(o$m)],
                        levels = c("low", "mid", "high"), ordered = TRUE)
  i <- d; i$m <- as.integer(i$m)
  expect_equal(ipwpoint(m, "ordinal", "logit", ~ 1, ~ l, o)$ipw.weights,
               ipwpoint(m, "ordinal", "logit", ~ 1, ~ l, i)$ipw.weights, tolerance = 1e-6)
})

test_that("missing values give an informative error", {
  d <- make_point()
  dna <- d; dna$l <- NA_real_
  expect_error(bin_w(dna), regexp = "Missing values")
  dnan <- d; dnan$l[3] <- NaN
  expect_error(bin_w(dnan), regexp = "Missing values")
  expect_error(ipwpoint(a, "binomial", "logit", ~ 1, ~ nope, d), regexp = "not found")
  expect_error(ipwpoint(a, "binomial", NULL, ~ 1, ~ l, d), regexp = "link")
})

test_that("infinite covariate values give an informative ipw error", {
  d <- make_point()
  dinf <- d; dinf$l[3] <- Inf
  expect_error(bin_w(dinf), regexp = "(infinite|non-finite)")
})

test_that("an exposure with a single observed value is rejected", {
  d <- make_point()
  dc <- d; dc$a <- 0
  expect_error(bin_w(dc), regexp = "(only one|single|constant|both 0 and 1)")
  expect_error(ipwpoint(a, "binomial", "logit", denominator = ~ l, data = dc, estimand = "ATT"),
               regexp = "(only one|single|constant|both 0 and 1)")
})

test_that("too few rows give an error", {
  d <- make_point()
  expect_error(bin_w(d[1, ]))
  expect_error(bin_w(d[0, ]))
})

test_that("ipwpoint leaves options and the RNG state alone and is reproducible", {
  d <- make_point()
  op <- options()
  invisible(bin_w(d))
  expect_identical(options(), op)
  set.seed(42); s0 <- .Random.seed
  invisible(bin_w(d))
  expect_identical(.Random.seed, s0)
  w1 <- ipwpoint(m, "multinomial", numerator = ~ 1, denominator = ~ l, data = d)$ipw.weights
  w2 <- ipwpoint(m, "multinomial", numerator = ~ 1, denominator = ~ l, data = d)$ipw.weights
  expect_identical(w1, w2)
})

test_that("simulation: stabilized binomial ATE weights recover the marginal effect", {
  skip_on_cran()
  set.seed(2026)
  est <- replicate(200, {
    n <- 800; l <- stats::rnorm(n); a <- stats::rbinom(n, 1, stats::plogis(-0.5 + l))
    d <- data.frame(l, a, y = 1 + 2 * a + l + stats::rnorm(n))
    w <- ipwpoint(a, "binomial", "logit", ~ 1, ~ l, d)$ipw.weights
    stats::coef(stats::lm(y ~ a, d, weights = w))[["a"]]
  })
  mcse <- stats::sd(est) / sqrt(length(est))
  expect_lt(abs(mean(est) - 2), 4 * mcse)
})
