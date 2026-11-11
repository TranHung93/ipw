# Second review: ipwpoint() against by-hand fits, and findings S-01 to S-03.

r2_point <- function(n = 400, seed = 41) {
  set.seed(seed)
  d <- data.frame(l = stats::rnorm(n), z = stats::rnorm(n))
  d$a <- stats::rbinom(n, 1, stats::plogis(-0.3 + 0.8 * d$l))
  d$off <- 0.7 * d$z
  k <- cut(d$l + stats::rnorm(n), c(-Inf, -0.5, 0.5, Inf), labels = FALSE)
  d$m <- c("low", "mid", "high")[k]
  d$mo <- factor(d$m, levels = c("low", "mid", "high"), ordered = TRUE)
  d$cnt <- stats::rpois(n, exp(0.5 + 0.3 * d$l))
  d$pos <- stats::rgamma(n, shape = 3, rate = 3 / exp(0.2 + 0.3 * d$l))
  d$cont <- 0.5 * d$l + stats::rnorm(n)
  d
}
ate_w <- function(a, p) ifelse(a == 1, 1 / p, 1 / (1 - p))

test_that("binomial weights equal by-hand glm for every link and estimand", {
  d <- r2_point()
  for (lk in c("logit", "probit", "cauchit", "cloglog")) {
    p <- stats::fitted(stats::glm(a ~ l, stats::binomial(lk), d))
    expect_equal(ipwpoint(a, "binomial", lk, NULL, ~ l, d)$ipw.weights, unname(ate_w(d$a, p)), tolerance = 1e-10)
    expect_equal(ipwpoint(a, "binomial", lk, NULL, ~ l, d, estimand = "ATT")$ipw.weights,
                 unname(ifelse(d$a == 1, 1, p / (1 - p))), tolerance = 1e-10)
    expect_equal(ipwpoint(a, "binomial", lk, NULL, ~ l, d, estimand = "ATO")$ipw.weights,
                 unname(ifelse(d$a == 1, 1 - p, p)), tolerance = 1e-10)
  }
})

test_that("count and continuous families equal by-hand density ratios", {
  d <- r2_point()
  g1 <- stats::glm(cnt ~ l, stats::poisson, d); g0 <- stats::glm(cnt ~ 1, stats::poisson, d)
  expect_equal(ipwpoint(cnt, "poisson", NULL, ~ 1, ~ l, d)$ipw.weights,
               unname(stats::dpois(d$cnt, stats::fitted(g0)) / stats::dpois(d$cnt, stats::fitted(g1))), tolerance = 1e-10)
  g1 <- stats::glm(pos ~ l, stats::Gamma("log"), d); g0 <- stats::glm(pos ~ 1, stats::Gamma("log"), d)
  sh <- function(g) 1 / summary(g)$dispersion
  expect_equal(ipwpoint(pos, "gamma", NULL, ~ 1, ~ l, d)$ipw.weights,
               unname(stats::dgamma(d$pos, sh(g0), sh(g0) / stats::fitted(g0)) /
                        stats::dgamma(d$pos, sh(g1), sh(g1) / stats::fitted(g1))), tolerance = 1e-10)
  g1 <- stats::lm(log(pos) ~ l, d); g0 <- stats::lm(log(pos) ~ 1, d)
  expect_equal(ipwpoint(pos, "lognormal", NULL, ~ 1, ~ l, d)$ipw.weights,
               unname(stats::dnorm(log(d$pos), stats::fitted(g0), stats::sd(stats::resid(g0))) /
                        stats::dnorm(log(d$pos), stats::fitted(g1), stats::sd(stats::resid(g1)))), tolerance = 1e-10)
  g1 <- stats::lm(cont ~ l, d); g0 <- stats::lm(cont ~ 1, d)
  expect_equal(ipwpoint(cont, "gaussian", NULL, ~ 1, ~ l, d)$ipw.weights,
               unname(stats::dnorm(d$cont, stats::fitted(g0), stats::sd(stats::resid(g0))) /
                        stats::dnorm(d$cont, stats::fitted(g1), stats::sd(stats::resid(g1)))), tolerance = 1e-10)
})

test_that("multinomial and ordinal without numerator equal by-hand fits (failed in 1.2.1.1)", {
  d <- r2_point()
  mm <- nnet::multinom(m ~ l, d, trace = FALSE); pp <- stats::fitted(mm)
  expect_equal(ipwpoint(m, "multinomial", NULL, NULL, ~ l, d)$ipw.weights,
               unname(1 / pp[cbind(seq_len(nrow(d)), match(d$m, colnames(pp)))]), tolerance = 1e-8)
  for (lk in c("logit", "probit", "cloglog", "cauchit")) {
    meth <- c(logit = "logistic", probit = "probit", cloglog = "cloglog", cauchit = "cauchit")[[lk]]
    pr <- MASS::polr(mo ~ l, d, method = meth, Hess = TRUE); pp <- stats::fitted(pr)
    expect_equal(ipwpoint(mo, "ordinal", lk, NULL, ~ l, d)$ipw.weights,
                 unname(1 / pp[cbind(seq_len(nrow(d)), match(as.character(d$mo), colnames(pp)))]), tolerance = 1e-8)
  }
})

test_that("S-01: an unordered factor is not silently used in alphabetical order as an ordinal exposure", {
  d <- r2_point()
  d$mf <- factor(d$m)   # default levels: high < low < mid (alphabetical)
  expect_warning(ipwpoint(mf, "ordinal", "logit", ~ 1, ~ l, d), regexp = "order")
})

test_that("S-02: offset() terms in the model formulas are kept", {
  d <- r2_point()
  p <- stats::fitted(stats::glm(a ~ l + offset(off), stats::binomial, d))
  expect_equal(ipwpoint(a, "binomial", "logit", NULL, ~ l + offset(off), d)$ipw.weights,
               unname(ate_w(d$a, p)), tolerance = 1e-10)
})

test_that("S-03: a dot formula uses all other columns, as in ipw 1.2.1.1", {
  d <- r2_point()[, c("a", "l", "z")]
  p <- stats::fitted(stats::glm(a ~ ., stats::binomial, d))
  expect_equal(ipwpoint(a, "binomial", "logit", NULL, ~ ., d)$ipw.weights,
               unname(ate_w(d$a, p)), tolerance = 1e-10)
})

test_that("S-10: a non-numeric exposure for a continuous family gives an informative error", {
  d <- r2_point()
  d$y <- factor(d$a)
  expect_error(ipwpoint(y, "gaussian", NULL, ~ 1, ~ l, d), regexp = "numeric")
})
