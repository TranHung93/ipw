test_that("ipwpoint correctly calculates stabilized and unstabilized weights for binomial exposure", {
  set.seed(123)
  n <- 200
  l <- rnorm(n, 10, 5)
  pa <- exp((l - 10) / 5) / (1 + exp((l - 10) / 5))
  a <- rbinom(n, 1, prob = pa)
  dat <- data.frame(l = l, a = a)

  # 1. Stabilized weights (numerator = ~ 1)
  res_stab <- ipwpoint(
    exposure = a,
    family = "binomial",
    link = "logit",
    numerator = ~ 1,
    denominator = ~ l,
    data = dat
  )

  expect_s3_class(res_stab, "ipwpoint")
  expect_named(res_stab, c("ipw.weights", "call", "num.mod", "den.mod"))

  # Math verification for first observation
  p_num <- predict(res_stab$num.mod, type = "response")[1]
  p_den <- predict(res_stab$den.mod, type = "response")[1]
  expected_w <- if (dat$a[1] == 1) p_num / p_den else (1 - p_num) / (1 - p_den)
  expect_equal(unname(res_stab$ipw.weights[1]), unname(expected_w))

  # Mean of stabilized weights should be approximately 1
  expect_equal(mean(res_stab$ipw.weights), 1.0, tolerance = 0.15)

  # 2. Unstabilized weights (numerator = NULL)
  res_unstab <- ipwpoint(
    exposure = a,
    family = "binomial",
    link = "logit",
    numerator = NULL,
    denominator = ~ l,
    data = dat
  )

  expect_null(res_unstab$num.mod)
  expect_named(res_unstab, c("ipw.weights", "call", "den.mod"))
  expected_unstab <- if (dat$a[1] == 1) 1 / p_den else 1 / (1 - p_den)
  expect_equal(unname(res_unstab$ipw.weights[1]), unname(expected_unstab))
})

test_that("ipwpoint handles multinomial exposure", {
  set.seed(123)
  n <- 150
  l <- rnorm(n)
  a <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  dat <- data.frame(a = a, l = l)

  res_multi <- ipwpoint(
    exposure = a,
    family = "multinomial",
    numerator = ~ 1,
    denominator = ~ l,
    data = dat
  )

  expect_s3_class(res_multi, "ipwpoint")
  expect_named(res_multi, c("ipw.weights", "call", "num.mod", "den.mod"))
  expect_true(all(res_multi$ipw.weights > 0))
})

test_that("ipwpoint handles ordinal exposure", {
  set.seed(123)
  n <- 150
  l <- rnorm(n)
  a <- factor(sample(1:3, n, replace = TRUE), ordered = TRUE)
  dat <- data.frame(a = a, l = l)

  res_ord <- ipwpoint(
    exposure = a,
    family = "ordinal",
    link = "logit",
    numerator = ~ 1,
    denominator = ~ l,
    data = dat
  )

  expect_s3_class(res_ord, "ipwpoint")
  expect_named(res_ord, c("ipw.weights", "call", "num.mod", "den.mod"))
  expect_true(all(res_ord$ipw.weights > 0))
})

test_that("ipwpoint handles gaussian exposure (continuous treatment)", {
  set.seed(123)
  n <- 150
  l <- rnorm(n)
  a <- 2 * l + rnorm(n)
  dat <- data.frame(a = a, l = l)

  # Gaussian family requires numerator formula
  res_gauss <- ipwpoint(
    exposure = a,
    family = "gaussian",
    numerator = ~ 1,
    denominator = ~ l,
    data = dat
  )

  expect_s3_class(res_gauss, "ipwpoint")
  expect_named(res_gauss, c("ipw.weights", "call", "num.mod", "den.mod"))
  expect_true(all(res_gauss$ipw.weights > 0))
})

test_that("ipwpoint correctly applies weight truncation", {
  set.seed(123)
  dat <- data.frame(a = rbinom(100, 1, 0.5), l = rnorm(100))

  res_trunc <- ipwpoint(
    exposure = a,
    family = "binomial",
    link = "logit",
    numerator = ~ 1,
    denominator = ~ l,
    data = dat,
    trunc = 0.05
  )

  expect_named(res_trunc, c("ipw.weights", "weights.trunc", "call", "num.mod", "den.mod"))

  w  <- res_trunc$ipw.weights
  wt <- res_trunc$weights.trunc
  q  <- unname(quantile(w, c(0.05, 0.95)))

  # Weights inside the bounds are untouched, those outside are set to the bound
  inside <- w >= q[1] & w <= q[2]
  expect_equal(wt[inside], w[inside])
  expect_equal(wt[w < q[1]], rep(q[1], sum(w < q[1])))
  expect_equal(wt[w > q[2]], rep(q[2], sum(w > q[2])))
  expect_equal(wt, pmin(pmax(w, q[1]), q[2]))
})

test_that("ipwpoint accepts column references as symbols and character strings", {
  set.seed(123)
  dat <- data.frame(treat = rbinom(50, 1, 0.5), covar = rnorm(50))

  res_sym <- ipwpoint(exposure = treat, family = "binomial", link = "logit", denominator = ~ covar, data = dat)
  res_str <- ipwpoint(exposure = "treat", family = "binomial", link = "logit", denominator = ~ covar, data = dat)

  expect_equal(res_sym$ipw.weights, res_str$ipw.weights)
})

test_that("ipwpoint rejects invalid input configurations", {
  dat <- data.frame(a = rbinom(50, 1, 0.5), l = rnorm(50))

  # Invalid family
  expect_error(ipwpoint(exposure = a, family = "invalid", denominator = ~ l, data = dat), "family")

  # Invalid link for binomial
  expect_error(ipwpoint(exposure = a, family = "binomial", link = "invalid_link", denominator = ~ l, data = dat), "link")

  # Missing denominator
  expect_error(ipwpoint(exposure = a, family = "binomial", link = "logit", data = dat), "denominator")

  # Missing numerator for gaussian family
  expect_error(ipwpoint(exposure = l, family = "gaussian", denominator = ~ a, data = dat), "numerator")

  # Invalid trunc values
  expect_error(ipwpoint(exposure = a, family = "binomial", link = "logit", denominator = ~ l, data = dat, trunc = 0.8), "trunc")
  expect_error(ipwpoint(exposure = a, family = "binomial", link = "logit", denominator = ~ l, data = dat, trunc = -0.1), "trunc")
})

test_that("ipwpoint handles a two-level multinomial exposure", {
  set.seed(1)
  n <- 200
  l <- rnorm(n)
  a <- factor(ifelse(runif(n) < plogis(l), "B", "A"))
  dat <- data.frame(a = a, l = l)

  res_multi <- ipwpoint(exposure = a, family = "multinomial",
                        numerator = ~ 1, denominator = ~ l, data = dat)

  # Must equal the binomial-logit weights for the same exposure
  dat$a01 <- as.numeric(dat$a == "B")
  res_bin <- ipwpoint(exposure = a01, family = "binomial", link = "logit",
                      numerator = ~ 1, denominator = ~ l, data = dat)
  expect_equal(unname(res_multi$ipw.weights), unname(res_bin$ipw.weights), tolerance = 1e-3)
})

test_that("ipwpoint models keep a readable call and support summary()", {
  set.seed(2)
  dat <- data.frame(l = rnorm(150))
  dat$a <- factor(cut(dat$l + rnorm(150), 3, labels = c("lo", "mid", "hi")), ordered = TRUE)

  res <- ipwpoint(exposure = a, family = "ordinal", link = "logit",
                  numerator = ~ 1, denominator = ~ l, data = dat)
  expect_identical(res$den.mod$call$data, quote(dat))
  expect_no_error(summary(res$den.mod))

  dat$b <- rbinom(150, 1, plogis(dat$l))
  res_b <- ipwpoint(exposure = b, family = "binomial", link = "probit",
                    denominator = ~ l, data = dat)
  expect_identical(res_b$den.mod$call$data, quote(dat))
  expect_identical(res_b$den.mod$call$family, quote(binomial(link = "probit")))
})

test_that("ipwpoint rejects missing values and non 0/1 binomial exposure", {
  dat <- data.frame(a = c(0, 1, 1, 0, 1, 0), l = c(1, 2, NA, 4, 5, 6))
  expect_error(ipwpoint(exposure = a, family = "binomial", link = "logit",
                        denominator = ~ l, data = dat), "Missing values")

  dat2 <- data.frame(a = c(0, 1, 2, 0, 1, 0), l = rnorm(6))
  expect_error(ipwpoint(exposure = a, family = "binomial", link = "logit",
                        denominator = ~ l, data = dat2), "0 or 1")

  expect_error(ipwpoint(exposure = notacolumn, family = "binomial", link = "logit",
                        denominator = ~ l, data = dat2), "not found")
})

test_that("ipwpoint gives unstabilized weights for multinomial and ordinal exposures (#3)", {
  # Reproducible example from GitHub issue #3
  set.seed(1)
  n <- 200
  d <- data.frame(l = rnorm(n))
  eta2 <- -0.5 + d$l
  eta3 <- -0.5 + 2 * d$l
  den <- 1 + exp(eta2) + exp(eta3)
  p1 <- 1 / den
  p2 <- exp(eta2) / den
  u <- runif(n)
  d$a <- factor(ifelse(u < p1, 1, ifelse(u < p1 + p2, 2, 3)))
  d$ao <- ordered(d$a)

  res_m <- ipwpoint(exposure = a, family = "multinomial", denominator = ~ l, data = d)
  expect_length(res_m$ipw.weights, n)
  expect_null(res_m$num.mod)
  probs_m <- predict(res_m$den.mod, type = "probs")
  p_obs_m <- probs_m[cbind(seq_len(n), match(as.character(d$a), colnames(probs_m)))]
  expect_equal(unname(res_m$ipw.weights), unname(1 / p_obs_m))
  # Mean of unstabilized weights is close to the number of categories
  expect_equal(mean(res_m$ipw.weights), 3, tolerance = 0.25)

  res_o <- ipwpoint(exposure = ao, family = "ordinal", link = "logit", denominator = ~ l, data = d)
  expect_length(res_o$ipw.weights, n)
  expect_null(res_o$num.mod)
  probs_o <- predict(res_o$den.mod, type = "probs")
  p_obs_o <- probs_o[cbind(seq_len(n), match(as.character(d$ao), colnames(probs_o)))]
  expect_equal(unname(res_o$ipw.weights), unname(1 / p_obs_o))
})
