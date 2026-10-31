library(testthat)
if (file.exists("ipwpoint_new.R")) source("ipwpoint_new.R")

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

  # Check truncation bounds
  p05 <- quantile(res_trunc$ipw.weights, 0.05)
  p95 <- quantile(res_trunc$ipw.weights, 0.95)
  expect_true(all(res_trunc$weights.trunc >= p05 - 1e-8))
  expect_true(all(res_trunc$weights.trunc <= p95 + 1e-8))
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
