library(testthat)
library(ipw)

test_that("ipwpoint correctly calculates weights for binomial exposure", {
  set.seed(123)
  n <- 200
  l <- rnorm(n, 10, 5)
  pa <- exp(l - 10) / (1 + exp(l - 10))
  a <- rbinom(n, 1, prob = pa)
  dat <- data.frame(l = l, a = a)

  # Run ipwpoint with numerator = ~ 1 (Stabilized weights)
  res <- ipwpoint(
    exposure = a,
    family = "binomial",
    link = "logit",
    numerator = ~ 1,
    denominator = ~ l,
    data = dat
  )

  # 1. Corrected Name Check: num.mod IS present because numerator was specified
  expect_named(res, c("ipw.weights", "call", "num.mod", "den.mod"))

  # 2. Corrected Value Check:
  # The function calculates stabilized weights: P(A)/P(A|L)
  p_num <- predict(res$num.mod, type = "response")[1]
  p_den <- predict(res$den.mod, type = "response")[1]

  # Calculate expected based on stabilized logic
  if(dat$a[1] == 1) {
    expected_w <- p_num / p_den
  } else {
    expected_w <- (1 - p_num) / (1 - p_den)
  }

  # Use unname() to avoid the 'names(expected) is a character vector' error
  expect_equal(unname(res$ipw.weights[1]), unname(expected_w))
})

test_that("ipwpoint truncation names match", {
  set.seed(123)
  dat <- data.frame(a = rbinom(100, 1, 0.5), l = rnorm(100))

  res_trunc <- ipwpoint(
    exposure = a, family = "binomial", link = "logit",
    numerator = ~ 1, denominator = ~ l, data = dat,
    trunc = 0.05
  )

  # Expecting 5 elements now
  expect_named(res_trunc, c("ipw.weights", "weights.trunc", "call", "num.mod", "den.mod"))
})
