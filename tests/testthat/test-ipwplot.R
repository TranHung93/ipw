library(testthat)
if (file.exists("ipwplot_new.R")) source("ipwplot_new.R")

test_that("ipwplot executes density plots for point treatment weights", {
  set.seed(123)
  w <- runif(100, 0.5, 5)

  # Default logscale = TRUE
  expect_no_error(ipwplot(weights = w))

  # Linear scale (logscale = FALSE)
  expect_no_error(ipwplot(weights = w, logscale = FALSE))

  # Without reference line
  expect_no_error(ipwplot(weights = w, ref = FALSE))

  # Custom main title and axis labels
  expect_no_error(ipwplot(weights = w, xlab = "Weights", ylab = "Density", main = "Density Plot"))
})

test_that("ipwplot executes boxplots for time-varying weights across follow-up strata", {
  set.seed(123)
  w <- runif(100, 0.5, 5)
  time <- rep(seq(0, 400, by = 100), length.out = 100)

  # Valid timevar and binwidth
  expect_no_error(ipwplot(weights = w, timevar = time, binwidth = 100))

  # Linear scale boxplots
  expect_no_error(ipwplot(weights = w, timevar = time, binwidth = 100, logscale = FALSE))

  # Custom graphical parameters
  expect_no_error(
    ipwplot(weights = w, timevar = time, binwidth = 100,
            xlab = "Follow-up Time (Days)", ylab = "log(IPW)", main = "Time-Varying IPW")
  )
})

test_that("ipwplot correctly handles complex expressions and label deparsing", {
  w_vec <- runif(50, 0.5, 2)

  # Complex expression in weights argument
  expect_no_error(ipwplot(weights = w_vec * 2 + 0.5, logscale = TRUE))
})

test_that("ipwplot rejects invalid arguments with informative errors", {
  w <- runif(20, 0.5, 3)
  time <- rep(1:5, 4)

  # Missing or non-numeric weights
  expect_error(ipwplot(), "weights")
  expect_error(ipwplot(weights = "invalid"), "weights")
  expect_error(ipwplot(weights = numeric(0)), "weights")

  # Timevar supplied without required binwidth
  expect_error(ipwplot(weights = w, timevar = time), "binwidth")

  # Invalid binwidth (<= 0 or non-numeric)
  expect_error(ipwplot(weights = w, timevar = time, binwidth = 0), "binwidth")
  expect_error(ipwplot(weights = w, timevar = time, binwidth = -10), "binwidth")

  # Mismatched length between weights and timevar
  expect_error(ipwplot(weights = w, timevar = time[1:10], binwidth = 2), "same length")
})
