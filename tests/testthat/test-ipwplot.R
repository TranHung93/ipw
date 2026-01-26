library(testthat)
library(ipw)

test_that("ipwplot executes without error for point treatment (density plot)", {
  # Generate dummy weights
  set.seed(123)
  w <- runif(100, 0.5, 5)

  # Check if it runs with default settings (logscale = TRUE)
  expect_error(ipwplot(weights = w), NA)

  # Check if it runs with logscale = FALSE
  expect_error(ipwplot(weights = w, logscale = FALSE), NA)

  # Check if it runs without a reference line
  expect_error(ipwplot(weights = w, ref = FALSE), NA)
})

test_that("ipwplot executes without error for time-varying weights (boxplots)", {
  # Generate dummy weights and time
  set.seed(123)
  w <- runif(100, 0.5, 5)
  time <- seq(1, 1000, length.out = 100)

  # Check if it runs with binwidth specified
  expect_error(ipwplot(weights = w, timevar = time, binwidth = 100), NA)

  # Check with custom labels
  expect_error(
    ipwplot(weights = w, timevar = time, binwidth = 100,
            xlab = "Days", ylab = "Weights", main = "Test Plot"),
    NA
  )
})

test_that("ipwplot correctly handles log(weights) and labels", {
  # This test ensures the internal label logic doesn't crash
  w_vec <- runif(10)

  # Testing the deparse(match.call()$weights) logic
  # Using a local capture to check if the function can handle complex expressions
  expect_error(ipwplot(weights = w_vec + 1, logscale = TRUE), NA)
})
