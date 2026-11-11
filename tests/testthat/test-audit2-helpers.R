# Second review: ipwprobtm(), tstartfun(), ipwplot().

test_that("ipwprobtm multiplies ratios in time order within factor ids", {
  w <- ipwprobtm(den = c(.5, .5, .25), id = factor(c("b", "a", "b")), timevar = c(2, 1, 1))$ipw.weights
  expect_equal(w, c(8, 2, 4))
})

test_that("tstartfun gives the previous time within subject on unsorted data", {
  d <- data.frame(id = c(2, 1, 2, 1, 2), t = c(5, 3, 1, 1, 9))
  expect_equal(tstartfun(id, t, d), c(1, 1, -1, -1, 5))
})

test_that("S-10: ipwplot gives informative errors for unusable input", {
  pdf(NULL); on.exit(dev.off())
  expect_error(ipwplot(c(1, 2, NA, 3), logscale = FALSE), regexp = "weights")
  expect_error(ipwplot(1:4, timevar = 1:4, binwidth = c(1, 2)), regexp = "binwidth")
})

test_that("ipwplot leaves open graphics devices alone", {
  pdf(NULL); pdf(NULL); on.exit(graphics.off())
  n_dev <- length(grDevices::dev.list())
  ipwplot(c(0.5, 1, 2, 3))
  expect_equal(length(grDevices::dev.list()), n_dev)
})
