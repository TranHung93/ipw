# Audit tests for the helper functions ipwprobtm(), tstartfun() and ipwplot().

test_that("ipwprobtm rejects values that cannot be probabilities", {
  expect_error(ipwprobtm(den = c(2, 3), id = c(1, 1), timevar = 1:2), regexp = "1")
  expect_error(ipwprobtm(den = c(0, .5), id = c(1, 1), timevar = 1:2), regexp = "positive")
  expect_error(ipwprobtm(den = c(.5, .5), id = c(1, 1), timevar = c(1, 1)), regexp = "only once")
})

test_that("ipwprobtm handles factor ids with unused levels and unsorted rows", {
  w <- ipwprobtm(den = c(.4, .5, .5), id = factor(c(2, 1, 1), levels = 1:3), timevar = c(1, 2, 1))
  expect_equal(w$ipw.weights, c(1 / .4, 1 / .25, 1 / .5))
})

test_that("tstartfun gives tstart < timevar, also when follow-up starts at or below 0", {
  d <- data.frame(id = c(1, 1, 2, 2), t = c(0, 3, 5, 9))
  expect_equal(tstartfun(id, t, d), c(-1, 0, -1, 5))
  # Documented: the first record gets -1. When follow-up starts at or below -1
  # that gives tstart >= timevar, an invalid interval; the user should be told.
  dn <- data.frame(id = c(1, 1), t = c(-5, 3))
  expect_warning(tstartfun(id, t, dn), regexp = "-1")
})

test_that("ipwplot does not change user-settable par() values", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  # usr, xaxp and yaxp are set by every high-level plot and are excluded.
  keep <- function(p) p[setdiff(names(p), c("usr", "xaxp", "yaxp"))]
  p0 <- keep(graphics::par(no.readonly = TRUE))
  ipwplot(c(0.5, 1, 2, 3))
  ipwplot(c(0.5, 1, 2, 3, 1, 2), timevar = c(1, 1, 2, 2, 3, 3), binwidth = 1)
  expect_identical(keep(graphics::par(no.readonly = TRUE)), p0)
})

test_that("ipwplot gives an informative error for weights that cannot be logged", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_error(ipwplot(c(0, 1, 2, 3)), regexp = "positive")
})
