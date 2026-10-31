# ipwtm(family = "survival", type = "all"): Andersen-Gill model for recurrent
# exposure events (Graffeo, Latouche, Geskus & Chevret 2018).

library(survival)

# Irregular start-stop data: each subject has 3-6 intervals with their own
# cut points; expo = 1 when an exposure event happens at the end of an interval.
make_recurrent <- function(seed = 21, n_id = 150) {
  set.seed(seed)
  rows <- lapply(seq_len(n_id), function(i) {
    k <- sample(3:6, 1)
    cuts <- sort(runif(k, 0, 1))
    data.frame(id = i, tstart = c(0, cuts[-k]), fuptime = cuts)
  })
  d <- do.call(rbind, rows)
  d$x <- rep(rbinom(n_id, 1, 0.5), times = as.vector(table(d$id)))
  d$xt <- rbinom(nrow(d), 1, 0.4)
  d$expo <- rbinom(nrow(d), 1, plogis(-1 + 0.5 * d$x - 1.0 * d$xt))
  d
}

# Probability of the observed status per row, written independently of the
# package code with stepfun(): P(no event) = exp(-(L0(stop) - L0(start)) * risk)
hand_weights <- function(fit, d) {
  bh <- basehaz(fit, centered = TRUE)
  L0 <- stepfun(bh$time, c(0, bh$hazard))
  h <- (L0(d$fuptime) - L0(pmax(d$tstart, 0))) * predict(fit, type = "risk")
  p <- ifelse(d$expo == 1, 1 - exp(-h), exp(-h))
  ave(p, d$id, FUN = cumprod)
}

test_that("survival type = 'all' matches a hand calculation on irregular intervals", {
  d <- make_recurrent()
  res <- ipwtm(exposure = expo, family = "survival", numerator = ~ x,
               denominator = ~ x + xt, id = id, tstart = tstart, timevar = fuptime,
               type = "all", data = d)

  expect_true(all(res$selvar == 1))
  den <- coxph(Surv(tstart, fuptime, expo) ~ x + xt, data = d)
  num <- coxph(Surv(tstart, fuptime, expo) ~ x, data = d)
  expected <- hand_weights(num, d) / hand_weights(den, d)
  expect_equal(unname(res$ipw.weights), unname(expected), tolerance = 1e-8)
  expect_equal(unname(coef(res$den.mod)), unname(coef(den)), tolerance = 1e-8)

  unstab <- ipwtm(exposure = expo, family = "survival", denominator = ~ x + xt,
                  id = id, tstart = tstart, timevar = fuptime, type = "all", data = d)
  expect_equal(unname(unstab$ipw.weights), unname(1 / hand_weights(den, d)), tolerance = 1e-8)
})

test_that("survival type = 'all' keeps the input row order", {
  d <- make_recurrent(seed = 22, n_id = 60)
  shuffled <- sample(nrow(d))
  a <- ipwtm(exposure = expo, family = "survival", numerator = ~ x, denominator = ~ x + xt,
             id = id, tstart = tstart, timevar = fuptime, type = "all", data = d)
  b <- ipwtm(exposure = expo, family = "survival", numerator = ~ x, denominator = ~ x + xt,
             id = id, tstart = tstart, timevar = fuptime, type = "all", data = d[shuffled, ])
  expect_equal(b$ipw.weights, a$ipw.weights[shuffled], tolerance = 1e-10)
})

test_that("survival type = 'all' equals type = 'first' when nobody has an event before their last row", {
  # Common time grid, exposure event only possible in a subject's last row:
  # then 'first' selects every row and both types must give the same weights.
  set.seed(23)
  n_id <- 200
  k <- sample(2:6, n_id, replace = TRUE)
  d <- data.frame(id = rep(seq_len(n_id), k), fuptime = unlist(lapply(k, seq_len)))
  d$tstart <- d$fuptime - 1
  d$x <- rep(rbinom(n_id, 1, 0.5), k)
  d$xt <- rbinom(nrow(d), 1, 0.4)
  last <- !duplicated(d$id, fromLast = TRUE)
  d$expo <- ifelse(last, rbinom(nrow(d), 1, 0.5), 0)

  first <- ipwtm(exposure = expo, family = "survival", numerator = ~ x, denominator = ~ x + xt,
                 id = id, tstart = tstart, timevar = fuptime, type = "first", data = d)
  all_t <- ipwtm(exposure = expo, family = "survival", numerator = ~ x, denominator = ~ x + xt,
                 id = id, tstart = tstart, timevar = fuptime, type = "all", data = d)
  expect_true(all(first$selvar == 1))
  expect_equal(all_t$ipw.weights, first$ipw.weights, tolerance = 1e-8)
})

test_that("survival type = 'all' uses rows after the first event, unlike type = 'first'", {
  d <- make_recurrent(seed = 24)
  first <- ipwtm(exposure = expo, family = "survival", numerator = ~ x, denominator = ~ x + xt,
                 id = id, tstart = tstart, timevar = fuptime, type = "first", data = d)
  all_t <- ipwtm(exposure = expo, family = "survival", numerator = ~ x, denominator = ~ x + xt,
                 id = id, tstart = tstart, timevar = fuptime, type = "all", data = d)
  expect_true(any(first$selvar == 0))
  expect_true(all(all_t$selvar == 1))
  expect_false(isTRUE(all.equal(all_t$ipw.weights, first$ipw.weights)))
  expect_true(all(is.finite(all_t$ipw.weights) & all_t$ipw.weights > 0))
  expect_equal(mean(all_t$ipw.weights), 1, tolerance = 0.2)
})

test_that("survival type = 'all' accepts tstart = -1 for the first interval", {
  d <- make_recurrent(seed = 25, n_id = 80)
  d2 <- d
  d2$tstart[d2$tstart == 0] <- -1
  a <- ipwtm(exposure = expo, family = "survival", denominator = ~ x + xt,
             id = id, tstart = tstart, timevar = fuptime, type = "all", data = d)
  b <- ipwtm(exposure = expo, family = "survival", denominator = ~ x + xt,
             id = id, tstart = tstart, timevar = fuptime, type = "all", data = d2)
  expect_equal(b$ipw.weights, a$ipw.weights, tolerance = 1e-8)
})
