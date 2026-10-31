test_that("ipwtm correctly calculates cumulative weights for survival family", {
  if (requireNamespace("ipw", quietly = TRUE)) {
    data(haartdat, package = "ipw")
    test_dat <- haartdat
  } else {
    set.seed(123)
    test_dat <- data.frame(
      patient = rep(1:20, each = 4),
      tstart = rep(c(-1, 30, 60, 90), 20),
      fuptime = rep(c(30, 60, 90, 120), 20),
      haartind = rbinom(80, 1, 0.3),
      sex = rep(rbinom(20, 1, 0.5), each = 4),
      age = rep(runif(20, 20, 50), each = 4),
      cd4.sqrt = rnorm(80, 300, 50)
    )
  }

  res <- ipwtm(
    exposure = haartind,
    family = "survival",
    numerator = ~ sex + age,
    denominator = ~ sex + age + cd4.sqrt,
    id = patient,
    tstart = tstart,
    timevar = fuptime,
    type = "first",
    data = test_dat
  )

  expect_s3_class(res, "ipwtm")
  expect_named(res, c("ipw.weights", "call", "selvar", "num.mod", "den.mod"))

  # Verify selvar logic for type = "first"
  first_p_id <- test_dat$patient[1]
  p_rows <- test_dat$patient == first_p_id
  p_exp <- test_dat$haartind[p_rows]
  p_sel <- res$selvar[p_rows]

  first_switch <- which(p_exp == 1)[1]
  if (!is.na(first_switch)) {
    expect_equal(p_sel[first_switch], 1)
    if (first_switch < length(p_sel)) {
      expect_equal(p_sel[(first_switch + 1):length(p_sel)], rep(0, length(p_sel) - first_switch))
    }
  }
})

test_that("ipwtm handles type = 'all' for binomial family", {
  set.seed(123)
  sim_long <- data.frame(
    id = rep(1:10, each = 3),
    time = rep(0:2, 10),
    a = rbinom(30, 1, 0.5),
    l = rnorm(30)
  )

  res_all <- ipwtm(
    exposure = a,
    family = "binomial",
    link = "logit",
    numerator = ~ 1,
    denominator = ~ l,
    id = id,
    timevar = time,
    type = "all",
    data = sim_long
  )

  expect_s3_class(res_all, "ipwtm")
  expect_true(all(res_all$selvar == 1))
  expect_true(all(res_all$ipw.weights > 0))
})

test_that("ipwtm correctly applies weight truncation in longitudinal setting", {
  set.seed(123)
  sim_long <- data.frame(
    id = rep(1:20, each = 3),
    time = rep(0:2, 20),
    a = rbinom(60, 1, 0.4),
    l = rnorm(60)
  )

  res_trunc <- ipwtm(
    exposure = a,
    family = "binomial",
    link = "logit",
    numerator = ~ 1,
    denominator = ~ l,
    id = id,
    timevar = time,
    type = "all",
    data = sim_long,
    trunc = 0.05
  )

  expect_named(res_trunc, c("ipw.weights", "weights.trunc", "call", "selvar", "num.mod", "den.mod"))
  w <- res_trunc$ipw.weights
  q <- unname(quantile(w, c(0.05, 0.95)))
  expect_equal(res_trunc$weights.trunc, pmin(pmax(w, q[1]), q[2]))
})

test_that("ipwtm preserves original row ordering on unsorted longitudinal datasets", {
  set.seed(123)
  sim_long <- data.frame(
    id = rep(1:10, each = 3),
    time = rep(0:2, 10),
    a = rbinom(30, 1, 0.5),
    l = rnorm(30)
  )
  shuffled_idx <- sample(1:30)
  shuffled_data <- sim_long[shuffled_idx, ]

  res_shuffled <- ipwtm(
    exposure = a,
    family = "binomial",
    link = "logit",
    numerator = ~ 1,
    denominator = ~ l,
    id = id,
    timevar = time,
    type = "all",
    data = shuffled_data
  )

  res_sorted <- ipwtm(
    exposure = a,
    family = "binomial",
    link = "logit",
    numerator = ~ 1,
    denominator = ~ l,
    id = id,
    timevar = time,
    type = "all",
    data = sim_long
  )

  expect_equal(length(res_shuffled$ipw.weights), nrow(shuffled_data))
  # Weights follow the rows of the input data, whatever their order
  expect_equal(res_shuffled$ipw.weights, res_sorted$ipw.weights[shuffled_idx])
  expect_equal(res_shuffled$selvar, res_sorted$selvar[shuffled_idx])
})

test_that("ipwtm rejects invalid family and argument combinations", {
  sim_long <- data.frame(id = rep(1:5, each = 2), time = rep(0:1, 5), a = rbinom(10, 1, 0.5), l = rnorm(10))

  # Missing tstart for survival family
  expect_error(
    ipwtm(exposure = a, family = "survival", denominator = ~ l, id = id, timevar = time, type = "first", data = sim_long),
    "tstart"
  )

  # Type 'all' not implemented for survival
  sim_long$tstart <- sim_long$time - 1
  expect_error(
    ipwtm(exposure = a, family = "survival", denominator = ~ l, id = id, tstart = tstart,
          timevar = time, type = "all", data = sim_long),
    "not yet implemented"
  )

  # strata() not supported for survival
  expect_error(
    ipwtm(exposure = a, family = "survival", denominator = ~ strata(l), id = id, tstart = tstart,
          timevar = time, type = "first", data = sim_long),
    "strata"
  )

  # Invalid family
  expect_error(
    ipwtm(exposure = a, family = "invalid", denominator = ~ l, id = id, timevar = time, type = "all", data = sim_long),
    "should be one of"
  )

  # Missing values and non 0/1 binomial exposure
  na_dat <- sim_long
  na_dat$l[3] <- NA
  expect_error(
    ipwtm(exposure = a, family = "binomial", link = "logit", denominator = ~ l, id = id,
          timevar = time, type = "all", data = na_dat),
    "Missing values"
  )
  bad_dat <- sim_long
  bad_dat$a[1] <- 2
  expect_error(
    ipwtm(exposure = a, family = "binomial", link = "logit", denominator = ~ l, id = id,
          timevar = time, type = "all", data = bad_dat),
    "0 or 1"
  )
})

test_that("ipwtm multinomial selection uses the global reference level", {
  # Subject 2 is at level 2 throughout: it has switched away from the
  # lowest level (1) at its first record, so only that record is selected.
  dat <- data.frame(
    id   = rep(1:3, each = 3),
    time = rep(1:3, 3),
    a    = c(1, 1, 2,   2, 2, 2,   1, 3, 3)
  )
  sel <- ipw:::.compute_selvar(dat$a, dat$id, type = "first", family = "multinomial")
  expect_equal(sel, c(1, 1, 1,   1, 0, 0,   1, 1, 0))
})

test_that("ipwtm handles multinomial and ordinal exposures, incl. two levels", {
  set.seed(42)
  n_id <- 150
  dat <- data.frame(id = rep(seq_len(n_id), each = 4), time = rep(1:4, n_id))
  dat$l <- rnorm(nrow(dat))
  dat$a <- 1
  for (i in seq_len(n_id)) {
    rows <- which(dat$id == i)
    for (r in rows[-1]) {
      prev <- dat$a[r - 1]
      dat$a[r] <- if (prev == 1) sample(1:3, 1, prob = c(0.7, 0.15, 0.15)) else prev
    }
  }

  res_m <- ipwtm(exposure = a, family = "multinomial", numerator = ~ 1, denominator = ~ l,
                 id = id, timevar = time, type = "first", data = dat)
  expect_true(all(is.finite(res_m$ipw.weights) & res_m$ipw.weights > 0))
  expect_identical(res_m$den.mod$call$data, quote(dat))

  res_o <- ipwtm(exposure = a, family = "ordinal", link = "logit", numerator = ~ 1,
                 denominator = ~ l, id = id, timevar = time, type = "first", data = dat)
  expect_true(all(is.finite(res_o$ipw.weights) & res_o$ipw.weights > 0))
  expect_no_error(summary(res_o$den.mod))

  # Two-level multinomial must match binomial "first"
  dat2 <- dat
  dat2$a <- as.numeric(dat2$a > 1)
  res_m2 <- ipwtm(exposure = a, family = "multinomial", numerator = ~ 1, denominator = ~ l,
                  id = id, timevar = time, type = "first", data = dat2)
  res_b2 <- ipwtm(exposure = a, family = "binomial", link = "logit", numerator = ~ 1,
                  denominator = ~ l, id = id, timevar = time, type = "first", data = dat2)
  expect_equal(unname(res_m2$ipw.weights), unname(res_b2$ipw.weights), tolerance = 1e-3)
})

test_that("ipwtm gaussian models keep a readable call", {
  skip_if_not_installed("geepack")
  set.seed(7)
  dat <- data.frame(id = rep(1:40, each = 3), time = rep(1:3, 40))
  dat$l <- rnorm(nrow(dat))
  dat$a <- 0.5 * dat$l + rnorm(nrow(dat))

  res <- ipwtm(exposure = a, family = "gaussian", numerator = ~ 1, denominator = ~ l,
               id = id, timevar = time, type = "all", corstr = "ar1", data = dat)
  expect_true(all(is.finite(res$ipw.weights) & res$ipw.weights > 0))
  expect_identical(res$den.mod$call$data, quote(dat))
  expect_identical(res$den.mod$call$id, quote(id))
})
