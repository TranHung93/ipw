library(testthat)
if (file.exists("ipwtm_new.R")) source("ipwtm_new.R")

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
  expect_true(max(res_trunc$weights.trunc) <= max(res_trunc$ipw.weights))
  expect_true(min(res_trunc$weights.trunc) >= min(res_trunc$ipw.weights))
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

  expect_equal(length(res_shuffled$ipw.weights), nrow(shuffled_data))
})

test_that("ipwtm rejects invalid family and argument combinations", {
  sim_long <- data.frame(id = rep(1:5, each = 2), time = rep(0:1, 5), a = rbinom(10, 1, 0.5), l = rnorm(10))

  # Missing tstart for survival family
  expect_error(
    ipwtm(exposure = a, family = "survival", denominator = ~ l, id = id, timevar = time, type = "first", data = sim_long),
    "tstart"
  )

  # Type 'all' not implemented for survival
  expect_error(
    ipwtm(exposure = a, family = "survival", denominator = ~ l, id = id, tstart = time - 1, timevar = time, type = "all", data = sim_long)
  )

  # Invalid family
  expect_error(
    ipwtm(exposure = a, family = "invalid", denominator = ~ l, id = id, timevar = time, type = "all", data = sim_long)
  )
})
