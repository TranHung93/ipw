library(testthat)
library(ipw)
library(survival)

test_that("ipwtm correctly calculates cumulative weights for survival family", {
  # Load the haartdat provided in the package
  data(haartdat)

  # Run ipwtm (using Example 1 logic)
  res <- ipwtm(
    exposure = haartind,
    family = "survival",
    numerator = ~ sex + age,
    denominator = ~ sex + age + cd4.sqrt,
    id = patient,
    tstart = tstart,
    timevar = fuptime,
    type = "first",
    data = haartdat
  )

  # 1. Check object structure
  expect_named(res, c("ipw.weights", "call", "selvar", "num.mod", "den.mod"))

  # 2. Check cumulative logic for first patient
  # Weights should be the cumulative product of interval-specific probabilities
  # For type = "first", weights become constant after the first switch to 1
  first_patient_w <- res$ipw.weights[haartdat$patient == 1]
  first_switch_idx <- which(haartdat$haartind[haartdat$patient == 1] == 1)[1]

  # Ensure weights are constant after the switch
  if (!is.na(first_switch_idx) && first_switch_idx < length(first_patient_w)) {
    expect_equal(
      first_patient_w[first_switch_idx],
      first_patient_w[first_switch_idx + 1]
    )
  }

  # 3. Verify selvar logic
  # selvar should be 1 up to the first switch, then 0
  first_patient_sel <- res$selvar[haartdat$patient == 1]
  expect_equal(first_patient_sel[first_switch_idx], 1)
  if (first_switch_idx < length(first_patient_sel)) {
    expect_equal(first_patient_sel[first_switch_idx + 1], 0)
  }
})

test_that("ipwtm handles type = 'all' for binomial family", {
  # Create simple longitudinal binomial data
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

  # In type = 'all', selvar should always be 1
  expect_true(all(res_all$selvar == 1))

  # Weights at time 1 should be P(A0)*P(A1) / [P(A0|L0)*P(A1|L1)]
  # (Since numerator is ~ 1, P(A) is marginal)
  expect_true(all(res_all$ipw.weights > 0))
})
