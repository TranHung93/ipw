# Compare current results with those of CRAN ipw 1.2.1.1 (the last version
# before the 1.3.x rewrite). The reference file is created once with
# data-raw/make-regression-fixtures.R; the test is skipped until it exists.
# Weights, truncated weights, selvar and tstart are compared at 1e-8, model
# coefficients at 1e-6 (see regression_tolerance()).

fixture <- test_path("fixtures", "ipw-1.2.1.1.rds")

test_that("results match ipw 1.2.1.1", {
  skip_if_not(file.exists(fixture), "Run data-raw/make-regression-fixtures.R to create the fixture")
  ref <- readRDS(fixture)
  cases <- regression_cases()

  # Every reference case must still exist (cases that need nlme may be absent)
  optional <- "tm_survival_tb_example"
  expect_true(all(setdiff(names(ref), optional) %in% names(cases)),
              label = "all fixture cases are still defined")

  for (nm in intersect(names(cases), names(ref))) {
    # glm warns "fitted probabilities numerically 0 or 1" for the probit and
    # cloglog links on simdat (steep propensity); 1.2.1.1 gives the same warning.
    new <- regression_extract(suppressWarnings(cases[[nm]](ipwpoint, ipwtm, tstartfun)))
    old <- ref[[nm]]
    expect_setequal(names(new), names(old))
    for (el in intersect(names(old), names(new))) {
      expect_equal(new[[el]], old[[el]], tolerance = regression_tolerance(el),
                   label = paste0(nm, "$", el, " (current)"),
                   expected.label = paste0(nm, "$", el, " (1.2.1.1)"))
    }
  }
})

test_that("every regression case has a 1.2.1.1 reference", {
  skip_if_not(file.exists(fixture), "Run data-raw/make-regression-fixtures.R to create the fixture")
  ref <- readRDS(fixture)
  expect_setequal(intersect(names(regression_cases()), names(ref)), names(regression_cases()))
})
