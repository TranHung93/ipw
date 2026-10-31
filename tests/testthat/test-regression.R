# Compare current results with those of CRAN ipw 1.2.1.1 (the last version
# before the 1.3.x rewrite). The reference file is created once with
# data-raw/make-regression-fixtures.R; the test is skipped until it exists.

fixture <- test_path("fixtures", "ipw-1.2.1.1.rds")

test_that("weights match ipw 1.2.1.1", {
  skip_if_not(file.exists(fixture), "Run data-raw/make-regression-fixtures.R to create the fixture")
  ref <- readRDS(fixture)
  cases <- regression_cases()

  for (nm in names(cases)) {
    if (is.null(ref[[nm]])) next
    new <- regression_extract(cases[[nm]](ipwpoint, ipwtm))
    old <- ref[[nm]]

    for (el in intersect(names(old), names(new))) {
      expect_equal(new[[el]], old[[el]], tolerance = 1e-6,
                   label = paste0(nm, "$", el, " (1.3.1)"),
                   expected.label = paste0(nm, "$", el, " (1.2.1.1)"))
    }
  }
})
