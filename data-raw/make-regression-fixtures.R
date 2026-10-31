# Create tests/testthat/fixtures/ipw-1.2.1.1.rds: reference weights from
# CRAN ipw 1.2.1.1, the last release before the 1.3.x rewrite.
#
# Run once from the package root:
#   Rscript data-raw/make-regression-fixtures.R
#
# ipw 1.2.1.1 is installed into a temporary library, so your installed
# version of ipw is not touched. Requires the 'remotes' package.

old_lib <- file.path(tempdir(), "ipw-old")
dir.create(old_lib, showWarnings = FALSE)
remotes::install_version("ipw", version = "1.2.1.1", lib = old_lib,
                         repos = "https://cloud.r-project.org",
                         upgrade = "never", dependencies = FALSE)

library(ipw, lib.loc = old_lib)
stopifnot(as.character(utils::packageVersion("ipw", lib.loc = old_lib)) == "1.2.1.1")

source(file.path("tests", "testthat", "helper-regression.R"))
cases <- regression_cases()

ref <- lapply(names(cases), function(nm) {
  res <- tryCatch(cases[[nm]](ipw::ipwpoint, ipw::ipwtm), error = function(e) {
    message("Case '", nm, "' failed in ipw 1.2.1.1: ", conditionMessage(e))
    NULL
  })
  if (is.null(res)) NULL else regression_extract(res)
})
names(ref) <- names(cases)
ref <- ref[!vapply(ref, is.null, logical(1))]

dir.create(file.path("tests", "testthat", "fixtures"), showWarnings = FALSE, recursive = TRUE)
saveRDS(ref, file.path("tests", "testthat", "fixtures", "ipw-1.2.1.1.rds"), version = 2)
message("Saved reference weights for: ", paste(names(ref), collapse = ", "))
