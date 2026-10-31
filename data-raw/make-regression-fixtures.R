# Create tests/testthat/fixtures/ipw-1.2.1.1.rds: reference weights from
# CRAN ipw 1.2.1.1, the last release before the 1.3.x rewrite.
#
# Run once from the package root:
#   Rscript data-raw/make-regression-fixtures.R
# or source() it from the package root.
#
# ipw 1.2.1.1 is installed into a temporary library and run in a separate R
# process (via callr), so this works even when ipw is loaded in your session,
# and your installed version of ipw is not touched.
# Needs: callr.

if (!requireNamespace("callr", quietly = TRUE)) install.packages("callr")

helper <- normalizePath(file.path("tests", "testthat", "helper-regression.R"), mustWork = TRUE)
old_lib <- file.path(tempdir(), "ipw-1.2.1.1")
dir.create(old_lib, showWarnings = FALSE)

# 1. Install ipw 1.2.1.1 in a fresh R process (ipw is not loaded there)
callr::r(function(lib) {
  utils::install.packages(
    "https://cloud.r-project.org/src/contrib/Archive/ipw/ipw_1.2.1.1.tar.gz",
    repos = NULL, type = "source", lib = lib, quiet = TRUE
  )
}, args = list(lib = old_lib), show = TRUE)
if (!dir.exists(file.path(old_lib, "ipw"))) stop("Installing ipw 1.2.1.1 failed.")

# 2. Compute the reference results with ipw 1.2.1.1 in another fresh process
ref <- callr::r(function(lib, helper) {
  suppressPackageStartupMessages({
    library(survival)
    library(ipw, lib.loc = lib)
  })
  ver <- read.dcf(file.path(getNamespaceInfo("ipw", "path"), "DESCRIPTION"),
                  fields = "Version")[1, 1]
  if (package_version(ver) != package_version("1.2.1.1")) {
    stop("Loaded ipw ", ver, " instead of 1.2.1.1")
  }

  source(helper, local = TRUE)
  cases <- regression_cases()
  out <- lapply(names(cases), function(nm) {
    res <- tryCatch(cases[[nm]](ipw::ipwpoint, ipw::ipwtm), error = function(e) {
      message("Case '", nm, "' failed in ipw 1.2.1.1: ", conditionMessage(e))
      NULL
    })
    if (is.null(res)) NULL else regression_extract(res)
  })
  names(out) <- names(cases)
  out[!vapply(out, is.null, logical(1))]
}, args = list(lib = old_lib, helper = helper), libpath = c(old_lib, .libPaths()), show = TRUE)

# 3. Save
dir.create(file.path("tests", "testthat", "fixtures"), showWarnings = FALSE, recursive = TRUE)
saveRDS(ref, file.path("tests", "testthat", "fixtures", "ipw-1.2.1.1.rds"), version = 2)
message("Saved reference weights for: ", paste(names(ref), collapse = ", "))
