## Update to ipw 1.3.1

This is a patch release (1.3.0 -> 1.3.1).

### Summary of changes
* Fixes incorrect truncated weights in `ipwpoint()` (`trunc` argument).
* Fixes probability-to-level matching for multinomial/ordinal exposures in
  `ipwtm()` and two-level multinomial exposures in `ipwpoint()`.
* Clearer input validation and error messages; internal refactoring.
* Expanded testthat test suite.

### Test environments
<!-- Fill in after running checks, e.g.: -->
<!-- * local: R x.y.z on ... -->
<!-- * win-builder: R-devel, R-release -->
<!-- * GitHub Actions (R-CMD-check): ubuntu (devel, release, oldrel-1), macOS, windows -->

### R CMD check results
<!-- Fill in after running `devtools::check(remote = TRUE, manual = TRUE)`. -->

### Reverse dependencies
<!-- Run revdepcheck::revdep_check(). Known reverse dependency: ipeval (Suggests). -->
