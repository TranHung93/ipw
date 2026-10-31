## Update to ipw 1.3.2

This is a patch release (1.3.1 -> 1.3.2), submitted shortly after 1.3.1
because 1.3.1 returns incorrect truncated weights from `ipwpoint()`.

### Summary of changes
* Fixes incorrect truncated weights in `ipwpoint()` (`trunc` argument): in
  1.3.1 every truncated weight was set to one of the two truncation bounds.
* `ipwtm()` and `ipwpoint()` now handle two-level multinomial exposures, and
  `ipwtm()` handles exposure levels absent among the selected rows; 1.3.1
  stopped with an error in these cases.
* `ipwtm()` with `type = "first"` for multinomial/ordinal exposures again uses
  the lowest level in the whole data set as reference, as in ipw <= 1.2.1.1.
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
