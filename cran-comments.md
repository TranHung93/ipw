## Submission: ipw 1.3.2

This is a patch release (1.3.1 -> 1.3.2). I am sorry to submit so soon after
1.3.1 (published 2026-10-30). Version 1.3.1 introduced a bug that returns
incorrect results: the truncated weights returned by `ipwpoint()` (argument
`trunc`) were wrong for almost all observations, because each weight was
replaced by either the lower or the upper truncation bound regardless of its
value. Versions 1.3.0 and earlier are not affected. Because users may get
incorrect analysis results, I would like to correct this as soon as possible.

### Changes in this version
* Fixes the incorrect truncated weights in `ipwpoint()` introduced in 1.3.1.
* `ipwtm()` and `ipwpoint()` now handle two-level multinomial exposures, and
  `ipwtm()` handles exposure levels absent among the selected rows; 1.3.1
  stopped with an error in these cases.
* `ipwtm()` with `type = "first"` for multinomial/ordinal exposures again uses
  the lowest level in the whole data set as reference, as in ipw <= 1.2.1.1.
* Clearer input validation and error messages.
* Expanded tests, including a regression test showing that results match
  ipw 1.2.1.1.

### Test environments
* local: Windows 11 x64, R 4.6.1
* GitHub Actions: ubuntu-latest (R-devel, R-release, R-oldrel-1),
  macos-latest (R-release), windows-latest (R-release)
* win-builder: R-devel  <!-- TODO: confirm after devtools::check_win_devel() -->

### R CMD check results
0 errors | 0 warnings | 0 notes

<!-- TODO: if win-builder reports a NOTE, list it here with a short explanation. -->

### Reverse dependencies
<!-- TODO: replace with the result of revdepcheck::revdep_check(), e.g.: -->
We checked 1 reverse dependency (ipeval), comparing R CMD check results
across the CRAN and development versions of this package.
* We saw 0 new problems.
* We failed to check 0 packages.
