## Submission: ipw 1.4.0

This is a minor release (1.3.2 -> 1.4.0). It comes soon after 1.3.2
(published 2026-11-01) because it fixes two bugs that give wrong weights
without any warning in 1.3.1 and 1.3.2:

* `offset()` terms in the exposure model formulas were silently dropped;
* a 0/1 exposure stored as a factor with levels `c("1", "0")` had the
  probabilities of the two exposure values swapped.

### Main changes

* Bug fixes: the two above, plus several smaller ones that now give clear
  errors instead of wrong results or obscure messages.
* `ipwtm()` survival weights (`type = "first"` or `"cens"`) now use each
  row's own interval, so data with individual visit times can be used
  directly. For data on a full time grid (as in the JSS paper and the
  package examples) the weights are unchanged.
* New features: ATT, ATC and overlap weights in `ipwpoint()` (also for
  categorical exposures); count and skewed exposure families; smooth terms
  via mgcv; `type = "all"` for more families in `ipwtm()`; new exported
  functions `ipwvcov()` (M-estimation standard errors, marked experimental)
  and `ipwprobtm()`. Existing calls give the same weights as before, checked
  by regression tests against ipw 1.2.1.1.
* License changed back to GPL (>= 2), the license of ipw before 1.3.0
  (1.3.0 to 1.3.2 were MIT). All authors agree to the change.
* A new vignette.

### Test environments

* local: Windows 11 x64, R 4.6.1 (2026-06-24 ucrt)
* GitHub Actions: ubuntu-latest (R-devel, R-release, R-oldrel-1),
  macos-latest (R-release), windows-latest (R-release)
* win-builder: R-devel  <!-- fill in from devtools::check_win_devel() -->
* R-hub  <!-- fill in from rhub::rhub_check() -->

### R CMD check results

0 errors | 0 warnings | 0 notes

* local: `devtools::check(remote = TRUE, manual = TRUE)` and
  `R CMD check --as-cran` (including the CRAN incoming checks).
* The package also passes `R CMD check` with only the packages in Depends
  and Imports installed (`_R_CHECK_DEPENDS_ONLY_=true`); all packages in
  Suggests are used conditionally.
<!-- Add here any NOTE from win-builder or R-hub (e.g. "Days since last
update: 12"; the reason is given at the top). -->

### Reverse dependencies

We checked 1 reverse dependency (ipeval 0.1.2), comparing R CMD check results
across the CRAN and development versions of this package
(`revdepcheck::revdep_check()`).
* We saw 0 new problems.
* We failed to check 0 packages.
