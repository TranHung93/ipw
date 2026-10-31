# Changelog

## ipw 1.3.2

### Bug fixes

- [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md):
  `weights.trunc` was wrong whenever `trunc` was set. The lower and
  upper quantiles were recycled across rows, so every truncated weight
  was replaced by one of the two bounds. Weights are now clipped at the
  `trunc` and `1 - trunc` quantiles, and weights between the bounds are
  left unchanged.
- [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) with
  `family = "multinomial"` or `"ordinal"`: predicted probabilities are
  now matched to the observed exposure by level name rather than by
  column position. 1.3.1 stopped with an error for a two-level exposure,
  or when an exposure level did not occur among the rows used to fit the
  models.
- [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  with `family = "multinomial"` now works for a two-level exposure;
  1.3.1 stopped with an error (“Observed exposure levels do not match
  fitted probability levels”).
- Added a test for unstabilized weights (no `numerator`) with
  `family = "multinomial"` and `"ordinal"` in
  [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md).
  These failed in 1.3.0 and were fixed by the 1.3.1 rewrite
  ([\#3](https://github.com/TranHung93/ipw/issues/3),
  [@adcascone](https://github.com/adcascone)).
- [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) with
  `type = "first"` and `family = "multinomial"` / `"ordinal"`: a switch
  is again defined as the first record at any level other than the
  lowest level in the whole data set, as in ipw \<= 1.2.1.1. In 1.3.1
  the reference level was chosen separately for each subject.

### Improvements

- [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) now
  checks that all model variables exist and contain no missing values,
  and that a binomial/survival exposure is coded 0/1, with clear error
  messages instead of errors from inside
  [`glm()`](https://rdrr.io/r/stats/glm.html) /
  [`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html).
- [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) with
  `family = "survival"` now stops with a clear message when a formula
  contains [`strata()`](https://rdrr.io/pkg/survival/man/strata.html),
  which is not supported.
- Fitted models returned in `num.mod` / `den.mod` no longer carry an
  inlined copy of the data in their `call`, so printing them is readable
  again.
- [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  and [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md)
  share one set of internal helpers for argument checks, formula
  construction, probability extraction and truncation.
- Removed blanket `import()` of whole packages; all calls use
  `pkg::fun`. `methods` is no longer imported.

### Tests

- Test suite expanded: exact checks of truncation, row-order invariance,
  two-level multinomial vs binomial equivalence, and error messages.
- New regression test against ipw 1.2.1.1 (reference weights created
  with `data-raw/make-regression-fixtures.R`).

## ipw 1.3.1

CRAN release: 2026-10-30

- **Core Refactoring (`ipwtm`, `ipwpoint`, `ipwplot`, `tstartfun`)**:
  - Improved defensive input validation, parameter checking, and error
    messages across all core functions.
  - Modernized formula construction
    ([`as.formula()`](https://rdrr.io/r/stats/formula.html)) to
    eliminate fragile `eval(parse(text = ...))` string manipulations.
  - Fixed variable symbol evaluation issues when passing variables as
    function arguments rather than raw names.
  - Modularized `ipwtm` into dedicated family modeling engines
    (`binomial`, `survival`, `multinomial`, `ordinal`, `gaussian`).
  - Replaced manual survival hazard step-differencing with continuous
    baseline hazard interpolation using
    [`survival::basehaz`](https://rdrr.io/pkg/survival/man/basehaz.html)
    and [`stats::approx`](https://rdrr.io/r/stats/approxfun.html).
  - Refactored `tstartfun` to use safe column evaluation and direct
    vector lag logic (`c(-1, t[-length(t)])`).
- **Test Suite Expansion**:
  - Upgraded `testthat` suite to 3rd edition specs, expanding coverage
    for `multinomial`, `ordinal`, `gaussian`, unstabilized models,
    weight truncation, and invalid parameter rejection.

## ipw 1.3.0

CRAN release: 2026-01-27

- New maintainer: Hung Thai Tran.
- License changed to MIT (earlier versions: GPL (\>= 2)).
