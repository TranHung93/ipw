# ipw (development version)

## New features

* `ipwtm()`: `type = "all"` now works for `family = "multinomial"` and
  `"ordinal"` (time-varying categorical exposures that switch back and forth).
* New `ipwprob()`: inverse probability weights from probabilities you
  estimated yourself (for example with machine learning or SuperLearner), for
  point treatments or cumulated over time within subjects.
* `ipwpoint()`: new families for count exposures (`"poisson"`, `"negbin"`) and
  skewed continuous exposures (`"gamma"`, `"lognormal"`).
* Smooth terms (`s()`, `te()`, `ti()`, `t2()`) in model formulas are fitted as
  generalized additive models with mgcv: in `ipwpoint()` for the binomial,
  gaussian and poisson families, and in `ipwtm()` for the binomial family.
* `ipwpoint()`: new `estimand` argument for a binomial exposure: `"ATE"`
  (default), `"ATT"` (effect in the treated) and `"ATO"` (overlap weights).

# ipw 1.3.2

## Bug fixes

* `ipwpoint()`: `weights.trunc` was wrong whenever `trunc` was set. The lower
  and upper quantiles were recycled across rows, so every truncated weight was
  replaced by one of the two bounds. Weights are now clipped at the `trunc` and
  `1 - trunc` quantiles, and weights between the bounds are left unchanged.
* `ipwtm()` with `family = "multinomial"` or `"ordinal"`: predicted
  probabilities are now matched to the observed exposure by level name rather
  than by column position. 1.3.1 stopped with an error for a two-level
  exposure, or when an exposure level did not occur among the rows used to fit
  the models.
* `ipwpoint()` with `family = "multinomial"` now works for a two-level
  exposure; 1.3.1 stopped with an error ("Observed exposure levels do not match
  fitted probability levels").
* Added a test for unstabilized weights (no `numerator`) with
  `family = "multinomial"` and `"ordinal"` in `ipwpoint()`. These failed in
  1.3.0 and were fixed by the 1.3.1 rewrite (#3, @adcascone).
* `ipwtm()` with `type = "first"` and `family = "multinomial"` / `"ordinal"`:
  a switch is again defined as the first record at any level other than the
  lowest level in the whole data set, as in ipw <= 1.2.1.1. In 1.3.1 the
  reference level was chosen separately for each subject.

## Improvements

* `ipwtm()` now checks that all model variables exist and contain no missing
  values, and that a binomial/survival exposure is coded 0/1, with clear error
  messages instead of errors from inside `glm()` / `coxph()`.
* `ipwtm()` with `family = "survival"` now stops with a clear message when a
  formula contains `strata()`, which is not supported.
* Fitted models returned in `num.mod` / `den.mod` no longer carry an inlined
  copy of the data in their `call`, so printing them is readable again.
* `ipwpoint()` and `ipwtm()` share one set of internal helpers for argument
  checks, formula construction, probability extraction and truncation.
* Removed blanket `import()` of whole packages; all calls use `pkg::fun`.
  `methods` is no longer imported.

## Tests

* Test suite expanded: exact checks of truncation, row-order invariance,
  two-level multinomial vs binomial equivalence, and error messages.
* New regression test against ipw 1.2.1.1 (reference weights created with
  `data-raw/make-regression-fixtures.R`).

# ipw 1.3.1

* **Core Refactoring (`ipwtm`, `ipwpoint`, `ipwplot`, `tstartfun`)**:
  * Improved defensive input validation, parameter checking, and error messages across all core functions.
  * Modernized formula construction (`as.formula()`) to eliminate fragile `eval(parse(text = ...))` string manipulations.
  * Fixed variable symbol evaluation issues when passing variables as function arguments rather than raw names.
  * Modularized `ipwtm` into dedicated family modeling engines (`binomial`, `survival`, `multinomial`, `ordinal`, `gaussian`).
  * Replaced manual survival hazard step-differencing with continuous baseline hazard interpolation using `survival::basehaz` and `stats::approx`.
  * Refactored `tstartfun` to use safe column evaluation and direct vector lag logic (`c(-1, t[-length(t)])`).

* **Test Suite Expansion**:
  * Upgraded `testthat` suite to 3rd edition specs, expanding coverage for `multinomial`, `ordinal`, `gaussian`, unstabilized models, weight truncation, and invalid parameter rejection.

# ipw 1.3.0

* New maintainer: Hung Thai Tran.
* License changed to MIT (earlier versions: GPL (>= 2)).
