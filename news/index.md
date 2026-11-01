# Changelog

## ipw (development version)

### New features

- [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md):
  `type = "all"` now works for `family = "multinomial"` and `"ordinal"`
  (time-varying categorical exposures that switch back and forth).

- [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md):
  `type = "all"` now works for `family = "survival"`, following Graffeo,
  Latouche, Geskus & Chevret (2018, *Biometrical Journal*): exposure
  events may recur, an Andersen-Gill Cox model is fitted on all rows,
  and each row’s probability uses the cumulative baseline hazard over
  its own interval.

- New
  [`ipwprob()`](https://tranhung93.github.io/ipw/reference/ipwprob.md):
  inverse probability weights from probabilities you estimated yourself
  (for example with machine learning or SuperLearner), for point
  treatments or cumulated over time within subjects.

- [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md):
  new families for count exposures (`"poisson"`, `"negbin"`) and skewed
  continuous exposures (`"gamma"`, `"lognormal"`).

- Smooth terms (`s()`, `te()`, `ti()`, `t2()`) in model formulas are
  fitted as generalized additive models with mgcv: in
  [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  for the binomial, gaussian and poisson families, and in
  [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) for
  the binomial family.

- [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md):
  new `estimand` argument for a binomial exposure: `"ATE"` (default),
  `"ATT"` (effect in the treated) and `"ATO"` (overlap weights).

- Existing calls are unchanged: without smooth terms the models are
  still fitted with [`glm()`](https://rdrr.io/r/stats/glm.html), and
  `estimand` defaults to `"ATE"`. A covariate named `s`, `te`, `ti` or
  `t2` is treated as an ordinary variable.

- New
  [`ipwboot()`](https://tranhung93.github.io/ipw/reference/ipwboot.md):
  bootstrap standard errors and percentile confidence intervals for a
  whole weighted analysis, re-estimating the weights in each bootstrap
  sample (Austin 2016). With `id`, whole subjects are resampled and each
  sampled copy gets a new identifier; optional `strata`; failed samples
  are recorded instead of stopping.
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`confint()`](https://rdrr.io/r/stats/confint.html) methods.

### Documentation

- New vignette, “Inverse probability weighting with ipw”
  ([`vignette("ipw")`](https://tranhung93.github.io/ipw/articles/ipw.md)),
  covering every function, with the new models and
  [`ipwboot()`](https://tranhung93.github.io/ipw/reference/ipwboot.md)
  highlighted.
- Bootstrap standard errors. The help pages of
  [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  and [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md)
  have a new “Standard errors” paragraph in Details and new examples: a
  bootstrap that re-estimates the weights and the marginal structural
  model in each bootstrap sample
  ([`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)),
  and a cluster bootstrap that resamples patients and re-runs both
  [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) calls
  for the HAART example
  ([`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md)).
  Robust standard errors
  ([`svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html),
  [`cluster()`](https://rdrr.io/pkg/survival/man/cluster.html)) treat
  the weights as known and tend to be conservative; Austin (2016,
  <doi:10.1002/sim.7084>) recommends the bootstrap.

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
