# Changelog

## ipw 1.4.0

### New features

- [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md):
  new `estimand` argument. For a binomial exposure: `"ATE"` (default),
  `"ATT"` (effect in the treated), `"ATC"` (effect in the untreated) and
  `"ATO"` (overlap weights). For `family = "multinomial"` and
  `"ordinal"`: `"ATT"`, with the new argument `focal` naming the level
  whose effect is targeted, and `"ATO"` (generalized overlap weights; Li
  & Li 2019). For the effect in the controls of a categorical exposure,
  use `"ATT"` with `focal` set to the control level.
- New
  [`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
  (**experimental**): the variance-covariance matrix of the coefficients
  of a marginal structural model, taking the estimation of the weights
  into account (M-estimation with stacked estimating equations).
  Supported:
  - [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
    weights for a binary exposure with a logistic model (any
    `estimand`), with an [`lm()`](https://rdrr.io/r/stats/lm.html),
    [`glm()`](https://rdrr.io/r/stats/glm.html) or
    [`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) model
    (right-censored data, one row per person);
  - [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md)
    weights with `family = "binomial"`, `link = "logit"` and any `type`,
    or a list of
    [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md)/[`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
    results whose weights multiply (e.g. treatment times censoring
    weights), with an [`lm()`](https://rdrr.io/r/stats/lm.html),
    [`glm()`](https://rdrr.io/r/stats/glm.html) or counting-process
    [`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) model on
    the long data. Scores are summed within subject.

  Other cases (other families or links, truncated weights, smooth terms)
  give an error; use the bootstrap there. Checked against numerical
  derivatives of the estimating equations,
  [`WeightIt::glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.html)
  and PSweight. With heavy-tailed weights no normal-theory interval is
  reliable, whichever method gives the standard error (see
  [`?ipwvcov`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)).
- New
  [`ipwprobtm()`](https://tranhung93.github.io/ipw/reference/ipwprobtm.md):
  time-varying weights from probabilities estimated elsewhere (for
  example with machine learning), multiplied over time within each
  subject. For a point treatment the weights are simply `num / den`.
- [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md):
  new families for count exposures (`"poisson"`, `"negbin"`) and skewed
  continuous exposures (`"gamma"`, `"lognormal"`).
- [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md):
  `type = "all"` (exposures that switch back and forth) now works for
  `family = "multinomial"`, `"ordinal"` and `"survival"` (recurrent
  exposure events, Graffeo et al. 2018).
- Smooth terms (`s()`, `te()`, `ti()`, `t2()`) in model formulas are
  fitted with [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html): in
  [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  for the binomial, gaussian and poisson families, and in
  [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) for
  the binomial family. Without smooth terms the models are still fitted
  with [`glm()`](https://rdrr.io/r/stats/glm.html).
- [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  results carry the attributes `"estimand"`, `"exposure"` and (for ATT
  and ATC) `"focal"`.
  [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md)
  results carry the attributes `"family"`, `"type"`, `"id"` and
  `"sort.order"`. The list elements and weights are unchanged.

### Changes

- **Survival weights and the time grid (caution for older versions).**
  [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) with
  `family = "survival"` and `type = "first"` or `"cens"` now computes
  each row’s probability from the baseline hazard over that row’s own
  interval `(tstart, timevar]`. ipw 1.3.2 and earlier used the increase
  in the baseline hazard since the previous time point in the whole data
  set, which is right only when every subject’s follow-up is split at
  all times that occur in the data (the “full time grid”, as in the JSS
  paper, `haartdat` and Example 2 of
  [`?ipwtm`](https://tranhung93.github.io/ipw/reference/ipwtm.md)).
  **For data on the full grid the weights are unchanged.** Data with
  individual visit times can now be used directly; if you used such data
  with an older version, re-run the analysis with this version. The
  article “Survival weights and the time grid” shows how to check
  whether your data were on the full grid.
- License changed back to GPL (\>= 2), the license of ipw before 1.3.0
  (1.3.0 to 1.3.2 were MIT).
- Two functions that were only in development versions on GitHub, never
  on CRAN, changed: `ipwprob()` was renamed
  [`ipwprobtm()`](https://tranhung93.github.io/ipw/reference/ipwprobtm.md)
  and now handles time-varying exposures only, and `ipwboot()` was
  removed (the bootstrap is shown in
  [`?ipwpoint`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  and
  [`vignette("ipw")`](https://tranhung93.github.io/ipw/articles/ipw.md)
  instead).
- Existing calls give the same weights as before (except for the
  time-grid change above), checked by regression tests against ipw
  1.2.1.1. A covariate named `s`, `te`, `ti` or `t2` is treated as an
  ordinary variable.

### Bug fixes

- [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  and [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md):
  [`offset()`](https://rdrr.io/r/stats/offset.html) terms in the
  `numerator` and `denominator` formulas were silently dropped, so **the
  weights were wrong**. This affected ipw 1.3.1 and 1.3.2; ipw 1.2.1.1
  was fine.
- [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  and [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md):
  a binomial (or survival) 0/1 exposure stored as a factor or character
  is now used by its labels. With factor levels in the order
  `c("1", "0")` the probabilities of the two exposure values were
  swapped, giving **wrong weights**.
- [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md):
  `denominator = ~ .` works again, as in ipw 1.2.1.1. In
  [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) a `.`
  is an error, because it would include the id, time and exposure
  columns.
- [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md),
  `family = "survival"`: `ties` passed through `...` to
  [`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) was ignored
  (Efron was always used); it is now passed on.
- [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md): a
  covariate called `selvar` was replaced by the internal selection
  indicator; it is now used as the user’s covariate.
- `family = "ordinal"`: a character exposure is now an error, and an
  unordered factor whose levels are not numbers gives a warning, because
  the levels would be ordered alphabetically. Use an ordered factor.
- Clearer errors instead of wrong results or obscure messages:
  duplicated `id`/`timevar` combinations in
  [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md),
  `tstart >= timevar`, an exposure that takes one value, infinite values
  in model variables, a non-numeric gaussian exposure, `den`/`num` above
  1 in
  [`ipwprobtm()`](https://tranhung93.github.io/ipw/reference/ipwprobtm.md),
  missing or non-positive weights in
  [`ipwplot()`](https://tranhung93.github.io/ipw/reference/ipwplot.md).
  [`tstartfun()`](https://tranhung93.github.io/ipw/reference/tstartfun.md)
  warns when follow-up times are \<= -1.
- The [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md)
  examples no longer call
  [`graphics.off()`](https://rdrr.io/r/grDevices/dev.html), which closed
  all of the user’s graphics devices.

### Documentation

- New vignette, “Inverse probability weighting with ipw”
  ([`vignette("ipw")`](https://tranhung93.github.io/ipw/articles/ipw.md)),
  covering every function, with sections on standard errors (bootstrap,
  [`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)),
  covariate balance with cobalt and the effective sample size.
- “Standard errors” paragraphs in
  [`?ipwpoint`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  and [`?ipwtm`](https://tranhung93.github.io/ipw/reference/ipwtm.md):
  when the bootstrap helps, when robust standard errors that treat the
  weights as known can be too small, and which cases
  [`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
  supports.
- New “ipw or WeightIt?” section in the README and the vignette.
- New articles on the package website
  (<https://tranhung93.github.io/ipw/>): “Methods” (formulas and code
  for the new features), “Simulation study”, “Bugs fixed after ipw
  1.3.2”, “Survival weights and the time grid” and “ipw, WeightIt and
  PSweight side by side”.

### Tests

- Weights and standard errors are compared with WeightIt and PSweight
  (both in Suggests), with numerical derivatives of the estimating
  equations, and with weights computed by hand; regression tests against
  ipw 1.2.1.1.

## ipw 1.3.2

CRAN release: 2026-11-01

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
