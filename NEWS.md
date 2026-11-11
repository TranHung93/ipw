# ipw (development version)

## Bug fixes

* `ipwtm()`, `family = "survival"`, `type = "first"` or `"cens"`: the
  probability of each row now uses the baseline hazard over the row's own
  interval `(tstart, timevar]`. Before, it used the increase in the baseline
  hazard since the previous time point in the whole data set, which is only
  right when all subjects share one time grid (as in `haartdat`). **Weights
  change for data in which subjects have their own follow-up times**; for data
  with a common time grid they are unchanged. `type = "all"` was not affected.
* `ipwpoint()` and `ipwtm()`: a binomial (or survival) 0/1 exposure stored as
  a factor or character is now used by its labels. With factor levels in the
  order `c("1", "0")` the probabilities of the two exposure values were
  swapped, giving wrong weights.
* `ipwpoint()` and `ipwtm()`, `family = "ordinal"`: a character exposure is
  now an error. It was ordered alphabetically (e.g. "high" < "low" < "mid").
  Use an ordered factor or numeric codes.
* `ipwtm()`: a covariate called `selvar` was replaced by the internal
  selection indicator; it is now used as the user's covariate.
* `ipwtm()`: duplicated combinations of `id` and `timevar` are now an error
  (as in `ipwprobtm()`).
* `ipwpoint()`: an exposure that takes only one value is now an error
  instead of returning meaningless weights.
* `ipwpoint()` and `ipwtm()`: infinite values in model variables now give an
  error that names the variables.
* `ipwprobtm()`: `den` and `num` larger than 1 are now an error.
* `tstartfun()`: warns when follow-up times are <= -1, because the starting
  time of -1 given to the first record is then not before the follow-up time.
* `ipwplot()`: zero or negative weights with `logscale = TRUE` are now an
  error; they were silently left out of the plot.
* `ipwtm()` examples: the cluster bootstrap moved to `vignette("ipw")`
  (section 5), and the part of Example 2 that needs nlme (in Suggests) only
  runs when nlme is installed.

## New features

* `ipwtm()`: `type = "all"` now works for `family = "multinomial"` and
  `"ordinal"` (time-varying categorical exposures that switch back and forth).
* `ipwtm()`: `type = "all"` now works for `family = "survival"`, following
  Graffeo, Latouche, Geskus & Chevret (2018, *Biometrical Journal*): exposure
  events may recur, an Andersen-Gill Cox model is fitted on all rows, and each
  row's probability uses the cumulative baseline hazard over its own interval.
* New `ipwprobtm()`: time-varying inverse probability weights from
  probabilities you estimated yourself (for example with machine learning or
  SuperLearner): the ratios are multiplied over time within each subject, in
  time order. For a point treatment no function is needed: the weights are
  `num / den`.
* `ipwpoint()`: new families for count exposures (`"poisson"`, `"negbin"`) and
  skewed continuous exposures (`"gamma"`, `"lognormal"`).
* Smooth terms (`s()`, `te()`, `ti()`, `t2()`) in model formulas are fitted as
  generalized additive models with mgcv: in `ipwpoint()` for the binomial,
  gaussian and poisson families, and in `ipwtm()` for the binomial family.
* `ipwpoint()`: new `estimand` argument for a binomial exposure: `"ATE"`
  (default), `"ATT"` (effect in the treated) and `"ATO"` (overlap weights).
* Existing calls are unchanged: without smooth terms the models are still
  fitted with `glm()`, and `estimand` defaults to `"ATE"`. A covariate named
  `s`, `te`, `ti` or `t2` is treated as an ordinary variable.

## Documentation

* New vignette, "Inverse probability weighting with ipw"
  (`vignette("ipw")`), covering every function, with the new models
  highlighted. Its standard-errors section describes, from the package's
  simulation study, when the bootstrap helps (well-behaved weights) and when
  no standard error is reliable (extreme weights).
* Bootstrap standard errors. The help pages of `ipwpoint()` and
  `ipwtm()` have a new "Standard errors" paragraph in Details and new examples:
  a bootstrap that re-estimates the weights and the marginal structural model
  in each bootstrap sample (`ipwpoint()`), and a cluster bootstrap that
  resamples patients and re-runs both `ipwtm()` calls for the HAART example
  (`ipwtm()`). Robust standard errors (`svyglm()`, `cluster()`) treat the
  weights as known and tend to be conservative; Austin (2016,
  doi:10.1002/sim.7084) recommends the bootstrap.

## License

* License changed back to GPL (>= 2), the license of ipw before version 1.3.0.

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
