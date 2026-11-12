# Changelog

## ipw (development version)

### Improvements

- [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md),
  `family = "survival"`, `type = "first"` or `"cens"` no longer needs
  the data on a full time grid. Up to ipw 1.3.2, the probability of each
  row used the increase in the baseline hazard since the previous time
  point in the whole data set, so the weights were correct only when all
  follow-up was split at all times that occur in the data, as in the
  examples of the JSS paper, `haartdat` and Example 2 of
  [`?ipwtm`](https://tranhung93.github.io/ipw/reference/ipwtm.md). Each
  row now uses the baseline hazard over its own interval
  `(tstart, timevar]`, so data with individual visit times can be used
  directly. **For data on the full grid the weights are unchanged**; for
  data that are not, they now differ from earlier versions, which gave
  incorrect weights for such data. `type = "all"` always used each row’s
  own interval.

### Bug fixes

- [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  and [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md):
  [`offset()`](https://rdrr.io/r/stats/offset.html) terms in the
  `numerator` and `denominator` formulas were silently dropped, so the
  models were fitted without the offset and **the weights were wrong**.
  This affected ipw 1.3.1 and 1.3.2 on CRAN; ipw 1.2.1.1 was fine.
  Formulas are now used as written.
- [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md):
  `denominator = ~ .` (all other columns of `data`) failed with
  “Variables not found in data: .” and works again, as in ipw 1.2.1.1.
  In [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) a
  `.` is now a clear error, because it would include the id, time and
  exposure columns; list the covariates instead.
- [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  and [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md),
  `family = "ordinal"`: an unordered factor exposure whose levels are
  not numbers now gives a warning, because its levels are used in
  alphabetical order unless set otherwise (e.g. “high” \< “low” \<
  “mid”). Use an ordered factor to make the order explicit.
- [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md),
  `family = "survival"`: a `ties` argument passed through `...` to
  [`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) was silently
  ignored (Efron was always used), and `method` gave the error “formal
  argument matched by multiple actual arguments”. Both are now passed
  on; the default is still `ties = "efron"`.
- [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md),
  `family = "survival"`: rows with `tstart >= timevar` now give a clear
  error instead of “missing values in object” or no error at all.
- Clearer errors:
  [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  with `family = "gaussian"` and a non-numeric exposure (was “NA/NaN/Inf
  in ‘y’”);
  [`ipwplot()`](https://tranhung93.github.io/ipw/reference/ipwplot.md)
  with missing weights (was “‘x’ contains missing values”) or a
  `binwidth` of length \> 1.
- [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  and [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md):
  a binomial (or survival) 0/1 exposure stored as a factor or character
  is now used by its labels. With factor levels in the order
  `c("1", "0")` the probabilities of the two exposure values were
  swapped, giving wrong weights.
- [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  and [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md),
  `family = "ordinal"`: a character exposure is now an error. It was
  ordered alphabetically (e.g. “high” \< “low” \< “mid”). Use an ordered
  factor or numeric codes.
- [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md): a
  covariate called `selvar` was replaced by the internal selection
  indicator; it is now used as the user’s covariate.
- [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md):
  duplicated combinations of `id` and `timevar` are now an error (as in
  [`ipwprobtm()`](https://tranhung93.github.io/ipw/reference/ipwprobtm.md)).
- [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md):
  an exposure that takes only one value is now an error instead of
  returning meaningless weights.
- [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  and [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md):
  infinite values in model variables now give an error that names the
  variables.
- [`ipwprobtm()`](https://tranhung93.github.io/ipw/reference/ipwprobtm.md):
  `den` and `num` larger than 1 are now an error.
- [`tstartfun()`](https://tranhung93.github.io/ipw/reference/tstartfun.md):
  warns when follow-up times are \<= -1, because the starting time of -1
  given to the first record is then not before the follow-up time.
- [`ipwplot()`](https://tranhung93.github.io/ipw/reference/ipwplot.md):
  zero or negative weights with `logscale = TRUE` are now an error; they
  were silently left out of the plot.
- [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md)
  examples: the cluster bootstrap moved to
  [`vignette("ipw")`](https://tranhung93.github.io/ipw/articles/ipw.md)
  (section 5), and the part of Example 2 that needs nlme (in Suggests)
  only runs when nlme is installed.

### Documentation

- [`vignette("ipw")`](https://tranhung93.github.io/ipw/articles/ipw.md),
  section 6: covariate balance with
  [`cobalt::bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.html)
  and `love.plot()` (ipw weights passed as a vector), for the binary
  point treatment and for the HAART example at two time points, and the
  Kish effective sample size overall and by exposure group.

- New section “ipw or WeightIt?” in the README (with a comparison table)
  and a short version near the start of
  [`vignette("ipw")`](https://tranhung93.github.io/ipw/articles/ipw.md):
  when to use which package.

- New article “ipw, WeightIt and PSweight side by side”: the same binary
  and categorical point treatments in all three packages (weights and
  standard errors), and longitudinal data with individual visit times in
  [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) and
  `weightitMSM()`. The README section “ipw or WeightIt?” also names
  PSweight as an option for point treatments.

- [`?ipwpoint`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  and
  [`vignette("ipw")`](https://tranhung93.github.io/ipw/articles/ipw.md),
  section 5:
  [`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
  is the recommended standard error for a binary point treatment; the
  bootstrap remains the advice otherwise. The methods article has new
  sections on ATC, categorical ATT/ATO weights and
  [`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md).

- New article “Bugs fixed after ipw 1.3.2” on the package website: for
  each bug above, an example run in ipw 1.2.1.1, ipw 1.3.2 and this
  version, with the weights compared with weights computed by hand.

- New article “Survival weights and the time grid”: why
  [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md)
  survival weights in ipw 1.3.2 and earlier need the data on a full time
  grid (as in the JSS paper), a worked example, `rows_affected()` to
  check whether your data were on the full grid, the results of every
  ipw version on CRAN (1.0 to 1.3.2), a simulation with a known answer,
  the old and the new code, and the tests.

### New features

- New
  [`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
  (**new, experimental; still under discussion**): the
  variance-covariance matrix of the coefficients of a marginal
  structural model fitted with [`lm()`](https://rdrr.io/r/stats/lm.html)
  or [`glm()`](https://rdrr.io/r/stats/glm.html) and
  [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  weights for a binary exposure, taking the estimation of the weights
  into account (M-estimation with stacked estimating equations,
  Lunceford & Davidian 2004). Works for `estimand = "ATE"` (stabilized
  or not), `"ATT"`, `"ATC"` and `"ATO"`, with logistic exposure models
  (`link = "logit"`). Not for other links, truncated weights, smooth
  terms, other families or
  [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md)
  weights: use the bootstrap there. Agrees with
  `WeightIt::glm_weightit(vcov = "asympt")` in the tests.

- [`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
  now also accepts a
  [`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) marginal
  structural model with right-censored data (`Surv(time, status)`, one
  row per person) and
  [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  weights: the weighted partial likelihood score is stacked with the
  exposure models (Binder 1992; Lin & Wei 1989). Tied event times need
  `ties = "breslow"`; counting-process data,
  [`strata()`](https://rdrr.io/pkg/survival/man/strata.html),
  [`cluster()`](https://rdrr.io/pkg/survival/man/cluster.html),
  [`frailty()`](https://rdrr.io/pkg/survival/man/frailty.html) and
  `tt()` are not supported. Checked against a numerical Jacobian of the
  stacked equations and a coverage simulation. Results for
  [`lm()`](https://rdrr.io/r/stats/lm.html) and
  [`glm()`](https://rdrr.io/r/stats/glm.html) models are unchanged.

- [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md):
  `estimand = "ATC"` (effect in the untreated) for a binomial exposure:
  treated get weight (1 - p)/p, untreated 1.

- [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md):
  `estimand = "ATT"` and `"ATO"` for `family = "multinomial"` and
  `"ordinal"`. The new argument `focal` gives the level whose effect the
  ATT targets: weight P(focal) / P(observed level), 1 at the focal
  level. `"ATO"` gives generalized overlap weights (Li & Li 2019). For
  two levels both equal the binomial weights. For the effect in the
  controls of a categorical exposure use `"ATT"` with `focal` set to the
  control level.

- [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  results now carry the attributes `"estimand"`, `"exposure"` and (for
  ATT and ATC) `"focal"`. The list elements are unchanged.

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
  [`ipwprobtm()`](https://tranhung93.github.io/ipw/reference/ipwprobtm.md):
  time-varying inverse probability weights from probabilities you
  estimated yourself (for example with machine learning or
  SuperLearner): the ratios are multiplied over time within each
  subject, in time order. For a point treatment no function is needed:
  the weights are `num / den`.

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

### Documentation

- New vignette, “Inverse probability weighting with ipw”
  ([`vignette("ipw")`](https://tranhung93.github.io/ipw/articles/ipw.md)),
  covering every function, with the new models highlighted. Its
  standard-errors section describes, from the package’s simulation
  study, when the bootstrap helps (well-behaved weights) and when no
  standard error is reliable (extreme weights).
- Bootstrap standard errors. The help pages of
  [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  and [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md)
  have a new “Standard errors” paragraph in Details.
  [`?ipwpoint`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  has a bootstrap example that re-estimates the weights and the marginal
  structural model in each bootstrap sample; the cluster bootstrap that
  resamples patients and re-runs both
  [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) calls
  for the HAART example is in section 5 of
  [`vignette("ipw")`](https://tranhung93.github.io/ipw/articles/ipw.md).
  Robust standard errors
  ([`svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html),
  [`cluster()`](https://rdrr.io/pkg/survival/man/cluster.html)) treat
  the weights as known. They were conservative for a binary point
  treatment with good overlap, but can be too small for a continuous
  exposure and for time-varying weights; Austin (2016,
  <doi:10.1002/sim.7084>) recommends the bootstrap.
- [`ipwplot()`](https://tranhung93.github.io/ipw/reference/ipwplot.md):
  the help page now documents the return value (`NULL`, invisibly).
- The [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md)
  examples no longer call
  [`graphics.off()`](https://rdrr.io/r/grDevices/dev.html), which closed
  all of the user’s open graphics devices.

### Tests

- New tests against PSweight (in Suggests) as a second independent
  check:
  [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  weights for a binary exposure (ATE, ATT, ATC, ATO) and a three-level
  multinomial exposure (ATE, ATT with a focal level, overlap weights),
  and
  [`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
  standard errors against PSweight’s sandwich standard errors of the
  difference in means.

### License

- License changed back to GPL (\>= 2), the license of ipw before version
  1.3.0.

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
