# Standard Errors that Account for Estimated Weights (M-estimation)

**New, experimental.** Variance-covariance matrix of the coefficients of
a marginal structural model fitted with
[`lm()`](https://rdrr.io/r/stats/lm.html),
[`glm()`](https://rdrr.io/r/stats/glm.html) or
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) and
weights from
[`ipwpoint`](https://tranhung93.github.io/ipw/reference/ipwpoint.md) for
a binary point treatment, or with time-varying weights from
[`ipwtm`](https://tranhung93.github.io/ipw/reference/ipwtm.md), taking
the estimation of the weights into account. The marginal structural
model and the exposure models are treated as one set of stacked
estimating equations (M-estimation; Lunceford & Davidian 2004), and the
sandwich variance of all parameters is computed; the block for the
coefficients of the marginal structural model is returned.

## Usage

``` r
ipwvcov(object, msm)
```

## Arguments

- object:

  Result of
  [`ipwpoint`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  with `family = "binomial"` and `link = "logit"`; or result of
  [`ipwtm`](https://tranhung93.github.io/ipw/reference/ipwtm.md) with
  `family = "binomial"` and `link = "logit"`; or a list of
  [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) (and
  [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md))
  results whose weights multiply to the weights of `msm`.

- msm:

  Marginal structural model: an `lm`, `glm` or `coxph` fit with weights
  equal to `object$ipw.weights` (or their product for a list), on the
  same rows in the same order as the data given to
  [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  or [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md).

## Value

Variance-covariance matrix of the coefficients of `msm`, with the
coefficient names as row and column names. It can be passed, for
example, to `lmtest::coeftest(msm, vcov. = ...)`.

## Details

Robust (sandwich) standard errors of a weighted model, such as those
from [`svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html) or
[`sandwich::vcovHC()`](https://zeileis.codeberg.page/sandwich/reference/vcovHC.html),
treat the weights as known. Here the estimating equations of the
denominator (and numerator) model are stacked under those of the
marginal structural model. With \\\theta = (\beta, \alpha, \gamma)\\ the
coefficients of the marginal structural model, the denominator model and
the numerator model, and \\\psi_i(\theta)\\ the contribution of person
\\i\\ to the stacked estimating equations, the variance is \\A^{-1} B
A^{-T}\\ with \\A = -\sum_i \partial \psi_i / \partial \theta\\ and \\B
= \sum_i \psi_i \psi_i^T\\, evaluated at the estimates. The derivative
of each weight with respect to the exposure-model coefficients is worked
out for each estimand (ATE, ATT, ATC, ATO). The exposure models must use
the logit link: the information matrices use the expected (Fisher)
information, which equals the observed information, and so gives the
exact sandwich, only for canonical links. For the marginal structural
model any [`glm()`](https://rdrr.io/r/stats/glm.html) link is accepted.
With a non-canonical link (for example a log link for risk ratios) and a
saturated model such as `y ~ a` for the binary exposure, the expected
and observed information are exactly equal in any sample: the fitted
mean in each exposure group is that group's weighted mean of `y`, so the
extra term of the observed information, a weighted sum of residuals
times a factor that is constant within each group, sums to zero within
each group. For other models with a non-canonical link the two are
asymptotically equivalent when the model is correctly specified. No
small-sample correction is applied. Proportional rescaling of the
weights does not change the result.

**Cox models.** For a
[`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) marginal
structural model the weighted partial likelihood score is used as the
estimating equation (Binder 1992; Lin & Wei 1989), with the Breslow form
for tied event times. Person \\i\\'s contribution is \\w_i\\ times the
score residual, \\U_i = w_i \[\delta_i (X_i - \bar X(t_i)) - \sum_k
\delta_k w_k Y_i(t_k) e^{\beta' X_i} (X_i - \bar X(t_k)) / S_0(t_k)\]\\,
with \\S_0(t) = \sum_j w_j Y_j(t) e^{\beta' X_j}\\ and \\\bar X(t)\\ the
weighted mean of \\X\\ in the risk set; `residuals(msm, type = "score")`
returns the part in brackets (without the weight). The bread for
\\\beta\\ is the weighted information matrix (the inverse of
`msm$naive.var`). The weights enter both each person's own term and the
risk-set sums; differentiating the total score with respect to one
weight gives exactly that person's score residual, so the cross block is
\\-\sum_i U_i d_i^T\\ with \\d_i = \partial \log w_i / \partial
\alpha\\, as for [`glm()`](https://rdrr.io/r/stats/glm.html). It is
computed analytically; the package tests check it against central finite
differences of the stacked estimating equations. Without the exposure
models, the sandwich reduces to
[`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)'s robust
variance (`msm$var`), which treats the weights as known.

Supported: `object` from `ipwpoint(family = "binomial", link = "logit")`
with denominator (and, if present, numerator) models fitted by
[`glm()`](https://rdrr.io/r/stats/glm.html) (no smooth terms) and any
`estimand`; `msm` of class `"lm"` or `"glm"`, fitted on the same rows,
in the same order, as the weights, with `weights = object$ipw.weights`;
or `msm` from [`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)
with right-censored data `Surv(time, status)`, one row per person in the
same order, and the same weights. For
[`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html), tied event
times need `ties = "breslow"`; without tied event times all tie methods
give the same result.

**Time-varying weights.** `object` can be an
[`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) result,
or a list of
[`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) (and
[`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md))
results whose `ipw.weights` multiply to the weights of `msm`, typically
treatment weights times censoring weights. Each model is a separate
parameter block. The unit of independence is the subject (the `id` of
[`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md)): all
scores are summed within subject before the meat \\B = \sum\_{subjects}
\psi_p \psi_p^T\\ is formed. The weight of row \\r\\ is the product,
over the rows \\s \le r\\ of the same subject, of numerator over
denominator factors \\f_s\\, as defined by
[`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md):
\\P(A_s = a_s)\\ for `type = "first"` and `"all"`, \\1 - p_s\\ for
`type = "cens"` (the censoring row included), and 1 after the first
switch. With logistic models, \\\partial \log f_s / \partial \alpha\\ is
\\(a_s - p_s) Z_s\\ (`"first"`, `"all"`) or \\-p_s Z_s\\ (`"cens"`), so
\\d_r = \partial \log w_r / \partial \alpha\\ is minus the cumulative
sum of these over the rows used by the model, within subject in time
order (plus, for the numerator). The cross block is \\-\sum_r U_r
d_r^T\\ for [`lm()`](https://rdrr.io/r/stats/lm.html),
[`glm()`](https://rdrr.io/r/stats/glm.html) and Cox models; for a Cox
model on counting-process data `Surv(tstart, time, event)` the identity
above holds row by row, because each row enters the risk sets with its
own weight. Without the exposure models the sandwich reduces to the
cluster-robust variance that treats the weights as known (for a Cox
model, [`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) with
`cluster(id)`).

Supported for time-varying weights:
`ipwtm(family = "binomial", link = "logit")` with any `type`, glm models
without smooth terms, no truncation, made with this version of ipw;
`msm` fitted with [`lm()`](https://rdrr.io/r/stats/lm.html),
[`glm()`](https://rdrr.io/r/stats/glm.html) or
[`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)
(counting-process data;
[`cluster()`](https://rdrr.io/pkg/survival/man/cluster.html) is allowed,
`ties = "breslow"` when event times are tied) on all rows of the long
data, in the same order, with weights equal to the product of the
weights in `object`.

**Heavy-tailed weights.** Cumulative weights can have a heavy right tail
when some people have a very low probability of their observed exposure
history (a near-violation of positivity). Normal-theory intervals are
then unreliable whatever method gives the standard error, the bootstrap
included. In simulations where treatment depended strongly on the
confounder and on past treatment, the tail index of the final cumulative
weights was about 2.3 (the weights have a finite variance, but not a
finite third or fourth moment). 95% intervals then covered 92 to 93% of
the time with `ipwvcov()` and 90 to 95% with the cluster bootstrap, with
no trend towards 95% between 300 and 4800 people. With a tail index of
about 3.9 both covered about 95%. Check the weight distribution
([`summary()`](https://rdrr.io/r/base/summary.html),
[`ipwplot`](https://tranhung93.github.io/ipw/reference/ipwplot.md) and
the effective sample size `sum(w)^2 / sum(w^2)`) before trusting any
standard error.

Not supported, and an error: other exposure families or links, smooth
terms (GAMs), truncated weights (truncation at estimated percentiles is
not a smooth function of the model coefficients),
[`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) weights
from Cox (`family = "survival"`), multinomial, ordinal or gaussian
models, and [`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)
models with [`strata()`](https://rdrr.io/pkg/survival/man/strata.html),
[`frailty()`](https://rdrr.io/pkg/survival/man/frailty.html) or `tt()`
terms, or Efron or exact ties with tied event times. With
[`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
weights the [`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)
model must have right-censored data (one row per person) and no
[`cluster()`](https://rdrr.io/pkg/survival/man/cluster.html); with
[`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) weights
it must have counting-process data. Use the bootstrap there,
re-estimating the weights in each bootstrap sample (see
[`ipwpoint`](https://tranhung93.github.io/ipw/reference/ipwpoint.md) and
[`vignette("ipw")`](https://tranhung93.github.io/ipw/articles/ipw.md);
for long data, the cluster bootstrap that resamples subjects).

This is the same kind of standard error as `vcov = "asympt"` in
[`WeightIt::glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.html);
in the package tests the two agree for binary treatments. For Cox
models, `ipwvcov()` agrees with a numerical Jacobian of the stacked
estimating equations to about \\10^{-7}\\.
`WeightIt::coxph_weightit(vcov = "asympt")` gives slightly different
values: for the ATE and overlap weights the relative difference was
below 0.1% and shrank with the sample size; for ATT and ATC weights its
standard error was 1 to 4% larger and the difference did not shrink. Its
`vcov = "HC0"` equals
[`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)'s robust
variance exactly.

## References

Lunceford, J.K. & Davidian, M. (2004). Stratification and weighting via
the propensity score in estimation of causal treatment effects: a
comparative study. *Statistics in Medicine*, **23**, 2937-2960.
[doi:10.1002/sim.1903](https://doi.org/10.1002/sim.1903) .

Binder, D.A. (1992). Fitting Cox's proportional hazards models from
survey data. *Biometrika*, **79**, 139-147.
[doi:10.1093/biomet/79.1.139](https://doi.org/10.1093/biomet/79.1.139) .

Lin, D.Y. & Wei, L.J. (1989). The robust inference for the Cox
proportional hazards model. *Journal of the American Statistical
Association*, **84**, 1074-1078.
[doi:10.1080/01621459.1989.10478874](https://doi.org/10.1080/01621459.1989.10478874)
.

Stefanski, L.A. & Boos, D.D. (2002). The calculus of M-estimation. *The
American Statistician*, **56**, 29-38.
[doi:10.1198/000313002753631330](https://doi.org/10.1198/000313002753631330)
.

Greifer, N. WeightIt: Weighting for Covariate Balance in Observational
Studies. R package, <https://ngreifer.github.io/WeightIt/>. (Its
`glm_weightit(vcov = "asympt")` is used in the package tests as an
independent comparison.)

## See also

[`ipwpoint`](https://tranhung93.github.io/ipw/reference/ipwpoint.md),
[`ipwtm`](https://tranhung93.github.io/ipw/reference/ipwtm.md)

## Examples

``` r
set.seed(16)
n <- 1000
simdat <- data.frame(l = rnorm(n, 10, 5))
simdat$a <- rbinom(n, 1, plogis(simdat$l - 10))
simdat$y <- 10 * simdat$a + 0.5 * simdat$l + rnorm(n, -10, 5)

w <- ipwpoint(exposure = a, family = "binomial", link = "logit",
              numerator = ~ 1, denominator = ~ l, data = simdat)
msm <- lm(y ~ a, data = simdat, weights = w$ipw.weights)

# standard errors that account for the estimated weights
sqrt(diag(ipwvcov(w, msm)))
#> (Intercept)           a 
#>    1.052322    1.086986 

# overlap weights
w.ato <- ipwpoint(exposure = a, family = "binomial", link = "logit",
                  denominator = ~ l, data = simdat, estimand = "ATO")
msm.ato <- lm(y ~ a, data = simdat, weights = w.ato$ipw.weights)
sqrt(diag(ipwvcov(w.ato, msm.ato)))
#> (Intercept)           a 
#>   0.4309341   0.6040691 

# Wald tests and confidence intervals
if (requireNamespace("lmtest", quietly = TRUE)) {
  print(lmtest::coeftest(msm, vcov. = ipwvcov(w, msm)))
  print(lmtest::coefci(msm, vcov. = ipwvcov(w, msm)))
}
#> 
#> t test of coefficients:
#> 
#>             Estimate Std. Error t value  Pr(>|t|)    
#> (Intercept)  -4.3755     1.0523 -4.1579 3.487e-05 ***
#> a            10.6466     1.0870  9.7946 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>                 2.5 %   97.5 %
#> (Intercept) -6.440496 -2.31046
#> a            8.513571 12.77965

# Cox marginal structural model for a survival outcome
library(survival)
simdat$time <- rexp(n, 0.1 * exp(-0.5 * simdat$a + 0.05 * simdat$l))
simdat$status <- as.numeric(simdat$time < 10)
simdat$time <- pmin(simdat$time, 10)
msm.cox <- coxph(Surv(time, status) ~ a, data = simdat, weights = w$ipw.weights)
cbind(robust = sqrt(diag(msm.cox$var)),          # weights treated as known
      ipwvcov = sqrt(diag(ipwvcov(w, msm.cox))))  # weights estimated
#>     robust   ipwvcov
#> a 0.184859 0.1646684

# Time-varying weights: treatment (type = "first") times censoring
# (type = "cens") weights from pooled logistic models, and a Cox MSM on
# the long data
data(haartdat)
w.haart <- ipwtm(exposure = haartind, family = "binomial", link = "logit",
                 numerator = ~ sex + age, denominator = ~ sex + age + cd4.sqrt,
                 id = patient, timevar = fuptime, type = "first", data = haartdat)
w.cens <- ipwtm(exposure = dropout, family = "binomial", link = "logit",
                numerator = ~ sex + age, denominator = ~ sex + age + cd4.sqrt,
                id = patient, timevar = fuptime, type = "cens", data = haartdat)
msm.haart <- coxph(Surv(tstart, fuptime, event) ~ haartind + cluster(patient),
                   data = haartdat, ties = "breslow",
                   weights = w.haart$ipw.weights * w.cens$ipw.weights)
cbind(robust.cluster = sqrt(diag(msm.haart$var)),
      ipwvcov = sqrt(diag(ipwvcov(list(w.haart, w.cens), msm.haart))))
#>          robust.cluster   ipwvcov
#> haartind      0.4633263 0.4585489
```
