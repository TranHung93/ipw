# Estimate Inverse Probability Weights (Point Treatment)

Estimates inverse probability weights to fit marginal structural models
in a point treatment situation. Exposure can be binomial, multinomial,
ordinal, a count (Poisson or negative binomial) or continuous (gaussian,
gamma or log-normal). Both stabilized and unstabilized weights can be
estimated, and for a binomial exposure also weights for the effect in
the treated (ATT), the effect in the untreated (ATC) and overlap weights
(ATO). For a multinomial or ordinal exposure, ATT weights (for a chosen
`focal` level) and generalized overlap weights are available as well.
Formulas may contain smooth terms such as `s(x)`, which are then fitted
as generalized additive models with mgcv.

## Usage

``` r
ipwpoint(
  exposure,
  family,
  link = NULL,
  numerator = NULL,
  denominator,
  data,
  trunc = NULL,
  estimand = c("ATE", "ATT", "ATC", "ATO"),
  focal = NULL,
  ...
)
```

## Arguments

- exposure:

  Name of the exposure column in `data` (unquoted or as a string).
  Binomial exposure variables should be coded as 0/1.

- family:

  Model family for exposure: `"binomial"`, `"multinomial"`, `"ordinal"`,
  `"gaussian"`, `"poisson"`, `"negbin"` (negative binomial, fitted with
  [`MASS::glm.nb`](https://rdrr.io/pkg/MASS/man/glm.nb.html)), `"gamma"`
  (gamma GLM with log link) or `"lognormal"` (normal model for the log
  of the exposure). Count exposures must be non-negative integers; gamma
  and log-normal exposures must be positive.

- link:

  Link function for binomial (`"logit"`, `"probit"`, `"cauchit"`,
  `"log"`, `"cloglog"`) or ordinal (`"logit"`, `"probit"`, `"cauchit"`,
  `"cloglog"`). Ignored for multinomial and gaussian.

- numerator:

  Formula specifying the right-hand side of the numerator model. If
  `NULL`, unstabilized weights (numerator = 1) are calculated. Required
  for the continuous families (`"gaussian"`, `"gamma"`, `"lognormal"`);
  must be `NULL` when `estimand` is `"ATT"`, `"ATC"` or `"ATO"`.

- denominator:

  Formula specifying the right-hand side of the denominator model.
  Smooth terms (`s()`, `te()`, `ti()`, `t2()`) are allowed for the
  binomial, gaussian and poisson families; the model is then fitted with
  [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html).

- data:

  Data frame containing exposure and model covariates.

- trunc:

  Optional truncation percentile between 0 and 0.5 (e.g., 0.01 for
  1st/99th percentiles).

- estimand:

  Target estimand. For a binomial exposure: `"ATE"` (default; inverse
  probability weights), `"ATT"` (effect in the treated: treated get
  weight 1, untreated \\p/(1-p)\\), `"ATC"` (effect in the untreated, or
  controls: treated get \\(1-p)/p\\, untreated 1) or `"ATO"` (overlap
  weights: treated \\1-p\\, untreated \\p\\), where \\p\\ is the fitted
  probability of exposure from the denominator model. For a multinomial
  or ordinal exposure: `"ATE"`, `"ATT"` (effect in the people at level
  `focal`: weight \\P_i(focal) / P_i(A_i)\\, 1 for people at the focal
  level) or `"ATO"` (generalized overlap weights, Li & Li 2019:
  \\(1/P_i(A_i)) / \sum_k 1/P_i(k)\\), where \\P_i(k)\\ is the fitted
  probability of level \\k\\ and \\A_i\\ the observed level. For two
  levels these equal the binomial ATT and ATO weights. `"ATC"` is only
  available for a binomial exposure; for a categorical exposure use
  `"ATT"` with `focal` set to the control level. Other families only
  allow `"ATE"`.

- focal:

  For `estimand = "ATT"` with `family = "multinomial"` or `"ordinal"`:
  the exposure level whose effect is targeted (the "treated" group),
  given as a level of the exposure. Required there and must be `NULL`
  otherwise; for a binomial exposure the focal group is always the
  exposed (1).

- ...:

  Additional arguments passed to model fitting functions (`glm`,
  `multinom`, `polr`, `glm.nb`, `gam`).

## Value

A list of class `"ipwpoint"` with attributes `"estimand"`, `"exposure"`
(the name of the exposure column) and, for ATT and ATC weights,
`"focal"` (the level whose effect is targeted: 1 for a binomial ATT, 0
for ATC), containing:

- ipw.weights:

  Numeric vector of inverse probability weights.

- weights.trunc:

  Numeric vector of truncated weights (if `trunc` is specified).

- call:

  The original function call.

- num.mod:

  Numerator model object (if `numerator` is specified).

- den.mod:

  Denominator model object.

## Details

Standard errors: `ipwpoint` only estimates the weights. Standard errors
for the marginal structural model have to be obtained when that model is
fitted.
[`ipwvcov`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
(experimental) gives standard errors that take the estimation of the
weights into account (M-estimation), and is the recommended option where
it applies: weights from `ipwpoint()` with `family = "binomial"` and
`link = "logit"` (any `estimand`, stabilized or not, no smooth terms, no
truncation) and a marginal structural model fitted with
[`lm()`](https://rdrr.io/r/stats/lm.html),
[`glm()`](https://rdrr.io/r/stats/glm.html) (any link) or
[`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) (right-censored
data, one row per person); and time-varying weights from
[`ipwtm`](https://tranhung93.github.io/ipw/reference/ipwtm.md) with
`family = "binomial"` and `link = "logit"` (`type = "first"`, `"cens"`
or `"all"`, alone or multiplied, e.g. treatment times censoring weights)
with an [`lm()`](https://rdrr.io/r/stats/lm.html),
[`glm()`](https://rdrr.io/r/stats/glm.html) or counting-process
[`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) model on the
long data. For all other cases use the bootstrap described below. Robust
(sandwich) standard errors, for example from
[`svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html) in the survey
package as in Van der Wal and Geskus (2011), treat the weights as known.
For a binary point treatment with good overlap they were conservative
(too large) in simulations, but they can be too small for a continuous
exposure (`family = "gaussian"`) and for time-varying weights; use the
bootstrap there. Austin (2016) compared naive, robust and bootstrap
standard errors in simulations and recommends the bootstrap, in which
the weights are re-estimated in each bootstrap sample. This takes the
uncertainty in the estimated weights into account; see the bootstrap
example below. About 200 bootstrap samples are usually enough to
estimate a standard error. In the package's simulation study (article
"Simulation study" on the package website), the bootstrap standard error
was accurate for a point treatment with good overlap, while the variance
from robust standard errors was 1.25 times the true variance (the
standard error 12% too large). When some weights are extreme (poor
overlap), no standard error was reliable, the bootstrap included. Check
the weights first with
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`ipwplot`](https://tranhung93.github.io/ipw/reference/ipwplot.md), and
consider truncation (`trunc`) or overlap weights (`estimand = "ATO"`).

## References

Cole, S.R. & Hernán, M.A. (2008). Constructing inverse probability
weights for marginal structural models. *American Journal of
Epidemiology*, **168**(6), 656-664.

Robins, J.M., Hernán, M.A. & Brumback, B.A. (2000). Marginal structural
models and causal inference in epidemiology. *Epidemiology*, **11**,
550-560.

Van der Wal W.M. & Geskus R.B. (2011). ipw: An R Package for Inverse
Probability Weighting. *Journal of Statistical Software*, **43**(13),
1-23. [doi:10.18637/jss.v043.i13](https://doi.org/10.18637/jss.v043.i13)
.

Li, F., Morgan, K.L. & Zaslavsky, A.M. (2018). Balancing covariates via
propensity score weighting. *Journal of the American Statistical
Association*, **113**, 390-400.
[doi:10.1080/01621459.2016.1260466](https://doi.org/10.1080/01621459.2016.1260466)
.

Greifer, N. WeightIt: Weighting for Covariate Balance in Observational
Studies. R package, <https://ngreifer.github.io/WeightIt/>. (Used in the
package tests as an independent comparison of the weights.)

Li, F. & Li, F. (2019). Propensity score weighting for causal inference
with multiple treatments. *Annals of Applied Statistics*, **13**(4),
2389-2415.
[doi:10.1214/19-AOAS1282](https://doi.org/10.1214/19-AOAS1282) .

Austin, P.C. (2016). Variance estimation when using inverse probability
of treatment weighting (IPTW) with survival analysis. *Statistics in
Medicine*, **35**, 5642-5655.
[doi:10.1002/sim.7084](https://doi.org/10.1002/sim.7084) .

## See also

[`ipwvcov`](https://tranhung93.github.io/ipw/reference/ipwvcov.md),
[`basdat`](https://tranhung93.github.io/ipw/reference/basdat.md),
[`haartdat`](https://tranhung93.github.io/ipw/reference/haartdat.md),
[`ipwplot`](https://tranhung93.github.io/ipw/reference/ipwplot.md),
`ipwpoint`,
[`ipwtm`](https://tranhung93.github.io/ipw/reference/ipwtm.md),
[`timedat`](https://tranhung93.github.io/ipw/reference/timedat.md),
[`tstartfun`](https://tranhung93.github.io/ipw/reference/tstartfun.md).

## Author

Willem M. van der Wal <willem@vanderwalresearch.com>, Ronald B. Geskus
<rgeskus@oucru.org>, Hung Thai Tran <hungtt@oucru.org>

## Examples

``` r
# Simulate data with continuous confounder and outcome, binomial exposure.
# Marginal causal effect of exposure on outcome: 10.
set.seed(16)
n <- 1000
simdat <- data.frame(l = rnorm(n, 10, 5))
a.lin <- simdat$l - 10
pa <- exp(a.lin)/(1 + exp(a.lin))
simdat$a <- rbinom(n, 1, prob = pa)
simdat$y <- 10*simdat$a + 0.5*simdat$l + rnorm(n, -10, 5)
simdat[1:5,]
#>           l a          y
#> 1 12.382067 1   6.635898
#> 2  9.373100 0 -11.722042
#> 3 15.481081 1  11.501612
#> 4  2.778855 0  -9.464074
#> 5 15.739146 1  15.085261

# Estimate ipw weights
temp <- ipwpoint(
  exposure = a,
  family = "binomial",
  link = "logit",
  numerator = ~ 1,
  denominator = ~ l,
  data = simdat)
summary(temp$ipw.weights)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.4810  0.5127  0.5285  0.9095  0.6318 74.6994 

# Plot inverse probability weights
ipwplot(weights = temp$ipw.weights, logscale = FALSE,
        main = "Stabilized weights", xlim = c(0, 8))


# Examine numerator and denominator models
summary(temp$num.mod)
#> 
#> Call:
#> (function (formula, family = gaussian, data, weights, subset, 
#>     na.action, start = NULL, etastart, mustart, offset, control = list(...), 
#>     model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, singular.ok = TRUE, 
#>     contrasts = NULL, ...) 
#> {
#>     cal <- match.call()
#>     if (is.character(family)) 
#>         family <- get(family, mode = "function", envir = parent.frame())
#>     if (is.function(family)) 
#>         family <- family()
#>     if (is.null(family$family)) {
#>         print(family)
#>         stop("'family' not recognized")
#>     }
#>     if (missing(data)) 
#>         data <- environment(formula)
#>     mf <- match.call(expand.dots = FALSE)
#>     m <- match(c("formula", "data", "subset", "weights", "na.action", 
#>         "etastart", "mustart", "offset"), names(mf), 0L)
#>     mf <- mf[c(1L, m)]
#>     mf$drop.unused.levels <- TRUE
#>     mf[[1L]] <- quote(stats::model.frame)
#>     mf <- eval(mf, parent.frame())
#>     if (identical(method, "model.frame")) 
#>         return(mf)
#>     if (!is.character(method) && !is.function(method)) 
#>         stop("invalid 'method' argument")
#>     if (identical(method, "glm.fit")) 
#>         control <- do.call("glm.control", control)
#>     mt <- attr(mf, "terms")
#>     Y <- model.response(mf, "any")
#>     if (length(dim(Y)) == 1L) {
#>         nm <- rownames(Y)
#>         dim(Y) <- NULL
#>         if (!is.null(nm)) 
#>             names(Y) <- nm
#>     }
#>     X <- if (!is.empty.model(mt)) 
#>         model.matrix(mt, mf, contrasts)
#>     else matrix(, NROW(Y), 0L)
#>     weights <- as.vector(model.weights(mf))
#>     if (!is.null(weights) && !is.numeric(weights)) 
#>         stop("'weights' must be a numeric vector")
#>     if (!is.null(weights) && any(weights < 0)) 
#>         stop("negative weights not allowed")
#>     offset <- as.vector(model.offset(mf))
#>     if (!is.null(offset)) {
#>         if (length(offset) != NROW(Y)) 
#>             stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
#>                 length(offset), NROW(Y)), domain = NA)
#>     }
#>     mustart <- model.extract(mf, "mustart")
#>     etastart <- model.extract(mf, "etastart")
#>     fit <- eval(call(if (is.function(method)) "method" else method, 
#>         x = X, y = Y, weights = weights, start = start, etastart = etastart, 
#>         mustart = mustart, offset = offset, family = family, 
#>         control = control, intercept = attr(mt, "intercept") > 
#>             0L, singular.ok = singular.ok))
#>     if (length(offset) && attr(mt, "intercept") > 0L) {
#>         fit2 <- eval(call(if (is.function(method)) "method" else method, 
#>             x = X[, "(Intercept)", drop = FALSE], y = Y, mustart = fit$fitted.values, 
#>             weights = weights, offset = offset, family = family, 
#>             control = control, intercept = TRUE))
#>         if (!fit2$converged) 
#>             warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?")
#>         fit$null.deviance <- fit2$deviance
#>     }
#>     if (model) 
#>         fit$model <- mf
#>     fit$na.action <- attr(mf, "na.action")
#>     if (x) 
#>         fit$x <- X
#>     if (!y) 
#>         fit$y <- NULL
#>     structure(c(fit, list(call = cal, formula = formula, terms = mt, 
#>         data = data, offset = offset, control = control, method = method, 
#>         contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, 
#>             mf))), class = c(fit$class, c("glm", "lm")))
#> })(formula = a ~ 1, family = binomial(link = "logit"), data = simdat, 
#>     na.action = function (object, ...) 
#>     UseMethod("na.fail"))
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)
#> (Intercept)  0.07604    0.06329   1.201     0.23
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 1384.9  on 999  degrees of freedom
#> Residual deviance: 1384.9  on 999  degrees of freedom
#> AIC: 1386.9
#> 
#> Number of Fisher Scoring iterations: 3
#> 
summary(temp$den.mod)
#> 
#> Call:
#> (function (formula, family = gaussian, data, weights, subset, 
#>     na.action, start = NULL, etastart, mustart, offset, control = list(...), 
#>     model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, singular.ok = TRUE, 
#>     contrasts = NULL, ...) 
#> {
#>     cal <- match.call()
#>     if (is.character(family)) 
#>         family <- get(family, mode = "function", envir = parent.frame())
#>     if (is.function(family)) 
#>         family <- family()
#>     if (is.null(family$family)) {
#>         print(family)
#>         stop("'family' not recognized")
#>     }
#>     if (missing(data)) 
#>         data <- environment(formula)
#>     mf <- match.call(expand.dots = FALSE)
#>     m <- match(c("formula", "data", "subset", "weights", "na.action", 
#>         "etastart", "mustart", "offset"), names(mf), 0L)
#>     mf <- mf[c(1L, m)]
#>     mf$drop.unused.levels <- TRUE
#>     mf[[1L]] <- quote(stats::model.frame)
#>     mf <- eval(mf, parent.frame())
#>     if (identical(method, "model.frame")) 
#>         return(mf)
#>     if (!is.character(method) && !is.function(method)) 
#>         stop("invalid 'method' argument")
#>     if (identical(method, "glm.fit")) 
#>         control <- do.call("glm.control", control)
#>     mt <- attr(mf, "terms")
#>     Y <- model.response(mf, "any")
#>     if (length(dim(Y)) == 1L) {
#>         nm <- rownames(Y)
#>         dim(Y) <- NULL
#>         if (!is.null(nm)) 
#>             names(Y) <- nm
#>     }
#>     X <- if (!is.empty.model(mt)) 
#>         model.matrix(mt, mf, contrasts)
#>     else matrix(, NROW(Y), 0L)
#>     weights <- as.vector(model.weights(mf))
#>     if (!is.null(weights) && !is.numeric(weights)) 
#>         stop("'weights' must be a numeric vector")
#>     if (!is.null(weights) && any(weights < 0)) 
#>         stop("negative weights not allowed")
#>     offset <- as.vector(model.offset(mf))
#>     if (!is.null(offset)) {
#>         if (length(offset) != NROW(Y)) 
#>             stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
#>                 length(offset), NROW(Y)), domain = NA)
#>     }
#>     mustart <- model.extract(mf, "mustart")
#>     etastart <- model.extract(mf, "etastart")
#>     fit <- eval(call(if (is.function(method)) "method" else method, 
#>         x = X, y = Y, weights = weights, start = start, etastart = etastart, 
#>         mustart = mustart, offset = offset, family = family, 
#>         control = control, intercept = attr(mt, "intercept") > 
#>             0L, singular.ok = singular.ok))
#>     if (length(offset) && attr(mt, "intercept") > 0L) {
#>         fit2 <- eval(call(if (is.function(method)) "method" else method, 
#>             x = X[, "(Intercept)", drop = FALSE], y = Y, mustart = fit$fitted.values, 
#>             weights = weights, offset = offset, family = family, 
#>             control = control, intercept = TRUE))
#>         if (!fit2$converged) 
#>             warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?")
#>         fit$null.deviance <- fit2$deviance
#>     }
#>     if (model) 
#>         fit$model <- mf
#>     fit$na.action <- attr(mf, "na.action")
#>     if (x) 
#>         fit$x <- X
#>     if (!y) 
#>         fit$y <- NULL
#>     structure(c(fit, list(call = cal, formula = formula, terms = mt, 
#>         data = data, offset = offset, control = control, method = method, 
#>         contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, 
#>             mf))), class = c(fit$class, c("glm", "lm")))
#> })(formula = a ~ l, family = binomial(link = "logit"), data = simdat, 
#>     na.action = function (object, ...) 
#>     UseMethod("na.fail"))
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -9.69809    0.66066  -14.68   <2e-16 ***
#> l            0.97272    0.06525   14.91   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 1384.85  on 999  degrees of freedom
#> Residual deviance:  496.86  on 998  degrees of freedom
#> AIC: 500.86
#> 
#> Number of Fisher Scoring iterations: 7
#> 

# Paste inverse probability weights
simdat$sw <- temp$ipw.weights

# Marginal structural model for the causal effect of a on y
# corrected for confounding by l using inverse probability weighting
if (requireNamespace("survey", quietly = TRUE)) {
  library(survey)
  msm <- svyglm(y ~ a,
                design = svydesign(~1, weights = ~temp$ipw.weights,
                data = simdat))
  summary(msm)
}
#> Loading required package: grid
#> Loading required package: Matrix
#> Loading required package: survival
#> 
#> Attaching package: ‘survey’
#> The following object is masked from ‘package:graphics’:
#> 
#>     dotchart
#> 
#> Call:
#> svyglm(formula = y ~ a, design = svydesign(~1, weights = ~temp$ipw.weights, 
#>     data = simdat))
#> 
#> Survey design:
#> svydesign(~1, weights = ~temp$ipw.weights, data = simdat)
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)   -4.375      1.142  -3.832 0.000135 ***
#> a             10.647      1.190   8.948  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for gaussian family taken to be 29.58889)
#> 
#> Number of Fisher Scoring iterations: 2
#> 

# Check the weights before trusting any standard error. These weights are
# extreme (poor overlap): ipwdiag(temp) shows a tail index below 2. For
# standard errors that take the estimation of the weights into account
# (ipwvcov(), M-estimation), see the example on well-behaved weights in
# ?ipwvcov.

# Bootstrap standard error and percentile confidence interval for the MSM.
# The robust standard error from svyglm() treats the weights as known and
# can be too large or too small (simulation study, scenarios 7a and 7b).
# Austin (2016) recommends the bootstrap, with the
# weights re-estimated in each bootstrap sample: each bootstrap sample below
# re-fits ipwpoint() and the MSM. (A weighted lm() gives the same point
# estimate as svyglm() and is faster.)
if (requireNamespace("boot", quietly = TRUE)) {
  boot.fun <- function(data, index) {
    d <- data[index, ]
    w <- ipwpoint(
      exposure = a,
      family = "binomial",
      link = "logit",
      numerator = ~ 1,
      denominator = ~ l,
      data = d)
    d$sw <- w$ipw.weights
    coef(lm(y ~ a, data = d, weights = sw))["a"]
  }
  set.seed(1)
  bootres <- boot::boot(data = simdat, statistic = boot.fun, R = 200)
  print(bootres)
  # bootstrap standard error
  print(sd(bootres$t[, 1]))
  # percentile 95 percent confidence interval
  print(boot::boot.ci(bootres, type = "perc"))
}
#> 
#> ORDINARY NONPARAMETRIC BOOTSTRAP
#> 
#> 
#> Call:
#> boot::boot(data = simdat, statistic = boot.fun, R = 200)
#> 
#> 
#> Bootstrap Statistics :
#>     original    bias    std. error
#> t1* 10.64661 0.2146065    1.115705
#> [1] 1.115705
#> BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
#> Based on 200 bootstrap replicates
#> 
#> CALL : 
#> boot::boot.ci(boot.out = bootres, type = "perc")
#> 
#> Intervals : 
#> Level     Percentile     
#> 95%   ( 8.43, 12.79 )  
#> Calculations and Intervals on Original Scale
#> Some percentile intervals may be unstable
```
