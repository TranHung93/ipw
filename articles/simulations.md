# Simulation study: do the weights recover known effects?

The new weights in ipw recover known causal effects when their model is
correct, and the naive or misspecified alternatives are clearly biased.
This article shows that with a simulation study in the spirit of the
simulated example in the JSS paper (van der Wal & Geskus 2011): each
scenario generates data in which the true marginal causal effect is
known, estimates it, and repeats this many times.

**How the results were produced.** The study is in
[`dev/simulation-study.R`](https://github.com/TranHung93/ipw/blob/main/dev/simulation-study.R)
in the repository; the code below shows the essential lines of each
scenario. It takes several minutes, so this page shows saved results
rather than running it. To reproduce them, run
`Rscript dev/simulation-study.R` in the package folder (seed 20261101),
or start the `simulations` workflow from the GitHub Actions tab.
Scenarios 1 to 4 use 500 replications, 5 and 6 use 250, and 7 uses 100
with 200 bootstrap samples each.

**Reading the tables.** *Bias* is the mean estimate minus the truth;
*empirical SD* is the spread of the estimates over replications; *RMSE*
combines both. A good method has bias close to 0 and a small RMSE.

![](simulations_files/figure-html/bias-plot-1.png)

## 1. ATE, ATT and overlap weights, and `ipwprob()`

The effect of a binary exposure varies with the confounder `l`, so the
three estimands differ. Their true values, computed from 2 million
simulated people, are ATE = 1, ATT = 1.276 and ATO = 0.861. Each
weighting scheme is compared with its own target;
[`ipwprob()`](https://tranhung93.github.io/ipw/reference/ipwprob.md) is
given the true propensity scores.

``` r

d <- data.frame(l = rnorm(1000))
d$a <- rbinom(1000, 1, plogis(0.5 + 0.8 * d$l))
d$y <- d$a * (1 + d$l) + d$l + rnorm(1000)
w.att <- ipwpoint(exposure = a, family = "binomial", link = "logit",
                  denominator = ~ l, data = d, estimand = "ATT")$ipw.weights
coef(lm(y ~ a, data = d, weights = w.att))["a"]
```

| Method | Truth | Mean estimate | Bias | Empirical SD | RMSE | Replications |
|:---|---:|---:|---:|---:|---:|---:|
| unweighted (target ATE) | 1.000 | 1.989 | 0.989 | 0.114 | 0.995 | 500 |
| ATE weights | 1.000 | 1.002 | 0.002 | 0.089 | 0.089 | 500 |
| ATT weights (target ATT) | 1.276 | 1.275 | -0.001 | 0.107 | 0.107 | 500 |
| ATO weights (target ATO) | 0.861 | 0.864 | 0.002 | 0.082 | 0.082 | 500 |
| ipwprob, true probabilities (ATE) | 1.000 | 1.006 | 0.007 | 0.144 | 0.144 | 500 |

All weighted estimates are unbiased for their own estimand; the
unweighted difference is nearly twice the ATE.

## 2. Count exposure: Poisson or negative binomial

The exposure is an overdispersed count (negative binomial, variance
larger than the mean); each extra unit raises the outcome by 0.5.

``` r

d <- data.frame(l = rnorm(1000))
d$v <- rnbinom(1000, size = 2, mu = exp(0.5 + 0.4 * d$l))
d$y <- 0.5 * d$v + d$l + rnorm(1000)
w.nb <- ipwpoint(exposure = v, family = "negbin", numerator = ~ 1,
                 denominator = ~ l, data = d)$ipw.weights
```

| Method | Truth | Mean estimate | Bias | Empirical SD | RMSE | Replications |
|:---|---:|---:|---:|---:|---:|---:|
| unweighted | 0.5 | 0.671 | 0.171 | 0.022 | 0.172 | 500 |
| poisson weights (misspecified) | 0.5 | 0.329 | -0.171 | 0.078 | 0.188 | 500 |
| negbin weights (correct) | 0.5 | 0.507 | 0.007 | 0.038 | 0.039 | 500 |

The negative binomial weights are unbiased. The Poisson model ignores
the overdispersion, and its weights over-correct: the bias is as large
as that of the unweighted analysis, in the other direction.

## 3. Skewed continuous exposure: gaussian, gamma or log-normal

The dose is log-normal given `l`; one unit of dose raises the outcome by
1.

``` r

d <- data.frame(l = rnorm(2000))
d$dose <- exp(0.5 + 0.5 * d$l + rnorm(2000, 0, 0.6))
d$y <- d$dose + 2 * d$l + rnorm(2000)
w.ln <- ipwpoint(exposure = dose, family = "lognormal", numerator = ~ 1,
                 denominator = ~ l, data = d)$ipw.weights
```

| Method | Truth | Mean estimate | Bias | Empirical SD | RMSE | Replications |
|:---|---:|---:|---:|---:|---:|---:|
| unweighted | 1 | 1.534 | 0.534 | 0.042 | 0.536 | 500 |
| gaussian weights (misspecified) | 1 | 1.130 | 0.130 | 0.085 | 0.156 | 491 |
| gamma weights (misspecified) | 1 | 0.983 | -0.017 | 0.165 | 0.166 | 500 |
| lognormal weights (correct) | 1 | 1.061 | 0.061 | 0.087 | 0.106 | 500 |
| lognormal weights, truncated 1% | 1 | 1.179 | 0.179 | 0.036 | 0.183 | 500 |

The log-normal weights remove most of the confounding, with a small
remaining bias (0.06) that is typical of weights for continuous
exposures: the density ratio has a heavy tail. The gaussian model is
clearly biased, and in 9 of 500 replications gave infinite weights,
which
[`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
refuses. The gamma model happens to have little bias here but twice the
spread. Truncation at the 1st and 99th percentiles reduces the spread a
lot but adds bias: the usual trade-off.

## 4. Curved propensity: smooth terms with mgcv

The probability of exposure depends on `l` through a sine curve; the
true effect is 2.

``` r

d <- data.frame(l = runif(2000, -2, 2))
d$a <- rbinom(2000, 1, plogis(2 * sin(1.5 * d$l)))
d$y <- 2 * d$a + 3 * sin(1.5 * d$l) + rnorm(2000)
w.gam <- ipwpoint(exposure = a, family = "binomial", link = "logit",
                  numerator = ~ 1, denominator = ~ s(l), data = d)$ipw.weights
```

| Method | Truth | Mean estimate | Bias | Empirical SD | RMSE | Replications |
|:---|---:|---:|---:|---:|---:|---:|
| unweighted | 2 | 4.547 | 2.547 | 0.093 | 2.548 | 500 |
| linear logistic weights (misspecified) | 2 | 2.659 | 0.659 | 0.086 | 0.665 | 500 |
| s(l) weights | 2 | 2.070 | 0.070 | 0.064 | 0.095 | 500 |

A straight-line logistic model leaves about a quarter of the confounding
bias in place; `s(l)` lets mgcv estimate the curve and removes almost
all of it.

## 5. Time-varying ordinal dose, `type = "all"`

A dose level (0, 1, 2) is chosen at each of 4 visits from a
proportional-odds model given the confounder `L` and the previous dose.
`L` is affected by the previous dose and by an unmeasured `U`, which
also affects the outcome. The outcome rises by 1 per unit of cumulative
dose. This is the situation in the JSS paper where ordinary adjustment
fails: `L` is both a confounder and affected by earlier treatment.

``` r

for (t in 1:4) {
  l <- U + 0.4 * prev + rnorm(n, 0, 0.5)
  dose <- findInterval(0.8 * l + 0.6 * prev + rlogis(n), c(0, 1.5))   # 0, 1, 2
  prev <- dose
}
y <- 1 * cumulative.dose + 2 * U + rnorm(n)

w <- ipwtm(exposure = dose, family = "ordinal", link = "logit",
           numerator = ~ dose.lag, denominator = ~ dose.lag + l,
           id = id, timevar = time, type = "all", data = long)
# weight at the last visit, then lm(y ~ cumulative.dose, weights = ...)
```

| Method | Truth | Mean estimate | Bias | Empirical SD | RMSE | Replications |
|:---|---:|---:|---:|---:|---:|---:|
| unweighted | 1 | 1.492 | 0.492 | 0.035 | 0.493 | 250 |
| regression adjustment for L | 1 | 0.855 | -0.145 | 0.028 | 0.147 | 250 |
| ordinal weights (correct) | 1 | 1.027 | 0.027 | 0.110 | 0.113 | 250 |
| multinomial weights | 1 | 1.010 | 0.010 | 0.116 | 0.117 | 250 |

The unweighted analysis is confounded, and adjusting for `L` at every
visit is biased in the other direction, because it blocks part of the
effect and opens a path through `U`. Both ordinal weights (the correct
model) and multinomial weights recover the effect.

## 6. Recurrent exposure events, `family = "survival"`, `type = "all"`

Exposure events (for example new prescriptions) happen repeatedly with a
known proportional hazard: baseline 2 per unit time, log hazard ratio
0.5 for a fixed covariate `x` and -1 for a time-varying covariate `xt`.
Because the true hazard is known, the true weights are known too, so the
check is on the Cox coefficients and on how close the estimated weights
are to the true ones. In 6a all people share a grid of 10 intervals (as
in `haartdat`); in 6b the events happen in continuous time, `xt` changes
at different times for each person, and each event ends an interval
(start-stop data cut at events, as in Graffeo et al. 2018).

``` r

w <- ipwtm(exposure = expo, family = "survival", denominator = ~ x + xt,
           id = id, tstart = tstart, timevar = fuptime, type = "all", data = d)
coef(w$den.mod)
median(abs(log(w$ipw.weights / true.weights)))
```

| Method | Truth | Mean estimate | Bias | Empirical SD | RMSE | Replications |
|:---|---:|---:|---:|---:|---:|---:|
| Cox coefficient of x | 0.5 | 0.499 | -0.001 | 0.084 | 0.083 | 250 |
| Cox coefficient of xt | -1.0 | -0.983 | 0.017 | 0.099 | 0.100 | 250 |
| median \|log(estimated / true weight)\| | 0.0 | 0.069 | 0.069 | 0.020 | 0.072 | 250 |
| Cox coefficient of x | 0.5 | 0.500 | 0.000 | 0.091 | 0.091 | 250 |
| Cox coefficient of xt | -1.0 | -0.996 | 0.004 | 0.093 | 0.093 | 250 |
| median \|log(estimated / true weight)\| | 0.0 | 0.091 | 0.091 | 0.024 | 0.094 | 250 |

The Cox coefficients are unbiased in both designs, and the estimated
weights are typically within 7 to 9% of the true weights, which is the
estimation noise of the baseline hazard.

## 7. Standard errors: `ipwboot()` versus robust standard errors

The marginal structural model for a binary exposure (true effect 10, n =
500) is fitted with stabilized weights; 95% confidence intervals come
from
[`ipwboot()`](https://tranhung93.github.io/ipw/reference/ipwboot.md)
(200 bootstrap samples, re-estimating the weights) or from robust
standard errors
([`survey::svyglm`](https://rdrr.io/pkg/survey/man/svyglm.html)). In 7a,
exposure probabilities are mostly between 0.1 and 0.9 (good overlap); 7b
is the simulated example of the JSS paper, where the logit of the
propensity has a standard deviation of 5, so many people are almost
certain to be exposed or unexposed.

``` r

msm.effect <- function(dd) {
  w <- ipwpoint(exposure = a, family = "binomial", link = "logit", numerator = ~ 1,
                denominator = ~ l, data = dd)$ipw.weights
  coef(lm(y ~ a, data = dd, weights = w))["a"]
}
ipwboot(d, msm.effect, R = 200)
```

| Scenario | Method | Mean estimate | Mean SE | Empirical SD | Coverage of 95% CI | Replications |
|:---|:---|---:|---:|---:|---:|---:|
| 7a. good overlap | ipwboot, percentile CI | 10.054 | 0.492 | 0.457 | 0.95 | 100 |
| 7a. good overlap | robust SE (svyglm), Wald CI | 10.054 | 0.539 | 0.457 | 0.95 | 100 |
| 7b. JSS example (poor overlap) | ipwboot, percentile CI | 11.323 | 1.072 | 1.632 | 0.50 | 100 |
| 7b. JSS example (poor overlap) | robust SE (svyglm), Wald CI | 11.323 | 1.060 | 1.632 | 0.53 | 100 |

With good overlap, both intervals reach their nominal 95% coverage; the
robust standard error is about 18% larger than the true spread of the
estimates and the bootstrap about 8%, in line with Austin (2016). In the
JSS example, the estimate itself is biased at this sample size (11.3
instead of 10) and both standard errors are far too small, so coverage
is only about 50%. Near-violations of positivity, not the choice of
standard error, are the problem there: check the weights
([`summary()`](https://rdrr.io/r/base/summary.html),
[`ipwplot()`](https://tranhung93.github.io/ipw/reference/ipwplot.md))
and consider overlap weights or truncation.

## References

- Austin PC (2016). Variance estimation when using inverse probability
  of treatment weighting (IPTW) with survival analysis. *Statistics in
  Medicine* 35, 5642-5655.
- Graffeo N, Latouche A, Geskus RB, Chevret S (2018). Modeling
  time-varying exposure using inverse probability of treatment weights.
  *Biometrical Journal* 60, 323-332.
- Li F, Morgan KL, Zaslavsky AM (2018). Balancing covariates via
  propensity score weighting. *Journal of the American Statistical
  Association* 113, 390-400.
- van der Wal WM, Geskus RB (2011). ipw: An R package for inverse
  probability weighting. *Journal of Statistical Software* 43(13), 1-23.
