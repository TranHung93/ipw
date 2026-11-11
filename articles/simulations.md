# Simulation study: do the weights recover known effects?

This article checks, by simulation, whether the new features of ipw do
what they should: each scenario generates data in which the true causal
effect is known, estimates it, and repeats this hundreds of times, in
the spirit of the simulated example in the JSS paper (van der Wal &
Geskus 2011).

**Main findings**

- **The new weights remove confounding when their model is correct.**
  Naive analyses and weights from misspecified models are biased by 13%
  to over 100%; the new weights reduce this to under 6% in every
  scenario but one (8% for smooth terms in
  [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md),
  scenario 5b).
- **The remaining small biases (1 to 6%) are properties of inverse
  probability weighting, not of the implementation.** Weights computed
  from the true model (*oracle* weights) have the same bias for
  log-normal and gamma doses and for counts, and for smooth terms and
  ordinal exposures the bias shrinks as the sample grows.
- **A bootstrap that re-estimates the weights gives accurate standard
  errors when the weights are well behaved**, where robust standard
  errors are 10 to 12% too large. It does **not** fix analyses with
  extreme weights: with poor overlap or long products of time-varying
  weights, no standard error reaches 95% coverage, and the bootstrap is
  then not better than the robust standard error.

**How the results were produced.** The study is in
[`dev/simulation-study.R`](https://github.com/TranHung93/ipw/blob/main/dev/simulation-study.R);
the code below shows the essential lines of each scenario. It takes
about 20 minutes on 4 cores, so this page shows saved results. To
reproduce them, run `Rscript dev/simulation-study.R` in the package
folder (seed 20261101), or start the `simulations` workflow from the
GitHub Actions tab. Each replication has its own L’Ecuyer random number
stream, so the results are the same for any number of parallel workers
and on Windows (where the replications run one after another); only
different versions of R or of the packages used could change them
slightly. Replications: 500 for scenarios 1 to 4, 250 for 5 and 6, 200
for the large-sample and oracle checks, 1000 for 7a and 7b, and 500 for
7c and 7d, each with 200 bootstrap samples.

**Reading the tables.** *Bias* is the mean estimate minus the truth. The
*MCSE* (Monte Carlo standard error) shows how precisely the simulation
estimates that bias: a bias within about 2 MCSE of 0 could be simulation
noise. *Empirical SD* is the spread of the estimates; *RMSE* combines
bias and spread.

*Oracle* weights are computed from the true model that generated the
exposure, which is known only in a simulation. They show what weighting
achieves without any error in estimating that model: if oracle weights
have the same bias as the estimated weights, the bias comes from
weighting itself, not from the model.

![](simulations_files/figure-html/bias-plot-1.png)

## 1. ATE, ATT and overlap weights

The effect of a binary exposure varies with the confounder `l`, so the
three estimands differ. Their true values, from 2 million simulated
people, are ATE = 1.00, ATT = 1.27 and ATO = 0.86. Each weighting scheme
is compared with its own target. Weights from the true propensity scores
(1 / probability of the observed exposure) are oracle weights.

``` r

d <- data.frame(l = rnorm(1000))
d$a <- rbinom(1000, 1, plogis(0.5 + 0.8 * d$l))
d$y <- d$a * (1 + d$l) + d$l + rnorm(1000)
w.att <- ipwpoint(exposure = a, family = "binomial", link = "logit",
                  denominator = ~ l, data = d, estimand = "ATT")$ipw.weights
coef(lm(y ~ a, data = d, weights = w.att))["a"]
```

| Method | Truth | Mean estimate | Bias | MCSE of bias | Empirical SD | RMSE | Replications |
|:---|---:|---:|---:|---:|---:|---:|---:|
| unweighted (target ATE) | 0.999 | 1.983 | 0.985 | 0.005 | 0.114 | 0.991 | 500 |
| ATE weights | 0.999 | 1.002 | 0.003 | 0.004 | 0.089 | 0.089 | 500 |
| ATT weights (target ATT) | 1.275 | 1.281 | 0.006 | 0.005 | 0.106 | 0.106 | 500 |
| ATO weights (target ATO) | 0.861 | 0.864 | 0.003 | 0.004 | 0.079 | 0.079 | 500 |
| true probabilities, 1 / den (ATE) | 0.999 | 1.007 | 0.008 | 0.007 | 0.148 | 0.148 | 500 |

All weighted estimates are within 1% of their own estimand; the
unweighted difference is twice the ATE.

## 2. Skewed continuous exposures: gaussian, gamma or log-normal

One unit of dose raises the outcome by 1. In 2a the dose is log-normal
given `l`, in 2b gamma. *Oracle* weights divide the true marginal
density by the true conditional density. For the log-normal dose the
marginal has a closed form; for the gamma dose it is the average of the
true conditional density over the sample’s values of `l`. In 2a they are
also computed for a 5 times larger sample.

``` r

d <- data.frame(l = rnorm(2000))
d$dose <- exp(0.5 + 0.5 * d$l + rnorm(2000, 0, 0.6))   # 2b: rgamma(2000, 3, 3 / exp(0.5 + 0.3 * d$l))
d$y <- d$dose + 2 * d$l + rnorm(2000)
w.ln <- ipwpoint(exposure = dose, family = "lognormal", numerator = ~ 1,
                 denominator = ~ l, data = d)$ipw.weights
oracle <- dlnorm(d$dose, 0.5, sqrt(0.61)) / dlnorm(d$dose, 0.5 + 0.5 * d$l, 0.6)
# 2b: f <- function(x, l) dgamma(x, shape = 3, rate = 3 / exp(0.5 + 0.3 * l))
#     oracle <- rowMeans(outer(d$dose, d$l, f)) / f(d$dose, d$l)
```

### 2a. Log-normal dose

| Method | Truth | Mean estimate | Bias | MCSE of bias | Empirical SD | RMSE | Replications |
|:---|---:|---:|---:|---:|---:|---:|---:|
| unweighted | 1 | 1.535 | 0.535 | 0.002 | 0.044 | 0.537 | 500 |
| gaussian weights (misspecified) | 1 | 1.130 | 0.130 | 0.006 | 0.134 | 0.186 | 496 |
| gamma weights (misspecified) | 1 | 0.981 | -0.019 | 0.007 | 0.154 | 0.155 | 500 |
| lognormal weights (correct), n = 2000 | 1 | 1.055 | 0.055 | 0.005 | 0.107 | 0.120 | 500 |
| lognormal weights (correct), n = 10000 | 1 | 1.022 | 0.022 | 0.007 | 0.092 | 0.095 | 200 |
| lognormal weights, truncated 1% | 1 | 1.181 | 0.181 | 0.002 | 0.038 | 0.185 | 500 |
| oracle weights (true densities), n = 2000 | 1 | 1.056 | 0.056 | 0.005 | 0.113 | 0.126 | 500 |
| oracle weights (true densities), n = 10000 | 1 | 1.022 | 0.022 | 0.006 | 0.091 | 0.094 | 200 |

Log-normal weights remove most of the confounding (unweighted bias 54%).
The 5% that remains is **the same for oracle weights** (5.5% vs 5.6% at
n = 2000, 2.2% vs 2.2% at n = 10000), so it does not come from
estimating the model: it is the finite-sample bias of weighting with
densities of a skewed exposure, whose density ratios have a heavy tail,
and it shrinks slowly as the sample grows. The gaussian model is clearly
worse, and in 4 of 500 replications gave infinite weights, which
[`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
refuses. Gamma weights, the wrong skewed family here, have a small bias
but about 1.5 times the spread. Truncation reduces the spread a lot but
adds bias.

### 2b. Gamma dose

| Method | Truth | Mean estimate | Bias | MCSE of bias | Empirical SD | RMSE | Replications |
|:---|---:|---:|---:|---:|---:|---:|---:|
| unweighted | 1 | 1.758 | 0.758 | 0.002 | 0.038 | 0.759 | 500 |
| gaussian weights (misspecified) | 1 | 1.148 | 0.148 | 0.003 | 0.075 | 0.166 | 500 |
| gamma weights (correct) | 1 | 1.043 | 0.043 | 0.004 | 0.098 | 0.107 | 500 |
| lognormal weights (misspecified) | 1 | 1.025 | 0.025 | 0.009 | 0.204 | 0.205 | 500 |
| oracle weights (true densities) | 1 | 1.049 | 0.049 | 0.005 | 0.117 | 0.127 | 500 |

Gamma weights reduce the unweighted bias of 76% to 4.3%; gaussian
weights leave 15%. As in 2a, oracle weights have a similar bias (4.9%),
so it comes from weighting with a skewed exposure, not from estimating
the gamma model. Log-normal weights, the wrong skewed family here, again
have a small bias but twice the spread.

## 3. Count exposures: Poisson or negative binomial

Each extra unit of a count exposure raises the outcome by 0.5. In 3a the
counts are overdispersed (negative binomial); in 3b they are Poisson.
Oracle weights are computed as for the gamma dose, with the true
probability of each count.

``` r

d <- data.frame(l = rnorm(1000))
d$v <- rnbinom(1000, size = 2, mu = exp(0.5 + 0.4 * d$l))   # 3b: rpois(1000, exp(0.5 + 0.4 * d$l))
d$y <- 0.5 * d$v + d$l + rnorm(1000)
w.nb <- ipwpoint(exposure = v, family = "negbin", numerator = ~ 1,
                 denominator = ~ l, data = d)$ipw.weights
f <- function(v, l) dnbinom(v, size = 2, mu = exp(0.5 + 0.4 * l))   # 3b: dpois(v, exp(0.5 + 0.4 * l))
oracle <- rowMeans(outer(d$v, d$l, f)) / f(d$v, d$l)
```

### 3a. Overdispersed counts

| Method | Truth | Mean estimate | Bias | MCSE of bias | Empirical SD | RMSE | Replications |
|:---|---:|---:|---:|---:|---:|---:|---:|
| unweighted | 0.5 | 0.670 | 0.170 | 0.001 | 0.022 | 0.171 | 500 |
| poisson weights (misspecified) | 0.5 | 0.328 | -0.172 | 0.003 | 0.077 | 0.189 | 500 |
| negbin weights (correct) | 0.5 | 0.506 | 0.006 | 0.002 | 0.035 | 0.036 | 500 |
| oracle weights (true probabilities) | 0.5 | 0.508 | 0.008 | 0.002 | 0.038 | 0.039 | 500 |

Negative binomial weights, the correct model, are within 1.5% of the
truth, and oracle weights have a similar small bias (1.6%). Poisson
weights over-correct as strongly as the unweighted analysis
under-corrects.

### 3b. Poisson counts

| Method | Truth | Mean estimate | Bias | MCSE of bias | Empirical SD | RMSE | Replications |
|:---|---:|---:|---:|---:|---:|---:|---:|
| unweighted | 0.5 | 0.806 | 0.306 | 0.001 | 0.028 | 0.307 | 500 |
| poisson weights (correct) | 0.5 | 0.505 | 0.005 | 0.003 | 0.062 | 0.062 | 500 |
| negbin weights (also valid: theta large) | 0.5 | 0.515 | 0.015 | 0.003 | 0.063 | 0.065 | 500 |
| oracle weights (true probabilities) | 0.5 | 0.511 | 0.011 | 0.003 | 0.070 | 0.071 | 500 |

Poisson weights, now the correct model, are within 1.5% of the truth.
Oracle weights do no better (bias 2.2%, with a larger spread): weights
from an estimated model are often more precise than weights from the
true model (Hirano, Imbens & Ridder 2003), so the estimated weights lose
nothing. The negative binomial model remains usable for Poisson data
(bias 3.0%), so it is the safer choice when in doubt.

## 4. Curved propensity: smooth terms with mgcv

The probability of exposure depends on `l` through a sine curve; the
true effect is 2. Oracle weights use the true propensity.

``` r

d <- data.frame(l = runif(2000, -2, 2))
e <- plogis(2 * sin(1.5 * d$l))                        # true propensity
d$a <- rbinom(2000, 1, e)
d$y <- 2 * d$a + 3 * sin(1.5 * d$l) + rnorm(2000)
w.gam <- ipwpoint(exposure = a, family = "binomial", link = "logit",
                  numerator = ~ 1, denominator = ~ s(l), data = d)$ipw.weights
oracle <- ifelse(d$a == 1, mean(d$a) / e, (1 - mean(d$a)) / (1 - e))   # stabilized
```

| Method | Truth | Mean estimate | Bias | MCSE of bias | Empirical SD | RMSE | Replications |
|:---|---:|---:|---:|---:|---:|---:|---:|
| unweighted | 2 | 4.542 | 2.542 | 0.004 | 0.086 | 2.543 | 500 |
| linear logistic weights (misspecified) | 2 | 2.657 | 0.657 | 0.004 | 0.083 | 0.662 | 500 |
| s(l) weights, n = 2000 | 2 | 2.066 | 0.066 | 0.003 | 0.063 | 0.091 | 500 |
| s(l) weights, n = 10000 | 2 | 2.032 | 0.032 | 0.002 | 0.027 | 0.042 | 200 |
| oracle weights (true propensity), n = 2000 | 2 | 2.001 | 0.001 | 0.007 | 0.148 | 0.148 | 500 |
| oracle weights (true propensity), n = 10000 | 2 | 1.997 | -0.003 | 0.004 | 0.061 | 0.061 | 200 |

A straight-line logistic model leaves a quarter of the confounding bias.
With `s(l)` the bias is 3.3% at n = 2000 and 1.6% at n = 10000: it
halves as the sample grows five-fold, as expected for a smoothing
estimate that approaches the true curve. Oracle weights are unbiased,
which confirms that the remaining bias comes from estimating the curve.

## 5. Time-varying exposure, `type = "all"`

### 5a. Ordinal dose

A dose level (0, 1, 2) is chosen at each of 4 visits from a
proportional-odds model given the confounder `L` and the previous dose.
`L` is affected by the previous dose and by an unmeasured `U` that also
affects the outcome; the outcome rises by 1 per unit of cumulative dose.
This is the situation in the JSS paper where ordinary adjustment fails.

``` r

for (t in 1:4) {
  l <- U + 0.4 * prev + rnorm(n, 0, 0.5)
  dose <- findInterval(0.8 * l + 0.6 * prev + rlogis(n), c(0, 1.5))   # 0, 1, 2
  prev <- dose
}
y <- cumulative.dose + 2 * U + rnorm(n)
w <- ipwtm(exposure = dose, family = "ordinal", link = "logit",
           numerator = ~ dose.lag, denominator = ~ dose.lag + l,
           id = id, timevar = time, type = "all", data = long)
# oracle: true probability of the observed dose, multiplied over visits
w.oracle <- ipwprobtm(den = long$p.obs, id = long$id, timevar = long$time)
# weight at the last visit, then lm(y ~ cumulative.dose, weights = ...)
```

| Method | Truth | Mean estimate | Bias | MCSE of bias | Empirical SD | RMSE | Replications |
|:---|---:|---:|---:|---:|---:|---:|---:|
| unweighted | 1 | 1.491 | 0.491 | 0.002 | 0.035 | 0.493 | 250 |
| regression adjustment for L | 1 | 0.856 | -0.144 | 0.002 | 0.030 | 0.147 | 250 |
| ordinal weights (correct), 500 people | 1 | 1.031 | 0.031 | 0.007 | 0.110 | 0.114 | 250 |
| ordinal weights (correct), 2500 people | 1 | 1.010 | 0.010 | 0.005 | 0.076 | 0.077 | 200 |
| multinomial weights, 500 people | 1 | 1.015 | 0.015 | 0.008 | 0.121 | 0.122 | 250 |
| oracle weights, ipwprobtm() with true probabilities, 500 people | 1 | 1.013 | 0.013 | 0.010 | 0.156 | 0.157 | 250 |
| oracle weights, ipwprobtm() with true probabilities, 2500 people | 1 | 1.007 | 0.007 | 0.006 | 0.089 | 0.090 | 200 |

Both naive analyses are biased (+49% and -14%). Ordinal and multinomial
weights recover the effect; the 3.1% bias of the ordinal weights at 500
people falls to 1.0% at 2500. Oracle weights, from
[`ipwprobtm()`](https://tranhung93.github.io/ipw/reference/ipwprobtm.md)
with the true probability of each observed dose, are within 1.5%, which
also checks that
[`ipwprobtm()`](https://tranhung93.github.io/ipw/reference/ipwprobtm.md)
multiplies the ratios correctly over time.

### 5b. Smooth terms in `ipwtm()`

A binary exposure is chosen at each visit with a curved effect of `L`
(logit = -0.3 + 2 sin(1.5 L) + 0.6 × previous exposure), in the same
design.

``` r

w <- ipwtm(exposure = a, family = "binomial", link = "logit", numerator = ~ a.lag,
           denominator = ~ a.lag + s(l), id = id, timevar = time, type = "all", data = long)
```

| Method | Truth | Mean estimate | Bias | MCSE of bias | Empirical SD | RMSE | Replications |
|:---|---:|---:|---:|---:|---:|---:|---:|
| unweighted | 1 | 1.662 | 0.662 | 0.004 | 0.061 | 0.665 | 250 |
| linear logistic weights (misspecified) | 1 | -0.111 | -1.111 | 0.031 | 0.489 | 1.213 | 250 |
| s(l) weights | 1 | 1.083 | 0.083 | 0.015 | 0.243 | 0.256 | 250 |

Here a straight-line model is worse than no weighting at all: its
weights turn the effect of +1 into a negative estimate. With `s(l)` the
bias is 8.3%, the largest of any correctly specified model in this
study; with a binary exposure at 4 visits and 500 people, part of the
curve rests on few observations.

## 6. Recurrent exposure events, `family = "survival"`, `type = "all"`

### 6a and 6b. Weights against the true weights

Exposure events (for example new prescriptions) happen repeatedly with a
known proportional hazard: baseline 2 per unit time, log hazard ratio
0.5 for a fixed covariate `x` and -1 for a time-varying covariate `xt`.
Because the true hazard is known, the true weights are known too. In 6a
all people share a grid of 10 intervals (as in `haartdat`); in 6b the
events happen in continuous time, `xt` changes at different times for
each person, and each event ends an interval (start-stop data cut at
events, as in Graffeo et al. 2018).

``` r

w <- ipwtm(exposure = expo, family = "survival", denominator = ~ x + xt,
           id = id, tstart = tstart, timevar = fuptime, type = "all", data = d)
coef(w$den.mod)
median(abs(log(w$ipw.weights / true.weights)))
```

| Method | Truth | Mean estimate | Bias | MCSE of bias | Empirical SD | RMSE | Replications |
|:---|---:|---:|---:|---:|---:|---:|---:|
| Cox coefficient of x | 0.5 | 0.499 | -0.001 | 0.005 | 0.084 | 0.084 | 250 |
| Cox coefficient of xt | -1.0 | -1.010 | -0.010 | 0.007 | 0.117 | 0.118 | 250 |
| median \|log(estimated / true weight)\| | 0.0 | 0.072 | 0.072 | 0.001 | 0.022 | 0.075 | 250 |

| Method | Truth | Mean estimate | Bias | MCSE of bias | Empirical SD | RMSE | Replications |
|:---|---:|---:|---:|---:|---:|---:|---:|
| Cox coefficient of x | 0.5 | 0.494 | -0.006 | 0.006 | 0.087 | 0.087 | 250 |
| Cox coefficient of xt | -1.0 | -0.998 | 0.002 | 0.007 | 0.104 | 0.103 | 250 |
| median \|log(estimated / true weight)\| | 0.0 | 0.090 | 0.090 | 0.001 | 0.021 | 0.093 | 250 |

The Cox coefficients are unbiased in both designs, and the estimated
weights are typically within 7 to 9% of the true weights, which is the
estimation noise of the baseline hazard.

### 6c. Causal effect of the number of exposure events

Exposure events happen in continuous time with hazard 1.5 exp(0.6 L).
`L` is updated at 4 visits, is lowered by earlier events (as treatment
raising CD4 lowers the need for more) and depends on an unmeasured `U`
that also affects the outcome. The outcome rises by 1 per event.

``` r

w <- ipwtm(exposure = expo, family = "survival", numerator = ~ 1, denominator = ~ l,
           id = id, tstart = tstart, timevar = fuptime, type = "all", data = d)
# weight at each person's last row, then lm(y ~ number.of.events, weights = ...)
```

| Method | Truth | Mean estimate | Bias | MCSE of bias | Empirical SD | RMSE | Replications |
|:---|---:|---:|---:|---:|---:|---:|---:|
| unweighted | 1 | 1.794 | 0.794 | 0.004 | 0.062 | 0.796 | 250 |
| regression adjustment for L | 1 | 1.214 | 0.214 | 0.003 | 0.045 | 0.219 | 250 |
| survival type all weights | 1 | 1.008 | 0.008 | 0.014 | 0.224 | 0.224 | 250 |

The weights reduce the bias from 79% (unweighted) and 21% (adjustment
for `L`) to 0.8%. The method of Graffeo et al. uses, for an interval
ending in an event, the probability of an event during that interval;
this simulation shows that the approximation works for estimating a
causal effect.

## 7. Standard errors: bootstrap, robust and naive

Each table shows the mean estimated standard error divided by the actual
spread of the estimates (1 is ideal) and the coverage of 95% confidence
intervals. The bootstrap re-estimates the weights and the MSM in each of
200 bootstrap samples (function `boot_sim()` in the script).

### 7a and 7b. Linear MSM for a binary exposure

True effect 10, n = 500, stabilized weights. In 7a, exposure
probabilities are mostly between 0.1 and 0.9 (good overlap); 7b is the
simulated example of the JSS paper, where many people are almost certain
to be exposed or unexposed.

``` r

msm.effect <- function(dd) {
  w <- ipwpoint(exposure = a, family = "binomial", link = "logit", numerator = ~ 1,
                denominator = ~ l, data = dd)$ipw.weights
  coef(lm(y ~ a, data = dd, weights = w))["a"]
}
boot.est <- replicate(200, msm.effect(d[sample(nrow(d), replace = TRUE), ]))
sd(boot.est)                          # bootstrap standard error
quantile(boot.est, c(0.025, 0.975))   # percentile interval
```

| Method | Mean estimate | Mean SE | Empirical SD | SE / empirical SD | Coverage of 95% CI | MCSE of coverage | Replications |
|:---|---:|---:|---:|---:|---:|---:|---:|
| bootstrap, percentile CI | 10.031 | 0.492 | 0.481 | 1.023 | 0.945 | 0.007 | 1000 |
| bootstrap SE, Wald CI | 10.031 | 0.492 | 0.481 | 1.023 | 0.955 | 0.007 | 1000 |
| robust SE (svyglm), Wald CI | 10.031 | 0.537 | 0.481 | 1.115 | 0.972 | 0.005 | 1000 |

| Method | Mean estimate | Mean SE | Empirical SD | SE / empirical SD | Coverage of 95% CI | MCSE of coverage | Replications |
|:---|---:|---:|---:|---:|---:|---:|---:|
| bootstrap, percentile CI | 11.055 | 1.135 | 1.848 | 0.614 | 0.575 | 0.016 | 1000 |
| bootstrap SE, Wald CI | 11.055 | 1.135 | 1.848 | 0.614 | 0.596 | 0.016 | 1000 |
| robust SE (svyglm), Wald CI | 11.055 | 1.144 | 1.848 | 0.619 | 0.600 | 0.015 | 1000 |

With good overlap the bootstrap standard error is accurate (ratio 1.02,
Wald coverage 95.5%), while the robust standard error is 12% too large
(coverage 97.2%), in line with Austin (2016); percentile intervals cover
slightly less than 95%. In the JSS example the estimate itself is biased
(11.1 instead of 10) and all standard errors are 38% too small, so
coverage is 60% or less whatever method is used.

### 7c. Weighted Cox MSM

Austin (2016) studied weighted Cox models. Here event times are
exponential with hazard 0.1 exp(-0.7 a + 0.8 l), censored at 5; n =
1000. The true marginal log hazard ratio, -0.619, is computed from both
potential outcomes of 200,000 people (it differs from -0.7 because
hazard ratios are not collapsible).

| Method | Mean estimate | Mean SE | Empirical SD | SE / empirical SD | Coverage of 95% CI | MCSE of coverage | Replications |
|:---|---:|---:|---:|---:|---:|---:|---:|
| naive SE (coxph), Wald CI | -0.617 | 0.117 | 0.11 | 1.065 | 0.958 | 0.009 | 500 |
| robust SE (coxph), Wald CI | -0.617 | 0.121 | 0.11 | 1.104 | 0.968 | 0.008 | 500 |
| bootstrap SE, Wald CI | -0.617 | 0.110 | 0.11 | 0.998 | 0.940 | 0.011 | 500 |
| bootstrap, percentile CI | -0.617 | 0.110 | 0.11 | 0.998 | 0.928 | 0.012 | 500 |

All three standard errors are close to the truth here. The bootstrap
standard error is accurate (ratio 1.00); the robust one is again too
large (+10%) and the naive one somewhat too large (+7%), as Austin
(2016) found. All Wald intervals cover the truth in 94 to 97% of
samples, so with stabilized weights and good overlap the choice matters
little in this setting.

### 7d. Time-varying exposure, bootstrap by person

The design of scenario 5 with a binary exposure; the bootstrap resamples
people with all their visits, and each sampled copy gets a new id.

``` r

boot_sim(long, tv.effect, id = "id", R = 200)
```

| Method | Mean estimate | Mean SE | Empirical SD | SE / empirical SD | Coverage of 95% CI | MCSE of coverage | Replications |
|:---|---:|---:|---:|---:|---:|---:|---:|
| bootstrap by person, percentile CI | 1.047 | 0.139 | 0.226 | 0.616 | 0.722 | 0.020 | 500 |
| bootstrap by person SE, Wald CI | 1.047 | 0.139 | 0.226 | 0.616 | 0.744 | 0.020 | 500 |
| robust SE (svyglm), Wald CI | 1.047 | 0.158 | 0.226 | 0.698 | 0.798 | 0.018 | 500 |

Neither method is adequate here: the bootstrap standard error is 38% too
small and the robust one 30% too small, with coverage of 74% and 80%.
The weights are products over 4 visits, so a few people get very large
weights and the estimate has a heavy-tailed distribution that 500 people
do not represent well, the same problem as in 7b. Truncation, a larger
sample or a better-balanced design would be needed; the bootstrap does
not repair it.

## Conclusions

The simulations support the new models: with the correct model, each of
them removes nearly all confounding, and the small biases that remain
are either shared by weights from the true model or shrink with the
sample size. They support a bootstrap that re-estimates the weights for
analyses with well-behaved weights, where it gives accurate standard
errors and robust standard errors were conservative in these scenarios
(binary point treatment, Cox model). They also show its limit: when
weights are extreme (poor overlap, long products of time-varying
weights), no standard error is reliable, and the first step is to
examine the weights ([`summary()`](https://rdrr.io/r/base/summary.html),
[`ipwplot()`](https://tranhung93.github.io/ipw/reference/ipwplot.md))
and consider overlap weights, truncation or a different model.

## References

- Austin PC (2016). Variance estimation when using inverse probability
  of treatment weighting (IPTW) with survival analysis. *Statistics in
  Medicine* 35, 5642-5655.
- Graffeo N, Latouche A, Geskus RB, Chevret S (2018). Modeling
  time-varying exposure using inverse probability of treatment weights.
  *Biometrical Journal* 60, 323-332.
- Hirano K, Imbens GW, Ridder G (2003). Efficient estimation of average
  treatment effects using the estimated propensity score. *Econometrica*
  71, 1161-1189.
- Li F, Morgan KL, Zaslavsky AM (2018). Balancing covariates via
  propensity score weighting. *Journal of the American Statistical
  Association* 113, 390-400.
- van der Wal WM, Geskus RB (2011). ipw: An R package for inverse
  probability weighting. *Journal of Statistical Software* 43(13), 1-23.
