# ipw, WeightIt and PSweight side by side

``` r

library(ipw)
library(survival)
```

[WeightIt](https://ngreifer.github.io/WeightIt/) (Greifer) and
[ipw](https://tranhung93.github.io/ipw/) (Tran, van der Wal & Geskus)
both estimate weights for marginal structural models;
[PSweight](https://CRAN.R-project.org/package=PSweight) (Zhou, Tong, Li,
Thomas & Li) estimates propensity score weighted treatment effects for
point treatments. This article runs the packages on the same data: a
binary and a categorical point treatment, where all three apply and
should give the same answer, and longitudinal data with individual visit
times, where ipw and WeightIt work differently. WeightIt and PSweight
are used only through their exported functions. This page was built with
ipw 1.4.0, WeightIt 2.1.0 and PSweight 2.1.2.

**What PSweight offers.** PSweight is designed for point treatments,
binary or with several levels. It offers inverse probability, ATT,
overlap, matching and entropy weights, propensity score trimming, and
augmented (doubly robust) estimators that add an outcome model. The
propensity score can be fitted by logistic or multinomial regression,
generalized boosted models or SuperLearner. Inference uses sandwich
(M-estimation) variances that account for the estimated propensity
score, or the bootstrap, for differences in means and, for binary or
count outcomes, risk ratios and odds ratios. It also has functions for
survey data with survey weights and for clustered data, both mainly for
binary treatments, and covariate balance summaries and plots. It does
not handle time-varying treatments.

## 1. Binary point treatment

A confounder `l`, a binary covariate `x`, a binary treatment `a` and an
outcome `y`; the effect of `a` on `y` is 2 for everyone, so all
estimands are 2 in truth.

``` r

set.seed(2026)
n <- 1000
d <- data.frame(l = rnorm(n), x = rbinom(n, 1, 0.4))
d$a <- rbinom(n, 1, plogis(-0.3 + 0.8 * d$l + 0.5 * d$x))
d$y <- 2 * d$a + d$l + d$x + rnorm(n)
stopifnot(all.equal(sum(d$y), stored$data_check))   # same data as the bootstrap
estimands <- c("ATE", "ATT", "ATC", "ATO")
```

### The weights

With a logistic model, both packages fit the same propensity score with
[`glm()`](https://rdrr.io/r/stats/glm.html) and apply the same formulas,
so the weights should agree up to rounding:

``` r

w_ipw <- lapply(setNames(estimands, estimands), function(est)
  ipwpoint(exposure = a, family = "binomial", link = "logit",
           denominator = ~ l + x, data = d, estimand = est))
w_wit <- lapply(setNames(estimands, estimands), function(est)
  WeightIt::weightit(a ~ l + x, data = d, method = "glm", estimand = est))
sapply(estimands, function(est)
  c(max.abs.difference = max(abs(w_ipw[[est]]$ipw.weights - w_wit[[est]]$weights)),
    max.weight = max(w_ipw[[est]]$ipw.weights)))
#>                             ATE          ATT          ATC          ATO
#> max.abs.difference 1.065814e-14 1.065814e-14 7.105427e-15 3.330669e-16
#> max.weight         1.180340e+01 8.832774e+00 1.080340e+01 9.152786e-01
```

The largest difference is at the level of floating-point rounding.
Stabilized ATE weights agree too: `numerator = ~ 1` in ipw and
`stabilize = TRUE` in WeightIt.

PSweight scales the weights to sum to 1 within each treatment group,
which changes no weighted group mean or effect estimate, so the weights
are compared after the same scaling. PSweight has no “ATC” option, but
its ATT weights with the untreated as the “treated” group
(`trtgrp = "0"`) are the ATC weights.

``` r

psw_args <- list(ATE = list(weight = "IPW"), ATT = list(weight = "treated"),
                 ATC = list(weight = "treated", trtgrp = "0"),
                 ATO = list(weight = "overlap"))
norm_by_group <- function(w, g) w / ave(w, g, FUN = sum)
w_ipw_bin <- lapply(setNames(estimands, estimands),
  function(est) ipwpoint(exposure = a, family = "binomial", link = "logit",
                         denominator = ~ l + x, data = d, estimand = est)$ipw.weights)
sapply(estimands, function(est) {
  ss <- do.call(PSweight::SumStat, c(list(ps.formula = a ~ l + x, data = d), psw_args[[est]]))
  max(abs(norm_by_group(w_ipw_bin[[est]], d$a) -
          norm_by_group(ss$ps.weights[[psw_args[[est]]$weight]], d$a)))
})
#>          ATE          ATT          ATC          ATO 
#> 0.000000e+00 4.336809e-19 3.469447e-18 8.673617e-19
```

### Standard errors

The marginal structural model is the weighted regression of `y` on `a`.
Five standard errors of the coefficient of `a`:

- **ipwvcov()**: ipw’s M-estimation standard error (new), which accounts
  for the estimation of the weights;
- **glm_weightit(vcov = “asympt”)**: WeightIt’s M-estimation standard
  error;
- **PSweight**: the sandwich standard error of the weighted difference
  in means from `PSweight()` and
  [`summary()`](https://rdrr.io/r/base/summary.html), which for the
  saturated model `y ~ a` is the same quantity;
- **robust (HC0)**: the sandwich standard error that treats the weights
  as known, as [`svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html)
  or `glm_weightit(vcov = "HC0")` give it;
- **bootstrap**: 2000 bootstrap samples, re-estimating the weights in
  each (computed once, see the setup chunk).

``` r

se_table <- t(sapply(estimands, function(est) {
  msm <- lm(y ~ a, data = d, weights = w_ipw[[est]]$ipw.weights)
  fit_m <- WeightIt::glm_weightit(y ~ a, data = d, weightit = w_wit[[est]], vcov = "asympt")
  fit_r <- WeightIt::glm_weightit(y ~ a, data = d, weightit = w_wit[[est]], vcov = "HC0")
  c(estimate = unname(coef(msm)["a"]),
    ipwvcov = sqrt(ipwvcov(w_ipw[[est]], msm)["a", "a"]),
    glm_weightit.asympt = sqrt(vcov(fit_m)["a", "a"]),
    PSweight = if (has_psweight) unname(summary(do.call(PSweight::PSweight,
      c(list(ps.formula = a ~ l + x, yname = "y", data = d), psw_args[[est]])))$estimates[1, "Std.Error"])
      else NA,
    robust.HC0 = sqrt(vcov(fit_r)["a", "a"]),
    bootstrap = unname(stored$boot_se[est]))
}))
round(se_table, 4)
#>     estimate ipwvcov glm_weightit.asympt PSweight robust.HC0 bootstrap
#> ATE   2.0259  0.0804              0.0804   0.0804     0.1234    0.0836
#> ATT   2.1543  0.0862              0.0862   0.0862     0.1134    0.0878
#> ATC   1.9160  0.1094              0.1094   0.1094     0.1380    0.1131
#> ATO   2.0807  0.0684              0.0684   0.0684     0.0962    0.0703
```

[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md),
`glm_weightit(vcov = "asympt")` and PSweight agree to many decimals:
they solve the same stacked estimating equations (Lunceford & Davidian
2004), and in the package tests they agree to about 1e-9 or better. All
three are close to the bootstrap. The robust standard error, which
ignores that the weights were estimated, is larger for all four
estimands here, so confidence intervals based on it would be wider than
needed. For the ATE this is the usual direction (estimating the weights
reduces the variance); for the ATT, ATC and ATO the robust standard
error is not guaranteed to be conservative and can also be too small.
The M-estimation or bootstrap standard error is therefore preferred.

### A Cox outcome model

For a time-to-event outcome the MSM is a weighted Cox model. We add a
survival time to the same data (hazard ratio exp(-0.5) for `a`,
administrative censoring at 5) and compare
[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
with [`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)’s robust
standard error, which treats the weights as known, and with
`coxph_weightit(vcov = "asympt")`:

``` r

set.seed(8)
t_ev <- rexp(nrow(d), 0.2 * exp(-0.5 * d$a + 0.5 * d$l + 0.3 * d$x))
d$status <- as.numeric(t_ev < 5)
d$time <- pmin(t_ev, 5)
cox_table <- t(sapply(estimands, function(est) {
  msm <- coxph(Surv(time, status) ~ a, data = d, weights = w_ipw[[est]]$ipw.weights)
  fit_w <- WeightIt::coxph_weightit(Surv(time, status) ~ a, data = d,
                                    weightit = w_wit[[est]], vcov = "asympt")
  c(log.hr = unname(coef(msm)),
    ipwvcov = sqrt(ipwvcov(w_ipw[[est]], msm)[1, 1]),
    coxph_weightit.asympt = sqrt(vcov(fit_w)[1, 1]),
    robust.coxph = sqrt(msm$var[1, 1]))
}))
round(cox_table, 4)
#>      log.hr ipwvcov coxph_weightit.asympt robust.coxph
#> ATE -0.4718  0.0833                0.0835       0.0941
#> ATT -0.4120  0.0826                0.0832       0.0884
#> ATC -0.5458  0.1058                0.1105       0.1145
#> ATO -0.4637  0.0844                0.0844       0.0896
```

[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
uses the weighted partial likelihood score as the estimating equation
and accounts for the weights in both each person’s own term and the
risk-set sums; the package tests check it against central finite
differences of the stacked estimating equations (agreement about 1e-7).
WeightIt’s values are close but not identical: its `vcov = "HC0"` equals
[`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)’s robust
variance exactly, so the difference lies in how the estimated weights
enter. In our checks the relative difference was below 0.1% for ATE and
overlap weights and shrank with the sample size; for ATT and ATC weights
WeightIt’s standard error was 1 to 4% larger and the difference did not
shrink (ATC in the table above).
[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
equals the finite-difference reference on these data too. The [methods
article](https://tranhung93.github.io/ipw/articles/methods.md) gives
details.

The table below summarises a simulation (500 people per data set,
treatment depending on a confounder, conditional hazard ratio exp(-0.5),
administrative censoring) of 95% Wald intervals for the marginal log
hazard ratio. The truth is the marginal log hazard ratio in the target
population, computed from both potential outcomes of a very large
sample. The bootstrap (200 samples, weights re-estimated in each) was
run on 300 of the 1000 data sets. Standard errors are judged on the
variance scale (Morris, White & Crowther 2019): the variance ratio is
the mean estimated variance, mean(SE²), divided by the empirical
variance of the estimates (1 is ideal), and the model SE,
sqrt(mean(SE²)), is compared with the empirical SE.

``` r

sim_table(cox_sim$summary, c(Estimand = "estimand", Truth = "truth"))
```

| estimand | truth | Method | Empirical SE | Model SE | Variance ratio | Relative error in model SE, % (MCSE) | Coverage | MCSE of coverage | Data sets |
|:---|---:|:---|---:|---:|---:|:---|---:|---:|---:|
| ATE | -0.4502 | ipwvcov | 0.1281 | 0.1280 | 1.00 | -0.1 (2.2) | 0.956 | 0.006 | 1000 |
| ATE | -0.4502 | robust | 0.1281 | 0.1404 | 1.20 | +9.6 (2.5) | 0.969 | 0.005 | 1000 |
| ATE | -0.4502 | bootstrap | 0.1281 | 0.1283 | 1.00 | +0.2 (2.3) | 0.950 | 0.013 | 300 |
| ATO | -0.4614 | ipwvcov | 0.1327 | 0.1325 | 1.00 | -0.2 (2.2) | 0.956 | 0.006 | 1000 |
| ATO | -0.4614 | robust | 0.1327 | 0.1402 | 1.12 | +5.6 (2.4) | 0.961 | 0.006 | 1000 |
| ATO | -0.4614 | bootstrap | 0.1327 | 0.1334 | 1.01 | +0.5 (2.3) | 0.950 | 0.013 | 300 |

The empirical variance of the estimates is what a variance estimator
should estimate, but with 1000 data sets it is itself uncertain by about
4.5% (relative Monte Carlo standard error sqrt(2 / (*R* − 1))). The
variance ratio of
[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md) is
1.00 (ATE) and 1.00 (ATO), and that of the bootstrap 1.00 and 1.01; the
relative errors in the model SE are all under 1% with Monte Carlo
standard errors above 2%: no detectable difference. The robust variance,
which treats the weights as known, is too large here: variance ratio
1.20 (ATE) and 1.12 (ATO), model SE +10% (MCSE 2) and +6% (MCSE 2). Its
intervals therefore cover more often than 95%: conservative, but wider
than needed. The coverages of
[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
(1000 data sets) and of the bootstrap (300 data sets) are both within
simulation error of 95% (the Monte Carlo standard error of a coverage
near 95% is about 0.7 percentage points with 1000 data sets and 1.3 with
300). Coverage does not depend on the empirical variance: it counts how
often the interval contains the known true value. The [methods
article](https://tranhung93.github.io/ipw/articles/methods.html#how-the-simulations-judge-a-standard-error)
explains these checks in more detail, including why variances rather
than standard errors are compared.

[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
covers binary exposures with logistic exposure models: point treatments
with [`lm()`](https://rdrr.io/r/stats/lm.html),
[`glm()`](https://rdrr.io/r/stats/glm.html) or
[`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) outcome models,
and time-varying weights from
[`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md)
(section 3). WeightIt’s outcome-model functions cover more:
multi-category and continuous treatments, many weighting methods, and
ordinal and multinomial outcome models.

## 2. Categorical point treatment

A treatment with three levels, A, B and C. With `estimand = "ATT"` and a
**focal** level, people at that level keep weight 1 and the others are
reweighted to resemble them; `estimand = "ATO"` gives generalized
overlap weights (Li & Li 2019).

``` r

set.seed(31)
n <- 1000
dc <- data.frame(l = rnorm(n), z = rnorm(n))
pr <- cbind(1, exp(0.5 + 0.8 * dc$l), exp(-0.3 - 0.6 * dc$l + 0.4 * dc$z))
pr <- pr / rowSums(pr)
dc$m <- factor(apply(pr, 1, function(p) sample(c("A", "B", "C"), 1, prob = p)))
table(dc$m)
#> 
#>   A   B   C 
#> 251 510 239
```

``` r

compare <- function(estimand, focal = NULL, ...) {
  wi <- ipwpoint(exposure = m, family = "multinomial", denominator = ~ l + z, data = dc,
                 estimand = estimand, focal = focal, ...)$ipw.weights
  ww <- WeightIt::weightit(m ~ l + z, data = dc, method = "glm",
                           estimand = estimand, focal = focal)$weights
  max(abs(wi / ww - 1))
}
rbind(
  default.multinom = c(ATT.focal.A = compare("ATT", "A"), ATO = compare("ATO")),
  tight.multinom = c(ATT.focal.A = compare("ATT", "A", reltol = 1e-12, maxit = 1000),
                     ATO = compare("ATO", reltol = 1e-12, maxit = 1000)))
#>                   ATT.focal.A          ATO
#> default.multinom 6.282293e-04 9.086836e-04
#> tight.multinom   2.844251e-08 3.185674e-08
```

The table shows the largest relative difference between the weights. ipw
fits the multinomial model with
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html), which
by default stops when the relative change in the fit criterion falls
below `1e-8`; that leaves relative differences of up to about 1e-3 from
WeightIt, whose own fitter converges further. With a tighter tolerance
(`reltol`, passed on through `...`), the weights agree to about 1e-8:
the formulas are the same.

PSweight fits the multinomial model with the same defaults as
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html), so its
weights (scaled to sum to 1 within each level) agree with ipw’s default
fit to rounding error; its `trtgrp` is ipw’s `focal`:

``` r

compare_psw <- function(estimand, weight, focal = NULL) {
  wi <- ipwpoint(exposure = m, family = "multinomial", denominator = ~ l + z, data = dc,
                 estimand = estimand, focal = focal)$ipw.weights
  wp <- PSweight::SumStat(ps.formula = m ~ l + z, data = dc, weight = weight,
                          trtgrp = focal)$ps.weights[[weight]]
  max(abs(wi / ave(wi, dc$m, FUN = sum) / (wp / ave(wp, dc$m, FUN = sum)) - 1))
}
c(ATE = compare_psw("ATE", "IPW"), ATT.focal.A = compare_psw("ATT", "treated", "A"),
  ATO = compare_psw("ATO", "overlap"))
#>          ATE  ATT.focal.A          ATO 
#> 0.000000e+00 0.000000e+00 2.220446e-16
```

Within each level, the ATT weights make the covariate means resemble
those of the focal group:

``` r

w_att <- ipwpoint(exposure = m, family = "multinomial", denominator = ~ l + z,
                  data = dc, estimand = "ATT", focal = "A")$ipw.weights
rbind(unweighted = tapply(dc$l, dc$m, mean),
      ATT.focal.A = sapply(levels(dc$m), function(k)
        weighted.mean(dc$l[dc$m == k], w_att[dc$m == k])))
#>                      A          B          C
#> unweighted  -0.2166951  0.4582985 -0.6187572
#> ATT.focal.A -0.2166951 -0.2541368 -0.2024986
```

For an ordered exposure, ipw uses
[`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html)
(`family = "ordinal"`); WeightIt also fits an ordinal model for an
ordered factor. In the package tests the two agree to 1e-4 or better.

## 3. Longitudinal data with individual visit times

Here the two packages are designed differently.

- [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) works
  on **long data**: one row per person per interval between visits, with
  `tstart` and `timevar`. The exposure model is fitted once, pooled over
  all rows (a Cox model, a pooled logistic regression, or a multinomial,
  ordinal or GEE model), so each person can have their own number of
  visits at their own times.
- `weightitMSM()` works on **wide data**: one row per person, and one
  treatment model per time point, given in `formula.list`. Its help page
  states that only wide data are supported. The time points are
  therefore the same for everyone. The weight is the product of the
  weights at the time points. This fits designs with fixed visits, such
  as yearly measurements in a cohort or a trial with scheduled visits,
  and it allows any WeightIt method at each time point.

### Data with individual visit times

We simulate 1000 people followed for 5 years. Each person is seen at
their own times, about twice a year on average with random gaps. At each
visit a confounder `l` is measured; it depends on its previous value and
on treatment. People not yet treated may start treatment at a visit,
more likely when `l` is low. Once started, treatment continues.

``` r

set.seed(11)
n <- 1000
sim_person <- function(i) {
  t <- 0; l <- rnorm(1); a <- 0; l0 <- l; rows <- list()
  repeat {
    t_next <- t + rexp(1, rate = 2)           # gap to the next visit
    if (t_next > 5) break
    l <- 0.8 * l + 0.4 * a + rnorm(1, 0, 0.5)  # confounder measured at the visit
    start <- if (a == 0) rbinom(1, 1, plogis(-2.5 - 1.0 * l)) else 0
    rows[[length(rows) + 1]] <- data.frame(id = i, tstart = t, time = t_next,
                                           l0 = l0, l = l, start = start)
    if (start == 1) a <- 1
    t <- t_next
  }
  do.call(rbind, rows)
}
long <- do.call(rbind, lapply(seq_len(n), sim_person))
long$a <- ave(long$start, long$id, FUN = cummax)   # 1 from the visit of starting on
head(long[, c("id", "tstart", "time", "l", "a")], 4)
#>   id     tstart       time           l a
#> 1  1 0.00000000 0.01060837 -1.57079786 1
#> 2  1 0.01060837 1.24229627 -1.13336715 1
#> 3  1 1.24229627 2.34270656  0.01376242 1
#> 4  1 2.34270656 2.41639228  0.38814845 1
c(people = n, rows = nrow(long), started = sum(long$start))
#>  people    rows started 
#>    1000   10040     576
summary(as.vector(table(long$id)))   # number of visits per person
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    2.00    8.00   10.00   10.04   12.00   22.00
```

### ipw: the long data as they are

[`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) fits a
Cox model for the time to starting treatment on the intervals between
each person’s own visits; with `type = "first"` the rows after the start
are not used for the model and keep the weight reached at the start.

``` r

w_tm <- ipwtm(exposure = a, family = "survival", numerator = ~ l0,
              denominator = ~ l0 + l, id = id, tstart = tstart, timevar = time,
              type = "first", data = long)
summary(w_tm$ipw.weights)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.1936  0.8119  0.9513  0.9980  1.0614  5.2554
```

The weights can be used directly in a weighted model on the same rows,
for example
`coxph(Surv(tstart, time, event) ~ a + l0, weights = w_tm$ipw.weights, ...)`
for an outcome event.

### WeightIt: a common time grid first

For `weightitMSM()`, the data must first be put on time points that are
the same for everyone, say the end of each year. At each time point `k`,
the analyst has to choose a value of `l` and of the treatment. A common
choice is the last value observed up to `k` (the baseline value if there
was no visit yet) and “started by `k`”:

``` r

grid <- 1:5
last_before <- function(k, column, default) {
  sub <- long[long$time <= k, ]
  last <- sub[!duplicated(sub$id, fromLast = TRUE), ]   # rows are in time order
  out <- default
  out[last$id] <- last[[column]]
  out
}
l0 <- long$l0[!duplicated(long$id)]
wide <- data.frame(id = seq_len(n), l0 = l0)
for (k in grid) {
  wide[[paste0("L", k)]] <- last_before(k, "l", l0)
  wide[[paste0("A", k)]] <- last_before(k, "a", rep(0, n))
}
head(wide, 3)
#>   id         l0          L1 A1         L2 A2         L3 A3         L4 A4
#> 1  1 -0.5910311 -1.57079786  1 -1.1333671  1  0.3881485  1  0.7743123  1
#> 2  2  0.8877716  0.07944407  0 -0.1577471  0 -0.1577471  0 -0.4809450  0
#> 3  3 -1.5774979 -1.10960567  0 -1.2076286  0 -0.7749272  0 -0.8303085  0
#>           L5 A5
#> 1  0.5324167  1
#> 2 -1.1703153  0
#> 3 -0.9115556  0
```

One model per time point, each with the treatment history and the
confounder history up to then (here the previous treatment and current
`l`), and a numerator model with the treatment history and baseline
`l0`:

``` r

formulas <- lapply(grid, function(k)
  as.formula(if (k == 1) "A1 ~ l0 + L1" else
               sprintf("A%d ~ A%d + l0 + L%d", k, k - 1, k)))
num_formulas <- lapply(grid, function(k)
  as.formula(if (k == 1) "~ l0" else sprintf("~ A%d + l0", k - 1)))
w_msm <- WeightIt::weightitMSM(formulas, data = wide, method = "glm",
                               stabilize = TRUE, num.formula = num_formulas)
summary(w_msm$weights)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.4485  0.9055  0.9797  1.0053  1.0731  2.5403
```

[`glm()`](https://rdrr.io/r/stats/glm.html) warns here (warnings hidden)
that fitted probabilities of 0 or 1 occurred: people who already started
have `A(k-1) = 1` and therefore `A(k) = 1`, so the model at each time
point also “predicts” a treatment that is no longer a choice; their
factor is about 1. This matches what `type = "first"` does in ipw by
leaving those rows out of the model.

### What the grid means for such data

The grid analysis uses the data in a coarser form. How much is lost
depends on how far the visit times are from the grid:

``` r

starts <- long[long$start == 1, ]
grid_k <- ceiling(starts$time)                     # grid point at which the start is seen
l_at_grid <- mapply(function(i, k) wide[i, paste0("L", k)], starts$id, grid_k)
last_visit <- sapply(grid, function(k) {
  sub <- long[long$time <= k, ]
  mean(k - tapply(sub$time, sub$id, max))
})
c(visits.per.person = nrow(long) / n,
  grid.points.per.person = length(grid),
  starts.l.differs.from.l.at.start = mean(abs(l_at_grid - starts$l) > 1e-12),
  mean.years.between.last.visit.and.grid.point = mean(last_visit))
#>                            visits.per.person 
#>                                   10.0400000 
#>                       grid.points.per.person 
#>                                    5.0000000 
#>             starts.l.differs.from.l.at.start 
#>                                    0.5677083 
#> mean.years.between.last.visit.and.grid.point 
#>                                    0.4598103
```

- **Timing.** About ten visits per person are reduced to five time
  points. A start at 1.1 and one at 1.9 years are both “started by year
  2”.
- **Which `l`.** For more than half of the people who started, the `l`
  used at the grid point is not the value measured at the visit where
  treatment was started. It comes from a later visit in the same year,
  measured *after* the start and so possibly affected by treatment,
  which the treatment model should not adjust for. Using the value at
  the start instead needs extra data preparation.
- **Carried-forward values.** On average the last visit was about half a
  year before the grid point, so for people who did not start, `l` is
  older than the time point it stands for.

A finer grid (for example monthly) reduces these problems, but means
many more time points with few or no visits, each with its own model and
many carried-forward values. Neither the grid analysis nor the long-data
analysis is wrong in general: with scheduled visits that everyone
attends, the grid is the design, and `weightitMSM()` offers the full
range of WeightIt methods at each time point. With visit times that
differ between people,
[`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) uses
each person’s own intervals directly.

### Standard errors on long data

With pooled logistic models
(`ipwtm(family = "binomial", link = "logit")`, any `type`),
[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
gives M-estimation standard errors for an MSM fitted on the long data
with [`lm()`](https://rdrr.io/r/stats/lm.html),
[`glm()`](https://rdrr.io/r/stats/glm.html) or a Cox model on
counting-process data, summing all scores within each person; a list of
weight objects, such as treatment and censoring weights, can be passed.
For the Cox exposure models used above (`family = "survival"`) it is not
available, and the cluster bootstrap remains the option.

When the visits are at the same times for everyone, the two packages can
be compared exactly: with time-specific columns in the pooled model the
weights equal those of `weightitMSM()`, and the package tests show that
[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
then reproduces `glm_weightit(vcov = "asympt")` for a person-level
outcome to about 1e-6.

The table summarises a simulation with individual visit times (300
people with 3 to 7 visits; a confounder that changes over time; a
treatment that switches back and forth and depends on the confounder and
on the previous treatment; an outcome at every visit with an effect of
current treatment of 1). In the “moderate” scenario treatment depends
moderately on the confounder and on the previous treatment; in the
“strong” scenario twice as strongly, which gives cumulative weights with
a heavy right tail (see the median of the largest weight). The weights
are `ipwtm(type = "all", numerator = ~ blag, denominator = ~ blag + l)`
and the MSM is `lm(y ~ b)` on all rows. The robust standard error is
clustered by person (HC0) and treats the weights as known; the cluster
bootstrap (200 samples of persons, weights re-estimated in each) was run
on a subset of the data sets.

``` r

sim_table(tm_sim$summary, c(Scenario = "scenario", `Median largest weight` = "median.max.weight"))
```

| scenario | median.max.weight | Method | Empirical SE | Model SE | Variance ratio | Relative error in model SE, % (MCSE) | Coverage | MCSE of coverage | Data sets |
|:---|---:|:---|---:|---:|---:|:---|---:|---:|---:|
| moderate | 4.152 | ipwvcov | 0.0574 | 0.0589 | 1.05 | +2.6 (2.3) | 0.958 | 0.006 | 1000 |
| moderate | 4.152 | robust | 0.0574 | 0.0826 | 2.07 | +43.8 (3.3) | 0.990 | 0.003 | 1000 |
| moderate | 4.152 | bootstrap | 0.0574 | 0.0586 | 1.04 | +2.0 (2.5) | 0.965 | 0.013 | 200 |
| strong | 12.013 | ipwvcov | 0.1012 | 0.0891 | 0.78 | -11.9 (2.6) | 0.930 | 0.008 | 1000 |
| strong | 12.013 | robust | 0.1012 | 0.1155 | 1.30 | +14.1 (3.2) | 0.973 | 0.005 | 1000 |
| strong | 12.013 | bootstrap | 0.1012 | 0.0932 | 0.85 | -7.9 (7.4) | 0.945 | 0.016 | 200 |

With moderate weights,
[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md) and
the cluster bootstrap agree with each other and with the empirical
variance of the estimates (variance ratio 1.05 and 1.04; model SE +3%
(MCSE 2) and +2% (MCSE 2)), and their intervals cover about 95%. The
robust variance, which treats the weights as known, is 2.07 times the
empirical variance here (model SE +44% (MCSE 3)), so its intervals are
much wider than needed. With strongly dependent treatment,
[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md) and
the bootstrap give similar standard errors, but they are too small: the
variance ratio of
[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md) is
0.78 (model SE -12% (MCSE 3)), and coverage drops to 93.0%. The
bootstrap’s variance ratio, 0.85 (model SE -8% (MCSE 7)), and coverage,
94.5%, rest on 200 data sets and are much less precise. The robust
standard error is conservative in both scenarios here, but that is not
guaranteed; the package’s simulation study found settings where it was
too small as well.

#### Heavy-tailed weights

Why do both methods fall short in the strong scenario? Cumulative
weights are products of many factors. When some people have a very low
probability of their observed treatment history (a near-violation of
positivity), the weights have a heavy right tail: P(W \> x) falls off
like $`x^{-\gamma}`$ for a tail index $`\gamma`$. The weights have a
finite variance only if $`\gamma > 2`$. The normal approximation also
relies on a finite third moment ($`\gamma > 3`$), and estimating the
variance well, by the sandwich or by the bootstrap, relies on a finite
fourth moment ($`\gamma > 4`$).

To check this explanation, the strong scenario and a milder version of
it (dependence on the confounder and on the previous treatment three
quarters as strong) were run with 300, 1200 and 4800 people. The weight
distribution refers to the final cumulative weight of each person: the
largest weight and the 99.9th percentile (medians over data sets), the
Kish effective sample size as a fraction of *n* (mean), and the Hill
estimate of the tail index from the top 5% (median). There were 1000,
1000 and 500 data sets for *n* = 300, 1200 and 4800, of which 200, 200
and 100 were bootstrapped (B = 200).

``` r

knitr::kable(setNames(tm_tail$weights,
                      c("Scenario", "n", "Largest weight", "99.9th percentile",
                        "ESS fraction", "Hill tail index")),
             digits = c(0, 0, 1, 1, 2, 2), row.names = FALSE)
```

| Scenario |    n | Largest weight | 99.9th percentile | ESS fraction | Hill tail index |
|:---------|-----:|---------------:|------------------:|-------------:|----------------:|
| strong   |  300 |            9.2 |               8.5 |         0.48 |            2.34 |
| strong   | 1200 |           17.2 |              11.5 |         0.44 |            2.29 |
| strong   | 4800 |           33.1 |              13.2 |         0.39 |            2.26 |
| milder   |  300 |            6.6 |               6.2 |         0.63 |            2.87 |
| milder   | 1200 |           10.7 |               7.6 |         0.61 |            2.83 |
| milder   | 4800 |           18.0 |               8.5 |         0.60 |            2.81 |

The standard errors, judged as above:

``` r

sim_table(tm_tail$summary, c(Scenario = "scenario", n = "n"))
```

| scenario | n | Method | Empirical SE | Model SE | Variance ratio | Relative error in model SE, % (MCSE) | Coverage | MCSE of coverage | Data sets |
|:---|---:|:---|---:|---:|---:|:---|---:|---:|---:|
| strong | 300 | ipwvcov | 0.0973 | 0.0891 | 0.84 | -8.4 (3.2) | 0.916 | 0.009 | 1000 |
| strong | 300 | bootstrap | 0.0973 | 0.0973 | 1.00 | -0.0 (11.9) | 0.950 | 0.015 | 200 |
| strong | 300 | robust | 0.0973 | 0.1150 | 1.40 | +18.2 (3.9) | 0.970 | 0.005 | 1000 |
| strong | 1200 | ipwvcov | 0.0499 | 0.0482 | 0.93 | -3.5 (2.9) | 0.932 | 0.008 | 1000 |
| strong | 1200 | bootstrap | 0.0499 | 0.0462 | 0.86 | -7.4 (4.5) | 0.945 | 0.016 | 200 |
| strong | 1200 | robust | 0.0499 | 0.0603 | 1.46 | +20.8 (3.3) | 0.975 | 0.005 | 1000 |
| strong | 4800 | ipwvcov | 0.0290 | 0.0275 | 0.89 | -5.4 (3.8) | 0.918 | 0.012 | 500 |
| strong | 4800 | bootstrap | 0.0290 | 0.0289 | 0.99 | -0.6 (7.3) | 0.900 | 0.030 | 100 |
| strong | 4800 | robust | 0.0290 | 0.0332 | 1.30 | +14.2 (4.2) | 0.980 | 0.006 | 500 |
| milder | 300 | ipwvcov | 0.0747 | 0.0707 | 0.90 | -5.4 (2.3) | 0.939 | 0.008 | 1000 |
| milder | 300 | bootstrap | 0.0747 | 0.0698 | 0.87 | -6.6 (2.5) | 0.935 | 0.017 | 200 |
| milder | 300 | robust | 0.0747 | 0.0961 | 1.65 | +28.6 (3.1) | 0.983 | 0.004 | 1000 |
| milder | 1200 | ipwvcov | 0.0404 | 0.0382 | 0.89 | -5.4 (2.8) | 0.934 | 0.008 | 1000 |
| milder | 1200 | bootstrap | 0.0404 | 0.0363 | 0.81 | -10.1 (2.5) | 0.945 | 0.016 | 200 |
| milder | 1200 | robust | 0.0404 | 0.0504 | 1.56 | +24.8 (3.3) | 0.987 | 0.004 | 1000 |
| milder | 4800 | ipwvcov | 0.0203 | 0.0193 | 0.90 | -5.1 (3.5) | 0.946 | 0.010 | 500 |
| milder | 4800 | bootstrap | 0.0203 | 0.0202 | 0.98 | -0.9 (7.1) | 0.930 | 0.026 | 100 |
| milder | 4800 | robust | 0.0203 | 0.0253 | 1.55 | +24.4 (4.3) | 0.996 | 0.003 | 500 |

The Hill estimate from one sample of 200,000 people per scenario, at
three tail fractions, with the largest weight and the effective sample
size fraction in that sample (“moderate” is the moderate scenario
above):

``` r

knitr::kable(tm_tail$hill_big, digits = 2)
```

|              | strong | milder | moderate |
|:-------------|-------:|-------:|---------:|
| hill_5pct    |   2.28 |   2.76 |     3.90 |
| hill_1pct    |   2.35 |   2.82 |     3.76 |
| hill_0.2pct  |   2.33 |   2.66 |     3.45 |
| max.weight   | 144.66 |  58.92 |    24.48 |
| ess.fraction |   0.40 |   0.59 |     0.80 |

The tail index is above 2 in every scenario, so the weights have a
finite variance and the intervals are valid in the limit. It is,
however, about 2.3 in the strong scenario and 2.8 in the milder one,
compared with 3.9 in the moderate scenario. The variance from
[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md) is
too small in both heavy-tailed scenarios, and by about the same amount
at every sample size: its ratio to the empirical variance is 0.84 to
0.93 in the strong scenario and 0.89 to 0.90 in the milder one (relative
error in model SE -8 to -3% and -5%, with MCSEs of 2 to 4 points),
against 1.05 in the moderate scenario. Coverage follows the tail index:
about 96% (moderate), 93.4 to 94.6% (milder) and 91.6 to 93.2% (strong).
Neither the variance ratio nor the coverage improves from 300 to 4800
people, while the largest weight roughly doubles with every fourfold
increase in *n* and the effective sample size fraction keeps falling.
The cluster bootstrap does no better: its variance ratio ranges from
0.81 to 1.00, but it rests on 100 to 200 data sets, so the MCSE of its
relative error in model SE is 2 to 12 points and that of its coverage
1.5 to 3.0 percentage points; it cannot be told apart from
[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md).
The robust variance is 1.30 to 1.65 times the empirical variance, so its
intervals over-cover (97.0 to 99.6%).

Heavy-tailed weights therefore make normal-theory intervals unreliable
for any method, the bootstrap included. The variance formula is not the
cause: it matches the numerical Jacobian of the estimating equations
exactly. Check the weight distribution
([`summary()`](https://rdrr.io/r/base/summary.html),
[`ipwplot()`](https://tranhung93.github.io/ipw/reference/ipwplot.md),
the effective sample size) before trusting any standard error. If a few
people carry much of the weight, consider truncation, a less extreme
weight model or a different estimand.

## 4. Which to use

- **WeightIt** for point treatments when you want machine-learning or
  balancing methods (GBM, SuperLearner, BART, CBPS, entropy balancing,
  energy balancing, …), many estimands, balance diagnostics through
  cobalt and M-estimation standard errors for outcome models
  (`glm_weightit()`, `coxph_weightit()`, …). Also for longitudinal
  treatments measured at fixed time points that are the same for
  everyone.
- **ipw** for exposures that change over time in long (counting-process)
  data with each person’s own visit times; Cox models for the time to
  the first switch or to censoring (`ipwtm(family = "survival")`);
  exposures that switch back and forth or recur (`type = "all"`); count
  and skewed continuous exposures (Poisson, negative binomial, gamma,
  log-normal); and analyses with few dependencies. With pooled logistic
  exposure models,
  [`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
  gives M-estimation standard errors on the long data.
- **PSweight** for point treatments when you want overlap, matching or
  entropy weights with built-in sandwich standard errors, augmented
  (doubly robust) estimators, ratio estimands for binary outcomes, or
  survey or clustered data.
- For a binary or categorical point treatment with a logistic model, the
  three give the same weights; for a binary one,
  [`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md),
  `glm_weightit(vcov = "asympt")` and PSweight give the same standard
  errors.

The
[README](https://tranhung93.github.io/ipw/index.html#ipw-or-weightit)
has a comparison table.

## References

- Greifer N. WeightIt: Weighting for Covariate Balance in Observational
  Studies. R package, <https://ngreifer.github.io/WeightIt/>.
- Li F, Li F (2019). Propensity score weighting for causal inference
  with multiple treatments. *Annals of Applied Statistics* 13,
  2389–2415.
- Lunceford JK, Davidian M (2004). Stratification and weighting via the
  propensity score in estimation of causal treatment effects: a
  comparative study. *Statistics in Medicine* 23, 2937–2960.
- Morris TP, White IR, Crowther MJ (2019). Using simulation studies to
  evaluate statistical methods. *Statistics in Medicine* 38, 2074–2102.
- Zhou T, Tong G, Li F, Thomas L, Li F. PSweight: Propensity Score
  Weighting for Causal Inference with Observational Studies and
  Randomized Trials. R package,
  <https://CRAN.R-project.org/package=PSweight>.
