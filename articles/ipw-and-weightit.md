# ipw, WeightIt and PSweight side by side

``` r

library(ipw)
library(survival)
```

[WeightIt](https://ngreifer.github.io/WeightIt/) (Greifer) and ipw both
estimate weights for marginal structural models;
[PSweight](https://CRAN.R-project.org/package=PSweight) (Zhou, Tong, Li,
Thomas & Li) estimates propensity score weighted treatment effects for
point treatments. This article runs the packages on the same data: a
binary and a categorical point treatment, where all three apply and
should give the same answer, and longitudinal data with individual visit
times, where ipw and WeightIt work differently. WeightIt and PSweight
are used only through their exported functions. This page was built with
ipw 1.3.2.9000, WeightIt 2.1.0 and PSweight 2.1.2.

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
#> max.abs.difference 1.065814e-14 2.220446e-15 8.881784e-15 2.220446e-16
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
#> 1.734723e-18 4.336809e-19 0.000000e+00 1.734723e-18
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
run on a subset of the data sets.

``` r

knitr::kable(cox_sim$summary, digits = 3, row.names = FALSE)
```

| estimand | truth | mean.estimate | empirical.sd | mean.se.ipwvcov | mean.se.robust | mean.se.bootstrap | coverage.ipwvcov | coverage.robust | coverage.bootstrap | replications | bootstrap.replications |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ATE | -0.450 | -0.447 | 0.128 | 0.128 | 0.14 | 0.128 | 0.956 | 0.969 | 0.95 | 1000 | 300 |
| ATO | -0.461 | -0.460 | 0.133 | 0.132 | 0.14 | 0.133 | 0.956 | 0.961 | 0.95 | 1000 | 300 |

The mean
[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
standard error is within 1% of the empirical standard deviation of the
estimates for both estimands, and so is the bootstrap’s. The robust
standard error, which treats the weights as known, is about 10% (ATE)
and 6% (ATO) too large here, so its intervals cover more often than 95%:
conservative, but wider than needed. The coverages of
[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
(1000 data sets) and of the bootstrap (300 data sets) are both within
simulation error of 95%.

[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
covers binary point treatments with
[`lm()`](https://rdrr.io/r/stats/lm.html),
[`glm()`](https://rdrr.io/r/stats/glm.html) or
[`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) outcome models.
WeightIt’s outcome-model functions cover more: multi-category and
continuous treatments, many weighting methods, ordinal and multinomial
outcome models, and longitudinal weights from `weightitMSM()`.

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
  log-normal); and analyses with few dependencies.
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
- Zhou T, Tong G, Li F, Thomas L, Li F. PSweight: Propensity Score
  Weighting for Causal Inference with Observational Studies and
  Randomized Trials. R package,
  <https://CRAN.R-project.org/package=PSweight>.
