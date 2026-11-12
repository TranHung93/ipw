# ipw and WeightIt side by side

``` r

library(ipw)
library(survival)
```

[WeightIt](https://ngreifer.github.io/WeightIt/) (Greifer) and ipw both
estimate weights for marginal structural models. This article runs the
two packages on the same data: a binary and a categorical point
treatment, where both apply and should give the same answer, and
longitudinal data with individual visit times, where they work
differently. WeightIt is used only through its exported functions. This
page was built with ipw 1.3.2.9000 and WeightIt 2.1.0.

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
```

### The weights

With a logistic model, both packages fit the same propensity score with
[`glm()`](https://rdrr.io/r/stats/glm.html) and apply the same formulas,
so the weights should agree up to rounding:

``` r

estimands <- c("ATE", "ATT", "ATC", "ATO")
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

### Standard errors

The marginal structural model is the weighted regression of `y` on `a`.
Four standard errors of the coefficient of `a`:

- **ipwvcov()**: ipw’s M-estimation standard error (new), which accounts
  for the estimation of the weights;
- **glm_weightit(vcov = “asympt”)**: WeightIt’s M-estimation standard
  error;
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
    robust.HC0 = sqrt(vcov(fit_r)["a", "a"]),
    bootstrap = unname(stored$boot_se[est]))
}))
round(se_table, 4)
#>     estimate ipwvcov glm_weightit.asympt robust.HC0 bootstrap
#> ATE   2.0259  0.0804              0.0804     0.1234    0.0836
#> ATT   2.1543  0.0862              0.0862     0.1134    0.0878
#> ATC   1.9160  0.1094              0.1094     0.1380    0.1131
#> ATO   2.0807  0.0684              0.0684     0.0962    0.0703
```

[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md) and
`glm_weightit(vcov = "asympt")` agree to many decimals: they solve the
same stacked estimating equations (Lunceford & Davidian 2004). Both are
close to the bootstrap. The robust standard error, which ignores that
the weights were estimated, is larger for all four estimands here, so
confidence intervals based on it would be wider than needed. For the ATE
this is the usual direction (estimating the weights reduces the
variance); for the ATT, ATC and ATO the robust standard error is not
guaranteed to be conservative and can also be too small. The
M-estimation or bootstrap standard error is therefore preferred.

[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
covers binary point treatments with
[`lm()`](https://rdrr.io/r/stats/lm.html) or
[`glm()`](https://rdrr.io/r/stats/glm.html) outcome models. WeightIt’s
outcome-model functions cover more: multi-category and continuous
treatments, many weighting methods, Cox, ordinal and multinomial outcome
models, and longitudinal weights from `weightitMSM()`.

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
the formulas are the same. Within each level, the ATT weights make the
covariate means resemble those of the focal group:

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
- For a binary or categorical point treatment with a logistic model, the
  two give the same weights; for a binary one,
  [`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
  and `glm_weightit(vcov = "asympt")` give the same standard errors.

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
