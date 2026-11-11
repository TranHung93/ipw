# Bugs fixed after ipw 1.3.2

This page describes the bugs that were found in a review of ipw after
version 1.3.2 and fixed in the next version. For each bug it shows a
small example and runs it in three versions:

- **ipw 1.2.1.1**, the last version before the 1.3 rewrite;
- **ipw 1.3.2**, the version on CRAN when the review was done;
- **this version**, the one this page was built with.

The weights of each version are compared with the *correct* weights,
which are computed by hand on this page with
[`glm()`](https://rdrr.io/r/stats/glm.html),
[`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html) or
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html),
without using ipw. A difference of about 1e-15 is rounding; a difference
of 0.1 means that some weights are off by 10%.

The results of 1.2.1.1 and 1.3.2 were made once in separate R sessions
with the script `data-raw/make-bugfix-article-data.R` and are stored
with this article; the correct weights and the results of this version
are computed every time the page is built.

## Who is affected?

You are affected if you used one of the following with an older version.
The sections below explain each case; the table at the end summarises
them.

- an [`offset()`](https://rdrr.io/r/stats/offset.html) term in a model
  formula (1.3.1 and 1.3.2);
- `denominator = ~ .` (1.3.1 and 1.3.2: an error, so no wrong results);
- a 0/1 exposure stored as a factor, with the level `"1"` first (all
  older versions);
- an ordinal exposure stored as text (1.3.1 and 1.3.2; 1.2.1.1 stopped
  with an unclear error), or as a factor whose levels are in
  alphabetical rather than dose order (all older versions);
- [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) with
  `family = "survival"` and `type = "first"` or `"cens"`, where subjects
  have their own visit times (all older versions);
- [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) with
  a covariate that is called `selvar` (1.3.1 and 1.3.2);
- [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) with
  `id`/`timevar` combinations that occur more than once (all older
  versions);
- [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) with
  `family = "survival"` and the `ties` argument (all older versions).

## 1. `offset()` terms were dropped (1.3.1 and 1.3.2)

An offset is a term with a fixed coefficient of 1, for example the log
of follow-up time when the chance of exposure grows with time at risk.
In ipw 1.3.1 and 1.3.2 the model formula was rebuilt from its term
labels, and R does not count
[`offset()`](https://rdrr.io/r/stats/offset.html) as a term label, so
the offset was silently left out: the model was fitted as `a ~ l`.

``` r

cases$offset$call
#> ipwpoint(exposure = a, family = "binomial", link = "logit", denominator = ~l + 
#>     offset(off), data = data)
```

The correct weights, by hand:

``` r

d <- cases$offset$data
p <- fitted(glm(a ~ l + offset(off), family = binomial, data = d))
ref_offset <- unname(1 / ifelse(d$a == 1, p, 1 - p))
```

| Version | Result | Largest relative difference from the correct weights |
|:---|:---|:---|
| 1.2.1.1 | Weights | 4.4e-16 |
| 1.3.2 (CRAN) | Weights | 1.1 |
| this version (1.3.2.9000) | Weights | 0 |

**Fix:** the model formula is now used exactly as written.

## 2. `denominator = ~ .` failed (1.3.1 and 1.3.2)

`~ .` means “all other columns of `data`”. In 1.3.1 and 1.3.2 the check
for missing variables looked for a column called `.` and stopped.

``` r

cases$dot$call
#> ipwpoint(exposure = a, family = "binomial", link = "logit", denominator = ~., 
#>     data = data)
names(cases$dot$data)
#> [1] "l1" "l2" "a"
```

``` r

d <- cases$dot$data
p <- fitted(glm(a ~ l1 + l2, family = binomial, data = d))
ref_dot <- unname(1 / ifelse(d$a == 1, p, 1 - p))
```

| Version | Result | Largest relative difference from the correct weights |
|:---|:---|:---|
| 1.2.1.1 | Weights | 8.9e-16 |
| 1.3.2 (CRAN) | Error: Variables not found in data: . |  |
| this version (1.3.2.9000) | Weights | 0 |

**Fix:** `~ .` works again in
[`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md).
In [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) it
is now a clear error, because the dot would also include the id, time
and exposure columns; list the covariates there.

## 3. A factor exposure with levels `c("1", "0")` swapped the probabilities

For a binomial exposure, ipw checked that the values were 0 and 1, but a
factor with those labels passed the check.
[`glm()`](https://rdrr.io/r/stats/glm.html) then models the probability
of the *second* level. With levels in the order `"1"`, `"0"` that is the
probability of being *unexposed*, so every exposed subject got the
weight meant for an unexposed one and vice versa. Factors with levels
`"0"`, `"1"` (the default order) were fine.

``` r

cases$factor$call
#> ipwpoint(exposure = a, family = "binomial", link = "logit", denominator = ~l1 + 
#>     l2, data = data)
levels(cases$factor$data$a)
#> [1] "1" "0"
```

``` r

d <- cases$factor$data
p <- fitted(glm(a_num ~ l1 + l2, family = binomial, data = d))
ref_factor <- unname(1 / ifelse(d$a_num == 1, p, 1 - p))
```

| Version | Result | Largest relative difference from the correct weights |
|:---|:---|:---|
| 1.2.1.1 | Weights | 8.7 |
| 1.3.2 (CRAN) | Weights | 8.7 |
| this version (1.3.2.9000) | Weights | 0 |

**Fix:** a 0/1 exposure stored as a factor or text is converted by its
labels, so the order of the levels no longer matters.

## 4. An ordinal exposure given as text was ordered alphabetically (1.3.1 and 1.3.2)

For `family = "ordinal"` the exposure is turned into a factor. Text is
then sorted alphabetically, so `"low"`, `"mid"`, `"high"` became
`"high" < "low" < "mid"`, and in 1.3.1 and 1.3.2 the ordinal model was
fitted with the wrong order without any message. ipw 1.2.1.1 stopped
with an error that did not say what was wrong. The same alphabetical
order is used in every version for an *unordered factor* with text
levels, unless its levels were set in dose order.

``` r

cases$ordinal$call
#> ipwpoint(exposure = dose, family = "ordinal", link = "logit", 
#>     denominator = ~l, data = data)
table(cases$ordinal$data$dose)
#> 
#> high  low  mid 
#>  146  193  161
```

The correct weights use the dose order:

``` r

d <- cases$ordinal$data
d$dose_ord <- factor(d$dose, levels = c("low", "mid", "high"), ordered = TRUE)
fit <- MASS::polr(dose_ord ~ l, data = d, method = "logistic")
pr <- predict(fit, type = "probs")
ref_ordinal <- unname(1 / pr[cbind(seq_len(nrow(d)), as.integer(d$dose_ord))])
```

| Version | Result | Largest relative difference from the correct weights |
|:---|:---|:---|
| 1.2.1.1 | Error: replacement has 0 rows, data has 500 |  |
| 1.3.2 (CRAN) | Weights | 1.8 |
| this version (1.3.2.9000) | Error: Ordinal exposure is a character vector, whose order is not defined (it would be sorted alphabetically). Use an… |  |

**Fix:** a text exposure is now an error that asks for an ordered factor
or numeric codes, and an unordered factor whose levels are not numbers
gives a warning. With an ordered factor the weights are correct in all
versions:

``` r

w <- ipwpoint(exposure = dose_ord, family = "ordinal", link = "logit",
              denominator = ~ l, data = d)
max_rel(w$ipw.weights, ref_ordinal)
#> [1] 0
```

## 5. Survival exposure: subjects with their own visit times

For `family = "survival"` with `type = "first"` (and the censoring
weights with `type = "cens"`), each row covers an interval
`(tstart, timevar]`, and its probability should use the baseline hazard
over that interval. Older versions used the increase in the baseline
hazard since the *previous time point in the whole data set*. That is
the same thing when all subjects are seen at the same times, as in
`haartdat`, but not when they have their own visit times: the hazard of
a row was then only counted from the last visit of *any* subject, so it
was too small.

``` r

cases$survival$call
#> ipwtm(exposure = a, family = "survival", denominator = ~l, id = id, 
#>     tstart = tstart, timevar = time, type = "first", data = data)
head(cases$survival$data[, c("id", "tstart", "time", "a", "l")], 4)
#>   id tstart time a          l
#> 1  1      0    5 0  1.3774531
#> 2  1      5   14 0  0.1963246
#> 3  1     14   28 1 -0.2319014
#> 4  2      0   13 0  0.9947626
```

The correct weights use the cumulative baseline hazard *H* at both ends
of each row’s interval:

``` r

d <- cases$survival$data
fit <- coxph(Surv(tstart, time, a) ~ l, data = d)
bh <- basehaz(fit, centered = TRUE)
H <- function(t) {
  h <- approx(bh$time, bh$hazard, xout = t, method = "constant", rule = 2)$y
  h[t < min(bh$time)] <- 0
  h
}
haz <- (H(d$time) - H(d$tstart)) * predict(fit, type = "risk")
p_obs <- ifelse(d$a == 1, 1 - exp(-haz), exp(-haz))
ref_survival <- unname(ave(1 / p_obs, d$id, FUN = cumprod))
```

(The package tests also check this against a baseline hazard computed
from the risk sets by hand, without
[`basehaz()`](https://rdrr.io/pkg/survival/man/basehaz.html).)

| Version | Result | Largest relative difference from the correct weights |
|:---|:---|:---|
| 1.2.1.1 | Weights | 25 |
| 1.3.2 (CRAN) | Weights | 25 |
| this version (1.3.2.9000) | Weights | 4.4e-16 |

With a common time grid, as in `haartdat`, nothing changes: all versions
give the same weights.

**Fix:** each row now uses the hazard over its own interval, as
`type = "all"` already did. The article [Survival weights with
individual visit
times](https://tranhung93.github.io/ipw/articles/survival-weights-bug.md)
describes this bug in detail: every version from ipw 1.0 to 1.3.2 is
affected, and it shows the old and the new code.

## 6. A covariate called `selvar` was overwritten (`ipwtm()`, 1.3.1 and 1.3.2)

From version 1.3.1,
[`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) adds a
column called `selvar` to the data to mark the rows that enter the
exposure model. A covariate of the user with that name was overwritten,
so the model used the indicator instead.

``` r

cases$selvar$call
#> ipwtm(exposure = a, family = "binomial", link = "logit", denominator = ~l + 
#>     selvar, id = id, timevar = time, type = "first", data = data)
```

By hand: the rows up to and including the first exposure enter the
model, and the weights are cumulative products within subject.

``` r

d <- cases$selvar$data
sel <- ave(d$a, d$id, FUN = function(x) cumsum(cumsum(x)) <= 1) == 1
p <- rep(1, nrow(d))
p[sel] <- fitted(glm(a ~ l + selvar, family = binomial, data = d[sel, ]))
p_obs <- ifelse(!sel, 1, ifelse(d$a == 1, p, 1 - p))
ref_selvar <- unname(ave(1 / p_obs, d$id, FUN = cumprod))
```

| Version | Result | Largest relative difference from the correct weights |
|:---|:---|:---|
| 1.2.1.1 | Weights | 1.2e-15 |
| 1.3.2 (CRAN) | Weights | 3.6 |
| this version (1.3.2.9000) | Weights | 4.4e-16 |

**Fix:** the internal column now has a name that does not clash
(`.ipw_selvar`); the `selvar` element of the result is unchanged.

## 7. Duplicated `id`/`timevar` rows were accepted (`ipwtm()`)

If the same subject and time occur twice, for example after a faulty
merge, the cumulative product is taken over both copies in an arbitrary
order and the weights are wrong.
[`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md)
accepted such data without a message
([`ipwprobtm()`](https://tranhung93.github.io/ipw/reference/ipwprobtm.md)
already refused it).

``` r

d <- cases$duplicates$data
sum(duplicated(d[, c("id", "time")]))
#> [1] 3
```

| Version | Result | Largest relative difference from the correct weights |
|:---|:---|:---|
| 1.2.1.1 | Weights |  |
| 1.3.2 (CRAN) | Weights |  |
| this version (1.3.2.9000) | Error: Each combination of ‘id’ and ‘timevar’ must occur only once. |  |

**Fix:** duplicated `id`/`timevar` combinations are now an error.

## 8. The `ties` argument was ignored (`ipwtm()`, survival)

Arguments in `...` are passed to the model function. For
`family = "survival"` the Cox model always used Efron’s method for tied
times, so `ties = "breslow"` had no effect, in 1.2.1.1 as well as in
1.3.x.

``` r

cases$ties$call
#> ipwtm(exposure = haartind, family = "survival", numerator = ~sex + 
#>     age, denominator = ~sex + age + cd4.sqrt, id = patient, tstart = tstart, 
#>     timevar = fuptime, type = "first", data = data, ties = "breslow")
```

By hand, with Breslow’s method for both models:

``` r

d <- haartdat
sel <- ave(d$haartind, d$patient, FUN = function(x) cumsum(cumsum(x)) <= 1) == 1
p_breslow <- function(rhs) {
  f <- update(rhs, Surv(tstart, fuptime, haartind) ~ .)
  fit <- coxph(f, data = d[sel, ], ties = "breslow")
  bh <- basehaz(fit, centered = TRUE)
  H <- function(t) {
    h <- approx(bh$time, bh$hazard, xout = t, method = "constant", rule = 2)$y
    h[t < min(bh$time)] <- 0
    h
  }
  p <- rep(1, nrow(d))
  haz <- (H(d$fuptime[sel]) - H(d$tstart[sel])) * predict(fit, type = "risk")
  p[sel] <- ifelse(d$haartind[sel] == 1, 1 - exp(-haz), exp(-haz))
  p
}
ref_ties <- unname(ave(p_breslow(~ sex + age) / p_breslow(~ sex + age + cd4.sqrt),
                       d$patient, FUN = cumprod))
```

| Version | Result | Largest relative difference from the correct weights |
|:---|:---|:---|
| 1.2.1.1 | Weights | 0.064 |
| 1.3.2 (CRAN) | Weights | 0.064 |
| this version (1.3.2.9000) | Weights | 1.1e-15 |

**Fix:** `ties` is passed on to
[`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html); the default is
still Efron.

## Summary

| Case | 1.2.1.1 | 1.3.2 (CRAN) | this version (1.3.2.9000) |
|:---|:---|:---|:---|
| 1\. offset() in a formula | correct | wrong (up to 1.1) | correct |
| 2\. denominator = ~ . | correct | error (should work) | correct |
| 3\. factor exposure, levels “1”, “0” | wrong (up to 8.7) | wrong (up to 8.7) | correct |
| 4\. ordinal exposure as text | error (unclear message) | wrong (up to 1.8) | error, as it should |
| 5\. survival, own visit times | wrong (up to 25) | wrong (up to 25) | correct |
| 6\. covariate called selvar | correct | wrong (up to 3.6) | correct |
| 7\. duplicated id/timevar rows | wrong (no error) | wrong (no error) | error, as it should |
| 8\. ties = “breslow” | wrong (up to 0.064) | wrong (up to 0.064) | correct |

All fixes are listed in the
[NEWS](https://tranhung93.github.io/ipw/news/index.md). The package
tests (`tests/testthat/test-audit*.R`, `test-regression*.R`) check each
case, and `test-regression.R` checks that, apart from these fixes, the
weights are the same as in ipw 1.2.1.1 for 44 settings.
