# Bootstrap Standard Errors for Inverse Probability Weighted Analyses

Bootstraps a whole inverse probability weighted analysis: in each
bootstrap sample the weights are re-estimated and the marginal
structural model is refitted, so the standard error accounts for the
uncertainty in the estimated weights. For data in long format (several
rows per subject), whole subjects are resampled and each sampled copy
gets a new identifier.

## Usage

``` r
ipwboot(data, statistic, id = NULL, R = 200, strata = NULL, conf = 0.95, ...)

# S3 method for class 'ipwboot'
print(x, digits = 4, ...)

# S3 method for class 'ipwboot'
confint(object, parm, level = object$conf, ...)
```

## Arguments

- data:

  Data frame with the original data.

- statistic:

  Function that takes a data frame (a bootstrap sample of `data`) as its
  first argument and returns a numeric vector of estimates, for example
  the coefficients of a weighted outcome model. It should compute the
  weights itself, e.g. by calling
  [`ipwpoint`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  or [`ipwtm`](https://tranhung93.github.io/ipw/reference/ipwtm.md).

- id:

  Optional name of the subject identifier column (unquoted or as a
  string). If given, subjects are resampled with all their rows, and in
  each bootstrap sample the column is replaced by new identifiers
  `1, 2, ...`, one per sampled copy, so that a subject drawn twice
  counts as two subjects. If `NULL`, rows are resampled.

- R:

  Number of bootstrap samples. About 200 are usually enough for a
  standard error; use more (e.g. 1000) for percentile confidence
  intervals.

- strata:

  Optional name of a column (unquoted or as a string) defining strata;
  resampling is then done within each stratum, keeping the number of
  subjects (or rows) per stratum fixed. With `id`, the column must be
  constant within each subject.

- conf:

  Confidence level of the percentile intervals.

- ...:

  Further arguments passed to `statistic`; for the methods, ignored.

- x, object:

  An object of class `"ipwboot"`.

- digits:

  Number of significant digits to print.

- parm:

  Estimates (names or positions) for which to give intervals; all by
  default.

- level:

  Confidence level.

## Value

A list of class `"ipwboot"` containing:

- t0:

  Estimates from the original data.

- t:

  Matrix of bootstrap estimates, one row per bootstrap sample (`NA` for
  failed samples).

- se:

  Bootstrap standard errors.

- ci:

  Percentile confidence intervals.

- R, conf:

  Number of bootstrap samples and confidence level.

- n.fail:

  Number of failed bootstrap samples.

- errors:

  Distinct error messages of failed samples.

- id, strata:

  Names of the identifier and strata columns, if used.

- call:

  The original function call.

## Details

Austin (2016) found in simulations that, for weighted Cox models, naive
standard errors were too small and robust (sandwich) standard errors too
large, while the bootstrap in which the weights are re-estimated in each
sample gave correct standard errors and confidence intervals. Austin
studied a point treatment with baseline confounders; for longitudinal
data the same idea applies with subjects as the resampling unit.

The package's simulation study (article "Simulation study" on the
package website) shows when the bootstrap helps. For a point treatment
with good overlap, the bootstrap standard error was accurate and robust
standard errors were about 10% too large. For a weighted Cox model with
good overlap, naive, robust and bootstrap standard errors were all close
to correct. When some weights were extreme, through poor overlap or
through products of time-varying weights, all standard errors were too
small, the bootstrap ones included, and 95% confidence intervals covered
the truth in only 60 to 80% of samples. The bootstrap does not repair
unstable weights: check the weights first with
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`ipwplot`](https://tranhung93.github.io/ipw/reference/ipwplot.md), and
consider truncation, overlap weights or a different exposure model.

Bootstrap samples in which `statistic` fails, or returns a result of a
different length, are left out of the summaries; their number is
reported. Use [`set.seed()`](https://rdrr.io/r/base/Random.html) before
calling `ipwboot` to make the results reproducible.

## References

Austin, P.C. (2016). Variance estimation when using inverse probability
of treatment weighting (IPTW) with survival analysis. *Statistics in
Medicine*, **35**, 5642-5655.
[doi:10.1002/sim.7084](https://doi.org/10.1002/sim.7084) .

## See also

[`ipwpoint`](https://tranhung93.github.io/ipw/reference/ipwpoint.md),
[`ipwtm`](https://tranhung93.github.io/ipw/reference/ipwtm.md).

## Author

Hung Thai Tran <hungtt@oucru.org>

## Examples

``` r
# Point treatment: true causal effect of a on y is 10
set.seed(16)
n <- 500
simdat <- data.frame(l = rnorm(n, 10, 5))
simdat$a <- rbinom(n, 1, plogis(simdat$l - 10))
simdat$y <- 10 * simdat$a + 0.5 * simdat$l + rnorm(n, -10, 5)

msm.effect <- function(d) {
  w <- ipwpoint(exposure = a, family = "binomial", link = "logit",
                numerator = ~ 1, denominator = ~ l, data = d)$ipw.weights
  coef(lm(y ~ a, data = d, weights = w))["a"]
}
set.seed(1)
res <- ipwboot(simdat, msm.effect, R = 50)  # use about 200 in practice
res
#> Bootstrap of an inverse probability weighted analysis
#> 50 bootstrap samples, resampling rows
#> 
#>   estimate     se 2.5 % 97.5 %
#> a    11.41 0.8103 10.19  12.94
confint(res)
#>      2.5 %   97.5 %
#> a 10.18977 12.93739

# \donttest{
# Time-varying treatment: resample patients, re-estimate both weights
data(haartdat)
haart.effect <- function(d) {
  w1 <- ipwtm(exposure = haartind, family = "survival", numerator = ~ sex + age,
              denominator = ~ sex + age + cd4.sqrt, id = patient,
              tstart = tstart, timevar = fuptime, type = "first", data = d)
  w2 <- ipwtm(exposure = dropout, family = "survival", numerator = ~ sex + age,
              denominator = ~ sex + age + cd4.sqrt, id = patient,
              tstart = tstart, timevar = fuptime, type = "cens", data = d)
  fit <- survival::coxph(survival::Surv(tstart, fuptime, event) ~ haartind,
                         data = d, weights = w1$ipw.weights * w2$ipw.weights)
  coef(fit)
}
set.seed(2)
ipwboot(haartdat, haart.effect, id = patient, R = 20)  # use about 200 in practice
#> Bootstrap of an inverse probability weighted analysis
#> 20 bootstrap samples, resampling subjects ('patient')
#> 
#>          estimate    se  2.5 %  97.5 %
#> haartind  -0.9363 0.362 -1.495 -0.2691
# }
```
