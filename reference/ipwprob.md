# Inverse Probability Weights From User-Supplied Probabilities

Computes inverse probability weights from probabilities (or densities)
of the observed exposure that you estimated yourself, for example with a
machine-learning model or SuperLearner. Point treatment weights are
`num / den`; with `id` and `timevar`, the ratios are multiplied
cumulatively over time within each subject, as in
[`ipwtm`](https://tranhung93.github.io/ipw/reference/ipwtm.md).

## Usage

``` r
ipwprob(den, num = NULL, id = NULL, timevar = NULL, trunc = NULL)
```

## Arguments

- den:

  Numeric vector: for each row, the probability (or density) of the
  observed exposure given the covariates of the denominator model.

- num:

  Optional numeric vector of the same length: the probability (or
  density) of the observed exposure given the covariates of the
  numerator model. If `NULL`, unstabilized weights (numerator = 1) are
  computed.

- id:

  Optional vector of subject identifiers, same length as `den`. When
  given, weights are cumulated over time within each subject.

- timevar:

  Numeric vector of follow-up times, same length as `den`; required when
  `id` is given. It sets the order in which the ratios are multiplied
  within a subject.

- trunc:

  Optional truncation percentile between 0 and 0.5 (e.g., 0.01 for
  1st/99th percentiles).

## Value

A list of class `"ipwprob"` containing:

- ipw.weights:

  Numeric vector of inverse probability weights, in the order of the
  input.

- weights.trunc:

  Numeric vector of truncated weights (if `trunc` is specified).

- call:

  The original function call.

## Details

For time-varying weights of `type = "first"` or `"cens"` (see
[`ipwtm`](https://tranhung93.github.io/ipw/reference/ipwtm.md)), set the
probabilities to 1 for the rows after the first switch, so that those
rows do not change the cumulative weight.

## See also

[`ipwpoint`](https://tranhung93.github.io/ipw/reference/ipwpoint.md),
[`ipwtm`](https://tranhung93.github.io/ipw/reference/ipwtm.md).

## Author

Hung Thai Tran <hungtt@oucru.org>

## Examples

``` r
# Point treatment: probabilities from any model, here a logistic regression
set.seed(1)
d <- data.frame(l = rnorm(200))
d$a <- rbinom(200, 1, plogis(d$l))
p <- fitted(glm(a ~ l, family = binomial, data = d))
w <- ipwprob(den = ifelse(d$a == 1, p, 1 - p),
             num = ifelse(d$a == 1, mean(d$a), 1 - mean(d$a)))
summary(w$ipw.weights)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.5892  0.7849  0.9089  0.9981  1.1440  2.6291 

# Time-varying: the ratios are multiplied within each subject over time
long <- data.frame(id = rep(1:3, each = 2), time = rep(1:2, 3),
                   p = c(0.5, 0.8, 0.4, 0.9, 0.6, 0.7))
ipwprob(den = long$p, id = long$id, timevar = long$time)$ipw.weights
#> [1] 2.000000 2.500000 2.500000 2.777778 1.666667 2.380952
```
