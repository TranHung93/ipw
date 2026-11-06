# Time-Varying Inverse Probability Weights From User-Supplied Probabilities

Computes time-varying inverse probability weights from probabilities (or
densities) of the observed exposure that you estimated yourself, for
example with a machine-learning model or SuperLearner. For each row the
ratio `num / den` is computed, and the ratios are multiplied
cumulatively over time within each subject, as in
[`ipwtm`](https://tranhung93.github.io/ipw/reference/ipwtm.md).

## Usage

``` r
ipwprobtm(den, num = NULL, id, timevar, trunc = NULL)
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

  Vector of subject identifiers, same length as `den`.

- timevar:

  Numeric vector of follow-up times, same length as `den`. It sets the
  order in which the ratios are multiplied within a subject, so the data
  do not need to be sorted.

- trunc:

  Optional truncation percentile between 0 and 0.5 (e.g., 0.01 for
  1st/99th percentiles).

## Value

A list of class `"ipwprobtm"` containing:

- ipw.weights:

  Numeric vector of inverse probability weights, in the order of the
  input.

- weights.trunc:

  Numeric vector of truncated weights (if `trunc` is specified).

- call:

  The original function call.

## Details

`ipwprobtm` does the last step of
[`ipwtm`](https://tranhung93.github.io/ipw/reference/ipwtm.md) only: the
exposure model is fitted by you, and `ipwprobtm` forms the ratios,
multiplies them within each subject in time order, returns the weights
in the order of the input and optionally truncates them.

For weights of `type = "first"` or `"cens"` (see
[`ipwtm`](https://tranhung93.github.io/ipw/reference/ipwtm.md)), set the
probabilities to 1 for the rows after the first switch, so that those
rows do not change the cumulative weight.

For a point treatment no function is needed: with your own probabilities
the weights are simply `num / den` (or `1 / den` for unstabilized
weights).

## See also

[`ipwtm`](https://tranhung93.github.io/ipw/reference/ipwtm.md),
[`ipwpoint`](https://tranhung93.github.io/ipw/reference/ipwpoint.md).

## Author

Hung Thai Tran <hungtt@oucru.org>

## Examples

``` r
# Probabilities from any model; the ratios are multiplied within each
# subject over time
long <- data.frame(id = rep(1:3, each = 2), time = rep(1:2, 3),
                   p = c(0.5, 0.8, 0.4, 0.9, 0.6, 0.7))
ipwprobtm(den = long$p, id = long$id, timevar = long$time)$ipw.weights
#> [1] 2.000000 2.500000 2.500000 2.777778 1.666667 2.380952

# Stabilized weights from a pooled logistic model fitted by hand; this
# reproduces ipwtm(..., type = "all")
set.seed(1)
n <- 200
long <- data.frame(id = rep(1:n, each = 3), time = rep(1:3, n))
long$l <- rnorm(3 * n)
long$a <- rbinom(3 * n, 1, plogis(-0.5 + long$l))
pd <- fitted(glm(a ~ l, family = binomial, data = long))
pn <- fitted(glm(a ~ 1, family = binomial, data = long))
w <- ipwprobtm(den = ifelse(long$a == 1, pd, 1 - pd),
               num = ifelse(long$a == 1, pn, 1 - pn),
               id = long$id, timevar = long$time)
summary(w$ipw.weights)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.2055  0.5656  0.7844  0.9958  1.1332 15.8338 
```
