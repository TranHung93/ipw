# IQ, Income and Health

A simulated dataset containing measurements of IQ, monthly income, and
health scores for 1000 individuals.

## Usage

``` r
data(healthdat)
```

## Format

A data frame with 1000 rows, with each row corresponding to a separate
individual. The following 4 variables are included:

- id:

  Individual ID.

- iq:

  IQ score.

- income:

  Gross monthly income (EUR).

- health:

  Health score (0-100).

## Details

In these simulated data, IQ acts as a confounder for the causal effect
of income on health. This dataset is primarily used to illustrate
point-treatment inverse probability weighting.

## References

Van der Wal W.M. & Geskus R.B. (2011). ipw: An R Package for Inverse
Probability Weighting. *Journal of Statistical Software*, **43**(13),
1-23. [doi:10.18637/jss.v043.i13](https://doi.org/10.18637/jss.v043.i13)

## See also

[`basdat`](https://tranhung93.github.io/ipw/reference/basdat.md),
[`haartdat`](https://tranhung93.github.io/ipw/reference/haartdat.md),
[`ipwplot`](https://tranhung93.github.io/ipw/reference/ipwplot.md),
[`ipwpoint`](https://tranhung93.github.io/ipw/reference/ipwpoint.md),
[`ipwtm`](https://tranhung93.github.io/ipw/reference/ipwtm.md),
[`timedat`](https://tranhung93.github.io/ipw/reference/timedat.md),
[`tstartfun`](https://tranhung93.github.io/ipw/reference/tstartfun.md)

## Author

Willem M. van der Wal <willem@vanderwalresearch.com>, Ronald B. Geskus
<rgeskus@oucru.org>

## Examples

``` r
# For an example of how to use this data with ipwpoint, see:
# ?ipwpoint
```
