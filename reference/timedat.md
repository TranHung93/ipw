# HIV: TB and Survival (Longitudinal Measurements)

A simulated dataset containing time-varying CD4 measurements for 386
HIV-positive individuals. Corresponding baseline data, including timing
of tuberculosis and death, are available in
[`basdat`](https://tranhung93.github.io/ipw/reference/basdat.md).

## Usage

``` r
data(timedat)
```

## Format

A data frame with 6291 observations on the following 3 variables:

- id:

  Patient ID.

- fuptime:

  Follow-up time (days since HIV seroconversion).

- cd4count:

  CD4 count measured at `fuptime`.

## Details

These simulated data are used together with
[`basdat`](https://tranhung93.github.io/ipw/reference/basdat.md) in a
detailed causal modeling example using inverse probability weighting
(IPW). See
[`ipwtm`](https://tranhung93.github.io/ipw/reference/ipwtm.md) for the
full example. Data were simulated using the algorithm described in Van
der Wal et al. (2009).

## References

Cole, S.R. & Hernán, M.A. (2008). Constructing inverse probability
weights for marginal structural models. *American Journal of
Epidemiology*, **168**(6), 656-664.

Robins, J.M., Hernán, M.A. & Brumback, B.A. (2000). Marginal structural
models and causal inference in epidemiology. *Epidemiology*, **11**,
550-560.

Van der Wal W.M. & Geskus R.B. (2011). ipw: An R Package for Inverse
Probability Weighting. *Journal of Statistical Software*, **43**(13),
1-23. [doi:10.18637/jss.v043.i13](https://doi.org/10.18637/jss.v043.i13)

Van der Wal W.M., Prins M., Lumbreras B. & Geskus R.B. (2009). A simple
G-computation algorithm to quantify the causal effect of a secondary
illness on the progression of a chronic disease. *Statistics in
Medicine*, **28**(18), 2325-2337.

## See also

[`basdat`](https://tranhung93.github.io/ipw/reference/basdat.md),
[`haartdat`](https://tranhung93.github.io/ipw/reference/haartdat.md),
[`ipwplot`](https://tranhung93.github.io/ipw/reference/ipwplot.md),
[`ipwpoint`](https://tranhung93.github.io/ipw/reference/ipwpoint.md),
[`ipwtm`](https://tranhung93.github.io/ipw/reference/ipwtm.md),
[`tstartfun`](https://tranhung93.github.io/ipw/reference/tstartfun.md)

## Author

Willem M. van der Wal <willem@vanderwalresearch.com>, Ronald B. Geskus
<rgeskus@oucru.org>

## Examples

``` r
# For an example of how to use these longitudinal measurements with basdat, see:
# ?ipwtm
```
