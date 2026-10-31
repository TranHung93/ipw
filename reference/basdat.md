# HIV: TB and Survival (Baseline Data)

Simulated dataset. Baseline data of 386 HIV positive individuals,
including time of first active tuberculosis, time of death, individual
end time. Time varying CD4 measurements of these patients are included
in dataset
[`timedat`](https://tranhung93.github.io/ipw/reference/timedat.md).

## Usage

``` r
data(basdat)
```

## Format

A data frame with 386 observations on the following 4 variables:

- id:

  Patient ID.

- Ttb:

  Time of first active tuberculosis, measured in days since HIV
  seroconversion.

- Tdeath:

  Time of death, measured in days since HIV seroconversion.

- Tend:

  Individual end time (either death or censoring), measured in days
  since HIV seroconversion.

## Details

These simulated data are used together with data in
[`timedat`](https://tranhung93.github.io/ipw/reference/timedat.md) in a
detailed causal modelling example using inverse probability weighting
(IPW). See
[`ipwtm`](https://tranhung93.github.io/ipw/reference/ipwtm.md) for the
example. Data were simulated using the algorithm described in Van der
Wal e.a. (2009).#' \[Image of a multistate model for HIV, TB, and
death\]

## References

Van der Wal W.M. & Geskus R.B. (2011). ipw: An R Package for Inverse
Probability Weighting. *Journal of Statistical Software*, **43**(13),
1-23. [doi:10.18637/jss.v043.i13](https://doi.org/10.18637/jss.v043.i13)
.

Van der Wal W.M., Prins M., Lumbreras B. & Geskus R.B. (2009). A simple
G-computation algorithm to quantify the causal effect of a secondary
illness on the progression of a chronic disease. *Statistics in
Medicine*, **28**(18), 2325-2337.

## See also

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
# Detailed examples can be found in the ipwtm documentation:
# ?ipwtm
```
