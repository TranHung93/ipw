# HAART and Survival in HIV Patients

Survival data measured in 1200 HIV positive patients. Start of follow-up
is HIV seroconversion. Each row corresponds to a 100 day interval of
follow-up time, using the counting process notation. Patients can
initiate HAART therapy. CD4 count is a confounder for the effect of
HAART on mortality.

## Usage

``` r
data(haartdat)
```

## Format

A data frame with 1200 patients and multiple observations per patient
(counting process notation) on the following 8 variables:

- patient:

  Patient ID.

- tstart:

  Starting time for each interval of follow-up, measured in days since
  HIV seroconversion.

- fuptime:

  End time for each interval of follow-up, measured in days since HIV
  seroconversion.

- haartind:

  Indicator for the initiation of HAART therapy at the end of the
  interval (0 = HAART not initiated / 1 = HAART initiated).

- event:

  Indicator for death at the end of the interval (0 = alive / 1 = died).

- sex:

  Sex (0 = male / 1 = female).

- age:

  Age at the start of follow-up (years).

- cd4.sqrt:

  Square root of CD4 count, measured at `fuptime`, before `haartind`.

- endtime:

  The final observed time point for the individual.

- dropout:

  Indicator for dropout/censoring at the end of the interval (0 = no, 1
  = yes).

## Details

These data were simulated to demonstrate Inverse Probability Weighting
(IPW). To allow for models predicting the initiation of HAART at
`fuptime = 0`, the starting time for the first interval of each patient
is set to -100.

## References

Van der Wal W.M. & Geskus R.B. (2011). ipw: An R Package for Inverse
Probability Weighting. *Journal of Statistical Software*, **43**(13),
1-23. [doi:10.18637/jss.v043.i13](https://doi.org/10.18637/jss.v043.i13)
.

## See also

[`basdat`](https://tranhung93.github.io/ipw/reference/basdat.md),
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
# For a full example of how to use this data with ipwtm, see:
# ?ipwtm
```
