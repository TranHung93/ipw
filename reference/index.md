# Package index

## Estimate weights

Fit the exposure (or censoring) models and compute inverse probability
weights.

- [`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md)
  : Estimate Inverse Probability Weights (Point Treatment)
- [`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) :
  Estimate Inverse Probability Weights (Time Varying)
- [`ipwprob()`](https://tranhung93.github.io/ipw/reference/ipwprob.md) :
  Inverse Probability Weights From User-Supplied Probabilities

## Standard errors

Bootstrap a weighted analysis, re-estimating the weights in each sample.

- [`ipwboot()`](https://tranhung93.github.io/ipw/reference/ipwboot.md)
  [`print(`*`<ipwboot>`*`)`](https://tranhung93.github.io/ipw/reference/ipwboot.md)
  [`confint(`*`<ipwboot>`*`)`](https://tranhung93.github.io/ipw/reference/ipwboot.md)
  : Bootstrap Standard Errors for Inverse Probability Weighted Analyses

## Diagnostics

- [`ipwplot()`](https://tranhung93.github.io/ipw/reference/ipwplot.md) :
  Plot Inverse Probability Weights

## Data preparation

- [`tstartfun()`](https://tranhung93.github.io/ipw/reference/tstartfun.md)
  : Compute Starting Time For Counting Process Notation

## Example data

- [`haartdat`](https://tranhung93.github.io/ipw/reference/haartdat.md) :
  HAART and Survival in HIV Patients
- [`basdat`](https://tranhung93.github.io/ipw/reference/basdat.md) :
  HIV: TB and Survival (Baseline Data)
- [`timedat`](https://tranhung93.github.io/ipw/reference/timedat.md) :
  HIV: TB and Survival (Longitudinal Measurements)
- [`healthdat`](https://tranhung93.github.io/ipw/reference/healthdat.md)
  : IQ, Income and Health
