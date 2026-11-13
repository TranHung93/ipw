# Performance of a standard error over simulated data sets, judged on the
# variance scale (Morris, White & Crowther 2019, Statistics in Medicine 38,
# 2074-2102, Table 6). Sourced by the simulation scripts in data-raw/ and by
# dev/simulation-study.R.
#
#   est      estimates from all data sets (gives the empirical variance)
#   se       standard errors, one per data set; NA for data sets without one
#            (e.g. the bootstrap, run on a subset)
#   truth    true value, for coverage of the 95% Wald interval est +- z * se
#   lo, hi   interval limits, if not Wald (e.g. percentile bootstrap)
#
# Returns one row:
#   mean.var  mean estimated variance, mean(se^2)
#   emp.var   empirical variance of the estimates, var(est)
#   var.ratio mean.var / emp.var (1 is ideal)
#   model.se  sqrt(mean.var), the model-based SE
#   emp.se    sqrt(emp.var), the empirical SE
#   rel.error relative % error in the model SE, 100 * (model.se / emp.se - 1)
#   mcse.rel.error  its Monte Carlo SE (Morris et al., Table 6):
#             100 * (model.se / emp.se) *
#               sqrt(var(se^2) / (4 * R_se * model.se^4) + 1 / (2 * (R - 1)))
#             where R_se is the number of data sets with a SE and R the
#             number used for emp.var. When R_se < R the first term uses
#             R_se and the second R; the two parts are treated as
#             independent, as in the formula.
#   coverage, mcse.coverage = sqrt(p (1 - p) / R_se), replications (R_se)
se_performance <- function(est, se, truth, lo = NULL, hi = NULL, z = stats::qnorm(0.975)) {
  ok_est <- is.finite(est)
  emp_var <- stats::var(est[ok_est])
  R <- sum(ok_est)
  ok <- ok_est & is.finite(se)
  if (!is.null(lo)) ok <- ok & is.finite(lo) & is.finite(hi)
  s2 <- se[ok]^2
  R_se <- sum(ok)
  if (is.null(lo)) {
    lo <- est - z * se
    hi <- est + z * se
  }
  cover <- mean(lo[ok] <= truth & hi[ok] >= truth)
  mean_var <- mean(s2)
  model_se <- sqrt(mean_var)
  emp_se <- sqrt(emp_var)
  data.frame(
    mean.var = mean_var, emp.var = emp_var, var.ratio = mean_var / emp_var,
    model.se = model_se, emp.se = emp_se,
    rel.error = 100 * (model_se / emp_se - 1),
    mcse.rel.error = 100 * (model_se / emp_se) *
      sqrt(stats::var(s2) / (4 * R_se * model_se^4) + 1 / (2 * (R - 1))),
    coverage = cover, mcse.coverage = sqrt(cover * (1 - cover) / R_se),
    replications = R_se)
}
