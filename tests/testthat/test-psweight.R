# Second independent check of the estimands and of ipwvcov(), against the
# PSweight package (used only through its exported functions).
#
# PSweight normalizes the weights within each exposure group (they sum to 1
# in each group), so weights are compared after the same normalization. A
# weighted mean within each group, and so the effect estimate, does not
# depend on that scaling.

skip_if_not_installed("PSweight")

norm_by_group <- function(w, g) w / stats::ave(w, g, FUN = sum)

psw_binary <- function(n = 800, seed = 5) {
  set.seed(seed)
  d <- data.frame(l = stats::rnorm(n), v = stats::rbinom(n, 1, 0.5))
  d$a <- stats::rbinom(n, 1, stats::plogis(-0.2 + 0.8 * d$l + 0.4 * d$v))
  d$y <- 2 * d$a + d$l + 0.5 * d$v + stats::rnorm(n)
  d
}

psw_multi <- function(n = 1000, seed = 31) {
  set.seed(seed)
  d <- data.frame(l = stats::rnorm(n), z = stats::rnorm(n))
  pr <- cbind(1, exp(0.5 + 0.8 * d$l), exp(-0.3 - 0.6 * d$l + 0.4 * d$z))
  pr <- pr / rowSums(pr)
  d$m <- factor(apply(pr, 1, function(p) sample(c("A", "B", "C"), 1, prob = p)))
  d
}

# ipw estimand -> PSweight weight name. PSweight has no "ATC" option, but
# its ATT weights with the untreated as "treated" group (trtgrp = "0") are
# the ATC weights.
psw_names <- c(ATE = "IPW", ATT = "treated", ATO = "overlap")
psw_binary_args <- list(ATE = list(weight = "IPW"), ATT = list(weight = "treated"),
                        ATC = list(weight = "treated", trtgrp = "0"),
                        ATO = list(weight = "overlap"))

test_that("binary weights agree with PSweight::SumStat (ATE, ATT, ATC, ATO)", {
  d <- psw_binary()
  expect_identical(PSweight::SumStat(ps.formula = a ~ l + v, data = d,
                                     weight = "treated")$trtgrp, "1")   # as in ipwpoint()
  for (est in names(psw_binary_args)) {
    ss <- do.call(PSweight::SumStat,
                  c(list(ps.formula = a ~ l + v, data = d), psw_binary_args[[est]]))
    w_ipw <- ipwpoint(a, "binomial", "logit", denominator = ~ l + v, data = d,
                      estimand = est)$ipw.weights
    w_psw <- ss$ps.weights[[psw_binary_args[[est]]$weight]]
    expect_equal(norm_by_group(w_ipw, d$a), norm_by_group(w_psw, d$a),
                 tolerance = 1e-10, label = est)
  }
})

test_that("three-level multinomial weights agree with PSweight::SumStat", {
  # PSweight fits the multinomial model with the same defaults as
  # nnet::multinom(), which ipwpoint() uses, so the fitted probabilities and
  # the weights agree to rounding error. (With a tighter reltol passed to
  # ipwpoint() they would differ by about 1e-3, the convergence tolerance of
  # the default fit, as for WeightIt in test-categorical-estimands.R.)
  d <- psw_multi()
  cases <- list(list("ATE", NULL), list("ATT", "A"), list("ATT", "C"), list("ATO", NULL))
  for (cs in cases) {
    est <- cs[[1]]
    foc <- cs[[2]]
    ss <- PSweight::SumStat(ps.formula = m ~ l + z, data = d, weight = psw_names[[est]],
                            trtgrp = foc)
    w_ipw <- ipwpoint(m, "multinomial", denominator = ~ l + z, data = d,
                      estimand = est, focal = foc)$ipw.weights
    w_psw <- ss$ps.weights[[psw_names[[est]]]]
    expect_equal(norm_by_group(w_ipw, d$m), norm_by_group(w_psw, d$m),
                 tolerance = 1e-8, label = paste(est, foc))
  }
})

test_that("ipwvcov SE equals PSweight's sandwich SE of the difference in means", {
  # Saturated MSM lm(y ~ a): its coefficient of a is the weighted difference
  # in means that PSweight estimates, and both standard errors are
  # M-estimation sandwiches that account for the logistic propensity model.
  for (seed in c(5, 6)) {
    d <- psw_binary(seed = seed)
    for (est in names(psw_binary_args)) {
      w <- ipwpoint(a, "binomial", "logit", denominator = ~ l + v, data = d, estimand = est)
      msm <- stats::lm(y ~ a, data = d, weights = w$ipw.weights)
      fit <- do.call(PSweight::PSweight,
                     c(list(ps.formula = a ~ l + v, yname = "y", data = d),
                       psw_binary_args[[est]]))
      res <- summary(fit)$estimates
      expect_equal(unname(stats::coef(msm)["a"]), unname(res[1, "Estimate"]),
                   tolerance = 1e-10, label = paste(seed, est, "estimate"))
      expect_equal(sqrt(ipwvcov(w, msm)["a", "a"]), unname(res[1, "Std.Error"]),
                   tolerance = 1e-8, label = paste(seed, est, "SE"))
    }
  }
})
