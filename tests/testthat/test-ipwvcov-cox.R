# ipwvcov() for a coxph() marginal structural model with ipwpoint() weights.

cox_data <- function(n = 400, seed = 7, round_times = FALSE) {
  set.seed(seed)
  d <- data.frame(l = stats::rnorm(n), v = stats::rbinom(n, 1, 0.5))
  d$a <- stats::rbinom(n, 1, stats::plogis(-0.2 + 0.8 * d$l + 0.4 * d$v))
  t_ev <- stats::rexp(n, 0.2 * exp(-0.5 * d$a + 0.6 * d$l + 0.3 * d$v))
  d$status <- as.numeric(t_ev < 4)
  d$time <- pmin(t_ev, 4)
  if (round_times) d$time <- ceiling(d$time * 4) / 4   # many tied event times
  d
}

# Reference: the stacked estimating equations written out directly, with
# brute-force risk sets (Breslow ties), the weights recomputed from the
# exposure-model coefficients, and the bread from central finite differences.
cox_numerical_vcov <- function(d, estimand, msm, den, num = NULL) {
  X <- stats::model.matrix(msm)
  Z <- stats::model.matrix(stats::formula(den), d)
  Zn <- if (!is.null(num)) stats::model.matrix(stats::formula(num), d)
  tm <- d$time
  st <- d$status
  a <- d$a
  kb <- ncol(X); ka <- ncol(Z); kg <- if (is.null(Zn)) 0 else ncol(Zn)
  at_risk <- outer(tm, tm, ">=")   # [j, k]: j at risk at time of k

  psi <- function(theta) {
    beta <- theta[seq_len(kb)]
    alpha <- theta[kb + seq_len(ka)]
    p <- stats::plogis(drop(Z %*% alpha))
    w <- switch(estimand,
      ATE = ifelse(a == 1, 1 / p, 1 / (1 - p)),
      ATT = ifelse(a == 1, 1, p / (1 - p)),
      ATC = ifelse(a == 1, (1 - p) / p, 1),
      ATO = ifelse(a == 1, 1 - p, p))
    out_num <- NULL
    if (kg > 0) {
      q <- stats::plogis(drop(Zn %*% theta[kb + ka + seq_len(kg)]))
      w <- w * ifelse(a == 1, q, 1 - q)
      out_num <- Zn * (a - q)
    }
    e <- exp(drop(X %*% beta))
    S0 <- colSums(at_risk * (w * e))
    Xbar <- crossprod(at_risk, X * (w * e)) / S0           # row k: Xbar(t_k)
    dk <- st * w / S0
    # w_i [ st_i (X_i - Xbar(t_i)) - e_i sum_k dk Y_i(t_k) (X_i - Xbar(t_k)) ]
    comp <- e * (X * drop(at_risk %*% dk) - at_risk %*% (Xbar * dk))
    U <- w * (st * (X - Xbar) - comp)
    cbind(U, Z * (a - p), out_num)
  }
  theta_hat <- c(stats::coef(msm), stats::coef(den), if (kg > 0) stats::coef(num))
  expect_lt(max(abs(colSums(psi(theta_hat)))) / nrow(d), 1e-6)
  k <- length(theta_hat)
  h <- 1e-6
  A <- matrix(0, k, k)
  for (j in seq_len(k)) {
    tp <- tm2 <- theta_hat
    tp[j] <- tp[j] + h
    tm2[j] <- tm2[j] - h
    A[, j] <- -(colSums(psi(tp)) - colSums(psi(tm2))) / (2 * h)
  }
  V <- solve(A) %*% crossprod(psi(theta_hat)) %*% t(solve(A))
  V[seq_len(kb), seq_len(kb), drop = FALSE]
}

test_that("Cox pieces: score residuals, information and the weights-known sandwich", {
  d <- cox_data()
  w <- ipwpoint(a, "binomial", "logit", denominator = ~ l + v, data = d)
  msm <- survival::coxph(survival::Surv(time, status) ~ a + v, data = d,
                         weights = w$ipw.weights)
  parts <- ipw:::.coxph_parts(msm, w$ipw.weights, stats::model.matrix(msm))
  # residuals(type = "score") is the unweighted-form residual: U_i = w_i * res_i
  expect_equal(unname(parts$res), unname(stats::residuals(msm, type = "score")),
               tolerance = 1e-10)
  expect_equal(unname(parts$A), unname(solve(msm$naive.var)), tolerance = 1e-8)
  # coxph's robust variance (weights treated as known) is the beta block of
  # the sandwich without the exposure model
  A_inv <- solve(parts$A)
  expect_equal(unname(A_inv %*% crossprod(parts$U) %*% A_inv), unname(msm$var),
               tolerance = 1e-8)
})

test_that("ipwvcov for coxph equals the numerical-Jacobian sandwich", {
  d <- cox_data()
  for (est in c("ATE", "ATT", "ATC", "ATO")) {
    w <- ipwpoint(a, "binomial", "logit", denominator = ~ l + v, data = d, estimand = est)
    msm <- survival::coxph(survival::Surv(time, status) ~ a, data = d,
                           weights = w$ipw.weights)
    ref <- cox_numerical_vcov(d, est, msm, w$den.mod)
    expect_equal(unname(ipwvcov(w, msm)), unname(ref), tolerance = 1e-5, label = est)
  }
})

test_that("coxph: stabilized weights, covariates in the MSM and tied times (Breslow)", {
  d <- cox_data()
  w <- ipwpoint(a, "binomial", "logit", numerator = ~ v, denominator = ~ l + v, data = d)
  msm <- survival::coxph(survival::Surv(time, status) ~ a + v, data = d,
                         weights = w$ipw.weights)
  V <- ipwvcov(w, msm)
  expect_equal(unname(V), unname(cox_numerical_vcov(d, "ATE", msm, w$den.mod, w$num.mod)),
               tolerance = 1e-5)
  expect_identical(dimnames(V), list(c("a", "v"), c("a", "v")))

  dt <- cox_data(round_times = TRUE)
  expect_true(anyDuplicated(dt$time[dt$status == 1]) > 0)
  wt <- ipwpoint(a, "binomial", "logit", denominator = ~ l + v, data = dt, estimand = "ATO")
  msm_t <- survival::coxph(survival::Surv(time, status) ~ a, data = dt,
                           weights = wt$ipw.weights, ties = "breslow")
  expect_equal(unname(ipwvcov(wt, msm_t)),
               unname(cox_numerical_vcov(dt, "ATO", msm_t, wt$den.mod)), tolerance = 1e-5)
})

test_that("without tied event times, Efron and Breslow give the same result", {
  d <- cox_data()
  expect_false(anyDuplicated(d$time[d$status == 1]) > 0)
  w <- ipwpoint(a, "binomial", "logit", denominator = ~ l + v, data = d)
  ve <- ipwvcov(w, survival::coxph(survival::Surv(time, status) ~ a, data = d,
                                   weights = w$ipw.weights, ties = "efron"))
  vb <- ipwvcov(w, survival::coxph(survival::Surv(time, status) ~ a, data = d,
                                   weights = w$ipw.weights, ties = "breslow"))
  expect_equal(ve, vb, tolerance = 1e-10)
})

test_that("ipwvcov gives clear errors for unsupported coxph models", {
  d <- cox_data(300)
  w <- ipwpoint(a, "binomial", "logit", denominator = ~ l + v, data = d)
  Surv <- survival::Surv
  fit <- function(f, ...) survival::coxph(f, data = d, weights = w$ipw.weights, ...)

  d$tstart <- 0
  expect_error(ipwvcov(w, fit(Surv(tstart, time, status) ~ a)), "counting-process data")
  expect_error(ipwvcov(w, fit(Surv(time, status) ~ a + strata(v))), "strata\\(\\)")
  d$id <- seq_len(nrow(d))
  expect_error(ipwvcov(w, fit(Surv(time, status) ~ a + cluster(id))), "cluster\\(\\)")
  expect_error(ipwvcov(w, survival::coxph(Surv(time, status) ~ a, data = d,
                                          weights = w$ipw.weights, cluster = id)),
               "cluster\\(\\)")
  expect_error(ipwvcov(w, survival::coxph(Surv(time, status) ~ a, data = d)),
               "must equal object\\$ipw.weights")
  expect_error(ipwvcov(w, survival::coxph(Surv(time, status) ~ a, data = d[-1, ],
                                          weights = w$ipw.weights[-1])),
               "fitted on 299 rows but there are 300 weights")
  d$grp <- factor(sample(1:20, nrow(d), replace = TRUE))
  expect_error(suppressWarnings(ipwvcov(w, fit(Surv(time, status) ~ a + frailty(grp)))),
               "penalized terms")

  dt <- cox_data(300, round_times = TRUE)
  wt <- ipwpoint(a, "binomial", "logit", denominator = ~ l + v, data = dt)
  expect_error(ipwvcov(wt, survival::coxph(Surv(time, status) ~ a, data = dt,
                                           weights = wt$ipw.weights)),
               "needs ties = \"breslow\" when event times are tied")
})

test_that("simulation: 95% intervals from ipwvcov for a Cox MSM cover about 95%", {
  skip_on_cran()
  set.seed(77)
  # true conditional log HR of a is -0.5; with no effect modification and a
  # covariate-free marginal model the marginal log HR differs (non-collapsible),
  # so the truth is computed from a large weighted fit of the same design
  gen <- function(n) {
    l <- stats::rnorm(n)
    a <- stats::rbinom(n, 1, stats::plogis(0.3 + 0.8 * l))
    t0 <- stats::rexp(n, 0.2 * exp(0.6 * l))      # untreated potential time
    t1 <- t0 * exp(0.5)                            # treated: hazard ratio exp(-0.5)
    tt <- ifelse(a == 1, t1, t0)
    data.frame(l, a, time = pmin(tt, 4), status = as.numeric(tt < 4), t0 = t0, t1 = t1)
  }
  big <- gen(2e5)
  truth <- unname(stats::coef(survival::coxph(
    survival::Surv(c(pmin(big$t0, 4), pmin(big$t1, 4)),
                   c(big$t0 < 4, big$t1 < 4)) ~ rep(0:1, each = nrow(big)))))
  covered <- replicate(300, {
    d <- gen(500)
    w <- ipwpoint(a, "binomial", "logit", denominator = ~ l, data = d)
    msm <- survival::coxph(survival::Surv(time, status) ~ a, data = d,
                           weights = w$ipw.weights)
    abs(stats::coef(msm) - truth) < stats::qnorm(0.975) * sqrt(ipwvcov(w, msm)[1, 1])
  })
  expect_gt(mean(covered), 0.91)
  expect_lt(mean(covered), 0.985)
})
