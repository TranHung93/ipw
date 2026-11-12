# ipwvcov(): M-estimation variance of the MSM coefficients, taking the
# estimation of the weights into account.

vc_data <- function(n = 800, seed = 5) {
  set.seed(seed)
  d <- data.frame(l = stats::rnorm(n), v = stats::rbinom(n, 1, 0.5))
  d$a <- stats::rbinom(n, 1, stats::plogis(-0.2 + 0.8 * d$l + 0.4 * d$v))
  d$y <- 2 * d$a + d$l + 0.5 * d$v + stats::rnorm(n)
  d$yb <- stats::rbinom(n, 1, stats::plogis(-0.5 + 0.7 * d$a + 0.5 * d$l))
  d
}

# Reference: the stacked estimating equations written out directly, with the
# derivative matrix A = -d sum(psi) / d theta obtained by central finite
# differences instead of formulas.
numerical_vcov <- function(d, outcome, estimand, link, msm_formula, msm_family,
                           den_formula, num_formula = NULL, msm, den, num = NULL) {
  X <- stats::model.matrix(msm_formula, d)
  Z <- stats::model.matrix(den_formula, d)
  Zn <- if (!is.null(num_formula)) stats::model.matrix(num_formula, d)
  y <- d[[outcome]]
  a <- d$a
  lk <- stats::binomial(link)
  kb <- ncol(X); ka <- ncol(Z); kg <- if (is.null(Zn)) 0 else ncol(Zn)

  psi <- function(theta) {
    beta <- theta[seq_len(kb)]
    alpha <- theta[kb + seq_len(ka)]
    p <- lk$linkinv(drop(Z %*% alpha))
    dp <- lk$mu.eta(drop(Z %*% alpha))
    w <- switch(estimand,
      ATE = ifelse(a == 1, 1 / p, 1 / (1 - p)),
      ATT = ifelse(a == 1, 1, p / (1 - p)),
      ATC = ifelse(a == 1, (1 - p) / p, 1),
      ATO = ifelse(a == 1, 1 - p, p))
    out_num <- NULL
    if (kg > 0) {
      gamma <- theta[kb + ka + seq_len(kg)]
      q <- stats::plogis(drop(Zn %*% gamma))
      w <- w * ifelse(a == 1, q, 1 - q)
      out_num <- Zn * (a - q)
    }
    eta <- drop(X %*% beta)
    mu <- if (msm_family == "gaussian") eta else stats::plogis(eta)
    cbind(X * (w * (y - mu)), Z * ((a - p) * dp / (p * (1 - p))), out_num)
  }

  theta_hat <- c(stats::coef(msm), stats::coef(den), if (kg > 0) stats::coef(num))
  # the estimates solve the estimating equations, so the reference reproduces
  # the package's weights and models
  expect_lt(max(abs(colSums(psi(theta_hat)))) / nrow(d), 1e-5)
  k <- length(theta_hat)
  A <- matrix(0, k, k)
  h <- 1e-6
  for (j in seq_len(k)) {
    tp <- tm <- theta_hat
    tp[j] <- tp[j] + h
    tm[j] <- tm[j] - h
    A[, j] <- -(colSums(psi(tp)) - colSums(psi(tm))) / (2 * h)
  }
  Psi <- psi(theta_hat)
  V <- solve(A) %*% crossprod(Psi) %*% t(solve(A))
  V[seq_len(kb), seq_len(kb)]
}

test_that("ipwvcov equals the sandwich of the stacked equations with a numerical Jacobian", {
  d <- vc_data()
  for (est in c("ATE", "ATT", "ATC", "ATO")) {
    w <- ipwpoint(a, "binomial", "logit", denominator = ~ l + v, data = d, estimand = est)
    msm <- stats::lm(y ~ a + l, data = d, weights = w$ipw.weights)
    ref <- numerical_vcov(d, "y", est, "logit", ~ a + l, "gaussian", ~ l + v,
                          msm = msm, den = w$den.mod)
    expect_equal(unname(ipwvcov(w, msm)), unname(ref), tolerance = 1e-5, label = est)

    msm_b <- suppressWarnings(stats::glm(yb ~ a, data = d, weights = w$ipw.weights,
                                         family = stats::binomial))
    ref_b <- numerical_vcov(d, "yb", est, "logit", ~ a, "binomial", ~ l + v,
                            msm = msm_b, den = w$den.mod)
    expect_equal(unname(ipwvcov(w, msm_b)), unname(ref_b), tolerance = 1e-5,
                 label = paste(est, "binomial MSM"))
  }
})

test_that("stabilized weights with a numerator model: numerical Jacobian agrees", {
  d <- vc_data()
  w <- ipwpoint(a, "binomial", "logit", numerator = ~ v, denominator = ~ l + v, data = d)
  msm <- stats::lm(y ~ a + v, data = d, weights = w$ipw.weights)
  ref <- numerical_vcov(d, "y", "ATE", "logit", ~ a + v, "gaussian", ~ l + v, ~ v,
                        msm = msm, den = w$den.mod, num = w$num.mod)
  V <- ipwvcov(w, msm)
  expect_equal(unname(V), unname(ref), tolerance = 1e-5)
  expect_identical(dimnames(V), list(c("(Intercept)", "a", "v"), c("(Intercept)", "a", "v")))
})

test_that("for a saturated MSM, stabilized and unstabilized ATE weights give the same SE", {
  d <- vc_data()
  ws <- ipwpoint(a, "binomial", "logit", numerator = ~ 1, denominator = ~ l + v, data = d)
  wu <- ipwpoint(a, "binomial", "logit", denominator = ~ l + v, data = d)
  vs <- ipwvcov(ws, stats::lm(y ~ a, data = d, weights = ws$ipw.weights))
  vu <- ipwvcov(wu, stats::lm(y ~ a, data = d, weights = wu$ipw.weights))
  expect_equal(vs, vu, tolerance = 1e-8)
})

test_that("non-logit exposure models are refused (expected information is not exact)", {
  d <- vc_data()
  for (lk in c("probit", "cloglog")) {
    w <- ipwpoint(a, "binomial", lk, denominator = ~ l + v, data = d)
    msm <- stats::lm(y ~ a, data = d, weights = w$ipw.weights)
    expect_error(ipwvcov(w, msm), "link = 'logit'", label = lk)
  }
  wn <- ipwpoint(a, "binomial", "probit", numerator = ~ v, denominator = ~ l + v, data = d)
  expect_error(ipwvcov(wn, stats::lm(y ~ a, data = d, weights = wn$ipw.weights)),
               "link = 'logit'")
})

test_that("ipwvcov agrees with WeightIt::glm_weightit(vcov = 'asympt')", {
  skip_if_not_installed("WeightIt")
  d <- vc_data()
  for (est in c("ATE", "ATT", "ATC", "ATO")) {
    w <- ipwpoint(a, "binomial", "logit", denominator = ~ l + v, data = d, estimand = est)
    W <- WeightIt::weightit(a ~ l + v, data = d, method = "glm", estimand = est)
    fit_g <- WeightIt::glm_weightit(y ~ a, data = d, weightit = W, vcov = "asympt")
    msm_g <- stats::lm(y ~ a, data = d, weights = w$ipw.weights)
    expect_equal(unname(ipwvcov(w, msm_g)), unname(stats::vcov(fit_g)), tolerance = 1e-6,
                 label = paste(est, "gaussian"))

    fit_b <- WeightIt::glm_weightit(yb ~ a, data = d, weightit = W, family = stats::binomial,
                                    vcov = "asympt")
    msm_b <- stats::glm(yb ~ a, data = d, weights = w$ipw.weights, family = stats::quasibinomial)
    expect_equal(unname(ipwvcov(w, msm_b)), unname(stats::vcov(fit_b)), tolerance = 1e-6,
                 label = paste(est, "binomial"))
  }
})

test_that("ipwvcov gives clear errors for unsupported input", {
  d <- vc_data(300)
  d$m <- factor(sample(c("x", "y", "z"), 300, replace = TRUE))
  w <- ipwpoint(a, "binomial", "logit", denominator = ~ l, data = d)
  msm <- stats::lm(y ~ a, data = d, weights = w$ipw.weights)

  wm <- ipwpoint(m, "multinomial", denominator = ~ l, data = d)
  expect_error(ipwvcov(wm, stats::lm(y ~ m, data = d, weights = wm$ipw.weights)),
               "only available for ipwpoint\\(\\) with family = 'binomial'")

  wt <- ipwpoint(a, "binomial", "logit", denominator = ~ l, data = d, trunc = 0.05)
  expect_error(ipwvcov(wt, stats::lm(y ~ a, data = d, weights = wt$weights.trunc)),
               "Truncated weights are not supported")
  # the untruncated weights of the same object are fine
  expect_no_error(ipwvcov(wt, stats::lm(y ~ a, data = d, weights = wt$ipw.weights)))

  expect_error(ipwvcov(w, stats::lm(y ~ a, data = d[-1, ], weights = w$ipw.weights[-1])),
               "fitted on 299 rows but there are 300 weights")
  expect_error(ipwvcov(w, stats::lm(y ~ a, data = d, weights = rev(w$ipw.weights))),
               "must equal object\\$ipw.weights")
  expect_error(ipwvcov(w, stats::lm(y ~ a, data = d)), "must equal object\\$ipw.weights")
  expect_error(ipwvcov(w, "not a model"), "fitted with lm\\(\\) or glm\\(\\)")
  expect_error(ipwvcov(list(ipw.weights = 1), msm), "must be the result of ipwpoint")

  long <- data.frame(id = rep(1:50, each = 3), time = rep(1:3, 50), l = stats::rnorm(150))
  long$x <- stats::rbinom(150, 1, stats::plogis(long$l))
  wtm <- ipwtm(exposure = x, family = "binomial", link = "logit", denominator = ~ l,
               id = id, timevar = time, type = "all", data = long)
  expect_error(ipwvcov(wtm, msm), "does not support time-varying weights from ipwtm")

  skip_if_not_installed("mgcv")
  wg <- ipwpoint(a, "binomial", "logit", denominator = ~ s(l), data = d)
  expect_error(ipwvcov(wg, stats::lm(y ~ a, data = d, weights = wg$ipw.weights)),
               "smooth terms \\(GAM\\)")
})

test_that("simulation: 95% intervals from ipwvcov cover the ATE about 95% of the time", {
  skip_on_cran()
  set.seed(2024)
  covered <- replicate(400, {
    n <- 500
    l <- stats::rnorm(n)
    a <- stats::rbinom(n, 1, stats::plogis(0.3 + 0.8 * l))
    y <- 2 * a + l + stats::rnorm(n)
    d <- data.frame(l, a, y)
    w <- ipwpoint(a, "binomial", "logit", denominator = ~ l, data = d)
    msm <- stats::lm(y ~ a, data = d, weights = w$ipw.weights)
    se <- sqrt(ipwvcov(w, msm)["a", "a"])
    abs(stats::coef(msm)["a"] - 2) < stats::qnorm(0.975) * se
  })
  expect_gt(mean(covered), 0.92)
  expect_lt(mean(covered), 0.98)
})
