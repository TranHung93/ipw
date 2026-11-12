# ATT (with a focal level) and generalized overlap weights (ATO) for
# multinomial and ordinal exposures in ipwpoint().

cat_data <- function(n = 600, seed = 31) {
  set.seed(seed)
  d <- data.frame(l = stats::rnorm(n), z = stats::rnorm(n))
  eta2 <- 0.5 + 0.8 * d$l
  eta3 <- -0.3 - 0.6 * d$l + 0.4 * d$z
  pr <- cbind(1, exp(eta2), exp(eta3))
  pr <- pr / rowSums(pr)
  d$m <- factor(apply(pr, 1, function(p) sample(c("A", "B", "C"), 1, prob = p)))
  d$o <- cut(0.9 * d$l + stats::rlogis(n), c(-Inf, -0.5, 0.8, Inf),
             labels = c("low", "mid", "high"), ordered_result = TRUE)
  d$a <- stats::rbinom(n, 1, stats::plogis(-0.2 + 0.7 * d$l))
  d$a2 <- factor(ifelse(d$a == 1, "treated", "control"), levels = c("control", "treated"))
  d
}

test_that("multinomial ATT and ATO weights equal by-hand nnet::multinom probabilities", {
  d <- cat_data()
  fit <- nnet::multinom(m ~ l + z, data = d, trace = FALSE)
  P <- stats::predict(fit, type = "probs")
  p_obs <- P[cbind(seq_len(nrow(d)), match(as.character(d$m), colnames(P)))]

  for (foc in levels(d$m)) {
    res <- ipwpoint(m, "multinomial", denominator = ~ l + z, data = d,
                    estimand = "ATT", focal = foc)
    expected <- ifelse(d$m == foc, 1, P[, foc] / p_obs)
    expect_equal(res$ipw.weights, unname(expected), tolerance = 1e-10, label = foc)
    expect_true(all(res$ipw.weights[d$m == foc] == 1))
    expect_identical(attr(res, "focal"), foc)
    expect_identical(attr(res, "estimand"), "ATT")
  }

  res <- ipwpoint(m, "multinomial", denominator = ~ l + z, data = d, estimand = "ATO")
  expect_equal(res$ipw.weights, unname((1 / p_obs) / rowSums(1 / P)), tolerance = 1e-10)
  expect_identical(attr(res, "estimand"), "ATO")
  expect_null(attr(res, "focal"))
  expect_identical(attr(res, "exposure"), "m")
})

test_that("ordinal ATT and ATO weights equal by-hand MASS::polr probabilities", {
  d <- cat_data()
  for (lk in c("logit", "probit")) {
    fit <- MASS::polr(o ~ l + z, data = d, method = ipw:::.polr_method(lk), Hess = TRUE)
    P <- stats::predict(fit, type = "probs")
    p_obs <- P[cbind(seq_len(nrow(d)), match(as.character(d$o), colnames(P)))]

    res <- ipwpoint(o, "ordinal", lk, denominator = ~ l + z, data = d,
                    estimand = "ATT", focal = "mid")
    expect_equal(res$ipw.weights, unname(ifelse(d$o == "mid", 1, P[, "mid"] / p_obs)),
                 tolerance = 1e-10, label = lk)

    res <- ipwpoint(o, "ordinal", lk, denominator = ~ l + z, data = d, estimand = "ATO")
    expect_equal(res$ipw.weights, unname((1 / p_obs) / rowSums(1 / P)),
                 tolerance = 1e-10, label = lk)
  }
})

test_that("numeric codes can be used as focal level", {
  d <- cat_data()
  d$k <- as.integer(d$o)
  a <- ipwpoint(k, "ordinal", "logit", denominator = ~ l, data = d, estimand = "ATT", focal = 2)
  b <- ipwpoint(o, "ordinal", "logit", denominator = ~ l, data = d, estimand = "ATT", focal = "mid")
  expect_equal(a$ipw.weights, b$ipw.weights)
  expect_identical(attr(a, "focal"), "2")
})

test_that("two-level multinomial ATT and ATO equal the binomial weights", {
  d <- cat_data()
  bin_att <- ipwpoint(a, "binomial", "logit", denominator = ~ l, data = d, estimand = "ATT")
  bin_ato <- ipwpoint(a, "binomial", "logit", denominator = ~ l, data = d, estimand = "ATO")
  mult_att <- ipwpoint(a2, "multinomial", denominator = ~ l, data = d,
                       estimand = "ATT", focal = "treated")
  mult_ato <- ipwpoint(a2, "multinomial", denominator = ~ l, data = d, estimand = "ATO")
  # multinom() and glm() converge to slightly different end points
  expect_equal(mult_att$ipw.weights, bin_att$ipw.weights, tolerance = 1e-5)
  expect_equal(mult_ato$ipw.weights, bin_ato$ipw.weights, tolerance = 1e-5)

  # ATT with the control level as focal equals the binomial ATC
  bin_atc <- ipwpoint(a, "binomial", "logit", denominator = ~ l, data = d, estimand = "ATC")
  mult_ctl <- ipwpoint(a2, "multinomial", denominator = ~ l, data = d,
                       estimand = "ATT", focal = "control")
  expect_equal(mult_ctl$ipw.weights, bin_atc$ipw.weights, tolerance = 1e-5)
})

test_that("binomial estimands record their focal group as attributes", {
  d <- cat_data()
  expect_identical(attr(ipwpoint(a, "binomial", "logit", denominator = ~ l, data = d,
                                 estimand = "ATT"), "focal"), 1)
  expect_identical(attr(ipwpoint(a, "binomial", "logit", denominator = ~ l, data = d,
                                 estimand = "ATC"), "focal"), 0)
  ate <- ipwpoint(a, "binomial", "logit", denominator = ~ l, data = d)
  expect_null(attr(ate, "focal"))
  expect_identical(attr(ate, "estimand"), "ATE")
  expect_identical(attr(ate, "exposure"), "a")
})

test_that("categorical ATT and ATO weights agree with WeightIt", {
  skip_if_not_installed("WeightIt")
  d <- cat_data(1000, 32)
  for (foc in levels(d$m)) {
    ipw_w <- ipwpoint(m, "multinomial", denominator = ~ l + z, data = d,
                      estimand = "ATT", focal = foc)$ipw.weights
    wit_w <- WeightIt::weightit(m ~ l + z, data = d, method = "glm",
                                estimand = "ATT", focal = foc)$weights
    # With its default reltol, nnet::multinom stops earlier than WeightIt's
    # fitter (relative differences up to a few times 1e-4); with a tight
    # reltol the two agree to 1e-6 or better.
    expect_equal(ipw_w, unname(wit_w), tolerance = 1e-3, label = foc)
    ipw_tight <- ipwpoint(m, "multinomial", denominator = ~ l + z, data = d,
                          estimand = "ATT", focal = foc, reltol = 1e-12, maxit = 1000)$ipw.weights
    expect_equal(ipw_tight, unname(wit_w), tolerance = 1e-6, label = foc)
  }
  ipw_w <- ipwpoint(m, "multinomial", denominator = ~ l + z, data = d, estimand = "ATO",
                    reltol = 1e-12, maxit = 1000)$ipw.weights
  wit_w <- WeightIt::weightit(m ~ l + z, data = d, method = "glm", estimand = "ATO")$weights
  expect_equal(ipw_w, unname(wit_w), tolerance = 1e-6)

  # ordered factor: ordinal logistic model (MASS::polr in ipw)
  for (est in c("ATT", "ATO")) {
    foc <- if (est == "ATT") "mid" else NULL
    ipw_o <- ipwpoint(o, "ordinal", "logit", denominator = ~ l + z, data = d,
                      estimand = est, focal = foc)$ipw.weights
    wit_o <- WeightIt::weightit(o ~ l + z, data = d, method = "glm",
                                estimand = est, focal = foc)$weights
    expect_equal(ipw_o, unname(wit_o), tolerance = 1e-4, label = paste("ordinal", est))
  }
})

test_that("estimand and focal errors are clear", {
  d <- cat_data(200)
  d$cnt <- stats::rpois(200, 2)
  expect_error(ipwpoint(m, "multinomial", denominator = ~ l, data = d, estimand = "ATT"),
               "needs 'focal'")
  expect_error(ipwpoint(m, "multinomial", denominator = ~ l, data = d, estimand = "ATT",
                        focal = "D"), "must be one level of the exposure: 'A', 'B', 'C'")
  expect_error(ipwpoint(m, "multinomial", denominator = ~ l, data = d, estimand = "ATT",
                        focal = c("A", "B")), "must be one level")
  expect_error(ipwpoint(m, "multinomial", denominator = ~ l, data = d, estimand = "ATC"),
               "only available for family = 'binomial'.*estimand = 'ATT' with 'focal' set to the control level")
  expect_error(ipwpoint(o, "ordinal", "logit", denominator = ~ l, data = d, estimand = "ATC"),
               "focal' set to the control level")
  expect_error(ipwpoint(m, "multinomial", numerator = ~ 1, denominator = ~ l, data = d,
                        estimand = "ATO"), "must be NULL")
  expect_error(ipwpoint(m, "multinomial", denominator = ~ l, data = d, estimand = "ATO",
                        focal = "A"), "'focal' is only used")
  expect_error(ipwpoint(m, "multinomial", denominator = ~ l, data = d, focal = "A"),
               "'focal' is only used")
  expect_error(ipwpoint(a, "binomial", "logit", denominator = ~ l, data = d,
                        estimand = "ATT", focal = 1), "focal group is always the exposed")
  expect_error(ipwpoint(cnt, "poisson", denominator = ~ l, data = d, estimand = "ATO"),
               "only available for family = 'binomial', 'multinomial', 'ordinal'")
})
