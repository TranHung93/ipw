# ipwvcov() for time-varying weights from ipwtm() (binomial, logistic models).

# Long data with individual visit times: a confounder l that changes over
# time, a treatment start A (type "first"), dropout C (type "cens"), an
# outcome event, a treatment B that switches back and forth (type "all") and
# a continuous outcome y measured at every visit (effect of B: 1).
tm_data <- function(n = 150, seed = 3) {
  set.seed(seed)
  rows <- lapply(seq_len(n), function(i) {
    k <- sample(3:6, 1)
    out <- vector("list", k)
    t <- 0; l <- stats::rnorm(1); a <- 0; b <- 0; x <- stats::rbinom(1, 1, 0.5)
    for (j in seq_len(k)) {
      t_new <- t + stats::runif(1, 0.5, 1.5)
      l <- 0.7 * l + stats::rnorm(1, 0, 0.7)
      b_lag <- b
      if (a == 0) a <- stats::rbinom(1, 1, stats::plogis(-1.3 - 0.8 * l + 0.3 * x))
      b <- stats::rbinom(1, 1, stats::plogis(-0.5 + 0.8 * l + 1.0 * b_lag))
      ev <- stats::rbinom(1, 1, stats::plogis(-2.5 - 0.5 * a + 0.5 * l))
      cens <- if (ev == 1) 0 else stats::rbinom(1, 1, stats::plogis(-2.2 + 0.4 * l))
      out[[j]] <- data.frame(id = i, tstart = t, time = t_new, x = x, l = l, A = a,
                             C = cens, event = ev, B = b, Blag = b_lag,
                             y = b + l + stats::rnorm(1))
      t <- t_new
      if (ev == 1 || cens == 1) break
    }
    do.call(rbind, out[!vapply(out, is.null, logical(1))])
  })
  do.call(rbind, rows)
}

# Reference: the stacked estimating equations written out directly, with the
# weights recomputed from the model coefficients following ipwtm()'s
# definition, brute-force risk sets for a Cox MSM, scores summed within
# person, and the bread from central finite differences.
#   specs: list of list(a, type, Zd, Zn) on the rows of d (sorted by id, time)
tm_numerical_vcov <- function(d, specs, X, msm_kind, theta_hat, y = NULL) {
  N <- nrow(d)
  id <- d$id
  sel_of <- function(a, type) {
    if (type == "all") return(rep(TRUE, N))
    unlist(lapply(split(a, id), function(x) {
      j <- match(1, x)
      if (is.na(j)) rep(TRUE, length(x)) else seq_along(x) <= j
    }), use.names = FALSE)
  }
  kb <- ncol(X)
  if (msm_kind == "cox") {
    at_risk <- outer(seq_len(N), seq_len(N), function(j, k)
      d$tstart[j] < d$time[k] & d$time[k] <= d$time[j])
  }
  psi <- function(theta) {
    pos <- kb
    log_w <- rep(0, N)
    scores <- list()
    for (sp in specs) {
      sel <- sel_of(sp$a, sp$type)
      for (part in c("den", "num")) {
        Zp <- if (part == "den") sp$Zd else sp$Zn
        if (is.null(Zp)) next
        coef_p <- theta[pos + seq_len(ncol(Zp))]
        pos <- pos + ncol(Zp)
        p <- stats::plogis(drop(Zp %*% coef_p))
        f <- if (sp$type == "cens") 1 - p else ifelse(sp$a == 1, p, 1 - p)
        f[!sel] <- 1
        lf <- unlist(lapply(split(log(f), id), cumsum), use.names = FALSE)
        log_w <- log_w + (if (part == "den") -lf else lf)
        scores <- c(scores, list(Zp * ((sp$a - p) * sel)))
      }
    }
    w <- exp(log_w)
    beta <- theta[seq_len(kb)]
    eta <- drop(X %*% beta)
    if (msm_kind == "gaussian") {
      U <- X * (w * (y - eta))
    } else {
      e <- exp(eta)
      S0 <- colSums(at_risk * (w * e))
      Xbar <- crossprod(at_risk, X * (w * e)) / S0
      dk <- d$event * w / S0
      dk[d$event == 0] <- 0
      Xbar[d$event == 0, ] <- 0
      comp <- e * (X * drop(at_risk %*% dk) - at_risk %*% (Xbar * dk))
      U <- w * (d$event * (X - Xbar) - comp)
    }
    list(rows = do.call(cbind, c(list(U), scores)), w = w)
  }
  k <- length(theta_hat)
  h <- 1e-6
  A <- matrix(0, k, k)
  for (j in seq_len(k)) {
    tp <- tm2 <- theta_hat
    tp[j] <- tp[j] + h
    tm2[j] <- tm2[j] - h
    A[, j] <- -(colSums(psi(tp)$rows) - colSums(psi(tm2)$rows)) / (2 * h)
  }
  at_hat <- psi(theta_hat)
  B <- crossprod(rowsum(at_hat$rows, id))
  V <- solve(A) %*% B %*% t(solve(A))
  list(V = V[seq_len(kb), seq_len(kb), drop = FALSE], w = at_hat$w,
       score_sum = colSums(at_hat$rows))
}

fit_tm <- function(dd) {
  list(
    first = ipwtm(exposure = A, family = "binomial", link = "logit", numerator = ~ x,
                  denominator = ~ x + l, id = id, tstart = tstart, timevar = time,
                  type = "first", data = dd),
    cens = ipwtm(exposure = C, family = "binomial", link = "logit", numerator = ~ x,
                 denominator = ~ x + l, id = id, tstart = tstart, timevar = time,
                 type = "cens", data = dd),
    all = ipwtm(exposure = B, family = "binomial", link = "logit", numerator = ~ Blag,
                denominator = ~ Blag + l, id = id, tstart = tstart, timevar = time,
                type = "all", data = dd))
}

test_that("glm MSM with type = 'all' and 'cens' weights equals the numerical-Jacobian sandwich", {
  d <- tm_data()
  set.seed(1)
  shuffled <- d[sample(nrow(d)), ]          # ipwtm() and the MSM use shuffled rows
  w <- fit_tm(shuffled)
  msm <- stats::lm(y ~ B, data = shuffled, weights = w$all$ipw.weights * w$cens$ipw.weights)
  V <- ipwvcov(list(w$all, w$cens), msm)

  specs <- list(
    list(a = d$B, type = "all", Zd = stats::model.matrix(~ Blag + l, d),
         Zn = stats::model.matrix(~ Blag, d)),
    list(a = d$C, type = "cens", Zd = stats::model.matrix(~ x + l, d),
         Zn = stats::model.matrix(~ x, d)))
  theta <- c(stats::coef(msm), stats::coef(w$all$den.mod), stats::coef(w$all$num.mod),
             stats::coef(w$cens$den.mod), stats::coef(w$cens$num.mod))
  ref <- tm_numerical_vcov(d, specs, stats::model.matrix(~ B, d), "gaussian", theta, y = d$y)
  # the reference reproduces ipwtm()'s weights and solves the estimating equations
  ord <- order(shuffled$id, shuffled$time)
  expect_equal(ref$w, (w$all$ipw.weights * w$cens$ipw.weights)[ord], tolerance = 1e-8)
  expect_lt(max(abs(ref$score_sum)) / nrow(d), 1e-6)
  expect_equal(unname(V), unname(ref$V), tolerance = 1e-5)
})

test_that("Cox MSM with type = 'first' x 'cens' weights equals the numerical-Jacobian sandwich", {
  d <- tm_data()
  w <- fit_tm(d)
  ww <- w$first$ipw.weights * w$cens$ipw.weights
  msm <- survival::coxph(survival::Surv(tstart, time, event) ~ A + x, data = d,
                         weights = ww, ties = "breslow")
  V <- ipwvcov(list(w$first, w$cens), msm)
  specs <- list(
    list(a = d$A, type = "first", Zd = stats::model.matrix(~ x + l, d),
         Zn = stats::model.matrix(~ x, d)),
    list(a = d$C, type = "cens", Zd = stats::model.matrix(~ x + l, d),
         Zn = stats::model.matrix(~ x, d)))
  theta <- c(stats::coef(msm), stats::coef(w$first$den.mod), stats::coef(w$first$num.mod),
             stats::coef(w$cens$den.mod), stats::coef(w$cens$num.mod))
  ref <- tm_numerical_vcov(d, specs, cbind(A = d$A, x = d$x), "cox", theta)
  expect_equal(ref$w, ww, tolerance = 1e-8)
  expect_lt(max(abs(ref$score_sum)) / nrow(d), 1e-6)
  expect_equal(unname(V), unname(ref$V), tolerance = 1e-5)
})

test_that("haartdat (subset): exposure x censoring weights and a Cox MSM", {
  data("haartdat", package = "ipw", envir = environment())
  # all 31 patients who died and 120 others, so that the brute-force risk
  # sets of the reference stay small (about 2500 rows)
  died <- unique(haartdat$patient[haartdat$event == 1])
  set.seed(4)
  keep <- c(died, sample(setdiff(unique(haartdat$patient), died), 120))
  d <- haartdat[haartdat$patient %in% keep, ]
  d <- d[order(d$patient, d$fuptime), ]
  w1 <- ipwtm(exposure = haartind, family = "binomial", link = "logit",
              numerator = ~ sex + age, denominator = ~ sex + age + cd4.sqrt,
              id = patient, tstart = tstart, timevar = fuptime, type = "first", data = d)
  w2 <- ipwtm(exposure = dropout, family = "binomial", link = "logit",
              numerator = ~ sex + age, denominator = ~ sex + age + cd4.sqrt,
              id = patient, tstart = tstart, timevar = fuptime, type = "cens", data = d)
  ww <- w1$ipw.weights * w2$ipw.weights
  msm <- survival::coxph(survival::Surv(tstart, fuptime, event) ~ haartind, data = d,
                         weights = ww, ties = "breslow")
  dd <- data.frame(id = d$patient, tstart = d$tstart, time = d$fuptime, event = d$event)
  specs <- list(
    list(a = d$haartind, type = "first", Zd = stats::model.matrix(~ sex + age + cd4.sqrt, d),
         Zn = stats::model.matrix(~ sex + age, d)),
    list(a = d$dropout, type = "cens", Zd = stats::model.matrix(~ sex + age + cd4.sqrt, d),
         Zn = stats::model.matrix(~ sex + age, d)))
  theta <- c(stats::coef(msm), stats::coef(w1$den.mod), stats::coef(w1$num.mod),
             stats::coef(w2$den.mod), stats::coef(w2$num.mod))
  ref <- tm_numerical_vcov(dd, specs, cbind(haartind = d$haartind), "cox", theta)
  expect_equal(ref$w, ww, tolerance = 1e-8)
  expect_equal(unname(ipwvcov(list(w1, w2), msm)), unname(ref$V), tolerance = 1e-5)
})

test_that("Cox pieces for counting-process data and the weights-known sandwich", {
  data("haartdat", package = "ipw", envir = environment())
  w1 <- ipwtm(exposure = haartind, family = "binomial", link = "logit",
              numerator = ~ sex + age, denominator = ~ sex + age + cd4.sqrt,
              id = patient, tstart = tstart, timevar = fuptime, type = "first",
              data = haartdat)
  msm <- survival::coxph(survival::Surv(tstart, fuptime, event) ~ haartind + cluster(patient),
                         data = haartdat, weights = w1$ipw.weights, ties = "breslow")
  parts <- ipw:::.coxph_parts_counting(msm, w1$ipw.weights, stats::model.matrix(msm))
  expect_equal(unname(parts$res[, 1]), unname(stats::residuals(msm, type = "score")),
               tolerance = 1e-10)
  expect_equal(c(parts$A), c(solve(msm$naive.var)), tolerance = 1e-8)
  # robust SE with cluster(patient) treats the weights as known: it is the
  # beta block of the sandwich without the exposure models
  Up <- rowsum(parts$U, haartdat$patient)
  expect_equal(c(solve(parts$A) %*% crossprod(Up) %*% solve(parts$A)), c(msm$var),
               tolerance = 1e-8)
  # the full ipwvcov() runs on all of haartdat and gives a positive variance
  V <- ipwvcov(w1, msm)
  expect_true(is.finite(V[1, 1]) && V[1, 1] > 0)
})

test_that("ipwvcov with ipwtm weights agrees with WeightIt (weightitMSM + glm_weightit)", {
  skip_if_not_installed("WeightIt")
  # Fixed visits 1..3 for everyone, so weightitMSM() (wide data, one model per
  # time point) applies; ipwtm() fits one pooled model whose time-specific
  # columns make it equal to the separate models.
  set.seed(12)
  n <- 400
  wide <- data.frame(id = seq_len(n), L1 = stats::rnorm(n))
  wide$A1 <- stats::rbinom(n, 1, stats::plogis(0.4 * wide$L1))
  wide$L2 <- 0.6 * wide$L1 + 0.3 * wide$A1 + stats::rnorm(n)
  wide$A2 <- stats::rbinom(n, 1, stats::plogis(-0.3 + 0.5 * wide$L2 + 0.8 * wide$A1))
  wide$L3 <- 0.6 * wide$L2 + 0.3 * wide$A2 + stats::rnorm(n)
  wide$A3 <- stats::rbinom(n, 1, stats::plogis(-0.2 + 0.5 * wide$L3 + 0.8 * wide$A2))
  wide$cumA <- wide$A1 + wide$A2 + wide$A3
  wide$y <- 0.5 * wide$cumA + wide$L3 + stats::rnorm(n)

  long <- data.frame(id = rep(wide$id, each = 3), time = rep(1:3, n))
  long$L <- c(t(wide[, c("L1", "L2", "L3")]))
  long$A <- c(t(wide[, c("A1", "A2", "A3")]))
  long$Alag <- c(t(cbind(0, wide$A1, wide$A2)))
  for (k in 1:3) long[[paste0("i", k)]] <- as.numeric(long$time == k)
  long$L_1 <- long$L * long$i1; long$L_2 <- long$L * long$i2; long$L_3 <- long$L * long$i3
  long$Alag_2 <- long$Alag * long$i2; long$Alag_3 <- long$Alag * long$i3
  long$is_last <- long$i3
  long$cumA <- rep(wide$cumA, each = 3) * long$is_last
  long$y <- rep(wide$y, each = 3) * long$is_last
  # one level per non-last row, one shared level for the last rows: the
  # non-last rows fit exactly (score 0), so the coefficients of 'last' and
  # cumA are the person-level MSM y ~ cumA with the final weights
  long$rowfac <- factor(ifelse(long$is_last == 1, "last", paste0("r", seq_len(nrow(long)))))

  w_tm <- ipwtm(exposure = A, family = "binomial", link = "logit",
                denominator = ~ i2 + i3 + L_1 + L_2 + L_3 + Alag_2 + Alag_3,
                id = id, timevar = time, type = "all", data = long)
  W <- WeightIt::weightitMSM(list(A1 ~ L1, A2 ~ L2 + A1, A3 ~ L3 + A2), data = wide,
                             method = "glm")
  expect_equal(w_tm$ipw.weights[long$is_last == 1], unname(W$weights), tolerance = 1e-8)

  msm_long <- stats::lm(y ~ 0 + rowfac + cumA, data = long, weights = w_tm$ipw.weights)
  fit_w <- WeightIt::glm_weightit(y ~ cumA, data = wide, weightit = W, vcov = "asympt")
  # the exposure models are fitted separately by the two packages, each to
  # glm()'s convergence tolerance, so the coefficients agree to about 1e-8
  expect_equal(unname(stats::coef(msm_long)[c("rowfaclast", "cumA")]),
               unname(stats::coef(fit_w)), tolerance = 1e-6)
  V <- ipwvcov(w_tm, msm_long)
  expect_equal(unname(V[c("rowfaclast", "cumA"), c("rowfaclast", "cumA")]),
               unname(stats::vcov(fit_w)), tolerance = 1e-6)
})

test_that("ipwvcov with ipwtm weights gives clear errors for unsupported input", {
  d <- tm_data(100)
  w <- fit_tm(d)
  msm <- stats::lm(y ~ B, data = d, weights = w$all$ipw.weights)

  ws <- ipwtm(exposure = A, family = "survival", numerator = ~ x, denominator = ~ x + l,
              id = id, tstart = tstart, timevar = time, type = "first", data = d)
  expect_error(ipwvcov(ws, stats::lm(y ~ A, data = d, weights = ws$ipw.weights)),
               "only for family = 'binomial'.*cluster bootstrap")
  d$m <- factor(ifelse(d$B == 1, "hi", ifelse(d$l > 0, "mid", "lo")))
  wm <- ipwtm(exposure = m, family = "multinomial", denominator = ~ l, id = id,
              timevar = time, type = "all", data = d)
  expect_error(ipwvcov(wm, stats::lm(y ~ m, data = d, weights = wm$ipw.weights)),
               "only for family = 'binomial'")
  wt <- ipwtm(exposure = B, family = "binomial", link = "logit", numerator = ~ Blag,
              denominator = ~ Blag + l, id = id, timevar = time, type = "all", data = d,
              trunc = 0.05)
  expect_error(ipwvcov(wt, stats::lm(y ~ B, data = d, weights = wt$weights.trunc)),
               "Truncated weights are not supported")
  expect_error(ipwvcov(w$all, stats::lm(y ~ B, data = d[-1, ], weights = w$all$ipw.weights[-1])),
               "fitted on .* rows but the long data have")
  expect_error(ipwvcov(list(w$all, w$cens), msm), "must equal the product")
  expect_error(ipwvcov(w$all, survival::coxph(survival::Surv(time, event) ~ B, data = d,
                                              weights = w$all$ipw.weights)),
               "counting-process data")
  expect_error(ipwvcov(w$all, survival::coxph(survival::Surv(tstart, time, event) ~ B +
                                                survival::strata(x), data = d,
                                              weights = w$all$ipw.weights, ties = "breslow")),
               "strata\\(\\)")
  wp <- ipwpoint(B, "binomial", "logit", denominator = ~ l, data = d)
  expect_error(ipwvcov(list(wp), stats::lm(y ~ B, data = d, weights = wp$ipw.weights)),
               "needs at least one ipwtm\\(\\) result")
  d2 <- d[d$id <= 80, ]
  w_other <- ipwtm(exposure = C, family = "binomial", link = "logit", denominator = ~ l,
                   id = id, timevar = time, type = "cens", data = d2)
  expect_error(ipwvcov(list(w$all, w_other), msm), "same long data")
  skip_if_not_installed("mgcv")
  wg <- ipwtm(exposure = B, family = "binomial", link = "logit", denominator = ~ Blag + s(l),
              id = id, timevar = time, type = "all", data = d)
  expect_error(ipwvcov(wg, stats::lm(y ~ B, data = d, weights = wg$ipw.weights)),
               "smooth terms \\(GAM\\)")
})

test_that("simulation: 95% intervals with ipwtm weights cover the effect about 95%", {
  skip_on_cran()
  set.seed(31)
  covered <- replicate(200, {
    d <- tm_data(300, seed = sample.int(1e6, 1))
    w <- ipwtm(exposure = B, family = "binomial", link = "logit", numerator = ~ Blag,
               denominator = ~ Blag + l, id = id, timevar = time, type = "all", data = d)
    msm <- stats::lm(y ~ B, data = d, weights = w$ipw.weights)
    abs(stats::coef(msm)["B"] - 1) < stats::qnorm(0.975) * sqrt(ipwvcov(w, msm)["B", "B"])
  })
  expect_gt(mean(covered), 0.91)
  expect_lt(mean(covered), 0.985)
})
