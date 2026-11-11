# Second review: ipwtm() against independent by-hand references, and findings.

# Survival data in which every subject has its own visit times
r2_irregular <- function(n = 150, seed = 1, ties = FALSE, start0 = 0) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n), function(i) {
    t <- cumsum(stats::runif(sample(2:6, 1), 0.5, 2))
    if (ties) t <- unique(ceiling(t))
    l <- stats::rnorm(length(t)); x <- stats::rbinom(1, 1, 0.5)
    h <- 1 - exp(-0.3 * exp(0.6 * l + 0.4 * x) * diff(c(start0, t)))
    data.frame(id = i, tstart = c(start0, utils::head(t, -1)), time = t, l = l, x = x,
               a = stats::rbinom(length(t), 1, h), cens = stats::rbinom(length(t), 1, 0.1))
  }))
}
r2_selvar <- function(x, id) {
  stats::ave(x, id, FUN = function(v) {
    k <- match(1, v); if (is.na(k)) rep(1, length(v)) else as.numeric(seq_along(v) <= k)
  })
}
# Efron cumulative baseline hazard from risk sets, by hand (no basehaz/survfit)
r2_H0 <- function(start, stop, event, r) {
  et <- sort(unique(stop[event == 1]))
  inc <- vapply(et, function(t) {
    S0 <- sum(r[start < t & stop >= t]); D <- which(stop == t & event == 1)
    d <- length(D); sum(1 / (S0 - (0:(d - 1)) / d * sum(r[D])))
  }, numeric(1))
  cum <- cumsum(inc)
  function(t) vapply(t, function(s) { k <- which(et <= s); if (length(k)) cum[max(k)] else 0 }, numeric(1))
}
# By-hand stabilised weights (numerator ~ x, denominator ~ x + l), data sorted by id, time
r2_surv_ref <- function(d, exposure, type) {
  d <- d[order(d$id, d$time), ]
  sel <- r2_selvar(d[[exposure]], d$id); act <- d[sel == 1, ]
  pfun <- function(covs) {
    f <- stats::as.formula(paste0("survival::Surv(tstart, time, ", exposure, ") ~ ", paste(covs, collapse = "+")))
    m <- survival::coxph(f, data = act, ties = "efron")
    X <- stats::model.matrix(stats::as.formula(paste("~", paste(covs, collapse = "+"))), d)[, -1, drop = FALSE]
    r <- as.numeric(exp(X %*% stats::coef(m)))
    H0 <- r2_H0(act$tstart, act$time, act[[exposure]], r[sel == 1])
    dH <- (H0(d$time) - H0(d$tstart)) * r
    p <- if (type == "first") ifelse(d[[exposure]] == 1, 1 - exp(-dH), exp(-dH)) else exp(-dH)
    p[sel == 0] <- 1
    p
  }
  list(w = stats::ave(pfun("x") / pfun(c("x", "l")), d$id, FUN = cumprod), key = paste(d$id, d$time))
}
r2_surv_w <- function(d, exposure, type) {
  r <- do.call(ipwtm, list(exposure = exposure, family = "survival", numerator = ~ x,
                           denominator = ~ x + l, id = "id", tstart = "tstart", timevar = "time",
                           type = type, data = d))
  r$ipw.weights
}

test_that("survival first/cens equal an independent by-hand Efron interval hazard", {
  sets <- list(irregular = r2_irregular(120, 1), ties = r2_irregular(120, 2, ties = TRUE),
               negative_start = r2_irregular(120, 3, start0 = -5))
  for (nm in names(sets)) for (tp in c("first", "cens")) {
    d <- sets[[nm]]
    d <- d[sample(nrow(d)), ]   # unsorted input
    e <- if (tp == "first") "a" else "cens"
    ref <- r2_surv_ref(d, e, tp)
    w <- r2_surv_w(d, e, tp)[match(ref$key, paste(d$id, d$time))]
    expect_equal(w, unname(ref$w), tolerance = 1e-10, label = paste(nm, tp))
  }
})

test_that("binomial type = 'cens' equals by-hand pooled logistic weights", {
  d <- r2_irregular(150, 4)
  d <- d[order(d$id, d$time), ]
  sel <- r2_selvar(d$cens, d$id)
  pf <- function(f) {
    p <- stats::predict(stats::glm(f, stats::binomial, d[sel == 1, ]), d, type = "response")
    p <- 1 - p; p[sel == 0] <- 1; p
  }
  ref <- stats::ave(pf(cens ~ x) / pf(cens ~ x + l), d$id, FUN = cumprod)
  w <- ipwtm(cens, "binomial", "logit", ~ x, ~ x + l, id = id, timevar = time, type = "cens", data = d)$ipw.weights
  expect_equal(w, unname(ref), tolerance = 1e-10)
})

test_that("multinomial and ordinal without numerator equal by-hand fits", {
  d <- r2_irregular(200, 5)
  d <- d[order(d$id, d$time), ]
  d$m <- pmin(d$a + stats::rbinom(nrow(d), 1, stats::plogis(d$l)), 2)
  d$o <- factor(d$m, levels = 0:2, ordered = TRUE)
  mm <- nnet::multinom(factor(m) ~ x + l, d, trace = FALSE); pp <- stats::fitted(mm)
  expect_equal(ipwtm(m, "multinomial", NULL, NULL, ~ x + l, id = id, timevar = time, type = "all", data = d)$ipw.weights,
               unname(stats::ave(1 / pp[cbind(seq_len(nrow(d)), d$m + 1)], d$id, FUN = cumprod)), tolerance = 1e-6)
  pr <- MASS::polr(o ~ x + l, d, Hess = TRUE); pp <- stats::fitted(pr)
  expect_equal(ipwtm(o, "ordinal", "logit", NULL, ~ x + l, id = id, timevar = time, type = "all", data = d)$ipw.weights,
               unname(stats::ave(1 / pp[cbind(seq_len(nrow(d)), as.integer(d$o))], d$id, FUN = cumprod)), tolerance = 1e-8)
})

test_that("S-02: offset() terms are kept in ipwtm models", {
  d <- r2_irregular(150, 6)
  d$off <- 0.5 * d$l
  rb <- ipwtm(a, "binomial", "logit", NULL, ~ x + offset(off), id = id, timevar = time, type = "first", data = d)
  rs <- ipwtm(a, "survival", NULL, NULL, ~ x + offset(off), id = id, tstart = tstart, timevar = time,
              type = "first", data = d)
  expect_true(any(grepl("offset(off)", deparse(stats::formula(rb$den.mod)), fixed = TRUE)))
  expect_true(any(grepl("offset(off)", deparse(stats::formula(rs$den.mod)), fixed = TRUE)))
})

test_that("S-05: ties = 'breslow' passed through ... is used by the Cox models", {
  d <- r2_irregular(200, 2, ties = TRUE)
  sel <- r2_selvar(d$a[order(d$id, d$time)], d$id[order(d$id, d$time)])
  ds <- d[order(d$id, d$time), ][sel == 1, ]
  mb <- survival::coxph(survival::Surv(tstart, time, a) ~ x + l, ds, ties = "breslow")
  r <- ipwtm(a, "survival", NULL, NULL, ~ x + l, id = id, tstart = tstart, timevar = time,
             type = "first", data = d, ties = "breslow")
  expect_equal(unname(stats::coef(r$den.mod)), unname(stats::coef(mb)), tolerance = 1e-8)
})

test_that("S-07: intervals with tstart >= timevar give an informative error", {
  d <- r2_irregular(100, 3)
  d <- d[order(d$id, d$time), ]
  d$tstart[1] <- d$time[1]
  expect_error(ipwtm(a, "survival", NULL, NULL, ~ l, id = id, tstart = tstart, timevar = time,
                     type = "first", data = d), regexp = "tstart")
})

test_that("S-01: an unordered factor is not silently used in alphabetical order as an ordinal exposure", {
  d <- r2_irregular(150, 7)
  d$m <- factor(c("low", "mid", "high")[pmin(d$a + stats::rbinom(nrow(d), 1, 0.4), 2) + 1])
  expect_warning(ipwtm(m, "ordinal", "logit", ~ 1, ~ l, id = id, timevar = time, type = "all", data = d),
                 regexp = "order")
})

test_that("simulation: time-varying feedback MSM is recovered (bias within 3 MCSE)", {
  skip_on_cran()
  est <- vapply(1:200, function(r) {
    set.seed(3000 + r); n <- 1000; U <- stats::rnorm(n); alag <- rep(0, n); rows <- vector("list", 3)
    for (k in 1:3) {
      L <- U - 0.7 * alag + stats::rnorm(n, 0, 0.5)
      A <- stats::rbinom(n, 1, stats::plogis(-0.5 + 0.8 * L + 0.5 * alag))
      rows[[k]] <- data.frame(id = seq_len(n), time = k, L = L, alag = alag, a = A); alag <- A
    }
    long <- do.call(rbind, rows)
    y <- 1 + tapply(long$a, long$id, sum) + U + stats::rnorm(n)
    w <- ipwtm(a, "binomial", "logit", ~ alag, ~ L + alag, id = id, timevar = time, type = "all", data = long)$ipw.weights
    last <- long$time == 3
    stats::coef(stats::lm(y ~ cumA, weights = w[last][order(long$id[last])],
                          data = data.frame(y = y, cumA = tapply(long$a, long$id, sum))))[["cumA"]]
  }, numeric(1))
  expect_lt(abs(mean(est) - 1), 3 * stats::sd(est) / sqrt(length(est)))
})
