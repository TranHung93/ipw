# Audit tests for ipwtm(): reference agreement, invariances, input checks.

make_long <- function(nid = 150, seed = 2) {
  set.seed(seed)
  d <- data.frame(id = rep(seq_len(nid), each = 5), time = rep(1:5, nid))
  d$l <- stats::rnorm(nrow(d))
  d$f <- factor(sample(c("a", "b"), nrow(d), TRUE))
  d$a <- stats::rbinom(nrow(d), 1, stats::plogis(-1 + d$l))
  d
}
tm_w <- function(d, ...) {
  ipwtm(exposure = a, family = "binomial", link = "logit", numerator = ~ 1,
        denominator = ~ l + f, id = id, timevar = time, type = "first", data = d, ...)$ipw.weights
}

# Survival data in which subjects have their own visit times
make_surv_irregular <- function(n = 300, seed = 3) {
  set.seed(seed)
  rows <- lapply(seq_len(n), function(i) {
    tt <- sort(unique(round(cumsum(stats::runif(6, 5, 15)))))
    data.frame(id = i, time = tt, l = stats::rnorm(length(tt)))
  })
  d <- do.call(rbind, rows)
  d$tstart <- stats::ave(d$time, d$id, FUN = function(t) c(0, utils::head(t, -1)))
  d$a <- stats::rbinom(nrow(d), 1, 1 - exp(-0.02 * (d$time - d$tstart) * exp(0.7 * d$l)))
  d[stats::ave(d$a, d$id, FUN = function(x) cumsum(cumsum(x)) <= 1) == 1, ]
}
# Correct weights: P(row) uses the baseline hazard over the row's own (tstart, time]
surv_first_reference <- function(d) {
  fit <- survival::coxph(survival::Surv(tstart, time, a) ~ l, data = d, method = "efron")
  bh <- survival::basehaz(fit, centered = TRUE)
  H <- function(t) {
    v <- stats::approx(bh$time, bh$hazard, t, method = "constant", rule = 2)$y
    v[t < min(bh$time)] <- 0
    v
  }
  h <- (H(d$time) - H(d$tstart)) * stats::predict(fit, type = "risk")
  p <- ifelse(d$a == 1, 1 - exp(-h), exp(-h))
  stats::ave(1 / p, d$id, FUN = cumprod)
}

test_that("binomial type = 'all' equals by-hand pooled logistic weights and ipwprobtm", {
  d <- make_long()
  w <- ipwtm(a, "binomial", "logit", ~ 1, ~ l, id = id, timevar = time, type = "all", data = d)
  pd <- stats::fitted(stats::glm(a ~ l, binomial, d))
  pn <- stats::fitted(stats::glm(a ~ 1, binomial, d))
  byhand <- stats::ave(ifelse(d$a == 1, pn / pd, (1 - pn) / (1 - pd)), d$id, FUN = cumprod)
  expect_equal(w$ipw.weights, unname(byhand))
  expect_equal(w$ipw.weights,
               ipwprobtm(den = ifelse(d$a == 1, pd, 1 - pd), num = ifelse(d$a == 1, pn, 1 - pn),
                         id = d$id, timevar = d$time)$ipw.weights)
})

test_that("survival type = 'first' uses the hazard over each row's own interval", {
  d <- make_surv_irregular()
  w <- ipwtm(a, "survival", denominator = ~ l, id = id, tstart = tstart, timevar = time,
             type = "first", data = d)
  expect_equal(w$ipw.weights, unname(surv_first_reference(d)), tolerance = 1e-6)
})

test_that("survival: type = 'first' and type = 'all' agree when every row is at risk once", {
  d <- make_surv_irregular()
  wf <- ipwtm(a, "survival", denominator = ~ l, id = id, tstart = tstart, timevar = time,
              type = "first", data = d)$ipw.weights
  wa <- ipwtm(a, "survival", denominator = ~ l, id = id, tstart = tstart, timevar = time,
              type = "all", data = d)$ipw.weights
  expect_equal(wf, wa, tolerance = 1e-6)
})

test_that("weights do not depend on data class, id type, row order or scale", {
  d <- make_long()
  ref <- tm_w(d)
  skip_if_not_installed("tibble")
  expect_equal(tm_w(tibble::as_tibble(d)), ref)
  d2 <- d; d2$id <- paste0("p", d2$id)
  expect_equal(tm_w(d2), ref)
  d3 <- d; d3$id <- factor(d3$id, levels = c(unique(d3$id), 9999))
  expect_equal(tm_w(d3), ref)
  o <- sample(nrow(d))
  expect_equal(tm_w(d[o, ])[order(o)], ref)
  d4 <- d; d4$l <- d4$l / 1000
  expect_equal(tm_w(d4), ref, tolerance = 1e-6)
  skip_if_not_installed("data.table")
  expect_equal(tm_w(data.table::as.data.table(d)), ref)
})

test_that("a factor 0/1 exposure gives the same weights whatever its level order", {
  d <- make_long()
  ref <- tm_w(d)
  d2 <- d; d2$a <- factor(d2$a, levels = c(1, 0))
  expect_equal(tm_w(d2), ref)
})

test_that("a covariate called 'selvar' is used as the user's covariate", {
  d <- make_long()
  set.seed(9); z <- stats::rnorm(nrow(d))
  w_z <- ipwtm(a, "binomial", "logit", ~ 1, ~ l + z, id = id, timevar = time,
               type = "first", data = transform(d, z = z))$ipw.weights
  w_s <- ipwtm(a, "binomial", "logit", ~ 1, ~ l + selvar, id = id, timevar = time,
               type = "first", data = transform(d, selvar = z))$ipw.weights
  expect_equal(w_s, w_z)
})

test_that("duplicated (id, timevar) rows are rejected", {
  d <- make_long()
  expect_error(tm_w(rbind(d, d[1:3, ])), regexp = "(only once|duplicate)")
})

test_that("infinite covariate values give an informative ipw error", {
  d <- make_long()
  dinf <- d; dinf$l[1] <- Inf
  expect_error(tm_w(dinf), regexp = "(infinite|non-finite)")
})

test_that("informative errors for missing values and unsupported combinations", {
  d <- make_long()
  dna <- d; dna$l <- NA
  expect_error(tm_w(dna), regexp = "Missing values")
  expect_error(ipwtm(a, "survival", denominator = ~ l, id = id, timevar = time, data = d),
               regexp = "tstart")
  expect_error(ipwtm(a, "gaussian", numerator = ~ 1, denominator = ~ l, id = id,
                     timevar = time, type = "first", data = d), regexp = "first")
})

test_that("ordinal exposure without a defined order is rejected", {
  d <- make_long()
  d$m <- sample(c("low", "mid", "high"), nrow(d), TRUE)
  expect_error(ipwtm(m, "ordinal", "logit", ~ 1, ~ l, id = id, timevar = time, type = "all",
                     data = d), regexp = "order")
})

test_that("survival type = 'all' uses the hazard over each row's own interval", {
  d <- make_surv_irregular()
  w <- ipwtm(a, "survival", denominator = ~ l, id = id, tstart = tstart, timevar = time,
             type = "all", data = d)
  expect_equal(w$ipw.weights, unname(surv_first_reference(d)), tolerance = 1e-6)
})

test_that("a subject with a single row is handled", {
  d <- make_long()
  d <- d[!(d$id == 1 & d$time > 1), ]
  w <- tm_w(d)
  expect_length(w, nrow(d))
  expect_true(all(is.finite(w)))
})

test_that("simulation: ipwtm recovers a two-period MSM with treatment-confounder feedback", {
  skip_on_cran()
  set.seed(2027)
  est <- t(replicate(150, {
    n <- 800; l0 <- stats::rnorm(n); a0 <- stats::rbinom(n, 1, stats::plogis(l0))
    l1 <- 0.5 * a0 + l0 + stats::rnorm(n)
    a1 <- stats::rbinom(n, 1, stats::plogis(-0.5 + l1 + 0.5 * a0))
    y <- 1 + a0 + a1 + l1 + stats::rnorm(n)
    long <- data.frame(id = rep(1:n, 2), time = rep(1:2, each = n), a = c(a0, a1),
                       l = c(l0, l1), alag = c(rep(0, n), a0))
    w <- ipwtm(a, "binomial", "logit", numerator = ~ factor(time) + alag,
               denominator = ~ factor(time) + l + alag, id = id, timevar = time,
               type = "all", data = long)$ipw.weights[long$time == 2]
    stats::coef(stats::lm(y ~ a0 + a1, weights = w))[c("a0", "a1")]
  }))
  mcse <- apply(est, 2, stats::sd) / sqrt(nrow(est))
  expect_lt(abs(mean(est[, "a0"]) - 1.5), 4 * mcse[["a0"]])
  expect_lt(abs(mean(est[, "a1"]) - 1.0), 4 * mcse[["a1"]])
})
