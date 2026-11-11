# Results that intentionally differ from ipw 1.2.1.1 (listed in NEWS.md) are
# checked against by-hand references here instead of the 1.2.1.1 fixture.

chg_irregular <- function() {
  set.seed(31)
  irr <- do.call(rbind, lapply(1:200, function(i) {
    t <- cumsum(stats::runif(sample(2:6, 1), 0.5, 2))
    data.frame(id = i, time = t, tstart = c(0, utils::head(t, -1)), l = stats::rnorm(length(t)))
  }))
  irr$a <- stats::rbinom(nrow(irr), 1, 1 - exp(-0.3 * exp(0.6 * irr$l) * (irr$time - irr$tstart)))
  irr$cens <- stats::rbinom(nrow(irr), 1, 0.1)
  irr
}
# By hand: Breslow-type risk sets with Efron ties; P(row) over its own (tstart, time]
chg_ref <- function(d, e, type) {
  sel <- stats::ave(d[[e]], d$id, FUN = function(v) { k <- match(1, v); if (is.na(k)) rep(1, length(v)) else as.numeric(seq_along(v) <= k) })
  act <- d[sel == 1, ]
  fit <- survival::coxph(stats::as.formula(paste0("survival::Surv(tstart, time, ", e, ") ~ l")), act, ties = "efron")
  r <- exp(stats::coef(fit) * d$l)
  et <- sort(unique(act$time[act[[e]] == 1]))
  inc <- vapply(et, function(t) {
    S0 <- sum(r[sel == 1][act$tstart < t & act$time >= t]); D <- which(act$time == t & act[[e]] == 1)
    k <- length(D); sum(1 / (S0 - (0:(k - 1)) / k * sum(r[sel == 1][D])))
  }, numeric(1))
  H <- function(s) vapply(s, function(x) sum(inc[et <= x]), numeric(1))
  dH <- (H(d$time) - H(d$tstart)) * r
  p <- if (type == "first") ifelse(d[[e]] == 1, 1 - exp(-dH), exp(-dH)) else exp(-dH)
  p[sel == 0] <- 1
  # numerator ~ 1: same with r = 1
  et0 <- et
  inc0 <- vapply(et0, function(t) { S0 <- sum(act$tstart < t & act$time >= t); k <- sum(act$time == t & act[[e]] == 1)
    sum(1 / (S0 - (0:(k - 1)) / k * k)) }, numeric(1))
  H0 <- function(s) vapply(s, function(x) sum(inc0[et0 <= x]), numeric(1))
  dH0 <- H0(d$time) - H0(d$tstart)
  p0 <- if (type == "first") ifelse(d[[e]] == 1, 1 - exp(-dH0), exp(-dH0)) else exp(-dH0)
  p0[sel == 0] <- 1
  stats::ave(p0 / p, d$id, FUN = cumprod)
}

test_that("survival first/cens on subject-specific visit times use each row's own interval (changed after 1.3.2)", {
  d <- chg_irregular()
  w1 <- ipwtm(a, "survival", NULL, ~ 1, ~ l, id = id, tstart = tstart, timevar = time, type = "first", data = d)$ipw.weights
  expect_equal(w1, unname(chg_ref(d, "a", "first")), tolerance = 1e-10)
  w2 <- ipwtm(cens, "survival", NULL, ~ 1, ~ l, id = id, tstart = tstart, timevar = time, type = "cens", data = d)$ipw.weights
  expect_equal(w2, unname(chg_ref(d, "cens", "cens")), tolerance = 1e-10)
})

test_that("a 0/1 factor exposure with levels c(1, 0) is used by its labels (changed after 1.3.2)", {
  set.seed(32)
  pt <- data.frame(l = stats::rnorm(500)); pt$a <- stats::rbinom(500, 1, stats::plogis(pt$l))
  pt$af <- factor(pt$a, levels = c(1, 0))
  p <- stats::fitted(stats::glm(a ~ l, stats::binomial, pt))
  expect_equal(ipwpoint(af, "binomial", "logit", ~ 1, ~ l, pt)$ipw.weights,
               unname(ifelse(pt$a == 1, mean(pt$a) / p, (1 - mean(pt$a)) / (1 - p))), tolerance = 1e-8)
})

test_that("a character ordinal exposure is an error (changed after 1.3.2)", {
  set.seed(32)
  pt <- data.frame(l = stats::rnorm(100)); pt$m <- sample(c("low", "mid", "high"), 100, TRUE)
  expect_error(ipwpoint(m, "ordinal", "logit", ~ 1, ~ l, pt), regexp = "order")
})
