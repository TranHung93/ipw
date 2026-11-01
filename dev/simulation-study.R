# Simulation study: do the weights of the new ipw features recover known
# causal effects, and are the bootstrap standard errors right?
#
# In the spirit of the simulated example in van der Wal & Geskus (2011, JSS),
# each scenario generates data in which the true marginal causal effect is
# known, estimates it, and repeats this many times. Results include Monte
# Carlo standard errors (MCSE) so that real bias can be told from simulation
# noise.
#
# Run from the package folder (ipw installed or loaded):
#   Rscript dev/simulation-study.R
# Environment variables (defaults in brackets):
#   IPW_SIM_REPS     replications for the bias scenarios [500]
#   IPW_SIM_REPS_BIG replications for the large-sample and oracle checks [200]
#   IPW_SIM_REPS_SE  replications for the standard-error scenarios 7a, 7b [1000]
#   IPW_SIM_REPS_SE2 replications for scenarios 7c, 7d [500]
#   IPW_SIM_CORES    parallel workers (Linux/macOS) [all cores]
# Results are written to vignettes/articles/simulation-results.csv and
# simulation-se.csv (read by the "Simulation study" article; commit them to
# update it) and reported on GitHub Actions when run there.

suppressPackageStartupMessages({
  library(ipw)
  library(survival)
  library(parallel)
})

env_int <- function(name, default) as.integer(Sys.getenv(name, as.character(default)))
S      <- env_int("IPW_SIM_REPS", 500)
S_BIG  <- env_int("IPW_SIM_REPS_BIG", 200)
S_SE   <- env_int("IPW_SIM_REPS_SE", 1000)
S_SE2  <- env_int("IPW_SIM_REPS_SE2", 500)
CORES  <- if (.Platform$OS.type == "windows") 1L else env_int("IPW_SIM_CORES", detectCores())
RNGkind("L'Ecuyer-CMRG")
set.seed(20261101)

results <- list()
coverage <- list()
t_start <- Sys.time()
timer <- function(label) {
  msg <- sprintf("%s: %.0f s", label, as.numeric(difftime(Sys.time(), t_start, units = "secs")))
  message(msg)
  if (nzchar(Sys.getenv("GITHUB_ACTIONS"))) cat(sprintf("::notice title=progress::%s\n", msg))
}

# Run f() n times in parallel; a replication that fails gives a row of NA and
# is left out (the number of failures is reported).
reps <- function(n, f) {
  out <- mclapply(seq_len(n), function(i) tryCatch(f(), error = function(e) conditionMessage(e)),
                  mc.cores = CORES, mc.set.seed = TRUE)
  bad <- vapply(out, function(x) is.character(x) || is.null(x), logical(1))
  if (any(bad)) message(sum(bad), " of ", n, " replications failed; first error: ", out[bad][[1]])
  if (all(bad)) stop("All replications failed: ", out[[1]], call. = FALSE)
  k <- max(lengths(out[!bad]))
  do.call(rbind, lapply(seq_along(out), function(i) if (bad[i]) rep(NA_real_, k) else out[[i]]))
}

add <- function(scenario, method, est, truth) {
  est <- est[is.finite(est)]
  results[[length(results) + 1]] <<- data.frame(
    scenario = scenario, method = method, truth = truth, mean = mean(est),
    bias = mean(est) - truth, mcse = sd(est) / sqrt(length(est)), emp.sd = sd(est),
    rmse = sqrt(mean((est - truth)^2)), reps = length(est))
}

# coverage rows: est, se and (optionally) interval limits per method
add_cov <- function(scenario, method, est, se, truth, lo = NULL, hi = NULL) {
  ok <- is.finite(est) & is.finite(se)
  if (!is.null(lo)) ok <- ok & is.finite(lo) & is.finite(hi)
  est <- est[ok]; se <- se[ok]
  if (is.null(lo)) { lo <- est - 1.96 * se; hi <- est + 1.96 * se } else { lo <- lo[ok]; hi <- hi[ok] }
  cov <- mean(lo <= truth & hi >= truth)
  coverage[[length(coverage) + 1]] <<- data.frame(
    scenario = scenario, method = method, truth = truth, mean.est = mean(est),
    mean.se = mean(se), emp.sd = sd(est), se.ratio = mean(se) / sd(est),
    coverage = cov, mcse.cov = sqrt(cov * (1 - cov) / length(est)), reps = length(est))
}

wlm <- function(formula, data, w = NULL) {
  # weights go in as a data column: lm() looks for 'weights' in the data and
  # the formula's environment, not in this function
  data$.w <- if (is.null(w)) 1 else w
  coef(lm(formula, data = data, weights = .w))[2]
}

# ---- 1. ATE, ATT and overlap (ATO) weights; ipwprob() ------------------------
e_fun <- function(l) plogis(0.5 + 0.8 * l)
truth1 <- local({
  l <- rnorm(2e6); e <- e_fun(l); tau <- 1 + l
  c(ATE = mean(tau), ATT = sum(e * tau) / sum(e), ATO = sum(e * (1 - e) * tau) / sum(e * (1 - e)))
})
r1 <- reps(S, function() {
  d <- data.frame(l = rnorm(1000))
  d$a <- rbinom(1000, 1, e_fun(d$l))
  d$y <- d$a * (1 + d$l) + d$l + rnorm(1000)
  fit <- function(est) ipwpoint(exposure = a, family = "binomial", link = "logit",
                                denominator = ~ l, data = d, estimand = est)$ipw.weights
  e_true <- e_fun(d$l)
  c(wlm(y ~ a, d),
    wlm(y ~ a, d, ipwpoint(exposure = a, family = "binomial", link = "logit",
                           numerator = ~ 1, denominator = ~ l, data = d)$ipw.weights),
    wlm(y ~ a, d, fit("ATT")), wlm(y ~ a, d, fit("ATO")),
    wlm(y ~ a, d, ipwprob(den = ifelse(d$a == 1, e_true, 1 - e_true))$ipw.weights))
})
sc <- "1. Binary exposure, effect varies"
add(sc, "unweighted (target ATE)", r1[, 1], truth1["ATE"])
add(sc, "ATE weights", r1[, 2], truth1["ATE"])
add(sc, "ATT weights (target ATT)", r1[, 3], truth1["ATT"])
add(sc, "ATO weights (target ATO)", r1[, 4], truth1["ATO"])
add(sc, "ipwprob, true probabilities (ATE)", r1[, 5], truth1["ATE"])
timer("1. ATE / ATT / ATO / ipwprob")

# ---- 2. Count exposures -----------------------------------------------------
count_reps <- function(gen_v) reps(S, function() {
  d <- data.frame(l = rnorm(1000))
  d$v <- gen_v(d$l)
  d$y <- 0.5 * d$v + d$l + rnorm(1000)
  w <- function(fam) ipwpoint(exposure = v, family = fam, numerator = ~ 1,
                              denominator = ~ l, data = d)$ipw.weights
  c(wlm(y ~ v, d), wlm(y ~ v, d, w("poisson")), wlm(y ~ v, d, w("negbin")))
})
r2 <- count_reps(function(l) rnbinom(length(l), size = 2, mu = exp(0.5 + 0.4 * l)))
sc <- "2a. Counts, overdispersed (negative binomial data)"
add(sc, "unweighted", r2[, 1], 0.5)
add(sc, "poisson weights (misspecified)", r2[, 2], 0.5)
add(sc, "negbin weights (correct)", r2[, 3], 0.5)
r2b <- count_reps(function(l) rpois(length(l), exp(0.5 + 0.4 * l)))
sc <- "2b. Counts, Poisson data"
add(sc, "unweighted", r2b[, 1], 0.5)
add(sc, "poisson weights (correct)", r2b[, 2], 0.5)
add(sc, "negbin weights (also valid: theta large)", r2b[, 3], 0.5)
timer("2. counts")

# ---- 3. Skewed continuous exposures -----------------------------------------
# 3a: log-normal dose. Oracle weights use the true densities: dose | l is
#     log-normal(0.5 + 0.5 l, 0.6) and, marginally, log-normal(0.5, sqrt(0.61)).
gen_lnorm <- function(n) {
  d <- data.frame(l = rnorm(n))
  d$dose <- exp(0.5 + 0.5 * d$l + rnorm(n, 0, 0.6))
  d$y <- d$dose + 2 * d$l + rnorm(n)
  d
}
w_dose <- function(d, fam, trunc = NULL) {
  r <- ipwpoint(exposure = dose, family = fam, numerator = ~ 1, denominator = ~ l,
                data = d, trunc = trunc)
  if (is.null(trunc)) r$ipw.weights else r$weights.trunc
}
oracle_lnorm <- function(d) dlnorm(d$dose, 0.5, sqrt(0.61)) / dlnorm(d$dose, 0.5 + 0.5 * d$l, 0.6)
r3 <- reps(S, function() {
  d <- gen_lnorm(2000)
  c(wlm(y ~ dose, d), wlm(y ~ dose, d, w_dose(d, "gaussian")), wlm(y ~ dose, d, w_dose(d, "gamma")),
    wlm(y ~ dose, d, w_dose(d, "lognormal")), wlm(y ~ dose, d, w_dose(d, "lognormal", 0.01)),
    wlm(y ~ dose, d, oracle_lnorm(d)))
})
r3big <- reps(S_BIG, function() {
  d <- gen_lnorm(10000)
  c(wlm(y ~ dose, d, w_dose(d, "lognormal")), wlm(y ~ dose, d, oracle_lnorm(d)))
})
sc <- "3a. Skewed exposure, log-normal dose"
add(sc, "unweighted", r3[, 1], 1)
add(sc, "gaussian weights (misspecified)", r3[, 2], 1)
add(sc, "gamma weights (misspecified)", r3[, 3], 1)
add(sc, "lognormal weights (correct), n = 2000", r3[, 4], 1)
add(sc, "lognormal weights (correct), n = 10000", r3big[, 1], 1)
add(sc, "lognormal weights, truncated 1%", r3[, 5], 1)
add(sc, "oracle weights (true densities), n = 2000", r3[, 6], 1)
add(sc, "oracle weights (true densities), n = 10000", r3big[, 2], 1)

# 3b: gamma dose (shape 3, mean exp(0.5 + 0.3 l)): gamma is the correct model
r3b <- reps(S, function() {
  d <- data.frame(l = rnorm(2000))
  d$dose <- rgamma(2000, shape = 3, rate = 3 / exp(0.5 + 0.3 * d$l))
  d$y <- d$dose + 2 * d$l + rnorm(2000)
  c(wlm(y ~ dose, d), wlm(y ~ dose, d, w_dose(d, "gaussian")), wlm(y ~ dose, d, w_dose(d, "gamma")),
    wlm(y ~ dose, d, w_dose(d, "lognormal")))
})
sc <- "3b. Skewed exposure, gamma dose"
add(sc, "unweighted", r3b[, 1], 1)
add(sc, "gaussian weights (misspecified)", r3b[, 2], 1)
add(sc, "gamma weights (correct)", r3b[, 3], 1)
add(sc, "lognormal weights (misspecified)", r3b[, 4], 1)
timer("3. skewed exposures")

# ---- 4. Curved propensity: smooth terms -------------------------------------
gen_curve <- function(n) {
  d <- data.frame(l = runif(n, -2, 2))
  d$e <- plogis(2 * sin(1.5 * d$l))
  d$a <- rbinom(n, 1, d$e)
  d$y <- 2 * d$a + 3 * sin(1.5 * d$l) + rnorm(n)
  d
}
w_curve <- function(d, f) ipwpoint(exposure = a, family = "binomial", link = "logit",
                                   numerator = ~ 1, denominator = f, data = d)$ipw.weights
oracle_curve <- function(d) ifelse(d$a == 1, mean(d$a) / d$e, (1 - mean(d$a)) / (1 - d$e))
if (requireNamespace("mgcv", quietly = TRUE)) {
  r4 <- reps(S, function() {
    d <- gen_curve(2000)
    c(wlm(y ~ a, d), wlm(y ~ a, d, w_curve(d, ~ l)), wlm(y ~ a, d, w_curve(d, ~ s(l))),
      wlm(y ~ a, d, oracle_curve(d)))
  })
  r4big <- reps(S_BIG, function() {
    d <- gen_curve(10000)
    c(wlm(y ~ a, d, w_curve(d, ~ s(l))), wlm(y ~ a, d, oracle_curve(d)))
  })
  sc <- "4. Curved propensity (mgcv)"
  add(sc, "unweighted", r4[, 1], 2)
  add(sc, "linear logistic weights (misspecified)", r4[, 2], 2)
  add(sc, "s(l) weights, n = 2000", r4[, 3], 2)
  add(sc, "s(l) weights, n = 10000", r4big[, 1], 2)
  add(sc, "oracle weights (true propensity), n = 2000", r4[, 4], 2)
  add(sc, "oracle weights (true propensity), n = 10000", r4big[, 2], 2)
}
timer("4. smooth terms")

# ---- 5. Time-varying exposure, type = "all" ----------------------------------
# U (unmeasured) affects L and Y; L_t is affected by the previous exposure and
# affects the current exposure. MSM: E[Y^abar] = b0 + 1 * sum(abar).
# Regression adjustment for L is biased (L is a collider on A -> L <- U -> Y).
sim_tv <- function(n = 500, n_t = 4, binary = FALSE, curved = FALSE) {
  U <- rnorm(n)
  long <- data.frame(id = rep(seq_len(n), each = n_t), time = rep(seq_len(n_t), n))
  long$l <- NA_real_; long$a <- NA_integer_; long$a.lag <- NA_integer_; long$p.obs <- NA_real_
  prev <- rep(0L, n)
  for (t in seq_len(n_t)) {
    rows <- long$time == t
    l <- U + 0.4 * prev + rnorm(n, 0, 0.5)
    if (binary) {
      lp <- if (curved) -0.3 + 2 * sin(1.5 * l) + 0.6 * prev else -0.5 + 0.8 * l + 0.6 * prev
      p1 <- plogis(lp)
      a <- rbinom(n, 1, p1)
      p_obs <- ifelse(a == 1, p1, 1 - p1)
    } else {
      eta <- 0.8 * l + 0.6 * prev            # proportional odds, cut points 0 and 1.5
      cum <- cbind(plogis(0 - eta), plogis(1.5 - eta), 1)
      a <- findInterval(eta + rlogis(n), c(0, 1.5))
      p_obs <- cum[cbind(seq_len(n), a + 1)] - cbind(0, cum[, 1:2])[cbind(seq_len(n), a + 1)]
    }
    long$l[rows] <- l; long$a.lag[rows] <- prev; long$a[rows] <- a; long$p.obs[rows] <- p_obs
    prev <- as.integer(a)
  }
  cum_a <- tapply(long$a, long$id, sum)
  long$cum <- cum_a[long$id]
  long$y <- (cum_a + 2 * U + rnorm(n))[long$id]
  long
}
last_rows <- function(long) long$time == max(long$time)
wide_l <- function(long) {
  w <- long[last_rows(long), ]
  lw <- reshape(long[, c("id", "time", "l")], idvar = "id", timevar = "time", direction = "wide")
  merge(w, lw, by = "id")
}
ord_fit <- function(n) {
  long <- sim_tv(n)
  long$a <- factor(long$a, levels = 0:2, ordered = TRUE)
  wide <- wide_l(long)
  last <- last_rows(long)
  w_ord <- ipwtm(exposure = a, family = "ordinal", link = "logit", numerator = ~ a.lag,
                 denominator = ~ a.lag + l, id = id, timevar = time, type = "all", data = long)
  w_mult <- ipwtm(exposure = a, family = "multinomial", numerator = ~ a.lag,
                  denominator = ~ a.lag + l, id = id, timevar = time, type = "all", data = long)
  # oracle: 1 / cumulative product of the true probabilities (unstabilized)
  c(wlm(y ~ cum, wide), wlm(y ~ cum + l.1 + l.2 + l.3 + l.4, wide),
    wlm(y ~ cum, wide, w_ord$ipw.weights[last]), wlm(y ~ cum, wide, w_mult$ipw.weights[last]),
    wlm(y ~ cum, wide, (ave(long$p.obs, long$id, FUN = cumprod)^-1)[last]))
}
r5 <- reps(max(100, round(S / 2)), function() ord_fit(500))
r5big <- reps(S_BIG, function() ord_fit(2500)[c(3, 5)])
sc <- "5a. Time-varying ordinal dose (type all)"
add(sc, "unweighted", r5[, 1], 1)
add(sc, "regression adjustment for L", r5[, 2], 1)
add(sc, "ordinal weights (correct), 500 people", r5[, 3], 1)
add(sc, "ordinal weights (correct), 2500 people", r5big[, 1], 1)
add(sc, "multinomial weights, 500 people", r5[, 4], 1)
add(sc, "oracle weights (true probabilities, unstabilized), 500 people", r5[, 5], 1)
add(sc, "oracle weights (true probabilities, unstabilized), 2500 people", r5big[, 2], 1)
timer("5a. time-varying ordinal")

# 5b: smooth terms in ipwtm(): binary exposure over time, curved effect of L
if (requireNamespace("mgcv", quietly = TRUE)) {
  r5b <- reps(max(100, round(S / 2)), function() {
    long <- sim_tv(500, binary = TRUE, curved = TRUE)
    wide <- wide_l(long)
    last <- last_rows(long)
    w <- function(f) ipwtm(exposure = a, family = "binomial", link = "logit", numerator = ~ a.lag,
                           denominator = f, id = id, timevar = time, type = "all",
                           data = long)$ipw.weights[last]
    c(wlm(y ~ cum, wide), wlm(y ~ cum, wide, w(~ a.lag + l)), wlm(y ~ cum, wide, w(~ a.lag + s(l))))
  })
  sc <- "5b. Time-varying binary exposure, curved effect of L (ipwtm + mgcv)"
  add(sc, "unweighted", r5b[, 1], 1)
  add(sc, "linear logistic weights (misspecified)", r5b[, 2], 1)
  add(sc, "s(l) weights", r5b[, 3], 1)
  timer("5b. ipwtm smooth terms")
}

# ---- 6. Recurrent exposure events, family = "survival", type = "all" --------
# 6a/6b: known proportional hazard (baseline 2, log HR 0.5 for x, -1 for xt);
# the true weights are known, so coefficients and weights are compared.
sim_rec_grid <- function(n = 300) {
  d <- data.frame(id = rep(seq_len(n), each = 10), fuptime = rep(seq(0.1, 1, by = 0.1), n))
  d$tstart <- d$fuptime - 0.1
  d$x <- rep(rbinom(n, 1, 0.5), each = 10)
  d$xt <- rbinom(nrow(d), 1, 0.4)
  p <- 1 - exp(-2 * 0.1 * exp(0.5 * d$x - 1.0 * d$xt))
  d$expo <- rbinom(nrow(d), 1, p)
  d$w.true <- 1 / ave(ifelse(d$expo == 1, p, 1 - p), d$id, FUN = cumprod)
  d
}
sim_rec_cont <- function(n = 300) {
  out <- vector("list", n)
  for (i in seq_len(n)) {
    x <- rbinom(1, 1, 0.5)
    bounds <- c(0, sort(runif(rpois(1, 3))), 1)
    rows <- list()
    for (k in seq_len(length(bounds) - 1)) {
      xt <- rbinom(1, 1, 0.4)
      h <- 2 * exp(0.5 * x - 1.0 * xt)
      t0 <- bounds[k]
      repeat {
        t1 <- t0 + rexp(1, h)
        if (t1 >= bounds[k + 1]) { rows[[length(rows) + 1]] <- c(t0, bounds[k + 1], 0, xt); break }
        rows[[length(rows) + 1]] <- c(t0, t1, 1, xt)
        t0 <- t1
      }
    }
    m <- do.call(rbind, rows)
    out[[i]] <- data.frame(id = i, tstart = m[, 1], fuptime = m[, 2], expo = m[, 3], x = x, xt = m[, 4])
  }
  d <- do.call(rbind, out)
  d <- d[d$fuptime > d$tstart, ]
  p <- 1 - exp(-2 * (d$fuptime - d$tstart) * exp(0.5 * d$x - 1.0 * d$xt))
  d$w.true <- 1 / ave(ifelse(d$expo == 1, p, 1 - p), d$id, FUN = cumprod)
  d
}
for (cont in c(FALSE, TRUE)) {
  sc <- if (cont) "6b. Recurrent exposure, continuous time, irregular intervals" else
    "6a. Recurrent exposure, common grid of intervals"
  r6 <- reps(max(100, round(S / 2)), function() {
    d <- if (cont) sim_rec_cont() else sim_rec_grid()
    fit <- ipwtm(exposure = expo, family = "survival", denominator = ~ x + xt, id = id,
                 tstart = tstart, timevar = fuptime, type = "all", data = d)
    c(coef(fit$den.mod), median(abs(log(fit$ipw.weights / d$w.true))))
  })
  add(sc, "Cox coefficient of x", r6[, 1], 0.5)
  add(sc, "Cox coefficient of xt", r6[, 2], -1.0)
  add(sc, "median |log(estimated / true weight)|", r6[, 3], 0)
}

# 6c: causal effect. Exposure events (e.g. doses) happen in continuous time
# with hazard 1.5 exp(0.6 L); L is updated at 4 visits, is lowered by earlier
# events (as treatment raising CD4 lowers the need for more) and depends on an
# unmeasured U that also affects Y. (A positive feedback would be explosive.) Y rises by 1
# per exposure event; MSM: E[Y] = b0 + 1 * (number of events).
sim_doses <- function(n = 400) {
  out <- vector("list", n); yv <- numeric(n); nv <- numeric(n)
  for (i in seq_len(n)) {
    U <- rnorm(1); nprev <- 0; rows <- list()
    for (k in 1:4) {
      l <- U - 0.3 * nprev + rnorm(1, 0, 0.5)
      h <- 1.5 * exp(0.6 * l)
      t0 <- (k - 1) / 4
      repeat {
        t1 <- t0 + rexp(1, h)
        if (t1 >= k / 4) { rows[[length(rows) + 1]] <- c(t0, k / 4, 0, l, k); break }
        rows[[length(rows) + 1]] <- c(t0, t1, 1, l, k)
        nprev <- nprev + 1
        t0 <- t1
      }
    }
    m <- do.call(rbind, rows)
    out[[i]] <- data.frame(id = i, tstart = m[, 1], fuptime = m[, 2], expo = m[, 3], l = m[, 4], visit = m[, 5])
    nv[i] <- nprev; yv[i] <- nprev + 2 * U + rnorm(1)
  }
  d <- do.call(rbind, out)
  d <- d[d$fuptime > d$tstart, ]
  list(long = d, wide = data.frame(id = seq_len(n), n.events = nv, y = yv))
}
r6c <- reps(max(100, round(S / 2)), function() {
  s <- sim_doses()
  d <- s$long; wide <- s$wide
  fit <- ipwtm(exposure = expo, family = "survival", numerator = ~ 1, denominator = ~ l, id = id,
               tstart = tstart, timevar = fuptime, type = "all", data = d)
  w_last <- fit$ipw.weights[!duplicated(d$id, fromLast = TRUE)]
  lv <- tapply(d$l, list(d$id, d$visit), function(z) z[1])
  wide <- cbind(wide, l1 = lv[, 1], l2 = lv[, 2], l3 = lv[, 3], l4 = lv[, 4])
  c(wlm(y ~ n.events, wide), wlm(y ~ n.events + l1 + l2 + l3 + l4, wide),
    wlm(y ~ n.events, wide, w_last))
})
sc <- "6c. Recurrent exposure events: causal effect (survival type all)"
add(sc, "unweighted", r6c[, 1], 1)
add(sc, "regression adjustment for L", r6c[, 2], 1)
add(sc, "survival type all weights", r6c[, 3], 1)
timer("6. recurrent exposure")

# ---- 7. Standard errors -----------------------------------------------------
has_survey <- requireNamespace("survey", quietly = TRUE)
robust_lm_se <- function(d, w) {
  d$.w <- w
  m <- survey::svyglm(y ~ a, design = survey::svydesign(~ 1, weights = ~ .w, data = d))
  sqrt(vcov(m)["a", "a"])
}
point_stat <- function(dd) {
  w <- ipwpoint(exposure = a, family = "binomial", link = "logit", numerator = ~ 1,
                denominator = ~ l, data = dd)$ipw.weights
  dd$.w <- w
  coef(lm(y ~ a, data = dd, weights = .w))["a"]
}
# 7a: good overlap; 7b: the JSS example (logit of the propensity has SD 5)
gen_overlap <- function() {
  d <- data.frame(l = rnorm(500))
  d$a <- rbinom(500, 1, plogis(0.8 * d$l))
  d$y <- 10 * d$a + 2 * d$l + rnorm(500, 0, 5)
  d
}
gen_jss <- function() {
  d <- data.frame(l = rnorm(500, 10, 5))
  d$a <- rbinom(500, 1, plogis(d$l - 10))
  d$y <- 10 * d$a + 0.5 * d$l + rnorm(500, -10, 5)
  d
}
for (sc in c("7a. Linear MSM, good overlap", "7b. Linear MSM, JSS example (poor overlap)")) {
  gen <- if (startsWith(sc, "7a")) gen_overlap else gen_jss
  r7 <- reps(S_SE, function() {
    d <- gen()
    b <- suppressWarnings(ipwboot(d, point_stat, R = 200))
    w <- ipwpoint(exposure = a, family = "binomial", link = "logit", numerator = ~ 1,
                  denominator = ~ l, data = d)$ipw.weights
    c(unname(b$t0), unname(b$se), b$ci[1, 1], b$ci[1, 2],
      if (has_survey) robust_lm_se(d, w) else NA)
  })
  add_cov(sc, "ipwboot, percentile CI", r7[, 1], r7[, 2], 10, r7[, 3], r7[, 4])
  add_cov(sc, "ipwboot SE, Wald CI", r7[, 1], r7[, 2], 10)
  if (has_survey) add_cov(sc, "robust SE (svyglm), Wald CI", r7[, 1], r7[, 5], 10)
  timer(sc)
}

# 7c: weighted Cox MSM (Austin 2016): exponential event times, hazard
# 0.1 exp(-0.7 a + 0.8 l), administrative censoring at 5. The true marginal
# log hazard ratio comes from both potential outcomes of 200000 people each.
cox_truth <- local({
  l <- rnorm(2e5)
  t0 <- rexp(2e5, 0.1 * exp(0.8 * l)); t1 <- rexp(2e5, 0.1 * exp(-0.7 + 0.8 * l))
  pop <- data.frame(a = rep(0:1, each = 2e5), t = pmin(c(t0, t1), 5), ev = as.integer(c(t0, t1) < 5))
  unname(coef(coxph(Surv(t, ev) ~ a, data = pop)))
})
gen_cox <- function(n = 1000) {
  d <- data.frame(l = rnorm(n))
  d$a <- rbinom(n, 1, plogis(-0.5 + 0.8 * d$l))
  tt <- rexp(n, 0.1 * exp(-0.7 * d$a + 0.8 * d$l))
  d$t <- pmin(tt, 5); d$ev <- as.integer(tt < 5)
  d
}
cox_stat <- function(dd) {
  dd$.w <- ipwpoint(exposure = a, family = "binomial", link = "logit", numerator = ~ 1,
                    denominator = ~ l, data = dd)$ipw.weights
  coef(coxph(Surv(t, ev) ~ a, data = dd, weights = .w))
}
r7c <- reps(S_SE2, function() {
  d <- gen_cox()
  d$.w <- ipwpoint(exposure = a, family = "binomial", link = "logit", numerator = ~ 1,
                   denominator = ~ l, data = d)$ipw.weights
  fit <- coxph(Surv(t, ev) ~ a, data = d, weights = .w, robust = TRUE)
  b <- suppressWarnings(ipwboot(d[, c("l", "a", "t", "ev")], cox_stat, R = 200))
  c(unname(coef(fit)), sqrt(fit$naive.var[1, 1]), sqrt(fit$var[1, 1]),
    unname(b$se), b$ci[1, 1], b$ci[1, 2])
})
sc <- "7c. Weighted Cox MSM (log hazard ratio)"
add_cov(sc, "naive SE (coxph), Wald CI", r7c[, 1], r7c[, 2], cox_truth)
add_cov(sc, "robust SE (coxph), Wald CI", r7c[, 1], r7c[, 3], cox_truth)
add_cov(sc, "ipwboot SE, Wald CI", r7c[, 1], r7c[, 4], cox_truth)
add_cov(sc, "ipwboot, percentile CI", r7c[, 1], r7c[, 4], cox_truth, r7c[, 5], r7c[, 6])
timer("7c. Cox MSM")

# 7d: time-varying binary exposure (scenario 5 design), bootstrap resampling
# people (id) versus robust SE of the weighted MSM on one row per person
tv_stat <- function(long) {
  last <- last_rows(long)
  w <- ipwtm(exposure = a, family = "binomial", link = "logit", numerator = ~ a.lag,
             denominator = ~ a.lag + l, id = id, timevar = time, type = "all",
             data = long)$ipw.weights[last]
  wide <- long[last, ]
  wide$.w <- w
  coef(lm(y ~ cum, data = wide, weights = .w))["cum"]
}
r7d <- reps(S_SE2, function() {
  long <- sim_tv(500, binary = TRUE)
  b <- suppressWarnings(ipwboot(long, tv_stat, id = id, R = 200))
  rob <- NA
  if (has_survey) {
    last <- last_rows(long)
    wide <- long[last, ]
    wide$.w <- ipwtm(exposure = a, family = "binomial", link = "logit", numerator = ~ a.lag,
                     denominator = ~ a.lag + l, id = id, timevar = time, type = "all",
                     data = long)$ipw.weights[last]
    m <- survey::svyglm(y ~ cum, design = survey::svydesign(~ 1, weights = ~ .w, data = wide))
    rob <- sqrt(vcov(m)["cum", "cum"])
  }
  c(unname(b$t0), unname(b$se), b$ci[1, 1], b$ci[1, 2], rob)
})
sc <- "7d. Time-varying exposure, bootstrap by person"
add_cov(sc, "ipwboot (id), percentile CI", r7d[, 1], r7d[, 2], 1, r7d[, 3], r7d[, 4])
add_cov(sc, "ipwboot (id) SE, Wald CI", r7d[, 1], r7d[, 2], 1)
if (has_survey) add_cov(sc, "robust SE (svyglm), Wald CI", r7d[, 1], r7d[, 5], 1)
timer("7d. time-varying bootstrap")

# ---- Report -----------------------------------------------------------------
res <- do.call(rbind, results)
cov_tab <- do.call(rbind, coverage)
round_df <- function(df, skip = "reps") {
  num <- vapply(df, is.numeric, logical(1)) & !(names(df) %in% skip)
  df[num] <- lapply(df[num], round, 3)
  df
}
res <- round_df(res); cov_tab <- round_df(cov_tab)

cat(sprintf("\nSimulation study: S = %d, S_BIG = %d, S_SE = %d, S_SE2 = %d, %d cores\n\n",
            S, S_BIG, S_SE, S_SE2, CORES))
print(res, row.names = FALSE)
cat("\nStandard errors\n")
print(cov_tab, row.names = FALSE)
cat(sprintf("\nTrue values: scenario 1 ATE %.3f, ATT %.3f, ATO %.3f; scenario 7c log HR %.4f\n",
            truth1["ATE"], truth1["ATT"], truth1["ATO"], cox_truth))

dir.create(file.path("vignettes", "articles"), showWarnings = FALSE, recursive = TRUE)
write.csv(res, file.path("vignettes", "articles", "simulation-results.csv"), row.names = FALSE)
write.csv(cov_tab, file.path("vignettes", "articles", "simulation-se.csv"), row.names = FALSE)

if (nzchar(Sys.getenv("GITHUB_ACTIONS"))) {
  md_table <- function(df) {
    c(paste("|", paste(names(df), collapse = " | "), "|"),
      paste("|", paste(rep("---", ncol(df)), collapse = " | "), "|"),
      apply(df, 1, function(r) paste("|", paste(trimws(format(r)), collapse = " | "), "|")))
  }
  summary_file <- Sys.getenv("GITHUB_STEP_SUMMARY")
  if (nzchar(summary_file)) {
    cat(c("## Simulation study", "", md_table(res), "", "### Standard errors", "", md_table(cov_tab)),
        file = summary_file, sep = "\n", append = TRUE)
  }
  # Annotations are truncated when long: emit the results as CSV in chunks
  esc <- function(x) gsub("\n", "%0A", gsub("%", "%25", x, fixed = TRUE), fixed = TRUE)
  emit <- function(df, title) {
    lines <- capture.output(write.csv(df, row.names = FALSE))
    chunks <- split(lines[-1], ceiling(seq_along(lines[-1]) / 9))
    for (i in seq_along(chunks)) {
      cat(sprintf("::notice title=%s part %d::%s\n", title, i, esc(paste(c(lines[1], chunks[[i]]), collapse = "\n"))))
    }
  }
  emit(res, "results")
  emit(cov_tab, "standard errors")
  cat(sprintf("::notice title=truth::ATE %.4f, ATT %.4f, ATO %.4f, Cox log HR %.4f\n",
              truth1["ATE"], truth1["ATT"], truth1["ATO"], cox_truth))
}
timer("finished")
