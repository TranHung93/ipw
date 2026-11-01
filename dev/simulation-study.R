# Simulation study: do the weights of the new ipw features recover known
# causal effects?
#
# In the spirit of the simulated example in van der Wal & Geskus (2011, JSS),
# each scenario generates data in which the true marginal causal effect is
# known, estimates it with the new weights, and repeats this S times. Bias,
# empirical SD and RMSE are compared with naive analyses or with weights from a
# misspecified model.
#
# Run from the package folder (ipw installed or loaded):
#   Rscript dev/simulation-study.R
# Number of replications: environment variable IPW_SIM_REPS (default 500;
# scenario 7 uses a fifth of it). Results are printed, written to
# dev/simulation-results.csv, and reported on GitHub Actions when run there.

suppressPackageStartupMessages({
  library(ipw)
  library(survival)
})

S <- as.integer(Sys.getenv("IPW_SIM_REPS", "500"))
set.seed(20261101)
results <- list()
coverage <- list()
t_start <- Sys.time()

add <- function(scenario, method, est, truth) {
  est <- est[is.finite(est)]
  results[[length(results) + 1]] <<- data.frame(
    scenario = scenario, method = method, truth = truth, mean = mean(est),
    bias = mean(est) - truth, emp.sd = sd(est), rmse = sqrt(mean((est - truth)^2)),
    reps = length(est))
}
wlm <- function(formula, data, w = NULL) {
  fit <- if (is.null(w)) lm(formula, data = data) else lm(formula, data = data, weights = w)
  coef(fit)[2]
}
# Run f() n times; a replication that fails (e.g. a model that does not
# converge) gives a row of NA and is left out of the summaries.
reps <- function(n, f) {
  out <- lapply(seq_len(n), function(i) tryCatch(f(), error = function(e) NULL))
  k <- max(lengths(out))
  do.call(rbind, lapply(out, function(x) if (is.null(x)) rep(NA_real_, k) else x))
}
timer <- function(label) message(sprintf("%-55s %5.0f s", label,
                                         as.numeric(difftime(Sys.time(), t_start, units = "secs"))))

# ---- 1. ATE, ATT and overlap (ATO) weights; ipwprob() ------------------------
# Effect varies with l: y = a (1 + l) + l + e, exposure P(a = 1 | l) = e(l).
e_fun <- function(l) plogis(0.5 + 0.8 * l)
truth1 <- local({
  l <- rnorm(2e6); e <- e_fun(l); tau <- 1 + l
  c(ATE = mean(tau), ATT = sum(e * tau) / sum(e), ATO = sum(e * (1 - e) * tau) / sum(e * (1 - e)))
})
r1 <- reps(S, function() {
  d <- data.frame(l = rnorm(1000))
  d$a <- rbinom(1000, 1, e_fun(d$l))
  d$y <- d$a * (1 + d$l) + d$l + rnorm(1000)
  fit <- function(est) {
    ipwpoint(exposure = a, family = "binomial", link = "logit", denominator = ~ l,
             data = d, estimand = est)$ipw.weights
  }
  e_true <- e_fun(d$l)
  c(unweighted = wlm(y ~ a, d),
    ATE = wlm(y ~ a, d, ipwpoint(exposure = a, family = "binomial", link = "logit",
                                 numerator = ~ 1, denominator = ~ l, data = d)$ipw.weights),
    ATT = wlm(y ~ a, d, fit("ATT")),
    ATO = wlm(y ~ a, d, fit("ATO")),
    ipwprob = wlm(y ~ a, d, ipwprob(den = ifelse(d$a == 1, e_true, 1 - e_true))$ipw.weights))
})
add("1. Binary exposure, effect varies", "unweighted (target ATE)", r1[, 1], truth1["ATE"])
add("1. Binary exposure, effect varies", "ATE weights", r1[, 2], truth1["ATE"])
add("1. Binary exposure, effect varies", "ATT weights (target ATT)", r1[, 3], truth1["ATT"])
add("1. Binary exposure, effect varies", "ATO weights (target ATO)", r1[, 4], truth1["ATO"])
add("1. Binary exposure, effect varies", "ipwprob, true probabilities (ATE)", r1[, 5], truth1["ATE"])
timer("1. ATE / ATT / ATO / ipwprob done")

# ---- 2. Count exposure: Poisson vs negative binomial ------------------------
r2 <- reps(S, function() {
  d <- data.frame(l = rnorm(1000))
  d$v <- rnbinom(1000, size = 2, mu = exp(0.5 + 0.4 * d$l))
  d$y <- 0.5 * d$v + d$l + rnorm(1000)
  w <- function(fam) ipwpoint(exposure = v, family = fam, numerator = ~ 1,
                              denominator = ~ l, data = d)$ipw.weights
  c(wlm(y ~ v, d), wlm(y ~ v, d, w("poisson")), wlm(y ~ v, d, w("negbin")))
})
add("2. Count exposure (overdispersed)", "unweighted", r2[, 1], 0.5)
add("2. Count exposure (overdispersed)", "poisson weights (misspecified)", r2[, 2], 0.5)
add("2. Count exposure (overdispersed)", "negbin weights (correct)", r2[, 3], 0.5)
timer("2. counts done")

# ---- 3. Skewed continuous exposure: gaussian vs gamma vs log-normal ---------
r3 <- reps(S, function() {
  d <- data.frame(l = rnorm(2000))
  d$dose <- exp(0.5 + 0.5 * d$l + rnorm(2000, 0, 0.6))
  d$y <- d$dose + 2 * d$l + rnorm(2000)
  w <- function(fam, trunc = NULL) {
    r <- ipwpoint(exposure = dose, family = fam, numerator = ~ 1, denominator = ~ l,
                  data = d, trunc = trunc)
    if (is.null(trunc)) r$ipw.weights else r$weights.trunc
  }
  c(wlm(y ~ dose, d), wlm(y ~ dose, d, w("gaussian")), wlm(y ~ dose, d, w("gamma")),
    wlm(y ~ dose, d, w("lognormal")), wlm(y ~ dose, d, w("lognormal", 0.01)))
})
add("3. Skewed exposure (log-normal dose)", "unweighted", r3[, 1], 1)
add("3. Skewed exposure (log-normal dose)", "gaussian weights (misspecified)", r3[, 2], 1)
add("3. Skewed exposure (log-normal dose)", "gamma weights (misspecified)", r3[, 3], 1)
add("3. Skewed exposure (log-normal dose)", "lognormal weights (correct)", r3[, 4], 1)
add("3. Skewed exposure (log-normal dose)", "lognormal weights, truncated 1%", r3[, 5], 1)
timer("3. skewed exposure done")

# ---- 4. Smooth terms: curved propensity -------------------------------------
if (requireNamespace("mgcv", quietly = TRUE)) {
  r4 <- reps(S, function() {
    d <- data.frame(l = runif(2000, -2, 2))
    d$a <- rbinom(2000, 1, plogis(2 * sin(1.5 * d$l)))
    d$y <- 2 * d$a + 3 * sin(1.5 * d$l) + rnorm(2000)
    w <- function(f) ipwpoint(exposure = a, family = "binomial", link = "logit",
                              numerator = ~ 1, denominator = f, data = d)$ipw.weights
    c(wlm(y ~ a, d), wlm(y ~ a, d, w(~ l)), wlm(y ~ a, d, w(~ s(l))))
  })
  add("4. Curved propensity (mgcv)", "unweighted", r4[, 1], 2)
  add("4. Curved propensity (mgcv)", "linear logistic weights (misspecified)", r4[, 2], 2)
  add("4. Curved propensity (mgcv)", "s(l) weights", r4[, 3], 2)
  timer("4. smooth terms done")
}

# ---- 5. Time-varying ordinal dose, type = "all" -----------------------------
# U (unmeasured) affects L and Y; L_t is affected by the previous dose and
# affects the current dose; dose follows a proportional-odds model. MSM:
# E[Y^abar] = b0 + theta * sum(abar), theta = 1. Regression adjustment for L
# is biased (L is a collider on A -> L <- U -> Y); IPW is not.
sim_tv <- function(n = 500, n_t = 4, theta = 1) {
  U <- rnorm(n)
  long <- data.frame(id = rep(seq_len(n), each = n_t), time = rep(seq_len(n_t), n))
  long$l <- NA_real_; long$dose <- NA_integer_; long$dose.lag <- NA_integer_
  prev <- rep(0L, n)
  for (t in seq_len(n_t)) {
    rows <- long$time == t
    l <- U + 0.4 * prev + rnorm(n, 0, 0.5)
    dose <- findInterval(0.8 * l + 0.6 * prev + rlogis(n), c(0, 1.5))   # 0, 1, 2
    long$l[rows] <- l; long$dose.lag[rows] <- prev; long$dose[rows] <- dose
    prev <- dose
  }
  cum <- tapply(long$dose, long$id, sum)
  y <- theta * cum + 2 * U + rnorm(n)
  long$cum <- cum[long$id]; long$y <- y[long$id]
  long$dose <- factor(long$dose, levels = 0:2, ordered = TRUE)
  long
}
r5 <- reps(max(100, round(S / 2)), function() {
  long <- sim_tv()
  last <- long$time == max(long$time)
  wide <- long[last, ]
  lw <- reshape(long[, c("id", "time", "l")], idvar = "id", timevar = "time", direction = "wide")
  wide <- merge(wide, lw, by = "id")
  w_ord <- ipwtm(exposure = dose, family = "ordinal", link = "logit", numerator = ~ dose.lag,
                 denominator = ~ dose.lag + l, id = id, timevar = time, type = "all",
                 data = long)$ipw.weights[last]
  w_mult <- ipwtm(exposure = dose, family = "multinomial", numerator = ~ dose.lag,
                  denominator = ~ dose.lag + l, id = id, timevar = time, type = "all",
                  data = long)$ipw.weights[last]
  c(wlm(y ~ cum, wide), wlm(y ~ cum + l.1 + l.2 + l.3 + l.4, wide),
    wlm(y ~ cum, wide, w_ord), wlm(y ~ cum, wide, w_mult))
})
add("5. Time-varying ordinal dose (type all)", "unweighted", r5[, 1], 1)
add("5. Time-varying ordinal dose (type all)", "regression adjustment for L", r5[, 2], 1)
add("5. Time-varying ordinal dose (type all)", "ordinal weights (correct)", r5[, 3], 1)
add("5. Time-varying ordinal dose (type all)", "multinomial weights", r5[, 4], 1)
timer("5. time-varying ordinal done")

# ---- 6. Recurrent exposure events, family = "survival", type = "all" --------
# Exposure events from a known proportional hazard: lambda0 = 2 per unit time,
# log HR 0.5 for x (fixed) and -1.0 for xt (time-varying). The probability of
# the observed status per interval, and so the true weight, is known exactly.
sim_rec <- function(n = 300, irregular = FALSE) {
  rows <- lapply(seq_len(n), function(i) {
    cuts <- if (irregular) sort(c(runif(sample(4:12, 1)), 1)) else seq(0.1, 1, by = 0.1)
    cuts <- unique(cuts)
    data.frame(id = i, tstart = c(0, cuts[-length(cuts)]), fuptime = cuts)
  })
  d <- do.call(rbind, rows)
  d$x <- rep(rbinom(n, 1, 0.5), times = as.vector(table(d$id)))
  d$xt <- rbinom(nrow(d), 1, 0.4)
  p <- 1 - exp(-2 * (d$fuptime - d$tstart) * exp(0.5 * d$x - 1.0 * d$xt))
  d$expo <- rbinom(nrow(d), 1, p)
  d$w.true <- 1 / ave(ifelse(d$expo == 1, p, 1 - p), d$id, FUN = cumprod)
  d
}
for (irr in c(FALSE, TRUE)) {
  lab <- if (irr) "6b. Recurrent exposure, irregular intervals" else "6a. Recurrent exposure, regular intervals"
  r6 <- reps(max(100, round(S / 2)), function() {
    d <- sim_rec(irregular = irr)
    fit <- ipwtm(exposure = expo, family = "survival", denominator = ~ x + xt, id = id,
                 tstart = tstart, timevar = fuptime, type = "all", data = d)
    c(coef(fit$den.mod), median(abs(log(fit$ipw.weights / d$w.true))))
  })
  add(lab, "Cox coefficient of x", r6[, 1], 0.5)
  add(lab, "Cox coefficient of xt", r6[, 2], -1.0)
  add(lab, "median |log(estimated / true weight)|", r6[, 3], 0)
}
timer("6. recurrent exposure done")

# ---- 7. Bootstrap vs robust standard errors (JSS-style example) -------------
S7 <- max(50, round(S / 5))
has_survey <- requireNamespace("survey", quietly = TRUE)
r7 <- reps(S7, function() {
  d <- data.frame(l = rnorm(500, 10, 5))
  d$a <- rbinom(500, 1, plogis(d$l - 10))
  d$y <- 10 * d$a + 0.5 * d$l + rnorm(500, -10, 5)
  stat <- function(dd) {
    w <- ipwpoint(exposure = a, family = "binomial", link = "logit", numerator = ~ 1,
                  denominator = ~ l, data = dd)$ipw.weights
    coef(lm(y ~ a, data = dd, weights = w))["a"]
  }
  b <- suppressWarnings(ipwboot(d, stat, R = 200))
  rob <- c(NA, NA)
  if (has_survey) {
    d$w <- ipwpoint(exposure = a, family = "binomial", link = "logit", numerator = ~ 1,
                    denominator = ~ l, data = d)$ipw.weights
    m <- survey::svyglm(y ~ a, design = survey::svydesign(~ 1, weights = ~ w, data = d))
    rob <- c(sqrt(vcov(m)["a", "a"]), NA)
  }
  c(est = unname(b$t0), boot.se = unname(b$se), boot.lo = b$ci[1, 1], boot.hi = b$ci[1, 2],
    rob.se = rob[1])
})
emp_sd <- sd(r7[, "est"])
cov_boot <- mean(r7[, "boot.lo"] <= 10 & r7[, "boot.hi"] >= 10)
coverage[[1]] <- data.frame(method = "ipwboot, percentile CI", mean.se = mean(r7[, "boot.se"]),
                            emp.sd = emp_sd, coverage = cov_boot, reps = S7)
if (has_survey) {
  cov_rob <- mean(abs(r7[, "est"] - 10) <= 1.96 * r7[, "rob.se"])
  coverage[[2]] <- data.frame(method = "robust SE (svyglm), Wald CI", mean.se = mean(r7[, "rob.se"]),
                              emp.sd = emp_sd, coverage = cov_rob, reps = S7)
}
timer("7. bootstrap coverage done")

# ---- Report -----------------------------------------------------------------
res <- do.call(rbind, results)
cov_tab <- do.call(rbind, coverage)
num <- vapply(res, is.numeric, logical(1))
res[num & names(res) != "reps"] <- lapply(res[num & names(res) != "reps"], round, 3)
cov_tab[2:4] <- lapply(cov_tab[2:4], round, 3)

cat("\nSimulation study, S =", S, "replications (scenario 5-6: S/2, scenario 7:", S7, ")\n\n")
print(res, row.names = FALSE)
cat("\n7. Standard errors for the ATE (true effect 10), n = 500, R = 200 bootstrap samples\n")
print(cov_tab, row.names = FALSE)
cat("\nTrue values, scenario 1: ATE", round(truth1["ATE"], 3), " ATT", round(truth1["ATT"], 3),
    " ATO", round(truth1["ATO"], 3), "\n")

write.csv(res, file.path("dev", "simulation-results.csv"), row.names = FALSE)

md_table <- function(df) {
  c(paste("|", paste(names(df), collapse = " | "), "|"),
    paste("|", paste(rep("---", ncol(df)), collapse = " | "), "|"),
    apply(df, 1, function(r) paste("|", paste(trimws(format(r)), collapse = " | "), "|")))
}
summary_file <- Sys.getenv("GITHUB_STEP_SUMMARY")
if (nzchar(Sys.getenv("GITHUB_ACTIONS"))) {
  md <- c(sprintf("## Simulation study (S = %d)", S), "", md_table(res), "",
          "### Standard errors (scenario 7, true effect 10)", "", md_table(cov_tab), "",
          sprintf("True values, scenario 1: ATE %.3f, ATT %.3f, ATO %.3f",
                  truth1["ATE"], truth1["ATT"], truth1["ATO"]))
  if (nzchar(summary_file)) cat(md, file = summary_file, sep = "\n", append = TRUE)
  esc <- function(x) gsub("\n", "%0A", gsub("%", "%25", x, fixed = TRUE), fixed = TRUE)
  txt <- paste(c(capture.output(print(res, row.names = FALSE)), "",
                 capture.output(print(cov_tab, row.names = FALSE)), "",
                 sprintf("Truth scenario 1: ATE %.3f, ATT %.3f, ATO %.3f",
                         truth1["ATE"], truth1["ATT"], truth1["ATO"])), collapse = "\n")
  cat(sprintf("::notice title=simulation results (S = %d)::%s\n", S, esc(txt)))
}
timer("finished")
