# Simulation for ipwvcov() with a coxph() marginal structural model
# (vignettes/articles/ipw-and-weightit.Rmd). Coverage of 95% Wald intervals
# for the marginal log hazard ratio with three standard errors:
#   ipwvcov()  - M-estimation, accounts for the estimated weights;
#   robust     - coxph's robust SE, treats the weights as known;
#   bootstrap  - weights and MSM re-estimated in each bootstrap sample.
# Rerun when the design changes. Needs: survival and the installed
# development ipw. Takes about 20-30 minutes.

library(ipw)
library(survival)
out_file <- file.path("vignettes", "articles", "ipwvcov-cox-simulation.rds")

# Treatment depends on l; the treated potential time is the untreated one
# stretched by exp(0.5), i.e. a conditional hazard ratio of exp(-0.5).
# Administrative censoring at 4.
gen <- function(n) {
  l <- stats::rnorm(n)
  ps <- stats::plogis(0.3 + 0.8 * l)
  a <- stats::rbinom(n, 1, ps)
  t0 <- stats::rexp(n, 0.2 * exp(0.6 * l))
  t1 <- t0 * exp(0.5)
  tt <- ifelse(a == 1, t1, t0)
  data.frame(l, a, ps, time = pmin(tt, 4), status = as.numeric(tt < 4), t0, t1)
}

# Marginal log hazard ratios: Cox model on both potential outcomes of a
# very large sample, weighted to the target population (ATE: everyone,
# ATO: weights e(1 - e)).
set.seed(1)
big <- gen(4e5)
pot <- data.frame(time = pmin(c(big$t0, big$t1), 4), status = as.numeric(c(big$t0, big$t1) < 4),
                  a = rep(0:1, each = nrow(big)), ps = rep(big$ps, 2))
truth <- c(
  ATE = unname(coef(coxph(Surv(time, status) ~ a, data = pot, ties = "breslow"))),
  ATO = unname(coef(coxph(Surv(time, status) ~ a, data = pot, ties = "breslow",
                          weights = ps * (1 - ps)))))

fit_one <- function(d, est) {
  w <- ipwpoint(a, "binomial", "logit", denominator = ~ l, data = d, estimand = est)
  msm <- coxph(Surv(time, status) ~ a, data = d, weights = w$ipw.weights, ties = "breslow")
  list(w = w, msm = msm)
}

n <- 500
R_sim <- 1000      # replications for ipwvcov and robust
R_boot_sim <- 300  # replications that also get a bootstrap
B <- 200           # bootstrap samples per replication
set.seed(2)
res <- lapply(seq_len(R_sim), function(r) {
  d <- gen(n)
  out <- lapply(c("ATE", "ATO"), function(est) {
    f <- fit_one(d, est)
    se_boot <- NA_real_
    if (r <= R_boot_sim) {
      bt <- replicate(B, {
        db <- d[sample(n, replace = TRUE), ]
        unname(coef(fit_one(db, est)$msm))
      })
      se_boot <- stats::sd(bt)
    }
    c(est = unname(coef(f$msm)), se_ipwvcov = sqrt(ipwvcov(f$w, f$msm)[1, 1]),
      se_robust = sqrt(f$msm$var[1, 1]), se_boot = se_boot)
  })
  names(out) <- c("ATE", "ATO")
  out
})

summarise <- function(est) {
  m <- do.call(rbind, lapply(res, `[[`, est))
  cover <- function(se, rows = seq_len(nrow(m))) {
    mean(abs(m[rows, "est"] - truth[[est]]) < stats::qnorm(0.975) * se[rows])
  }
  boot_rows <- which(!is.na(m[, "se_boot"]))
  data.frame(
    estimand = est, truth = truth[[est]], mean.estimate = mean(m[, "est"]),
    empirical.sd = stats::sd(m[, "est"]),
    mean.se.ipwvcov = mean(m[, "se_ipwvcov"]), mean.se.robust = mean(m[, "se_robust"]),
    mean.se.bootstrap = mean(m[boot_rows, "se_boot"]),
    coverage.ipwvcov = cover(m[, "se_ipwvcov"]), coverage.robust = cover(m[, "se_robust"]),
    coverage.bootstrap = cover(m[, "se_boot"], boot_rows),
    replications = nrow(m), bootstrap.replications = length(boot_rows))
}
saveRDS(list(summary = rbind(summarise("ATE"), summarise("ATO")),
             n = n, B = B, ipw_version = as.character(utils::packageVersion("ipw")),
             date = Sys.Date()),
        out_file)
cat("Saved", out_file, "\n")
print(readRDS(out_file)$summary)
