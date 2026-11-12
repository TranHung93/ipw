# Simulation for ipwvcov() with time-varying weights from ipwtm()
# (vignettes/articles/ipw-and-weightit.Rmd). Long data with individual visit
# times (3 to 7 visits); a confounder l that changes over time; a treatment b
# that switches back and forth and depends on l and on the previous
# treatment; an outcome measured at every visit with a known effect of
# current treatment of 1. Two scenarios: "moderate" and "strong" dependence of
# treatment on l and on the previous treatment (the latter gives cumulative
# weights with a heavy right tail).
# Weights: ipwtm(type = "all", numerator ~ blag, denominator ~ blag + l).
# MSM: lm(y ~ b) on all rows. Coverage of 95% Wald intervals for the effect:
#   ipwvcov()  - M-estimation, scores summed within person;
#   robust     - cluster-robust (HC0, by person) SE that treats the weights
#                as known (sandwich::vcovCL);
#   bootstrap  - persons resampled, ipwtm() and the MSM re-fitted.
# Needs: sandwich, parallel and the installed development ipw.

stopifnot(requireNamespace("sandwich"))
out_file <- file.path("vignettes", "articles", "ipwvcov-tm-simulation.rds")

scenarios <- list(moderate = c(b_l = 0.4, b_lag = 0.5, b_int = -0.3),
                  strong   = c(b_l = 0.8, b_lag = 1.0, b_int = -0.5))
n <- 300
R_sim <- 1000
R_boot_sim <- 200
B <- 200
truth <- 1

one_rep <- function(r, sc, n, B, R_boot_sim) {
  gen <- function(n) {
    rows <- lapply(seq_len(n), function(i) {
      k <- sample(3:7, 1)
      t <- cumsum(stats::runif(k, 0.5, 1.5))
      l <- numeric(k); b <- numeric(k); blag <- numeric(k)
      l_prev <- stats::rnorm(1); b_prev <- 0
      for (j in seq_len(k)) {
        l[j] <- 0.7 * l_prev + stats::rnorm(1, 0, 0.7)
        blag[j] <- b_prev
        b[j] <- stats::rbinom(1, 1, stats::plogis(sc[["b_int"]] + sc[["b_l"]] * l[j] +
                                                   sc[["b_lag"]] * b_prev))
        l_prev <- l[j]; b_prev <- b[j]
      }
      data.frame(id = i, time = t, l = l, b = b, blag = blag, y = b + l + stats::rnorm(k))
    })
    do.call(rbind, rows)
  }
  fit_one <- function(d) {
    w <- ipw::ipwtm(exposure = b, family = "binomial", link = "logit", numerator = ~ blag,
                    denominator = ~ blag + l, id = id, timevar = time, type = "all", data = d)
    list(w = w, msm = stats::lm(y ~ b, data = d, weights = w$ipw.weights))
  }
  cluster_sample <- function(d) {
    rows <- split(seq_len(nrow(d)), d$id)
    pick <- sample(length(rows), replace = TRUE)
    out <- d[unlist(rows[pick]), ]
    out$id <- rep(seq_along(pick), lengths(rows[pick]))
    out
  }
  d <- gen(n)
  f <- fit_one(d)
  se_boot <- NA_real_
  if (r <= R_boot_sim) {
    bt <- replicate(B, unname(stats::coef(fit_one(cluster_sample(d))$msm)["b"]))
    se_boot <- stats::sd(bt)
  }
  c(est = unname(stats::coef(f$msm)["b"]),
    se_ipwvcov = sqrt(ipw::ipwvcov(f$w, f$msm)["b", "b"]),
    se_robust = sqrt(sandwich::vcovCL(f$msm, cluster = d$id, type = "HC0",
                                      cadjust = FALSE)["b", "b"]),
    se_boot = se_boot, max_w = max(f$w$ipw.weights))
}

cl <- parallel::makeCluster(max(1, parallel::detectCores() - 2))
parallel::clusterSetRNGStream(cl, 2026)
summary <- do.call(rbind, lapply(names(scenarios), function(nm) {
  res <- do.call(rbind, parallel::parLapply(cl, seq_len(R_sim), one_rep,
                                            sc = scenarios[[nm]], n = n, B = B,
                                            R_boot_sim = R_boot_sim))
  cover <- function(se, rows = seq_len(nrow(res))) {
    mean(abs(res[rows, "est"] - truth) < stats::qnorm(0.975) * se[rows])
  }
  boot_rows <- which(!is.na(res[, "se_boot"]))
  data.frame(
    scenario = nm, median.max.weight = stats::median(res[, "max_w"]),
    mean.estimate = mean(res[, "est"]), empirical.sd = stats::sd(res[, "est"]),
    mean.se.ipwvcov = mean(res[, "se_ipwvcov"]), mean.se.robust = mean(res[, "se_robust"]),
    mean.se.bootstrap = mean(res[boot_rows, "se_boot"]),
    coverage.ipwvcov = cover(res[, "se_ipwvcov"]), coverage.robust = cover(res[, "se_robust"]),
    coverage.bootstrap = cover(res[, "se_boot"], boot_rows),
    replications = nrow(res), bootstrap.replications = length(boot_rows))
}))
parallel::stopCluster(cl)

saveRDS(list(summary = summary, n = n, B = B, truth = truth, scenarios = scenarios,
             ipw_version = as.character(utils::packageVersion("ipw")), date = Sys.Date()),
        out_file)
cat("Saved", out_file, "\n")
print(summary)
