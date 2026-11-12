# Why do ipwvcov() and the cluster bootstrap under-cover in the "strong"
# scenario of make-ipwvcov-tm-simulation.R? Test whether the cumulative
# weights are heavy-tailed: for n = 300, 1200 and 4800 people, the weight
# distribution (max, 99.9th percentile, Kish effective sample size / n), the
# Hill estimate of the tail index of the final cumulative weights (top 5%),
# and coverage of ipwvcov(), the cluster bootstrap and the robust cluster SE.
# Also a milder version of the strong scenario. A tail index below 2 means
# the weights have infinite variance.
# Needs: sandwich, parallel and the installed development ipw.
# Results: vignettes/articles/ipwvcov-tm-tail.rds

stopifnot(requireNamespace("sandwich"))
out_file <- file.path("vignettes", "articles", "ipwvcov-tm-tail.rds")

scenarios <- list(strong = c(b_l = 0.8, b_lag = 1.0, b_int = -0.5),
                  milder = c(b_l = 0.6, b_lag = 0.75, b_int = -0.4))
design <- data.frame(n = c(300, 1200, 4800), R = c(1000, 1000, 500), R_boot = c(200, 200, 100))
B <- 200
truth <- 1

gen <- function(n, sc) {
  k <- sample(3:7, n, replace = TRUE)
  id <- rep(seq_len(n), k)
  N <- length(id)
  l <- b <- blag <- numeric(N)
  # visit by visit, vectorised over people
  pos <- stats::ave(seq_len(N), id, FUN = seq_along)
  lp <- stats::rnorm(n); bp <- numeric(n)
  for (j in seq_len(max(k))) {
    rows <- which(pos == j)
    who <- id[rows]
    l[rows] <- 0.7 * lp[who] + stats::rnorm(length(rows), 0, 0.7)
    blag[rows] <- bp[who]
    b[rows] <- stats::rbinom(length(rows), 1, stats::plogis(sc[["b_int"]] + sc[["b_l"]] * l[rows] +
                                                             sc[["b_lag"]] * bp[who]))
    lp[who] <- l[rows]; bp[who] <- b[rows]
  }
  time <- stats::ave(stats::runif(N, 0.5, 1.5), id, FUN = cumsum)
  data.frame(id = id, time = time, l = l, b = b, blag = blag, y = b + l + stats::rnorm(N))
}

hill <- function(x, frac = 0.05) {
  x <- sort(x, decreasing = TRUE)
  k <- max(2, floor(frac * length(x)))
  1 / mean(log(x[seq_len(k)] / x[k + 1]))
}

one_rep <- function(r, sc, n, R_boot, B) {
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
  d <- gen(n, sc)
  f <- fit_one(d)
  w_final <- f$w$ipw.weights[!duplicated(d$id, fromLast = TRUE)]   # last row per person
  se_boot <- NA_real_
  if (r <= R_boot) {
    bt <- replicate(B, unname(stats::coef(fit_one(cluster_sample(d))$msm)["b"]))
    se_boot <- stats::sd(bt)
  }
  c(est = unname(stats::coef(f$msm)["b"]),
    se_ipwvcov = sqrt(ipw::ipwvcov(f$w, f$msm)["b", "b"]),
    se_robust = sqrt(sandwich::vcovCL(f$msm, cluster = d$id, type = "HC0",
                                      cadjust = FALSE)["b", "b"]),
    se_boot = se_boot,
    w_max = max(w_final), w_q999 = unname(stats::quantile(w_final, 0.999)),
    ess_frac = sum(w_final)^2 / sum(w_final^2) / n,
    hill = hill(w_final))
}

cl <- parallel::makeCluster(max(1, parallel::detectCores() - 2))
parallel::clusterExport(cl, c("gen", "hill"))
parallel::clusterSetRNGStream(cl, 2027)
rows <- list()
for (nm in names(scenarios)) {
  for (i in seq_len(nrow(design))) {
    n <- design$n[i]
    res <- do.call(rbind, parallel::parLapply(cl, seq_len(design$R[i]), one_rep,
                                              sc = scenarios[[nm]], n = n,
                                              R_boot = design$R_boot[i], B = B))
    cover <- function(se, keep = rep(TRUE, nrow(res))) {
      mean(abs(res[keep, "est"] - truth) < stats::qnorm(0.975) * se[keep])
    }
    bk <- !is.na(res[, "se_boot"])
    rows[[length(rows) + 1]] <- data.frame(
      scenario = nm, n = n,
      max.weight = stats::median(res[, "w_max"]), q999.weight = stats::median(res[, "w_q999"]),
      ess.fraction = mean(res[, "ess_frac"]), hill.tail.index = stats::median(res[, "hill"]),
      empirical.sd = stats::sd(res[, "est"]), mean.se.ipwvcov = mean(res[, "se_ipwvcov"]),
      mean.se.bootstrap = mean(res[bk, "se_boot"]),
      coverage.ipwvcov = cover(res[, "se_ipwvcov"]), coverage.bootstrap = cover(res[, "se_boot"], bk),
      coverage.robust = cover(res[, "se_robust"]),
      replications = nrow(res), bootstrap.replications = sum(bk))
    cat(nm, n, "done\n")
  }
}
parallel::stopCluster(cl)

# Tail index from one very large sample per scenario (less noisy than per replicate)
# ("moderate" is the moderate scenario of make-ipwvcov-tm-simulation.R; it is
# last so that the other two get the same draws as without it.)
tail_scenarios <- c(scenarios, list(moderate = c(b_l = 0.4, b_lag = 0.5, b_int = -0.3)))
RNGkind("L'Ecuyer-CMRG")   # as set by clusterSetRNGStream() above
set.seed(99)
hill_big <- vapply(tail_scenarios, function(sc) {
  d <- gen(2e5, sc)
  w <- ipw::ipwtm(exposure = b, family = "binomial", link = "logit", numerator = ~ blag,
                  denominator = ~ blag + l, id = id, timevar = time, type = "all", data = d)
  w_final <- w$ipw.weights[!duplicated(d$id, fromLast = TRUE)]
  c(hill_5pct = hill(w_final, 0.05), hill_1pct = hill(w_final, 0.01),
    hill_0.2pct = hill(w_final, 0.002), max.weight = max(w_final),
    ess.fraction = sum(w_final)^2 / sum(w_final^2) / length(w_final))
}, numeric(5))

saveRDS(list(summary = do.call(rbind, rows), hill_big = hill_big, scenarios = scenarios,
             B = B, ipw_version = as.character(utils::packageVersion("ipw")), date = Sys.Date()),
        out_file)
cat("Saved", out_file, "\n")
print(do.call(rbind, rows))
print(hill_big)
