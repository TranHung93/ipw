# Simulation for vignettes/articles/survival-weights-bug.Rmd, section 6.
#
# Target: the 1-year risk of death if nobody started treatment. Follow-up is
# censored when a patient starts treatment; censoring weights from ipwtm()
# (family = "survival", type = "cens") correct for it, and a weighted
# Kaplan-Meier estimate gives the risk. A time-varying covariate l raises both
# the hazard of starting treatment and the hazard of death, so patients who
# stay untreated are healthier and the unweighted estimate is too low.
#
# Every patient has their own visit times; l changes at visits. Two scenarios:
# visits 3 to 20 days apart ("long rows") and 1 to 5 days apart ("short
# rows"). ipwtm() gives each row the weight at its end, which is a small
# approximation within a row; it shrinks as rows get shorter. Treatment start and death happen at exact times, so the Cox models
# that ipwtm() fits are correctly specified. The true risk is computed by
# simulating the same process with treatment start switched off.
#
# The same data sets are analysed with ipw 1.2.1.1 (in a separate R process
# with callr) and with the installed development version.
# Needs: callr, survival, and the installed development ipw.

stopifnot(requireNamespace("callr"))
suppressPackageStartupMessages(library(survival))
out_file <- file.path("vignettes", "articles", "survival-weights-bug-sim.rds")

gen <- function(n, with_treatment = TRUE, gap_range = c(3, 20)) {
  lam_c <- 0.003; lam_d <- 0.001   # baseline hazards per day: start of treatment, death
  b_c <- 0.8; b_d <- 0.8           # log hazard ratios for l
  tau <- 365                       # one year of follow-up
  t <- rep(0, n); l <- stats::rnorm(n); active <- rep(TRUE, n)
  res <- list(); k <- 0
  while (any(active)) {
    idx <- which(active); m <- length(idx)
    gap <- stats::runif(m, gap_range[1], gap_range[2])
    tc <- if (with_treatment) stats::rexp(m, lam_c * exp(b_c * l[idx])) else rep(Inf, m)
    td <- stats::rexp(m, lam_d * exp(b_d * l[idx]))
    first <- pmin(t[idx] + gap, t[idx] + tc, t[idx] + td, tau)
    st <- as.numeric(first == t[idx] + tc)
    de <- as.numeric(first == t[idx] + td & st == 0)
    # times to 0.001 day; rows at least 0.001 day long (survival rejects rows
    # of effectively zero length)
    end <- pmin(pmax(round(first, 3), t[idx] + 0.001), tau)
    k <- k + 1
    res[[k]] <- data.frame(id = idx, tstart = t[idx], time = end, l = l[idx],
                           start_trt = st, death = de)
    stop_now <- st == 1 | de == 1 | end >= tau
    t[idx] <- end
    l[idx] <- 0.8 * l[idx] + stats::rnorm(m, 0, 0.6)
    active[idx[stop_now]] <- FALSE
  }
  d <- do.call(rbind, res)
  d <- d[order(d$id, d$tstart), ]
  rownames(d) <- NULL
  d
}

scenarios <- list(long = c(3, 20), short = c(1, 5))
R_SIM <- 500
N_SIM <- 1000
seeds <- 1000 + seq_len(R_SIM)

sim_fun <- function(seeds, n, gen, gap_range) {
  suppressPackageStartupMessages({ library(ipw); library(survival) })
  risk_at <- function(fit) 1 - summary(fit, times = 365, extend = TRUE)$surv
  out <- t(vapply(seeds, function(s) {
    set.seed(s)
    d <- gen(n, gap_range = gap_range)
    naive <- risk_at(survfit(Surv(tstart, time, death) ~ 1, data = d))
    w <- tryCatch(
      ipwtm(exposure = start_trt, family = "survival", denominator = ~ l, id = id,
            tstart = tstart, timevar = time, type = "cens", data = d)$ipw.weights,
      error = function(e) rep(NA_real_, nrow(d)))
    weighted <- if (anyNA(w)) NA_real_ else
      risk_at(survfit(Surv(tstart, time, death) ~ 1, data = d, weights = w))
    c(naive = naive, weighted = weighted, max_w = max(w), mean_w = mean(w),
      rows = nrow(d), starts = sum(d$start_trt), deaths = sum(d$death))
  }, numeric(7)))
  list(version = as.character(utils::packageVersion("ipw")), results = out)
}

old_lib <- file.path(tempdir(), "ipw-1.2.1.1")
dir.create(old_lib, showWarnings = FALSE)
utils::install.packages("https://cloud.r-project.org/src/contrib/Archive/ipw/ipw_1.2.1.1.tar.gz",
                        repos = NULL, type = "source", lib = old_lib, quiet = TRUE)

results <- lapply(names(scenarios), function(sc) {
  gr <- scenarios[[sc]]
  set.seed(1)
  d0 <- gen(400000, with_treatment = FALSE, gap_range = gr)
  truth <- mean(tapply(d0$death, d0$id, max))
  rm(d0)
  cat("Scenario", sc, ": true 1-year risk without treatment", round(truth, 4), "\n")
  old <- callr::r(sim_fun, args = list(seeds, N_SIM, gen, gr), libpath = c(old_lib, .libPaths()))
  dev <- callr::r(sim_fun, args = list(seeds, N_SIM, gen, gr))
  stopifnot(old$version == "1.2.1.1")
  list(gap_range = gr, truth = truth, truth_mcse = sqrt(truth * (1 - truth) / 400000),
       old = old, dev = dev)
})
names(results) <- names(scenarios)

# One data set (long rows) with its old weights, for the plot in the article
set.seed(seeds[1])
example_data <- gen(N_SIM, gap_range = scenarios$long)
example_old_w <- callr::r(function(d) {
  suppressPackageStartupMessages(library(ipw))
  ipwtm(exposure = start_trt, family = "survival", denominator = ~ l, id = id,
        tstart = tstart, timevar = time, type = "cens", data = d)$ipw.weights
}, args = list(example_data), libpath = c(old_lib, .libPaths()))

saveRDS(list(scenarios = results, R = R_SIM, n = N_SIM, seeds = seeds, gen = gen,
             example = list(data = example_data, old_w = example_old_w, scenario = "long"),
             made = format(Sys.Date()), R_version = R.version.string),
        out_file, compress = "xz")

summ <- function(x) sprintf("mean %.4f | bias %+.4f (MCSE %.4f, bias/MCSE %+.1f) | SD %.4f | NA %d",
  mean(x, na.rm = TRUE), mean(x, na.rm = TRUE) - truth, sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))),
  (mean(x, na.rm = TRUE) - truth) / (sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))), sd(x, na.rm = TRUE), sum(is.na(x)))
for (sc in names(results)) {
  r <- results[[sc]]; truth <- r$truth
  cat("\n== Scenario", sc, "(visits", r$gap_range[1], "to", r$gap_range[2], "days apart ) truth", round(truth, 4), "\n")
  cat("No weights      :", summ(r$dev$results[, "naive"]), "\n")
  cat("ipw 1.2.1.1     :", summ(r$old$results[, "weighted"]), "\n")
  cat("this version    :", summ(r$dev$results[, "weighted"]), " (", r$dev$version, ")\n")
  cat("mean of max weight: old", round(mean(r$old$results[, "max_w"]), 3),
      " new", round(mean(r$dev$results[, "max_w"]), 3), "\n")
  cat("rows per data set:", round(mean(r$dev$results[, "rows"])), " treatment starts:",
      round(mean(r$dev$results[, "starts"])), " deaths:", round(mean(r$dev$results[, "deaths"])), "\n")
}
cat("Saved", out_file, "(", file.size(out_file), "bytes )\n")
