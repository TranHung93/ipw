# Regression cases shared by tests/testthat/test-regression.R and
# data-raw/make-regression-fixtures.R.
#
# Each case is a function taking `ipwpoint`, `ipwtm` and `tstartfun` (from
# whichever ipw version is being run) and returning either an ipwpoint/ipwtm
# result or a tstart vector. The same inputs are used for the old (CRAN
# 1.2.1.1) and the current implementation. Only settings whose results should
# not change between the versions belong here; intended changes (listed in
# NEWS.md) are tested against by-hand references in test-regression-changes.R.

regression_data <- function() {
  env <- new.env()
  utils::data("haartdat", "healthdat", "basdat", "timedat", package = "ipw", envir = env)

  set.seed(16)
  n <- 1000
  simdat <- data.frame(l = stats::rnorm(n, 10, 5))
  a.lin <- simdat$l - 10
  pa <- exp(a.lin) / (1 + exp(a.lin))
  simdat$a <- stats::rbinom(n, 1, prob = pa)
  simdat$ac <- 0.3 * simdat$l + stats::rnorm(n)
  simdat$am <- cut(simdat$l + stats::rnorm(n, 0, 3), 3, labels = c("1", "2", "3"))

  set.seed(17)
  n_id <- 200
  long <- data.frame(id = rep(seq_len(n_id), each = 4), time = rep(1:4, n_id))
  long$l <- stats::rnorm(nrow(long))
  long$a <- stats::rbinom(nrow(long), 1, stats::plogis(-0.5 + long$l))

  # Wider set (added in the second review)
  set.seed(18)
  long$x <- rep(stats::rbinom(n_id, 1, 0.5), each = 4)
  long$m <- pmin(long$a + stats::rbinom(nrow(long), 1, stats::plogis(long$l)), 2)
  long$g <- 0.5 * long$l + 0.3 * long$x + stats::rnorm(nrow(long))
  long$c <- stats::rbinom(nrow(long), 1, stats::plogis(-2 + 0.5 * long$l))
  long$tstart <- long$time - 1
  shuffle <- sample(nrow(long))

  health <- env$healthdat
  health$hi <- as.numeric(health$health > stats::median(health$health))

  # Example 2 of ?ipwtm (TB and mortality); tstart computed without tstartfun.
  # Needs nlme (in Suggests); the case is left out when nlme is not installed.
  startstop <- NULL
  if (requireNamespace("nlme", quietly = TRUE)) {
    basdat <- env$basdat
    timedat <- env$timedat
    timedat$cd4.sqrt <- sqrt(timedat$cd4count)
    timedat <- merge(timedat, basdat[, c("id", "Ttb")], by = "id", all.x = TRUE)
    timedat$tb.lag <- ifelse(with(timedat, !is.na(Ttb) & fuptime > Ttb), 1, 0)
    cd4.lme <- nlme::lme(cd4.sqrt ~ fuptime + tb.lag, random = ~ fuptime | id, data = timedat)
    times <- sort(unique(c(basdat$Ttb, basdat$Tend)))
    startstop <- data.frame(id = rep(basdat$id, each = length(times)),
                            fuptime = rep(times, nrow(basdat)))
    startstop <- merge(startstop, basdat, by = "id", all.x = TRUE)
    startstop <- startstop[with(startstop, fuptime <= Tend), ]
    startstop <- startstop[order(startstop$id, startstop$fuptime), ]
    startstop$tstart <- stats::ave(startstop$fuptime, startstop$id,
                                   FUN = function(t) c(-1, utils::head(t, -1)))
    startstop$tb <- ifelse(with(startstop, !is.na(Ttb) & fuptime >= Ttb), 1, 0)
    startstop$tb.lag <- ifelse(with(startstop, !is.na(Ttb) & fuptime > Ttb), 1, 0)
    startstop$cd4.sqrt <- stats::predict(cd4.lme, newdata = data.frame(
      id = startstop$id, fuptime = startstop$fuptime, tb.lag = startstop$tb.lag))
    rownames(startstop) <- NULL
  }

  list(haartdat = env$haartdat, simdat = simdat, long = long, shuffle = shuffle,
       health = health, startstop = startstop)
}

regression_cases <- function() {
  D <- regression_data()
  haartdat <- D$haartdat; simdat <- D$simdat; long <- D$long; health <- D$health
  startstop <- D$startstop
  long_shuf <- long[D$shuffle, ]
  haart_shuf <- haartdat[rev(seq_len(nrow(haartdat))), ]

  cases <- list(
    # ---- original six cases -------------------------------------------------
    point_binomial = function(ipwpoint, ipwtm, tstartfun) {
      ipwpoint(exposure = a, family = "binomial", link = "logit",
               numerator = ~ 1, denominator = ~ l, data = simdat, trunc = 0.05)
    },
    point_gaussian = function(ipwpoint, ipwtm, tstartfun) {
      ipwpoint(exposure = ac, family = "gaussian",
               numerator = ~ 1, denominator = ~ l, data = simdat)
    },
    point_ordinal = function(ipwpoint, ipwtm, tstartfun) {
      ipwpoint(exposure = am, family = "ordinal", link = "logit",
               numerator = ~ 1, denominator = ~ l, data = simdat)
    },
    tm_survival_first = function(ipwpoint, ipwtm, tstartfun) {
      ipwtm(exposure = haartind, family = "survival",
            numerator = ~ sex + age, denominator = ~ sex + age + cd4.sqrt,
            id = patient, tstart = tstart, timevar = fuptime, type = "first",
            data = haartdat, trunc = 0.01)
    },
    tm_survival_cens = function(ipwpoint, ipwtm, tstartfun) {
      ipwtm(exposure = dropout, family = "survival",
            numerator = ~ sex + age, denominator = ~ sex + age + cd4.sqrt,
            id = patient, tstart = tstart, timevar = fuptime, type = "cens",
            data = haartdat)
    },
    tm_binomial_all = function(ipwpoint, ipwtm, tstartfun) {
      ipwtm(exposure = a, family = "binomial", link = "logit",
            numerator = ~ 1, denominator = ~ l, id = id, timevar = time,
            type = "all", data = long)
    },

    # ---- ipwpoint, wider ----------------------------------------------------
    point_binomial_unstab = function(ipwpoint, ipwtm, tstartfun) {
      ipwpoint(exposure = a, family = "binomial", link = "logit",
               denominator = ~ l, data = simdat)
    },
    point_binomial_shuffled = function(ipwpoint, ipwtm, tstartfun) {
      ipwpoint(exposure = a, family = "binomial", link = "logit",
               numerator = ~ 1, denominator = ~ l, data = simdat[rev(seq_len(nrow(simdat))), ])
    },
    point_multinomial = function(ipwpoint, ipwtm, tstartfun) {
      ipwpoint(exposure = am, family = "multinomial",
               numerator = ~ 1, denominator = ~ l, data = simdat)
    },
    point_multinomial_trunc = function(ipwpoint, ipwtm, tstartfun) {
      ipwpoint(exposure = am, family = "multinomial", numerator = ~ 1,
               denominator = ~ l, data = simdat, trunc = 0.02)
    },
    point_gaussian_health = function(ipwpoint, ipwtm, tstartfun) {
      ipwpoint(exposure = income, family = "gaussian",
               numerator = ~ 1, denominator = ~ iq, data = health, trunc = 0.01)
    },
    point_binomial_health = function(ipwpoint, ipwtm, tstartfun) {
      ipwpoint(exposure = hi, family = "binomial", link = "logit",
               numerator = ~ 1, denominator = ~ iq, data = health)
    },

    # ---- ipwtm survival, wider ----------------------------------------------
    tm_survival_first_unstab = function(ipwpoint, ipwtm, tstartfun) {
      ipwtm(exposure = haartind, family = "survival",
            denominator = ~ sex + age + cd4.sqrt,
            id = patient, tstart = tstart, timevar = fuptime, type = "first",
            data = haartdat)
    },
    tm_survival_cens_trunc = function(ipwpoint, ipwtm, tstartfun) {
      ipwtm(exposure = dropout, family = "survival",
            numerator = ~ sex + age, denominator = ~ sex + age + cd4.sqrt,
            id = patient, tstart = tstart, timevar = fuptime, type = "cens",
            data = haartdat, trunc = 0.02)
    },
    tm_survival_first_unsorted = function(ipwpoint, ipwtm, tstartfun) {
      ipwtm(exposure = haartind, family = "survival",
            numerator = ~ sex + age, denominator = ~ sex + age + cd4.sqrt,
            id = patient, tstart = tstart, timevar = fuptime, type = "first",
            data = haart_shuf)
    },
    tm_survival_common_grid_first = function(ipwpoint, ipwtm, tstartfun) {
      ipwtm(exposure = a, family = "survival", numerator = ~ x,
            denominator = ~ x + l, id = id, tstart = tstart, timevar = time,
            type = "first", data = long)
    },
    tm_survival_common_grid_cens = function(ipwpoint, ipwtm, tstartfun) {
      ipwtm(exposure = c, family = "survival", numerator = ~ x,
            denominator = ~ x + l, id = id, tstart = tstart, timevar = time,
            type = "cens", data = long_shuf)
    },

    # ---- ipwtm binomial, multinomial, ordinal, gaussian ---------------------
    tm_binomial_all_unsorted = function(ipwpoint, ipwtm, tstartfun) {
      ipwtm(exposure = a, family = "binomial", link = "logit",
            numerator = ~ x, denominator = ~ x + l, id = id, timevar = time,
            type = "all", data = long_shuf, trunc = 0.05)
    },
    tm_multinomial_first = function(ipwpoint, ipwtm, tstartfun) {
      ipwtm(exposure = m, family = "multinomial", numerator = ~ x,
            denominator = ~ x + l, id = id, timevar = time, type = "first", data = long)
    },
    tm_multinomial_first_unstab = function(ipwpoint, ipwtm, tstartfun) {
      ipwtm(exposure = m, family = "multinomial",
            denominator = ~ x + l, id = id, timevar = time, type = "first", data = long_shuf)
    },
    tm_gaussian_ar1 = function(ipwpoint, ipwtm, tstartfun) {
      ipwtm(exposure = g, family = "gaussian", numerator = ~ x,
            denominator = ~ x + l, id = id, timevar = time, type = "all",
            corstr = "ar1", data = long)
    },
    tm_gaussian_exch_unsorted = function(ipwpoint, ipwtm, tstartfun) {
      ipwtm(exposure = g, family = "gaussian", numerator = ~ x,
            denominator = ~ x + l, id = id, timevar = time, type = "all",
            corstr = "exchangeable", data = long_shuf, trunc = 0.01)
    },

    # ---- tstartfun ------------------------------------------------------------
    tstartfun_haart = function(ipwpoint, ipwtm, tstartfun) {
      tstartfun(patient, fuptime, haartdat)
    },
    tstartfun_unsorted = function(ipwpoint, ipwtm, tstartfun) {
      tstartfun(patient, fuptime, haart_shuf)
    },
    tstartfun_irregular_unsorted = function(ipwpoint, ipwtm, tstartfun) {
      d <- long_shuf
      d$t2 <- d$time + d$l^2   # own visit times per subject
      tstartfun(id, t2, d)
    }
  )

  if (!is.null(startstop)) {
    cases$tm_survival_tb_example <- function(ipwpoint, ipwtm, tstartfun) {
      ipwtm(exposure = tb, family = "survival", numerator = ~ 1,
            denominator = ~ cd4.sqrt, id = id, tstart = tstart, timevar = fuptime,
            type = "first", data = startstop)
    }
  }

  # Every link of the binomial and ordinal families. ipw 1.2.1.1 reads 'link'
  # and 'type' unevaluated from the call, so the values are inlined (bquote).
  here <- environment()
  for (lk in c("probit", "cauchit", "cloglog")) {
    cases[[paste0("point_binomial_", lk)]] <- eval(bquote(function(ipwpoint, ipwtm, tstartfun) {
      ipwpoint(exposure = a, family = "binomial", link = .(lk),
               numerator = ~ 1, denominator = ~ l, data = simdat)
    }), here)
    cases[[paste0("point_ordinal_", lk)]] <- eval(bquote(function(ipwpoint, ipwtm, tstartfun) {
      ipwpoint(exposure = am, family = "ordinal", link = .(lk),
               numerator = ~ 1, denominator = ~ l, data = simdat)
    }), here)
  }
  for (tp in c("first", "cens", "all")) for (lk in c("logit", "probit", "cloglog")) {
    cases[[paste("tm_binomial", tp, lk, sep = "_")]] <- eval(bquote(function(ipwpoint, ipwtm, tstartfun) {
      ipwtm(exposure = a, family = "binomial", link = .(lk), numerator = ~ x,
            denominator = ~ x + l, id = id, timevar = time, type = .(tp), data = long)
    }), here)
  }
  for (lk in c("logit", "probit", "cloglog")) {
    cases[[paste0("tm_ordinal_first_", lk)]] <- eval(bquote(function(ipwpoint, ipwtm, tstartfun) {
      ipwtm(exposure = m, family = "ordinal", link = .(lk), numerator = ~ x,
            denominator = ~ x + l, id = id, timevar = time, type = "first", data = long)
    }), here)
  }
  cases
}

# Coefficients of a fitted model (polr: slopes then intercepts)
regression_coef <- function(m) {
  if (is.null(m)) return(NULL)
  if (inherits(m, "polr")) return(unname(c(m$coefficients, m$zeta)))
  unname(as.vector(stats::coef(m)))
}

# Keep only the numeric outputs that should not change between versions
regression_extract <- function(res) {
  if (is.numeric(res)) return(list(tstart = unname(as.numeric(res))))
  out <- list(ipw.weights = unname(as.numeric(res$ipw.weights)))
  if (!is.null(res$weights.trunc)) out$weights.trunc <- unname(as.numeric(res$weights.trunc))
  if (!is.null(res$selvar)) out$selvar <- unname(as.numeric(res$selvar))
  if (!is.null(res$num.mod)) out$num.coef <- regression_coef(res$num.mod)
  if (!is.null(res$den.mod)) out$den.coef <- regression_coef(res$den.mod)
  out
}

# Tolerance used for each element: weights 1e-8, model coefficients 1e-6
regression_tolerance <- function(element, case = "") {
  if (element %in% c("num.coef", "den.coef")) return(1e-6)
  # Ordinal models are fitted by MASS::polr(), an iterative optimiser whose end
  # point differs in about the 8th digit between platforms (macOS arm64 vs the
  # Windows machine that made the fixture). A real regression shows up as a
  # difference of order 1, so 1e-6 still catches it.
  if (grepl("ordinal", case) && element %in% c("ipw.weights", "weights.trunc")) return(1e-6)
  1e-8
}
