# Regression cases shared by tests/testthat/test-regression.R and
# data-raw/make-regression-fixtures.R.
#
# Each case is a function taking `ipwpoint` and `ipwtm` (from whichever ipw
# version is being run) and returning a named list of weight vectors. The same
# inputs are used for the old (CRAN 1.2.1.1) and the current implementation.

regression_cases <- function() {
  env <- new.env()
  utils::data("haartdat", package = "ipw", envir = env)
  haartdat <- env$haartdat

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

  list(
    point_binomial = function(ipwpoint, ipwtm) {
      ipwpoint(exposure = a, family = "binomial", link = "logit",
               numerator = ~ 1, denominator = ~ l, data = simdat, trunc = 0.05)
    },
    point_gaussian = function(ipwpoint, ipwtm) {
      ipwpoint(exposure = ac, family = "gaussian",
               numerator = ~ 1, denominator = ~ l, data = simdat)
    },
    point_ordinal = function(ipwpoint, ipwtm) {
      ipwpoint(exposure = am, family = "ordinal", link = "logit",
               numerator = ~ 1, denominator = ~ l, data = simdat)
    },
    tm_survival_first = function(ipwpoint, ipwtm) {
      ipwtm(exposure = haartind, family = "survival",
            numerator = ~ sex + age, denominator = ~ sex + age + cd4.sqrt,
            id = patient, tstart = tstart, timevar = fuptime, type = "first",
            data = haartdat, trunc = 0.01)
    },
    tm_survival_cens = function(ipwpoint, ipwtm) {
      ipwtm(exposure = dropout, family = "survival",
            numerator = ~ sex + age, denominator = ~ sex + age + cd4.sqrt,
            id = patient, tstart = tstart, timevar = fuptime, type = "cens",
            data = haartdat)
    },
    tm_binomial_all = function(ipwpoint, ipwtm) {
      ipwtm(exposure = a, family = "binomial", link = "logit",
            numerator = ~ 1, denominator = ~ l, id = id, timevar = time,
            type = "all", data = long)
    }
  )
}

# Keep only the numeric outputs that should not change between versions
regression_extract <- function(res) {
  out <- list(ipw.weights = unname(as.numeric(res$ipw.weights)))
  if (!is.null(res$weights.trunc)) out$weights.trunc <- unname(as.numeric(res$weights.trunc))
  if (!is.null(res$selvar)) out$selvar <- unname(as.numeric(res$selvar))
  out
}
