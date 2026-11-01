#' Estimate Inverse Probability Weights (Point Treatment)
#'
#' @description
#' Estimates inverse probability weights to fit marginal structural models in a
#' point treatment situation. Exposure can be binomial, multinomial, ordinal,
#' a count (Poisson or negative binomial) or continuous (gaussian, gamma or
#' log-normal). Both stabilized and unstabilized weights can be estimated, and
#' for a binomial exposure also weights for the effect in the treated (ATT) and
#' overlap weights (ATO). Formulas may contain smooth terms such as \code{s(x)},
#' which are then fitted as generalized additive models with \pkg{mgcv}.
#'
#' @details
#' Standard errors: \code{ipwpoint} only estimates the weights. Standard errors
#' for the marginal structural model have to be obtained when that model is fitted.
#' Robust (sandwich) standard errors, for example from \code{svyglm()} in the survey
#' package as in Van der Wal and Geskus (2011), treat the weights as known and tend to
#' be conservative (too large). Austin (2016) compared naive, robust and bootstrap
#' standard errors in simulations and recommends the bootstrap, in which the weights
#' are re-estimated in each bootstrap sample. This takes the uncertainty in the
#' estimated weights into account; see the bootstrap example below, or use
#' \code{\link{ipwboot}}, which does this for you. About 200 bootstrap samples are
#' usually enough to estimate a standard error.
#'
#' @param exposure Name of the exposure column in \code{data} (unquoted or as a string).
#'   Binomial exposure variables should be coded as 0/1.
#' @param family Model family for exposure: \code{"binomial"}, \code{"multinomial"},
#'   \code{"ordinal"}, \code{"gaussian"}, \code{"poisson"}, \code{"negbin"}
#'   (negative binomial, fitted with \code{MASS::glm.nb}), \code{"gamma"}
#'   (gamma GLM with log link) or \code{"lognormal"} (normal model for the log
#'   of the exposure). Count exposures must be non-negative integers; gamma and
#'   log-normal exposures must be positive.
#' @param link Link function for binomial (\code{"logit"}, \code{"probit"}, \code{"cauchit"},
#'   \code{"log"}, \code{"cloglog"}) or ordinal (\code{"logit"}, \code{"probit"},
#'   \code{"cauchit"}, \code{"cloglog"}). Ignored for multinomial and gaussian.
#' @param numerator Formula specifying the right-hand side of the numerator model.
#'   If \code{NULL}, unstabilized weights (numerator = 1) are calculated. Required
#'   for the continuous families (\code{"gaussian"}, \code{"gamma"},
#'   \code{"lognormal"}); must be \code{NULL} when \code{estimand} is
#'   \code{"ATT"} or \code{"ATO"}.
#' @param denominator Formula specifying the right-hand side of the denominator model.
#'   Smooth terms (\code{s()}, \code{te()}, \code{ti()}, \code{t2()}) are allowed
#'   for the binomial, gaussian and poisson families; the model is then fitted with
#'   \code{mgcv::gam}.
#' @param data Data frame containing exposure and model covariates.
#' @param trunc Optional truncation percentile between 0 and 0.5 (e.g., 0.01 for 1st/99th percentiles).
#' @param estimand Target estimand for a binomial exposure: \code{"ATE"} (default;
#'   inverse probability weights), \code{"ATT"} (treated get weight 1, untreated
#'   \eqn{p/(1-p)}) or \code{"ATO"} (overlap weights: treated \eqn{1-p}, untreated
#'   \eqn{p}), where \eqn{p} is the fitted probability of exposure from the
#'   denominator model.
#' @param ... Additional arguments passed to model fitting functions (\code{glm},
#'   \code{multinom}, \code{polr}, \code{glm.nb}, \code{gam}).
#'
#' @return A list of class \code{"ipwpoint"} containing:
#' \item{ipw.weights}{Numeric vector of inverse probability weights.}
#' \item{weights.trunc}{Numeric vector of truncated weights (if \code{trunc} is specified).}
#' \item{call}{The original function call.}
#' \item{num.mod}{Numerator model object (if \code{numerator} is specified).}
#' \item{den.mod}{Denominator model object.}
#' @author Willem M. van der Wal \email{willem@vanderwalresearch.com}, Ronald B. Geskus \email{rgeskus@oucru.org},
#'   Hung Thai Tran \email{hungtt@oucru.org}
#' @references
#' Cole, S.R. & Hernán, M.A. (2008). Constructing inverse probability weights for marginal structural models. \emph{American Journal of Epidemiology}, \bold{168}(6), 656-664.
#'
#' Robins, J.M., Hernán, M.A. & Brumback, B.A. (2000). Marginal structural models and causal inference in epidemiology. \emph{Epidemiology}, \bold{11}, 550-560.
#'
#' Van der Wal W.M. & Geskus R.B. (2011). ipw: An R Package for Inverse Probability Weighting. \emph{Journal of Statistical Software}, \bold{43}(13), 1-23. \doi{10.18637/jss.v043.i13}.
#'
#' Austin, P.C. (2016). Variance estimation when using inverse probability of treatment weighting (IPTW) with survival analysis. \emph{Statistics in Medicine}, \bold{35}, 5642-5655. \doi{10.1002/sim.7084}.
#' @export
#' @seealso \code{\link{basdat}}, \code{\link{haartdat}}, \code{\link{ipwplot}}, \code{\link{ipwpoint}}, \code{\link{ipwtm}}, \code{\link{timedat}}, \code{\link{tstartfun}}.
#'
#' @examples
#' # Simulate data with continuous confounder and outcome, binomial exposure.
#' # Marginal causal effect of exposure on outcome: 10.
#' set.seed(16)
#' n <- 1000
#' simdat <- data.frame(l = rnorm(n, 10, 5))
#' a.lin <- simdat$l - 10
#' pa <- exp(a.lin)/(1 + exp(a.lin))
#' simdat$a <- rbinom(n, 1, prob = pa)
#' simdat$y <- 10*simdat$a + 0.5*simdat$l + rnorm(n, -10, 5)
#' simdat[1:5,]
#'
#' # Estimate ipw weights
#' temp <- ipwpoint(
#'   exposure = a,
#'   family = "binomial",
#'   link = "logit",
#'   numerator = ~ 1,
#'   denominator = ~ l,
#'   data = simdat)
#' summary(temp$ipw.weights)
#'
#' # Plot inverse probability weights
#' ipwplot(weights = temp$ipw.weights, logscale = FALSE,
#'         main = "Stabilized weights", xlim = c(0, 8))
#'
#' # Examine numerator and denominator models
#' summary(temp$num.mod)
#' summary(temp$den.mod)
#'
#' # Paste inverse probability weights
#' simdat$sw <- temp$ipw.weights
#'
#' # Marginal structural model for the causal effect of a on y
#' # corrected for confounding by l using inverse probability weighting
#' if (requireNamespace("survey", quietly = TRUE)) {
#'   library(survey)
#'   msm <- svyglm(y ~ a,
#'                 design = svydesign(~1, weights = ~temp$ipw.weights,
#'                 data = simdat))
#'   summary(msm)
#' }
#'
#' # Bootstrap standard error and percentile confidence interval for the MSM.
#' # The robust standard error from svyglm() treats the weights as known and
#' # tends to be conservative. Austin (2016) recommends the bootstrap, with the
#' # weights re-estimated in each bootstrap sample: each bootstrap sample below
#' # re-fits ipwpoint() and the MSM. (A weighted lm() gives the same point
#' # estimate as svyglm() and is faster.)
#' if (requireNamespace("boot", quietly = TRUE)) {
#'   boot.fun <- function(data, index) {
#'     d <- data[index, ]
#'     w <- ipwpoint(
#'       exposure = a,
#'       family = "binomial",
#'       link = "logit",
#'       numerator = ~ 1,
#'       denominator = ~ l,
#'       data = d)
#'     d$sw <- w$ipw.weights
#'     coef(lm(y ~ a, data = d, weights = sw))["a"]
#'   }
#'   set.seed(1)
#'   bootres <- boot::boot(data = simdat, statistic = boot.fun, R = 200)
#'   print(bootres)
#'   # bootstrap standard error
#'   print(sd(bootres$t[, 1]))
#'   # percentile 95 percent confidence interval
#'   print(boot::boot.ci(bootres, type = "perc"))
#' }
ipwpoint <- function(
    exposure,
    family,
    link = NULL,
    numerator = NULL,
    denominator,
    data,
    trunc = NULL,
    estimand = c("ATE", "ATT", "ATO"),
    ...)
{
  call_obj <- match.call()
  dots <- list(...)
  estimand <- match.arg(estimand)

  # 1. Validate arguments ---------------------------------------------------
  valid_families <- c("binomial", "multinomial", "ordinal", "gaussian",
                      "poisson", "negbin", "gamma", "lognormal")
  if (missing(family) || length(family) != 1 || !(family %in% valid_families)) {
    stop("Invalid 'family' specified. Must be one of: ",
         paste(sprintf('"%s"', valid_families), collapse = ", "), call. = FALSE)
  }
  .check_link(family, link)

  if (missing(denominator) || !inherits(denominator, "formula")) {
    stop("A valid 'denominator' formula must be specified.", call. = FALSE)
  }
  if (!is.null(numerator) && !inherits(numerator, "formula")) {
    stop("When specified, 'numerator' must be a valid formula.", call. = FALSE)
  }
  continuous <- c("gaussian", "gamma", "lognormal")
  if (family %in% continuous && is.null(numerator)) {
    stop("'numerator' formula is required when family = '", family, "'.", call. = FALSE)
  }
  if (estimand != "ATE") {
    if (family != "binomial") {
      stop("estimand = '", estimand, "' is only available for family = 'binomial'.",
           call. = FALSE)
    }
    if (!is.null(numerator)) {
      stop("'numerator' must be NULL when estimand = '", estimand, "'.", call. = FALSE)
    }
  }
  smooth_ok <- c("binomial", "gaussian", "poisson")
  if ((.uses_smooths(denominator) || .uses_smooths(numerator)) && !(family %in% smooth_ok)) {
    stop("Smooth terms (s(), te(), ti(), t2()) are only supported for family = ",
         paste(sprintf("'%s'", smooth_ok), collapse = ", "), ".", call. = FALSE)
  }
  if (missing(data) || !is.data.frame(data)) {
    stop("A valid 'data' frame must be provided.", call. = FALSE)
  }
  .check_trunc(trunc)

  # 2. Exposure column ------------------------------------------------------
  if (missing(exposure)) stop("No 'exposure' variable specified.", call. = FALSE)
  exp_name <- .get_col_name(substitute(exposure))
  if (!(exp_name %in% names(data))) {
    stop("Exposure variable '", exp_name, "' not found in 'data'.", call. = FALSE)
  }
  exp_vec <- data[[exp_name]]

  # 3. Variables present and complete ---------------------------------------
  model_vars <- unique(c(exp_name, all.vars(denominator),
                         if (!is.null(numerator)) all.vars(numerator)))
  missing_vars <- setdiff(model_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  if (anyNA(data[, model_vars, drop = FALSE])) {
    stop("Missing values detected in exposure or covariates. 'ipwpoint' requires complete cases.",
         call. = FALSE)
  }
  if (family == "binomial" && !all(exp_vec %in% c(0, 1))) {
    stop("Binomial exposure must be coded as 0 or 1.", call. = FALSE)
  }
  if (family %in% c("poisson", "negbin") &&
      (!is.numeric(exp_vec) || any(exp_vec < 0) || any(exp_vec != round(exp_vec)))) {
    stop("Count exposure must contain non-negative integers for family = '", family, "'.",
         call. = FALSE)
  }
  if (family %in% c("gamma", "lognormal") && (!is.numeric(exp_vec) || any(exp_vec <= 0))) {
    stop("Exposure must be positive for family = '", family, "'.", call. = FALSE)
  }

  lhs <- switch(family,
                ordinal   = sprintf("as.factor(%s)", exp_name),
                lognormal = sprintf("log(%s)", exp_name),
                exp_name)

  # 4. One fitting routine per family: returns model + P(observed exposure) --
  fit_one <- function(rhs) {
    form <- .build_model_formula(lhs, rhs)
    base <- list(formula = form, data = data, na.action = stats::na.fail)

    p <- NULL
    if (family == "binomial") {
      mod <- .fit_glm_or_gam(c(base, list(family = stats::binomial(link = link)), dots))
      p <- as.numeric(stats::predict(mod, type = "response"))
      dens <- ifelse(exp_vec == 1, p, 1 - p)
    } else if (family == "poisson") {
      mod <- .fit_glm_or_gam(c(base, list(family = stats::poisson()), dots))
      dens <- stats::dpois(exp_vec, lambda = as.numeric(stats::predict(mod, type = "response")))
    } else if (family == "negbin") {
      mod <- do.call(MASS::glm.nb, c(base, dots))
      dens <- stats::dnbinom(exp_vec, size = mod$theta,
                             mu = as.numeric(stats::predict(mod, type = "response")))
    } else if (family == "gamma") {
      mod <- do.call(stats::glm, c(base, list(family = stats::Gamma(link = "log")), dots))
      shape <- 1 / summary(mod)$dispersion
      mu <- as.numeric(stats::predict(mod, type = "response"))
      dens <- stats::dgamma(exp_vec, shape = shape, rate = shape / mu)
    } else if (family == "lognormal") {
      mod <- do.call(stats::glm, c(base, list(family = stats::gaussian()), dots))
      dens <- stats::dlnorm(exp_vec, meanlog = as.numeric(stats::predict(mod)),
                            sdlog = stats::sd(stats::residuals(mod)))
    } else if (family == "multinomial") {
      mod <- do.call(nnet::multinom,
                     .with_defaults(base, dots, list(trace = FALSE, Hess = TRUE)))
      dens <- .prob_of_observed(stats::predict(mod, type = "probs"), exp_vec, mod$lev)
    } else if (family == "ordinal") {
      mod <- do.call(MASS::polr,
                     .with_defaults(c(base, list(method = .polr_method(link))), dots,
                                    list(Hess = TRUE)))
      dens <- .prob_of_observed(stats::predict(mod, type = "probs"), exp_vec, mod$lev)
    } else {
      mod <- .fit_glm_or_gam(c(base, list(family = stats::gaussian()), dots))
      sigma <- if (inherits(mod, "gam")) sqrt(mod$sig2) else stats::sd(stats::residuals(mod))
      dens <- stats::dnorm(exp_vec, mean = as.numeric(stats::predict(mod)), sd = sigma)
    }
    mod <- .tidy_model_call(mod, list(
      data   = call_obj$data,
      family = switch(family,
                      binomial  = call("binomial", link = link),
                      gaussian  = quote(gaussian()),
                      lognormal = quote(gaussian()),
                      poisson   = quote(poisson()),
                      gamma     = quote(Gamma(link = "log")),
                      NULL)
    ))
    list(model = mod, dens = unname(dens), p = p)
  }

  den <- fit_one(denominator)
  num <- if (!is.null(numerator)) fit_one(numerator) else list(model = NULL, dens = 1)

  # 5. Weights --------------------------------------------------------------
  ipw_weights <- switch(estimand,
    ATE = num$dens / den$dens,
    ATT = ifelse(exp_vec == 1, 1, den$p / (1 - den$p)),
    ATO = ifelse(exp_vec == 1, 1 - den$p, den$p)
  )
  if (anyNA(ipw_weights) || any(is.infinite(ipw_weights))) {
    stop("Invalid IPW weights calculated (NA or Inf detected). ",
         "Check for near-zero denominator probabilities.", call. = FALSE)
  }

  res <- list(
    ipw.weights   = ipw_weights,
    weights.trunc = if (!is.null(trunc)) .truncate_weights(ipw_weights, trunc) else NULL,
    call          = call_obj,
    num.mod       = num$model,
    den.mod       = den$model
  )
  res <- res[!vapply(res, is.null, logical(1))]
  class(res) <- c("ipwpoint", "list")
  res
}
