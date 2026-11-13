#' Estimate Inverse Probability Weights (Point Treatment)
#'
#' @description
#' Estimates inverse probability weights to fit marginal structural models in a
#' point treatment situation. Exposure can be binomial, multinomial, ordinal,
#' a count (Poisson or negative binomial) or continuous (gaussian, gamma or
#' log-normal). Both stabilized and unstabilized weights can be estimated, and
#' for a binomial exposure also weights for the effect in the treated (ATT), the
#' effect in the untreated (ATC) and overlap weights (ATO). For a multinomial
#' or ordinal exposure, ATT weights (for a chosen \code{focal} level) and
#' generalized overlap weights are available as well. Formulas may
#' contain smooth terms such as \code{s(x)}, which are then fitted as
#' generalized additive models with \pkg{mgcv}.
#'
#' @details
#' Standard errors: \code{ipwpoint} only estimates the weights. Standard errors
#' for the marginal structural model have to be obtained when that model is fitted.
#' \code{\link{ipwvcov}} (experimental) gives standard errors that take the
#' estimation of the weights into account (M-estimation), and is the recommended
#' option where it applies: weights from \code{ipwpoint()} with
#' \code{family = "binomial"} and \code{link = "logit"} (any \code{estimand},
#' stabilized or not, no smooth terms, no truncation) and a marginal structural
#' model fitted with \code{lm()}, \code{glm()} (any link) or \code{coxph()}
#' (right-censored data, one row per person); and time-varying weights from
#' \code{\link{ipwtm}} with \code{family = "binomial"} and \code{link = "logit"}
#' (\code{type = "first"}, \code{"cens"} or \code{"all"}, alone or multiplied,
#' e.g. treatment times censoring weights) with an \code{lm()}, \code{glm()} or
#' counting-process \code{coxph()} model on the long data. For all other cases
#' use the bootstrap described below.
#' Robust (sandwich) standard errors, for example from \code{svyglm()} in the survey
#' package as in Van der Wal and Geskus (2011), treat the weights as known. For a
#' binary point treatment with good overlap they were conservative (too large) in
#' simulations, but they
#' can be too small for a continuous exposure (\code{family = "gaussian"}) and for
#' time-varying weights; use the bootstrap there. Austin (2016) compared naive, robust and bootstrap
#' standard errors in simulations and recommends the bootstrap, in which the weights
#' are re-estimated in each bootstrap sample. This takes the uncertainty in the
#' estimated weights into account; see the bootstrap example below. About 200
#' bootstrap samples are usually enough to estimate a standard error. In the package's simulation study
#' (article "Simulation study" on the package website), the bootstrap standard
#' error was accurate for a point treatment with good overlap, while the
#' variance from robust standard errors was 1.25 times the true variance (the
#' standard error 12\% too large). When some weights are extreme
#' (poor overlap), no standard error was reliable, the bootstrap included. Check
#' the weights first with \code{summary()} and \code{\link{ipwplot}}, and consider
#' truncation (\code{trunc}) or overlap weights (\code{estimand = "ATO"}).
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
#'   \code{"ATT"}, \code{"ATC"} or \code{"ATO"}.
#' @param denominator Formula specifying the right-hand side of the denominator model.
#'   Smooth terms (\code{s()}, \code{te()}, \code{ti()}, \code{t2()}) are allowed
#'   for the binomial, gaussian and poisson families; the model is then fitted with
#'   \code{mgcv::gam}.
#' @param data Data frame containing exposure and model covariates.
#' @param trunc Optional truncation percentile between 0 and 0.5 (e.g., 0.01 for 1st/99th percentiles).
#' @param estimand Target estimand. For a binomial exposure: \code{"ATE"} (default;
#'   inverse probability weights), \code{"ATT"} (effect in the treated: treated
#'   get weight 1, untreated \eqn{p/(1-p)}), \code{"ATC"} (effect in the
#'   untreated, or controls: treated get \eqn{(1-p)/p}, untreated 1) or
#'   \code{"ATO"} (overlap weights: treated \eqn{1-p}, untreated \eqn{p}), where
#'   \eqn{p} is the fitted probability of exposure from the denominator model.
#'   For a multinomial or ordinal exposure: \code{"ATE"}, \code{"ATT"} (effect
#'   in the people at level \code{focal}: weight \eqn{P_i(focal) / P_i(A_i)},
#'   1 for people at the focal level) or \code{"ATO"} (generalized overlap
#'   weights, Li & Li 2019: \eqn{(1/P_i(A_i)) / \sum_k 1/P_i(k)}), where
#'   \eqn{P_i(k)} is the fitted probability of level \eqn{k} and \eqn{A_i} the
#'   observed level. For two levels these equal the binomial ATT and ATO
#'   weights. \code{"ATC"} is only available for a binomial exposure; for a
#'   categorical exposure use \code{"ATT"} with \code{focal} set to the control
#'   level. Other families only allow \code{"ATE"}.
#' @param focal For \code{estimand = "ATT"} with \code{family = "multinomial"}
#'   or \code{"ordinal"}: the exposure level whose effect is targeted (the
#'   "treated" group), given as a level of the exposure. Required there and
#'   must be \code{NULL} otherwise; for a binomial exposure the focal group is
#'   always the exposed (1).
#' @param ... Additional arguments passed to model fitting functions (\code{glm},
#'   \code{multinom}, \code{polr}, \code{glm.nb}, \code{gam}).
#'
#' @return A list of class \code{"ipwpoint"} with attributes \code{"estimand"},
#'   \code{"exposure"} (the name of the exposure column) and, for ATT and ATC
#'   weights, \code{"focal"} (the level whose effect is targeted: 1 for a
#'   binomial ATT, 0 for ATC), containing:
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
#' Li, F., Morgan, K.L. & Zaslavsky, A.M. (2018). Balancing covariates via propensity score weighting. \emph{Journal of the American Statistical Association}, \bold{113}, 390-400. \doi{10.1080/01621459.2016.1260466}.
#'
#' Greifer, N. WeightIt: Weighting for Covariate Balance in Observational Studies. R package, \url{https://ngreifer.github.io/WeightIt/}. (Used in the package tests as an independent comparison of the weights.)
#'
#' Li, F. & Li, F. (2019). Propensity score weighting for causal inference with multiple treatments. \emph{Annals of Applied Statistics}, \bold{13}(4), 2389-2415. \doi{10.1214/19-AOAS1282}.
#'
#' Austin, P.C. (2016). Variance estimation when using inverse probability of treatment weighting (IPTW) with survival analysis. \emph{Statistics in Medicine}, \bold{35}, 5642-5655. \doi{10.1002/sim.7084}.
#' @export
#' @seealso \code{\link{ipwvcov}}, \code{\link{basdat}}, \code{\link{haartdat}}, \code{\link{ipwplot}}, \code{\link{ipwpoint}}, \code{\link{ipwtm}}, \code{\link{timedat}}, \code{\link{tstartfun}}.
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
#' # Standard errors that take the estimation of the weights into account
#' # (M-estimation; binomial exposure, lm() or glm() MSM)
#' msm.lm <- lm(y ~ a, data = simdat, weights = temp$ipw.weights)
#' sqrt(diag(ipwvcov(temp, msm.lm)))
#'
#' # Bootstrap standard error and percentile confidence interval for the MSM.
#' # The robust standard error from svyglm() treats the weights as known and
#' # can be too large or too small (simulation study, scenarios 7a and 7b).
#' # Austin (2016) recommends the bootstrap, with the
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
    estimand = c("ATE", "ATT", "ATC", "ATO"),
    focal = NULL,
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
  categorical <- c("multinomial", "ordinal")
  if (estimand != "ATE") {
    estimand_families <- if (estimand == "ATC") "binomial" else c("binomial", categorical)
    if (!(family %in% estimand_families)) {
      stop("estimand = '", estimand, "' is only available for family = ",
           paste(sprintf("'%s'", estimand_families), collapse = ", "), ".",
           if (estimand == "ATC" && family %in% categorical)
             " For a categorical exposure, use estimand = 'ATT' with 'focal' set to the control level.",
           call. = FALSE)
    }
    if (!is.null(numerator)) {
      stop("'numerator' must be NULL when estimand = '", estimand, "'.", call. = FALSE)
    }
  }
  categorical_att <- estimand == "ATT" && family %in% categorical
  if (!is.null(focal) && !categorical_att) {
    stop("'focal' is only used with estimand = 'ATT' and family = 'multinomial' or ",
         "'ordinal'. For a binomial exposure the focal group is always the exposed (1).",
         call. = FALSE)
  }
  if (categorical_att && is.null(focal)) {
    stop("estimand = 'ATT' with family = '", family, "' needs 'focal': the exposure ",
         "level whose effect is targeted.", call. = FALSE)
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
  model_vars <- setdiff(unique(c(exp_name, all.vars(denominator),
                                 if (!is.null(numerator)) all.vars(numerator))), ".")
  missing_vars <- setdiff(model_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  if (anyNA(data[, model_vars, drop = FALSE])) {
    stop("Missing values detected in exposure or covariates. 'ipwpoint' requires complete cases.",
         call. = FALSE)
  }
  .check_infinite(data, model_vars)
  if (family == "binomial" && (is.factor(exp_vec) || is.character(exp_vec))) {
    exp_vec <- .as_binary01(exp_vec)
    data <- as.data.frame(data)
    data[[exp_name]] <- exp_vec
  }
  if (family == "binomial" &&
      (!(is.numeric(exp_vec) || is.logical(exp_vec)) || !all(exp_vec %in% c(0, 1)))) {
    stop("Binomial exposure must be coded as 0 or 1.", call. = FALSE)
  }
  if (family %in% c("binomial", "multinomial", "ordinal") && length(unique(exp_vec)) < 2) {
    stop("The exposure takes only one value; weights cannot be estimated.", call. = FALSE)
  }
  if (family == "ordinal") .check_ordinal_exposure(exp_vec)
  if (categorical_att) {
    exp_levels <- if (is.factor(exp_vec)) levels(droplevels(exp_vec)) else
      as.character(sort(unique(exp_vec)))
    if (length(focal) != 1 || is.na(focal) || !(as.character(focal) %in% exp_levels)) {
      stop("'focal' must be one level of the exposure: ",
           paste(sprintf("'%s'", exp_levels), collapse = ", "), ".", call. = FALSE)
    }
    focal <- as.character(focal)
  }
  if (family %in% c("poisson", "negbin") &&
      (!is.numeric(exp_vec) || any(exp_vec < 0) || any(exp_vec != round(exp_vec)))) {
    stop("Count exposure must contain non-negative integers for family = '", family, "'.",
         call. = FALSE)
  }
  if (family == "gaussian" && !is.numeric(exp_vec)) {
    stop("Exposure must be numeric for family = 'gaussian'.", call. = FALSE)
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
    pmat <- NULL
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
      pmat <- .prob_matrix(stats::predict(mod, type = "probs"), mod$lev, length(exp_vec))
      dens <- .prob_of_observed(pmat, exp_vec, mod$lev)
    } else if (family == "ordinal") {
      mod <- do.call(MASS::polr,
                     .with_defaults(c(base, list(method = .polr_method(link))), dots,
                                    list(Hess = TRUE)))
      pmat <- .prob_matrix(stats::predict(mod, type = "probs"), mod$lev, length(exp_vec))
      dens <- .prob_of_observed(pmat, exp_vec, mod$lev)
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
    list(model = mod, dens = unname(dens), p = p, pmat = pmat)
  }

  den <- fit_one(denominator)
  num <- if (!is.null(numerator)) fit_one(numerator) else list(model = NULL, dens = 1)

  # 5. Weights --------------------------------------------------------------
  if (family %in% categorical && estimand != "ATE") {
    # den$pmat holds P_i(k) for every level k; den$dens is P_i(observed level)
    ipw_weights <- unname(switch(estimand,
      ATT = ifelse(as.character(exp_vec) == focal, 1, den$pmat[, focal] / den$dens),
      ATO = (1 / den$dens) / rowSums(1 / den$pmat)
    ))
  } else {
    ipw_weights <- switch(estimand,
      ATE = num$dens / den$dens,
      ATT = ifelse(exp_vec == 1, 1, den$p / (1 - den$p)),
      ATC = ifelse(exp_vec == 1, (1 - den$p) / den$p, 1),
      ATO = ifelse(exp_vec == 1, 1 - den$p, den$p)
    )
  }
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
  attr(res, "estimand") <- estimand
  attr(res, "exposure") <- exp_name
  attr(res, "focal") <- if (categorical_att) focal else
    switch(estimand, ATT = 1, ATC = 0, NULL)
  res
}
