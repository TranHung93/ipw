#' Estimate Inverse Probability Weights (Point Treatment)
#'
#' @description
#' Estimates inverse probability weights to fit marginal structural models in a
#' point treatment situation. Exposure can be binomial, multinomial, ordinal,
#' or continuous (gaussian). Both stabilized and unstabilized weights can be estimated.
#'
#' @param exposure Exposure variable of interest (unquoted symbol or string).
#'   Binomial exposure variables should be coded as 0/1.
#' @param family Model family for exposure: \code{"binomial"}, \code{"multinomial"},
#'   \code{"ordinal"}, or \code{"gaussian"}.
#' @param link Link function for binomial (\code{"logit"}, \code{"probit"}, \code{"cauchit"},
#'   \code{"log"}, \code{"cloglog"}) or ordinal (\code{"logit"}, \code{"probit"},
#'   \code{"cauchit"}, \code{"cloglog"}). Ignored for multinomial and gaussian.
#' @param numerator Formula specifying the right-hand side of the numerator model.
#'   If \code{NULL}, unstabilized weights (numerator = 1) are calculated.
#' @param denominator Formula specifying the right-hand side of the denominator model.
#' @param data Data frame containing exposure and model covariates.
#' @param trunc Optional truncation percentile between 0 and 0.5 (e.g., 0.01 for 1st/99th percentiles).
#' @param ... Additional arguments passed to model fitting functions (\code{glm}, \code{multinom}, \code{polr}).
#'
#' @return A list of class \code{"ipwpoint"} containing:
#' \item{ipw.weights}{Numeric vector of inverse probability weights.}
#' \item{weights.trunc}{Numeric vector of truncated weights (if \code{trunc} is specified).}
#' \item{call}{The original function call.}
#' \item{num.mod}{Numerator model object (if \code{numerator} is specified).}
#' \item{den.mod}{Denominator model object.}
#' @author Willem M. van der Wal \email{willem@vanderwalresearch.com}, Ronald B. Geskus \email{rgeskus@oucru.org}
#' @references
#' Cole, S.R. & Hernán, M.A. (2008). Constructing inverse probability weights for marginal structural models. \emph{American Journal of Epidemiology}, \bold{168}(6), 656-664.
#'
#' Robins, J.M., Hernán, M.A. & Brumback, B.A. (2000). Marginal structural models and causal inference in epidemiology. \emph{Epidemiology}, \bold{11}, 550-560.
#'
#' Van der Wal W.M. & Geskus R.B. (2011). ipw: An R Package for Inverse Probability Weighting. \emph{Journal of Statistical Software}, \bold{43}(13), 1-23. \doi{10.18637/jss.v043.i13}.
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
ipwpoint <- function(
    exposure,
    family,
    link = NULL,
    numerator = NULL,
    denominator,
    data,
    trunc = NULL,
    ...)
{
  # 1. Retain original function call for return object
  call_obj <- match.call()

  # 2. Validate family
  valid_families <- c("binomial", "multinomial", "ordinal", "gaussian")
  if (missing(family) || !(family %in% valid_families)) {
    stop("Invalid 'family' specified. Must be one of: ",
         paste(sprintf('"%s"', valid_families), collapse = ", "))
  }

  # 3. Validate link function based on family
  if (family == "binomial") {
    valid_links <- c("logit", "probit", "cauchit", "log", "cloglog")
    if (is.null(link) || !(link %in% valid_links)) {
      stop("Invalid 'link' specified for family = 'binomial'. Must be one of: ",
           paste(sprintf('"%s"', valid_links), collapse = ", "))
    }
  } else if (family == "ordinal") {
    valid_links <- c("logit", "probit", "cauchit", "cloglog")
    if (is.null(link) || !(link %in% valid_links)) {
      stop("Invalid 'link' specified for family = 'ordinal'. Must be one of: ",
           paste(sprintf('"%s"', valid_links), collapse = ", "))
    }
  }

  # 4. Validate formulas, data, and truncation arguments
  if (missing(denominator) || !inherits(denominator, "formula")) {
    stop("A valid 'denominator' formula must be specified.")
  }
  if (!is.null(numerator) && !inherits(numerator, "formula")) {
    stop("When specified, 'numerator' must be a valid formula.")
  }
  if (family == "gaussian" && is.null(numerator)) {
    stop("'numerator' formula is required when family = 'gaussian'.")
  }
  if (missing(data) || !is.data.frame(data)) {
    stop("A valid 'data' frame must be provided.")
  }
  if (!is.null(trunc)) {
    if (!is.numeric(trunc) || length(trunc) != 1 || trunc < 0 || trunc > 0.5) {
      stop("'trunc' must be a numeric percentile between 0 and 0.5.")
    }
  }

  # 5. Extract exposure column name and vector safely
  exp_sub <- substitute(exposure)
  if (is.name(exp_sub) || is.character(exp_sub)) {
    exp_name <- as.character(exp_sub)
  } else {
    exp_name <- deparse(exp_sub)
  }

  if (exp_name %in% names(data)) {
    exp_vec <- data[[exp_name]]
  } else {
    exp_vec <- tryCatch(eval(exp_sub, data, parent.frame()), error = function(e) NULL)
    if (is.null(exp_vec)) {
      stop("Exposure variable '", exp_name, "' not found in 'data'.")
    }
  }

  # 6. Helper function: Construct model formulas programmatically
  build_model_formula <- function(exp_var_name, rhs_formula) {
    rhs_terms <- attr(terms(rhs_formula), "term.labels")
    intercept <- attr(terms(rhs_formula), "intercept") == 1
    lhs <- if (family == "ordinal") sprintf("as.factor(%s)", exp_var_name) else exp_var_name

    if (length(rhs_terms) == 0) {
      f_str <- paste(lhs, if (intercept) "~ 1" else "~ 0")
    } else {
      f_str <- paste(lhs, "~", paste(rhs_terms, collapse = " + "))
      if (!intercept) f_str <- paste(f_str, "- 1")
    }
    as.formula(f_str, env = environment(rhs_formula))
  }

  # 7. Check for missing values in required variables
  model_vars <- unique(c(exp_name, all.vars(denominator), if (!is.null(numerator)) all.vars(numerator)))
  missing_vars <- setdiff(model_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }
  if (anyNA(data[, model_vars, drop = FALSE])) {
    stop("Missing values detected in exposure or covariates. 'ipwpoint' requires complete cases.")
  }

  w_num <- NULL
  w_den <- NULL
  mod_num <- NULL
  mod_den <- NULL

  # ------------------------------------------------------------------
  # FIT MODELS & ESTIMATE WEIGHTS BY FAMILY
  # ------------------------------------------------------------------

  # --- BINOMIAL FAMILY ---
  if (family == "binomial") {
    link_family <- switch(link,
                          "logit"   = stats::binomial(link = "logit"),
                          "probit"  = stats::binomial(link = "probit"),
                          "cauchit" = stats::binomial(link = "cauchit"),
                          "log"     = stats::binomial(link = "log"),
                          "cloglog" = stats::binomial(link = "cloglog")
    )

    if (!all(exp_vec %in% c(0, 1))) {
      stop("Binomial exposure must be coded as 0 or 1.")
    }

    form_den <- build_model_formula(exp_name, denominator)
    mod_den  <- stats::glm(form_den, family = link_family, data = data, na.action = stats::na.fail, ...)
    p_den    <- stats::predict(mod_den, type = "response")
    w_den    <- ifelse(exp_vec == 1, p_den, 1 - p_den)

    if (is.null(numerator)) {
      w_num <- rep(1, nrow(data))
    } else {
      form_num <- build_model_formula(exp_name, numerator)
      mod_num  <- stats::glm(form_num, family = link_family, data = data, na.action = stats::na.fail, ...)
      p_num    <- stats::predict(mod_num, type = "response")
      w_num    <- ifelse(exp_vec == 1, p_num, 1 - p_num)
    }

    # --- MULTINOMIAL FAMILY ---
  } else if (family == "multinomial") {
    if (!requireNamespace("nnet", quietly = TRUE)) {
      stop("Package 'nnet' is required for family = 'multinomial'.")
    }

    extract_multinomial_probs <- function(model, obs_values) {
      probs <- as.matrix(stats::predict(model, type = "probs"))
      obs_char <- as.character(obs_values)
      col_idx <- match(obs_char, colnames(probs))
      if (any(is.na(col_idx))) {
        stop("Observed exposure levels do not match fitted probability levels.")
      }
      probs[cbind(seq_len(nrow(probs)), col_idx)]
    }

    form_den <- build_model_formula(exp_name, denominator)
    mod_den  <- nnet::multinom(form_den, data = data, na.action = stats::na.fail, trace = FALSE, ...)
    w_den    <- extract_multinomial_probs(mod_den, exp_vec)

    if (is.null(numerator)) {
      w_num <- rep(1, nrow(data))
    } else {
      form_num <- build_model_formula(exp_name, numerator)
      mod_num  <- nnet::multinom(form_num, data = data, na.action = stats::na.fail, trace = FALSE, ...)
      w_num    <- extract_multinomial_probs(mod_num, exp_vec)
    }

    # --- ORDINAL FAMILY ---
  } else if (family == "ordinal") {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package 'MASS' is required for family = 'ordinal'.")
    }

    method_map <- c(logit = "logistic", probit = "probit", cloglog = "cloglog", cauchit = "cauchit")
    polr_method <- method_map[link]

    extract_ordinal_probs <- function(model, obs_values) {
      probs <- as.matrix(stats::predict(model, type = "probs"))
      obs_char <- as.character(obs_values)
      col_idx <- match(obs_char, colnames(probs))
      if (any(is.na(col_idx))) {
        stop("Observed exposure levels do not match fitted probability levels.")
      }
      probs[cbind(seq_len(nrow(probs)), col_idx)]
    }

    form_den <- build_model_formula(exp_name, denominator)
    mod_den  <- MASS::polr(form_den, data = data, method = polr_method, na.action = stats::na.fail, ...)
    w_den    <- extract_ordinal_probs(mod_den, exp_vec)

    if (is.null(numerator)) {
      w_num <- rep(1, nrow(data))
    } else {
      form_num <- build_model_formula(exp_name, numerator)
      mod_num  <- MASS::polr(form_num, data = data, method = polr_method, na.action = stats::na.fail, ...)
      w_num    <- extract_ordinal_probs(mod_num, exp_vec)
    }

    # --- GAUSSIAN FAMILY ---
  } else if (family == "gaussian") {
    form_den <- build_model_formula(exp_name, denominator)
    mod_den  <- stats::glm(form_den, family = stats::gaussian(), data = data, na.action = stats::na.fail, ...)
    res_den  <- stats::residuals(mod_den)
    w_den    <- stats::dnorm(exp_vec, mean = stats::predict(mod_den), sd = stats::sd(res_den))

    form_num <- build_model_formula(exp_name, numerator)
    mod_num  <- stats::glm(form_num, family = stats::gaussian(), data = data, na.action = stats::na.fail, ...)
    res_num  <- stats::residuals(mod_num)
    w_num    <- stats::dnorm(exp_vec, mean = stats::predict(mod_num), sd = stats::sd(res_num))
  }

  # 8. Calculate IPW weights and check for numerical issues
  ipw_weights <- w_num / w_den
  if (anyNA(ipw_weights) || any(is.infinite(ipw_weights))) {
    stop("Invalid IPW weights calculated (NA or Inf detected). Check for near-zero denominator probabilities.")
  }

  # 9. Perform truncation if requested
  weights_trunc <- NULL
  if (!is.null(trunc)) {
    q_bounds <- stats::quantile(ipw_weights, probs = c(trunc, 1 - trunc), na.rm = TRUE)
    weights_trunc <- pmin(pmax(ipw_weights, q_bounds), q_bounds)
  }

  # 10. Assemble return list
  res <- list(
    ipw.weights = ipw_weights,
    weights.trunc = weights_trunc,
    call = call_obj,
    num.mod = mod_num,
    den.mod = mod_den
  )

  res <- res[!sapply(res, is.null)]
  class(res) <- c("ipwpoint", "list")
  return(res)
}
