#' Standard Errors that Account for Estimated Weights (M-estimation)
#'
#' @description
#' \strong{New, experimental.} Variance-covariance matrix of the coefficients
#' of a marginal structural model fitted with \code{lm()} or \code{glm()} and
#' weights from \code{\link{ipwpoint}} for a binary point treatment, taking the
#' estimation of the weights into account. The marginal structural model and
#' the exposure models are treated as one set of stacked estimating equations
#' (M-estimation; Lunceford & Davidian 2004), and the sandwich variance of all
#' parameters is computed; the block for the coefficients of the marginal
#' structural model is returned.
#'
#' @details
#' Robust (sandwich) standard errors of a weighted model, such as those from
#' \code{svyglm()} or \code{sandwich::vcovHC()}, treat the weights as known.
#' Here the estimating equations of the denominator (and numerator) model are
#' stacked under those of the marginal structural model. With \eqn{\theta =
#' (\beta, \alpha, \gamma)} the coefficients of the marginal structural model,
#' the denominator model and the numerator model, and \eqn{\psi_i(\theta)}
#' the contribution of person \eqn{i} to the stacked estimating equations,
#' the variance is \eqn{A^{-1} B A^{-T}} with \eqn{A = -\sum_i \partial
#' \psi_i / \partial \theta} and \eqn{B = \sum_i \psi_i \psi_i^T}, evaluated
#' at the estimates. The derivative of each weight with respect to the
#' exposure-model coefficients is worked out for each estimand (ATE, ATT,
#' ATC, ATO). The exposure models must use the logit link: the information
#' matrices use the expected (Fisher) information, which equals the observed
#' information, and so gives the exact sandwich, only for canonical links. For
#' the marginal structural model any \code{glm()} link is accepted. With a
#' non-canonical link (for example a log link for risk ratios) and a
#' saturated model such as \code{y ~ a} for the binary exposure, the expected
#' and observed information are exactly equal in any sample: the fitted mean
#' in each exposure group is that group's weighted mean of \code{y}, so the
#' extra term of the observed information, a weighted sum of residuals times
#' a factor that is constant within each group, sums to zero within each
#' group. For other models with a non-canonical link the two are
#' asymptotically equivalent when the model is correctly specified.
#' No small-sample correction is applied. Proportional rescaling of the
#' weights does not change the result.
#'
#' Supported: \code{object} from \code{ipwpoint(family = "binomial",
#' link = "logit")} with denominator (and, if present, numerator) models
#' fitted by \code{glm()} (no smooth terms) and any \code{estimand}; \code{msm} of class
#' \code{"lm"} or \code{"glm"}, fitted on the same rows, in the same order, as
#' the weights, with \code{weights = object$ipw.weights}.
#'
#' Not supported, and an error: other exposure families or links, smooth terms
#' (GAMs), truncated weights (truncation at estimated percentiles is not a
#' smooth function of the model coefficients) and time-varying weights from
#' \code{\link{ipwtm}}. Use the bootstrap there, re-estimating the weights in
#' each bootstrap sample (see \code{\link{ipwpoint}} and \code{vignette("ipw")}).
#'
#' This is the same kind of standard error as \code{vcov = "asympt"} in
#' \code{WeightIt::glm_weightit()}; in the package tests the two agree for
#' binary treatments.
#'
#' @param object Result of \code{\link{ipwpoint}} with \code{family = "binomial"}
#'   and \code{link = "logit"}.
#' @param msm Marginal structural model: an \code{lm} or \code{glm} fit with
#'   prior weights equal to \code{object$ipw.weights}, on the same rows in the
#'   same order as the data given to \code{ipwpoint()}.
#'
#' @return Variance-covariance matrix of the coefficients of \code{msm}, with
#'   the coefficient names as row and column names. It can be passed, for
#'   example, to \code{lmtest::coeftest(msm, vcov. = ...)}.
#'
#' @references
#' Lunceford, J.K. & Davidian, M. (2004). Stratification and weighting via the
#' propensity score in estimation of causal treatment effects: a comparative
#' study. \emph{Statistics in Medicine}, \bold{23}, 2937-2960. \doi{10.1002/sim.1903}.
#'
#' Stefanski, L.A. & Boos, D.D. (2002). The calculus of M-estimation.
#' \emph{The American Statistician}, \bold{56}, 29-38. \doi{10.1198/000313002753631330}.
#'
#' Greifer, N. WeightIt: Weighting for Covariate Balance in Observational
#' Studies. R package, \url{https://ngreifer.github.io/WeightIt/}. (Its
#' \code{glm_weightit(vcov = "asympt")} is used in the package tests as an
#' independent comparison.)
#'
#' @seealso \code{\link{ipwpoint}}
#' @export
#'
#' @examples
#' set.seed(16)
#' n <- 1000
#' simdat <- data.frame(l = rnorm(n, 10, 5))
#' simdat$a <- rbinom(n, 1, plogis(simdat$l - 10))
#' simdat$y <- 10 * simdat$a + 0.5 * simdat$l + rnorm(n, -10, 5)
#'
#' w <- ipwpoint(exposure = a, family = "binomial", link = "logit",
#'               numerator = ~ 1, denominator = ~ l, data = simdat)
#' msm <- lm(y ~ a, data = simdat, weights = w$ipw.weights)
#'
#' # standard errors that account for the estimated weights
#' sqrt(diag(ipwvcov(w, msm)))
#'
#' # overlap weights
#' w.ato <- ipwpoint(exposure = a, family = "binomial", link = "logit",
#'                   denominator = ~ l, data = simdat, estimand = "ATO")
#' msm.ato <- lm(y ~ a, data = simdat, weights = w.ato$ipw.weights)
#' sqrt(diag(ipwvcov(w.ato, msm.ato)))
#'
#' # Wald tests and confidence intervals
#' if (requireNamespace("lmtest", quietly = TRUE)) {
#'   print(lmtest::coeftest(msm, vcov. = ipwvcov(w, msm)))
#'   print(lmtest::coefci(msm, vcov. = ipwvcov(w, msm)))
#' }
ipwvcov <- function(object, msm) {
  # 1. Checks ----------------------------------------------------------------
  if (inherits(object, "ipwtm")) {
    stop("ipwvcov() does not support time-varying weights from ipwtm(). ",
         "Use the cluster bootstrap that re-estimates the weights ",
         "(vignette(\"ipw\"), section 5).", call. = FALSE)
  }
  if (!inherits(object, "ipwpoint")) {
    stop("'object' must be the result of ipwpoint().", call. = FALSE)
  }
  estimand <- attr(object, "estimand")
  if (is.null(estimand)) {
    stop("'object' has no 'estimand' attribute; re-run ipwpoint() with this ",
         "version of ipw.", call. = FALSE)
  }
  .check_vcov_model(object$den.mod, "denominator")
  if (!is.null(object$num.mod)) .check_vcov_model(object$num.mod, "numerator")

  if (missing(msm) || !inherits(msm, "lm") || inherits(msm, c("mlm", "gam"))) {
    stop("'msm' must be a model fitted with lm() or glm().", call. = FALSE)
  }
  if (anyNA(stats::coef(msm))) {
    stop("'msm' has aliased (NA) coefficients; remove them from the model.", call. = FALSE)
  }
  w <- object$ipw.weights
  n <- length(w)
  X <- stats::model.matrix(msm)
  if (nrow(X) != n) {
    stop("'msm' was fitted on ", nrow(X), " rows but there are ", n, " weights. ",
         "Fit it on the same rows, in the same order, as ipwpoint().", call. = FALSE)
  }
  msm_w <- if (inherits(msm, "glm")) msm$prior.weights else msm$weights
  if (is.null(msm_w)) msm_w <- rep(1, n)
  msm_w <- as.numeric(msm_w)
  if (!isTRUE(all.equal(msm_w, as.numeric(w)))) {
    if (!is.null(object$weights.trunc) &&
        isTRUE(all.equal(msm_w, as.numeric(object$weights.trunc)))) {
      stop("Truncated weights are not supported: truncation at estimated percentiles ",
           "is not a smooth function of the model coefficients. Use the bootstrap, ",
           "re-estimating and truncating the weights in each sample.", call. = FALSE)
    }
    stop("The weights of 'msm' must equal object$ipw.weights, on the same rows in ",
         "the same order (fit it with weights = object$ipw.weights).", call. = FALSE)
  }

  # 2. Marginal structural model: scores U_i and information A_bb --------------
  if (inherits(msm, "glm")) {
    mu <- msm$fitted.values
    resid <- msm$y - mu
    mu_eta <- msm$family$mu.eta(msm$linear.predictors)
    var_mu <- msm$family$variance(mu)
  } else {
    resid <- msm$residuals
    mu_eta <- rep(1, n)
    var_mu <- rep(1, n)
  }
  U <- X * (w * resid * mu_eta / var_mu)
  A_bb <- crossprod(X * (w * mu_eta^2 / var_mu), X)

  # 3. Denominator model: scores S_i, information, d log w / d alpha ----------
  den <- .binomial_parts(object$den.mod)
  a <- den$a
  p <- den$p
  g <- switch(estimand,
    ATE = ifelse(a == 1, -1 / p, 1 / (1 - p)),
    ATT = ifelse(a == 1, 0, 1 / (p * (1 - p))),
    ATC = ifelse(a == 1, -1 / (p * (1 - p)), 0),
    ATO = ifelse(a == 1, -1 / (1 - p), 1 / p)
  )
  D <- den$Z * (g * den$dp)
  blocks_score <- list(U, den$S)
  blocks_info <- list(den$A)
  blocks_cross <- list(-crossprod(U, D))

  # 4. Numerator model (stabilized ATE weights) --------------------------------
  if (!is.null(object$num.mod)) {
    num <- .binomial_parts(object$num.mod)
    q <- num$p
    E <- num$Z * (ifelse(a == 1, 1 / q, -1 / (1 - q)) * num$dp)
    blocks_score <- c(blocks_score, list(num$S))
    blocks_info <- c(blocks_info, list(num$A))
    blocks_cross <- c(blocks_cross, list(-crossprod(U, E)))
  }

  # 5. Sandwich: A is block upper-triangular ------------------------------------
  Psi <- do.call(cbind, blocks_score)
  k <- ncol(X)
  sizes <- vapply(blocks_info, ncol, integer(1))
  A <- matrix(0, ncol(Psi), ncol(Psi))
  A[seq_len(k), seq_len(k)] <- A_bb
  start <- k
  for (j in seq_along(blocks_info)) {
    idx <- start + seq_len(sizes[j])
    A[idx, idx] <- blocks_info[[j]]
    A[seq_len(k), idx] <- blocks_cross[[j]]
    start <- start + sizes[j]
  }
  B <- crossprod(Psi)
  A_inv <- solve(A)
  V <- A_inv %*% B %*% t(A_inv)
  V <- V[seq_len(k), seq_len(k), drop = FALSE]
  dimnames(V) <- list(names(stats::coef(msm)), names(stats::coef(msm)))
  V
}

# An exposure model that ipwvcov() can handle: a logistic glm (not a GAM),
# unweighted, without aliased coefficients.
.check_vcov_model <- function(mod, which) {
  if (inherits(mod, "gam")) {
    stop("ipwvcov() does not support smooth terms (GAM) in the ", which,
         " model. Use the bootstrap.", call. = FALSE)
  }
  if (!inherits(mod, "glm") || !identical(mod$family$family, "binomial")) {
    stop("ipwvcov() is only available for ipwpoint() with family = 'binomial'. ",
         "Use the bootstrap for other exposures.", call. = FALSE)
  }
  if (!identical(mod$family$link, "logit")) {
    stop("ipwvcov() needs link = 'logit' for the ", which, " model (it has '",
         mod$family$link, "'). With other links the result is only approximate; ",
         "use the bootstrap.", call. = FALSE)
  }
  if (any(mod$prior.weights != 1)) {
    stop("ipwvcov() does not support prior weights in the ", which, " model.",
         call. = FALSE)
  }
  if (anyNA(stats::coef(mod))) {
    stop("The ", which, " model has aliased (NA) coefficients.", call. = FALSE)
  }
  invisible(TRUE)
}

# Pieces of a binomial glm used in the stacked estimating equations:
# design matrix Z, exposure a, fitted p, dp/deta, scores S_i and the
# (expected) information A.
.binomial_parts <- function(mod) {
  Z <- stats::model.matrix(mod)
  p <- as.numeric(mod$fitted.values)
  a <- as.numeric(mod$y)
  dp <- as.numeric(mod$family$mu.eta(mod$linear.predictors))
  v <- p * (1 - p)
  list(Z = Z, a = a, p = p, dp = dp,
       S = Z * ((a - p) * dp / v),
       A = crossprod(Z * (dp^2 / v), Z))
}
