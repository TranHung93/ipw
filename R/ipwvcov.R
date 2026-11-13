#' Standard Errors that Account for Estimated Weights (M-estimation)
#'
#' @description
#' \strong{Experimental:} this function is new in ipw 1.4.0. Its interface
#' and output may change in future versions, and the methods are still being
#' evaluated. The bootstrap, re-estimating the weights in
#' each bootstrap sample, remains the established alternative, especially
#' with heavy-tailed weights.
#'
#' Variance-covariance matrix of the coefficients
#' of a marginal structural model fitted with \code{lm()}, \code{glm()} or
#' \code{survival::coxph()} and
#' weights from \code{\link{ipwpoint}} for a binary point treatment, or with
#' time-varying weights from \code{\link{ipwtm}}, taking the
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
#' \strong{Cox models.} For a \code{coxph()} marginal structural model the
#' weighted partial likelihood score is used as the estimating equation
#' (Binder 1992; Lin & Wei 1989), with the Breslow form for tied event
#' times. Person \eqn{i}'s contribution is \eqn{w_i} times the score
#' residual, \eqn{U_i = w_i [\delta_i (X_i - \bar X(t_i)) - \sum_k \delta_k
#' w_k Y_i(t_k) e^{\beta' X_i} (X_i - \bar X(t_k)) / S_0(t_k)]}, with
#' \eqn{S_0(t) = \sum_j w_j Y_j(t) e^{\beta' X_j}} and \eqn{\bar X(t)} the
#' weighted mean of \eqn{X} in the risk set; \code{residuals(msm, type =
#' "score")} returns the part in brackets (without the weight). The bread
#' for \eqn{\beta} is the weighted information matrix (the inverse of
#' \code{msm$naive.var}). The weights enter both each person's own term and
#' the risk-set sums; differentiating the total score with respect to one
#' weight gives exactly that person's score residual, so the cross block is
#' \eqn{-\sum_i U_i d_i^T} with \eqn{d_i = \partial \log w_i / \partial
#' \alpha}, as for \code{glm()}. It is computed analytically; the package
#' tests check it against central finite differences of the stacked
#' estimating equations. Without the exposure models, the sandwich reduces
#' to \code{coxph()}'s robust variance (\code{msm$var}), which treats the
#' weights as known.
#'
#' Supported: \code{object} from \code{ipwpoint(family = "binomial",
#' link = "logit")} with denominator (and, if present, numerator) models
#' fitted by \code{glm()} (no smooth terms) and any \code{estimand}; \code{msm} of class
#' \code{"lm"} or \code{"glm"}, fitted on the same rows, in the same order, as
#' the weights, with \code{weights = object$ipw.weights}; or \code{msm} from
#' \code{coxph()} with right-censored data \code{Surv(time, status)}, one row
#' per person in the same order, and the same weights. For \code{coxph()},
#' tied event times need \code{ties = "breslow"}; without tied event times
#' all tie methods give the same result.
#'
#' \strong{Time-varying weights.} \code{object} can be an \code{ipwtm()}
#' result, or a list of \code{ipwtm()} (and \code{ipwpoint()}) results whose
#' \code{ipw.weights} multiply to the weights of \code{msm}, typically
#' treatment weights times censoring weights. Each model is a separate
#' parameter block. The unit of independence is the subject (the \code{id}
#' of \code{ipwtm()}): all scores are summed within subject before the meat
#' \eqn{B = \sum_{subjects} \psi_p \psi_p^T} is formed. The weight of row
#' \eqn{r} is the product, over the rows \eqn{s \le r} of the same subject,
#' of numerator over denominator factors \eqn{f_s}, as defined by
#' \code{ipwtm()}: \eqn{P(A_s = a_s)} for \code{type = "first"} and
#' \code{"all"}, \eqn{1 - p_s} for \code{type = "cens"} (the censoring row
#' included), and 1 after the first switch. With logistic models,
#' \eqn{\partial \log f_s / \partial \alpha} is \eqn{(a_s - p_s) Z_s}
#' (\code{"first"}, \code{"all"}) or \eqn{-p_s Z_s} (\code{"cens"}), so
#' \eqn{d_r = \partial \log w_r / \partial \alpha} is minus the cumulative
#' sum of these over the rows used by the model, within subject in time
#' order (plus, for the numerator). The cross block is
#' \eqn{-\sum_r U_r d_r^T} for \code{lm()}, \code{glm()} and Cox models; for a
#' Cox model on counting-process data \code{Surv(tstart, time, event)} the
#' identity above holds row by row, because each row enters the risk sets
#' with its own weight. Without the exposure models the sandwich reduces to
#' the cluster-robust variance that treats the weights as known (for a Cox
#' model, \code{coxph()} with \code{cluster(id)}).
#'
#' Supported for time-varying weights: \code{ipwtm(family = "binomial",
#' link = "logit")} with any \code{type}, glm models without smooth terms, no
#' truncation, made with this version of ipw; \code{msm} fitted with
#' \code{lm()}, \code{glm()} or \code{coxph()} (counting-process data;
#' \code{cluster()} is allowed, \code{ties = "breslow"} when event times are
#' tied) on all rows of the long data, in the same order, with weights equal
#' to the product of the weights in \code{object}.
#'
#' \strong{Heavy-tailed weights.} Cumulative weights can have a heavy right
#' tail when some people have a very low probability of their observed
#' exposure history (a near-violation of positivity). Normal-theory intervals
#' are then unreliable whatever method gives the standard error, the
#' bootstrap included. In simulations where treatment depended strongly on the
#' confounder and on past treatment, the tail index of the final cumulative
#' weights was about 2.3 (the weights have a finite variance, but not a
#' finite third or fourth moment). The mean variance from \code{ipwvcov()}
#' was then 0.84 to 0.93 times the true variance of the estimates (the
#' standard error 3 to 8\% too small), and 95\% intervals covered 92 to 93\%
#' of the time, with no trend towards 95\% between 300 and 4800 people. The
#' cluster bootstrap did no better (variance ratio 0.86 to 1.00 on fewer data
#' sets, coverage 90 to 95\%). With a tail index of about 3.9 the variance
#' ratios were 1.05 (\code{ipwvcov()}) and 1.04 (bootstrap), and both
#' covered 95.8\% and 96.5\% of the time. Check the weight distribution
#' (\code{summary()}, \code{\link{ipwplot}} and the effective sample size
#' \code{sum(w)^2 / sum(w^2)}) before trusting any standard error.
#'
#' Not supported, and an error: other exposure families or links, smooth terms
#' (GAMs), truncated weights (truncation at estimated percentiles is not a
#' smooth function of the model coefficients), \code{ipwtm()} weights from Cox
#' (\code{family = "survival"}), multinomial, ordinal or gaussian models, and
#' \code{coxph()} models with \code{strata()}, \code{frailty()} or \code{tt()}
#' terms, or Efron or exact ties with tied event times. With \code{ipwpoint()}
#' weights the \code{coxph()} model must have right-censored data (one row per
#' person) and no \code{cluster()}; with \code{ipwtm()} weights it must have
#' counting-process data. Use the bootstrap there, re-estimating the weights in
#' each bootstrap sample (see \code{\link{ipwpoint}} and \code{vignette("ipw")};
#' for long data, the cluster bootstrap that resamples subjects).
#'
#' This is the same kind of standard error as \code{vcov = "asympt"} in
#' \code{WeightIt::glm_weightit()}; in the package tests the two agree for
#' binary treatments. For Cox models, \code{ipwvcov()} agrees with a
#' numerical Jacobian of the stacked estimating equations to about
#' \eqn{10^{-7}}. \code{WeightIt::coxph_weightit(vcov = "asympt")} gives
#' slightly different values: for the ATE and overlap weights the relative
#' difference was below 0.1\% and shrank with the sample size; for ATT and
#' ATC weights its standard error was 1 to 4\% larger and the difference did
#' not shrink. Its \code{vcov = "HC0"} equals \code{coxph()}'s robust
#' variance exactly.
#'
#' @param object Result of \code{\link{ipwpoint}} with \code{family = "binomial"}
#'   and \code{link = "logit"}; or result of \code{\link{ipwtm}} with
#'   \code{family = "binomial"} and \code{link = "logit"}; or a list of
#'   \code{ipwtm()} (and \code{ipwpoint()}) results whose weights multiply to
#'   the weights of \code{msm}.
#' @param msm Marginal structural model: an \code{lm}, \code{glm} or
#'   \code{coxph} fit with weights equal to \code{object$ipw.weights} (or their
#'   product for a list), on the same rows in the same order as the data given
#'   to \code{ipwpoint()} or \code{ipwtm()}.
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
#' Binder, D.A. (1992). Fitting Cox's proportional hazards models from
#' survey data. \emph{Biometrika}, \bold{79}, 139-147. \doi{10.1093/biomet/79.1.139}.
#'
#' Lin, D.Y. & Wei, L.J. (1989). The robust inference for the Cox proportional
#' hazards model. \emph{Journal of the American Statistical Association},
#' \bold{84}, 1074-1078. \doi{10.1080/01621459.1989.10478874}.
#'
#' Stefanski, L.A. & Boos, D.D. (2002). The calculus of M-estimation.
#' \emph{The American Statistician}, \bold{56}, 29-38. \doi{10.1198/000313002753631330}.
#'
#' Greifer, N. WeightIt: Weighting for Covariate Balance in Observational
#' Studies. R package, \url{https://ngreifer.github.io/WeightIt/}. (Its
#' \code{glm_weightit(vcov = "asympt")} is used in the package tests as an
#' independent comparison.)
#'
#' @seealso \code{\link{ipwpoint}}, \code{\link{ipwtm}}
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
#'
#' # Cox marginal structural model for a survival outcome
#' library(survival)
#' simdat$time <- rexp(n, 0.1 * exp(-0.5 * simdat$a + 0.05 * simdat$l))
#' simdat$status <- as.numeric(simdat$time < 10)
#' simdat$time <- pmin(simdat$time, 10)
#' msm.cox <- coxph(Surv(time, status) ~ a, data = simdat, weights = w$ipw.weights)
#' cbind(robust = sqrt(diag(msm.cox$var)),          # weights treated as known
#'       ipwvcov = sqrt(diag(ipwvcov(w, msm.cox))))  # weights estimated
#'
#' # Time-varying weights: treatment (type = "first") times censoring
#' # (type = "cens") weights from pooled logistic models, and a Cox MSM on
#' # the long data
#' data(haartdat)
#' w.haart <- ipwtm(exposure = haartind, family = "binomial", link = "logit",
#'                  numerator = ~ sex + age, denominator = ~ sex + age + cd4.sqrt,
#'                  id = patient, timevar = fuptime, type = "first", data = haartdat)
#' w.cens <- ipwtm(exposure = dropout, family = "binomial", link = "logit",
#'                 numerator = ~ sex + age, denominator = ~ sex + age + cd4.sqrt,
#'                 id = patient, timevar = fuptime, type = "cens", data = haartdat)
#' msm.haart <- coxph(Surv(tstart, fuptime, event) ~ haartind + cluster(patient),
#'                    data = haartdat, ties = "breslow",
#'                    weights = w.haart$ipw.weights * w.cens$ipw.weights)
#' cbind(robust.cluster = sqrt(diag(msm.haart$var)),
#'       ipwvcov = sqrt(diag(ipwvcov(list(w.haart, w.cens), msm.haart))))
ipwvcov <- function(object, msm) {
  # 1. Checks ----------------------------------------------------------------
  if (inherits(object, "ipwtm") || .is_weight_list(object)) {
    return(.ipwvcov_tm(object, msm))
  }
  if (!inherits(object, "ipwpoint")) {
    stop("'object' must be the result of ipwpoint() or ipwtm(), or a list of such ",
         "results.", call. = FALSE)
  }
  estimand <- attr(object, "estimand")
  if (is.null(estimand)) {
    stop("'object' has no 'estimand' attribute; re-run ipwpoint() with this ",
         "version of ipw.", call. = FALSE)
  }
  .check_vcov_model(object$den.mod, "denominator")
  if (!is.null(object$num.mod)) .check_vcov_model(object$num.mod, "numerator")

  is_cox <- !missing(msm) && inherits(msm, "coxph")
  if (missing(msm) || !(is_cox || inherits(msm, "lm")) || inherits(msm, c("mlm", "gam"))) {
    stop("'msm' must be a model fitted with lm() or glm(), or with coxph() for ",
         "right-censored survival data.", call. = FALSE)
  }
  if (anyNA(stats::coef(msm))) {
    stop("'msm' has aliased (NA) coefficients; remove them from the model.", call. = FALSE)
  }
  if (is_cox) .check_coxph_msm(msm)
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
  if (is_cox) {
    cox <- .coxph_parts(msm, w, X)
    U <- cox$U
    A_bb <- cox$A
  } else if (inherits(msm, "glm")) {
    mu <- msm$fitted.values
    resid <- msm$y - mu
    mu_eta <- msm$family$mu.eta(msm$linear.predictors)
    var_mu <- msm$family$variance(mu)
  } else {
    resid <- msm$residuals
    mu_eta <- rep(1, n)
    var_mu <- rep(1, n)
  }
  if (!is_cox) {
    U <- X * (w * resid * mu_eta / var_mu)
    A_bb <- crossprod(X * (w * mu_eta^2 / var_mu), X)
  }

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

# A coxph() marginal structural model that ipwvcov() can handle: right-censored
# Surv(time, status), no strata, clusters, frailties or time-transforms, and
# Breslow ties whenever event times are tied (without tied event times all
# tie methods give the same partial likelihood).
.check_coxph_msm <- function(msm) {
  if (!identical(class(msm), "coxph")) {
    stop("ipwvcov() supports coxph() models without penalized terms (frailty(), ",
         "pspline()) and with at least one covariate.", call. = FALSE)
  }
  if (!identical(attr(msm$y, "type"), "right")) {
    stop("ipwvcov() supports a coxph() model with right-censored data, ",
         "Surv(time, status), one row per person; counting-process data ",
         "Surv(tstart, time, status) are not supported for ipwpoint() weights.",
         call. = FALSE)
  }
  specials <- attr(stats::terms(msm), "specials")
  used <- names(Filter(Negate(is.null), specials[c("strata", "cluster", "frailty", "tt")]))
  if (length(used) > 0 || !is.null(msm$strata) || !is.null(msm$call$cluster) ||
      !is.null(msm$call$tt)) {
    stop("ipwvcov() does not support strata(), cluster(), frailty() or tt() terms ",
         "in a coxph() model.", call. = FALSE)
  }
  event_times <- msm$y[msm$y[, 2] == 1, 1]
  if (anyDuplicated(event_times) && !identical(msm$method, "breslow")) {
    stop("The coxph() model has tied event times and uses ties = '", msm$method,
         "'. ipwvcov() needs ties = \"breslow\" when event times are tied.",
         call. = FALSE)
  }
  invisible(TRUE)
}

# Weighted Cox partial likelihood (Breslow) as an estimating equation.
# With r_j = w_j exp(eta_j), S0(t) = sum_{time_j >= t} r_j and
# S1(t) = sum_{time_j >= t} r_j X_j, person i's score is w_i times the
# (unweighted-form) score residual
#   res_i = status_i (X_i - Xbar(t_i)) - exp(eta_i) (X_i H0(t_i) - H1(t_i)),
#   H0(t) = sum_{event times u <= t} dW(u) / S0(u),
#   H1(t) = sum_{event times u <= t} dW(u) Xbar(u) / S0(u),
# with dW(u) the weighted number of events at u; res_i equals
# residuals(msm, type = "score"). The information is
#   sum_u dW(u) [S2(u) / S0(u) - Xbar(u) Xbar(u)'].
# Sums over risk sets use cumulative sums over the distinct times, so
# memory is linear in the number of rows.
.coxph_parts <- function(msm, w, X) {
  time <- msm$y[, 1]
  status <- msm$y[, 2]
  X <- matrix(X, nrow = nrow(X))   # drop model.matrix attributes
  k <- ncol(X)
  eta <- msm$linear.predictors   # centred; the constant cancels in all ratios
  r <- w * exp(eta)
  ut <- sort(unique(time))
  idx <- match(time, ut)
  m <- length(ut)
  by_time <- function(M) rowsum(as.matrix(M), idx, reorder = TRUE)
  rev_cumsum <- function(M) {
    rows <- rev(seq_len(nrow(M)))
    apply(M[rows, , drop = FALSE], 2, cumsum)[rows, , drop = FALSE]
  }
  rev_cumsum_m <- function(M) matrix(rev_cumsum(M), nrow = m)
  a_idx <- rep(seq_len(k), times = k)
  b_idx <- rep(seq_len(k), each = k)

  S0 <- rev_cumsum_m(by_time(r))[, 1]
  S1 <- rev_cumsum_m(by_time(X * r))
  S2 <- rev_cumsum_m(by_time(X[, a_idx, drop = FALSE] * X[, b_idx, drop = FALSE] * r))
  dW <- by_time(w * status)[, 1]
  Xbar <- S1 / S0

  info <- colSums(dW * (S2 / S0 - Xbar[, a_idx, drop = FALSE] * Xbar[, b_idx, drop = FALSE]))
  info <- matrix(info, k, k)
  H0 <- cumsum(dW / S0)
  H1 <- matrix(apply(dW * Xbar / S0, 2, cumsum), nrow = m)
  res <- status * (X - Xbar[idx, , drop = FALSE]) -
    exp(eta) * (X * H0[idx] - H1[idx, , drop = FALSE])
  list(U = res * w, A = info, res = res)
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
