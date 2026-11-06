#' Time-Varying Inverse Probability Weights From User-Supplied Probabilities
#'
#' @description
#' Computes time-varying inverse probability weights from probabilities (or
#' densities) of the observed exposure that you estimated yourself, for
#' example with a machine-learning model or SuperLearner. For each row the
#' ratio \code{num / den} is computed, and the ratios are multiplied
#' cumulatively over time within each subject, as in \code{\link{ipwtm}}.
#'
#' @param den Numeric vector: for each row, the probability (or density) of the
#'   observed exposure given the covariates of the denominator model.
#' @param num Optional numeric vector of the same length: the probability (or
#'   density) of the observed exposure given the covariates of the numerator
#'   model. If \code{NULL}, unstabilized weights (numerator = 1) are computed.
#' @param id Vector of subject identifiers, same length as \code{den}.
#' @param timevar Numeric vector of follow-up times, same length as \code{den}.
#'   It sets the order in which the ratios are multiplied within a subject, so
#'   the data do not need to be sorted.
#' @param trunc Optional truncation percentile between 0 and 0.5 (e.g., 0.01 for
#'   1st/99th percentiles).
#'
#' @details
#' \code{ipwprobtm} does the last step of \code{\link{ipwtm}} only: the
#' exposure model is fitted by you, and \code{ipwprobtm} forms the ratios,
#' multiplies them within each subject in time order, returns the weights in
#' the order of the input and optionally truncates them.
#'
#' For weights of \code{type = "first"} or \code{"cens"} (see
#' \code{\link{ipwtm}}), set the probabilities to 1 for the rows after the first
#' switch, so that those rows do not change the cumulative weight.
#'
#' For a point treatment no function is needed: with your own probabilities
#' the weights are simply \code{num / den} (or \code{1 / den} for
#' unstabilized weights).
#'
#' @return A list of class \code{"ipwprobtm"} containing:
#' \item{ipw.weights}{Numeric vector of inverse probability weights, in the
#'   order of the input.}
#' \item{weights.trunc}{Numeric vector of truncated weights (if \code{trunc} is
#'   specified).}
#' \item{call}{The original function call.}
#' @author Hung Thai Tran \email{hungtt@oucru.org}
#' @seealso \code{\link{ipwtm}}, \code{\link{ipwpoint}}.
#' @export
#'
#' @examples
#' # Probabilities from any model; the ratios are multiplied within each
#' # subject over time
#' long <- data.frame(id = rep(1:3, each = 2), time = rep(1:2, 3),
#'                    p = c(0.5, 0.8, 0.4, 0.9, 0.6, 0.7))
#' ipwprobtm(den = long$p, id = long$id, timevar = long$time)$ipw.weights
#'
#' # Stabilized weights from a pooled logistic model fitted by hand; this
#' # reproduces ipwtm(..., type = "all")
#' set.seed(1)
#' n <- 200
#' long <- data.frame(id = rep(1:n, each = 3), time = rep(1:3, n))
#' long$l <- rnorm(3 * n)
#' long$a <- rbinom(3 * n, 1, plogis(-0.5 + long$l))
#' pd <- fitted(glm(a ~ l, family = binomial, data = long))
#' pn <- fitted(glm(a ~ 1, family = binomial, data = long))
#' w <- ipwprobtm(den = ifelse(long$a == 1, pd, 1 - pd),
#'                num = ifelse(long$a == 1, pn, 1 - pn),
#'                id = long$id, timevar = long$time)
#' summary(w$ipw.weights)
ipwprobtm <- function(den, num = NULL, id, timevar, trunc = NULL) {
  call_obj <- match.call()

  .check_prob <- function(x, name) {
    if (!is.numeric(x) || length(x) == 0) {
      stop("'", name, "' must be a non-empty numeric vector.", call. = FALSE)
    }
    if (anyNA(x) || any(!is.finite(x)) || any(x <= 0)) {
      stop("'", name, "' must contain positive, finite values without missing values.",
           call. = FALSE)
    }
  }
  .check_prob(den, "den")
  n <- length(den)
  if (is.null(num)) {
    num <- rep(1, n)
  } else {
    .check_prob(num, "num")
    if (length(num) != n) stop("'num' and 'den' must have the same length.", call. = FALSE)
  }
  if (missing(id) || is.null(id)) {
    stop("'id' is required. For a point treatment the weights are simply num / den.",
         call. = FALSE)
  }
  if (missing(timevar) || is.null(timevar)) {
    stop("'timevar' is required.", call. = FALSE)
  }
  .check_trunc(trunc)

  if (length(id) != n) stop("'id' and 'den' must have the same length.", call. = FALSE)
  if (!is.numeric(timevar) || length(timevar) != n) {
    stop("'timevar' must be numeric and have the same length as 'den'.", call. = FALSE)
  }
  if (anyNA(id) || anyNA(timevar)) {
    stop("'id' and 'timevar' must not contain missing values.", call. = FALSE)
  }
  if (anyDuplicated(data.frame(id, timevar))) {
    stop("Each combination of 'id' and 'timevar' must occur only once.", call. = FALSE)
  }

  ratio <- num / den
  ord <- order(id, timevar)
  cum <- unsplit(lapply(split(ratio[ord], id[ord]), cumprod), id[ord])
  weights <- numeric(n)
  weights[ord] <- cum

  res <- list(
    ipw.weights   = weights,
    weights.trunc = if (!is.null(trunc)) .truncate_weights(weights, trunc) else NULL,
    call          = call_obj
  )
  res <- res[!vapply(res, is.null, logical(1))]
  class(res) <- c("ipwprobtm", "list")
  res
}
