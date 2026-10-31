#' Inverse Probability Weights From User-Supplied Probabilities
#'
#' @description
#' Computes inverse probability weights from probabilities (or densities) of
#' the observed exposure that you estimated yourself, for example with a
#' machine-learning model or SuperLearner. Point treatment weights are
#' \code{num / den}; with \code{id} and \code{timevar}, the ratios are
#' multiplied cumulatively over time within each subject, as in
#' \code{\link{ipwtm}}.
#'
#' @param den Numeric vector: for each row, the probability (or density) of the
#'   observed exposure given the covariates of the denominator model.
#' @param num Optional numeric vector of the same length: the probability (or
#'   density) of the observed exposure given the covariates of the numerator
#'   model. If \code{NULL}, unstabilized weights (numerator = 1) are computed.
#' @param id Optional vector of subject identifiers, same length as \code{den}.
#'   When given, weights are cumulated over time within each subject.
#' @param timevar Numeric vector of follow-up times, same length as \code{den};
#'   required when \code{id} is given. It sets the order in which the ratios are
#'   multiplied within a subject.
#' @param trunc Optional truncation percentile between 0 and 0.5 (e.g., 0.01 for
#'   1st/99th percentiles).
#'
#' @details
#' For time-varying weights of \code{type = "first"} or \code{"cens"} (see
#' \code{\link{ipwtm}}), set the probabilities to 1 for the rows after the first
#' switch, so that those rows do not change the cumulative weight.
#'
#' @return A list of class \code{"ipwprob"} containing:
#' \item{ipw.weights}{Numeric vector of inverse probability weights, in the
#'   order of the input.}
#' \item{weights.trunc}{Numeric vector of truncated weights (if \code{trunc} is
#'   specified).}
#' \item{call}{The original function call.}
#' @author Hung Thai Tran \email{hungtt@oucru.org}
#' @seealso \code{\link{ipwpoint}}, \code{\link{ipwtm}}.
#' @export
#'
#' @examples
#' # Point treatment: probabilities from any model, here a logistic regression
#' set.seed(1)
#' d <- data.frame(l = rnorm(200))
#' d$a <- rbinom(200, 1, plogis(d$l))
#' p <- fitted(glm(a ~ l, family = binomial, data = d))
#' w <- ipwprob(den = ifelse(d$a == 1, p, 1 - p),
#'              num = ifelse(d$a == 1, mean(d$a), 1 - mean(d$a)))
#' summary(w$ipw.weights)
#'
#' # Time-varying: the ratios are multiplied within each subject over time
#' long <- data.frame(id = rep(1:3, each = 2), time = rep(1:2, 3),
#'                    p = c(0.5, 0.8, 0.4, 0.9, 0.6, 0.7))
#' ipwprob(den = long$p, id = long$id, timevar = long$time)$ipw.weights
ipwprob <- function(den, num = NULL, id = NULL, timevar = NULL, trunc = NULL) {
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
  .check_trunc(trunc)

  ratio <- num / den

  if (!is.null(id)) {
    if (length(id) != n) stop("'id' and 'den' must have the same length.", call. = FALSE)
    if (is.null(timevar)) stop("'timevar' is required when 'id' is given.", call. = FALSE)
    if (!is.numeric(timevar) || length(timevar) != n) {
      stop("'timevar' must be numeric and have the same length as 'den'.", call. = FALSE)
    }
    if (anyNA(id) || anyNA(timevar)) {
      stop("'id' and 'timevar' must not contain missing values.", call. = FALSE)
    }
    if (anyDuplicated(data.frame(id, timevar))) {
      stop("Each combination of 'id' and 'timevar' must occur only once.", call. = FALSE)
    }
    ord <- order(id, timevar)
    cum <- unsplit(lapply(split(ratio[ord], id[ord]), cumprod), id[ord])
    weights <- numeric(n)
    weights[ord] <- cum
  } else {
    if (!is.null(timevar)) stop("'id' is required when 'timevar' is given.", call. = FALSE)
    weights <- ratio
  }

  res <- list(
    ipw.weights   = weights,
    weights.trunc = if (!is.null(trunc)) .truncate_weights(weights, trunc) else NULL,
    call          = call_obj
  )
  res <- res[!vapply(res, is.null, logical(1))]
  class(res) <- c("ipwprob", "list")
  res
}
