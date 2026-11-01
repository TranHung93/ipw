#' Bootstrap Standard Errors for Inverse Probability Weighted Analyses
#'
#' @description
#' Bootstraps a whole inverse probability weighted analysis: in each bootstrap
#' sample the weights are re-estimated and the marginal structural model is
#' refitted, so the standard error accounts for the uncertainty in the
#' estimated weights. For data in long format (several rows per subject), whole
#' subjects are resampled and each sampled copy gets a new identifier.
#'
#' @param data Data frame with the original data.
#' @param statistic Function that takes a data frame (a bootstrap sample of
#'   \code{data}) as its first argument and returns a numeric vector of
#'   estimates, for example the coefficients of a weighted outcome model. It
#'   should compute the weights itself, e.g. by calling \code{\link{ipwpoint}}
#'   or \code{\link{ipwtm}}.
#' @param id Optional name of the subject identifier column (unquoted or as a
#'   string). If given, subjects are resampled with all their rows, and in each
#'   bootstrap sample the column is replaced by new identifiers \code{1, 2, ...},
#'   one per sampled copy, so that a subject drawn twice counts as two subjects.
#'   If \code{NULL}, rows are resampled.
#' @param R Number of bootstrap samples. About 200 are usually enough for a
#'   standard error; use more (e.g. 1000) for percentile confidence intervals.
#' @param strata Optional name of a column (unquoted or as a string) defining
#'   strata; resampling is then done within each stratum, keeping the number of
#'   subjects (or rows) per stratum fixed. With \code{id}, the column must be
#'   constant within each subject.
#' @param conf Confidence level of the percentile intervals.
#' @param ... Further arguments passed to \code{statistic}; for the methods,
#'   ignored.
#' @param x,object An object of class \code{"ipwboot"}.
#' @param digits Number of significant digits to print.
#' @param parm Estimates (names or positions) for which to give intervals;
#'   all by default.
#' @param level Confidence level.
#'
#' @details
#' Austin (2016) found in simulations that, for weighted Cox models, naive
#' standard errors were too small and robust (sandwich) standard errors too
#' large, while the bootstrap in which the weights are re-estimated in each
#' sample gave correct standard errors and confidence intervals. Austin studied
#' a point treatment with baseline confounders; for longitudinal data the same
#' idea applies with subjects as the resampling unit.
#'
#' Bootstrap samples in which \code{statistic} fails, or returns a result of a
#' different length, are left out of the summaries; their number is reported.
#' Use \code{set.seed()} before calling \code{ipwboot} to make the results
#' reproducible.
#'
#' @return A list of class \code{"ipwboot"} containing:
#' \item{t0}{Estimates from the original data.}
#' \item{t}{Matrix of bootstrap estimates, one row per bootstrap sample
#'   (\code{NA} for failed samples).}
#' \item{se}{Bootstrap standard errors.}
#' \item{ci}{Percentile confidence intervals.}
#' \item{R, conf}{Number of bootstrap samples and confidence level.}
#' \item{n.fail}{Number of failed bootstrap samples.}
#' \item{errors}{Distinct error messages of failed samples.}
#' \item{id, strata}{Names of the identifier and strata columns, if used.}
#' \item{call}{The original function call.}
#' @references
#' Austin, P.C. (2016). Variance estimation when using inverse probability of
#' treatment weighting (IPTW) with survival analysis. \emph{Statistics in
#' Medicine}, \bold{35}, 5642-5655. \doi{10.1002/sim.7084}.
#' @author Hung Thai Tran \email{hungtt@oucru.org}
#' @seealso \code{\link{ipwpoint}}, \code{\link{ipwtm}}.
#' @export
#'
#' @examples
#' # Point treatment: true causal effect of a on y is 10
#' set.seed(16)
#' n <- 500
#' simdat <- data.frame(l = rnorm(n, 10, 5))
#' simdat$a <- rbinom(n, 1, plogis(simdat$l - 10))
#' simdat$y <- 10 * simdat$a + 0.5 * simdat$l + rnorm(n, -10, 5)
#'
#' msm.effect <- function(d) {
#'   w <- ipwpoint(exposure = a, family = "binomial", link = "logit",
#'                 numerator = ~ 1, denominator = ~ l, data = d)$ipw.weights
#'   coef(lm(y ~ a, data = d, weights = w))["a"]
#' }
#' set.seed(1)
#' res <- ipwboot(simdat, msm.effect, R = 50)  # use about 200 in practice
#' res
#' confint(res)
#'
#' \donttest{
#' # Time-varying treatment: resample patients, re-estimate both weights
#' data(haartdat)
#' haart.effect <- function(d) {
#'   w1 <- ipwtm(exposure = haartind, family = "survival", numerator = ~ sex + age,
#'               denominator = ~ sex + age + cd4.sqrt, id = patient,
#'               tstart = tstart, timevar = fuptime, type = "first", data = d)
#'   w2 <- ipwtm(exposure = dropout, family = "survival", numerator = ~ sex + age,
#'               denominator = ~ sex + age + cd4.sqrt, id = patient,
#'               tstart = tstart, timevar = fuptime, type = "cens", data = d)
#'   fit <- survival::coxph(survival::Surv(tstart, fuptime, event) ~ haartind,
#'                          data = d, weights = w1$ipw.weights * w2$ipw.weights)
#'   coef(fit)
#' }
#' set.seed(2)
#' ipwboot(haartdat, haart.effect, id = patient, R = 20)  # use about 200 in practice
#' }
ipwboot <- function(data, statistic, id = NULL, R = 200, strata = NULL,
                    conf = 0.95, ...) {
  call_obj <- match.call()

  # 1. Validate arguments ---------------------------------------------------
  if (missing(data) || !is.data.frame(data) || nrow(data) == 0) {
    stop("A non-empty 'data' frame must be provided.", call. = FALSE)
  }
  if (missing(statistic) || !is.function(statistic)) {
    stop("'statistic' must be a function that takes a data frame.", call. = FALSE)
  }
  if (!is.numeric(R) || length(R) != 1 || is.na(R) || R < 2 || R != round(R)) {
    stop("'R' must be a whole number of at least 2.", call. = FALSE)
  }
  if (!is.numeric(conf) || length(conf) != 1 || is.na(conf) || conf <= 0 || conf >= 1) {
    stop("'conf' must be a number between 0 and 1.", call. = FALSE)
  }

  id_expr <- substitute(id)
  strata_expr <- substitute(strata)
  id_name <- if (is.null(id_expr)) NULL else .get_col_name(id_expr)
  strata_name <- if (is.null(strata_expr)) NULL else .get_col_name(strata_expr)
  for (v in c(id_name, strata_name)) {
    if (!(v %in% names(data))) stop("Column '", v, "' not found in 'data'.", call. = FALSE)
    if (anyNA(data[[v]])) stop("Column '", v, "' must not contain missing values.", call. = FALSE)
  }

  # 2. Resampling units: subjects (with their rows) or rows -----------------
  if (!is.null(id_name)) {
    unit_ids <- unique(data[[id_name]])
    unit_rows <- split(seq_len(nrow(data)), factor(data[[id_name]], levels = unit_ids))
    names(unit_rows) <- NULL
  } else {
    unit_rows <- NULL
  }
  n_units <- if (is.null(unit_rows)) nrow(data) else length(unit_rows)

  unit_strata <- NULL
  if (!is.null(strata_name)) {
    if (is.null(unit_rows)) {
      unit_strata <- data[[strata_name]]
    } else {
      per_unit <- lapply(unit_rows, function(r) unique(data[[strata_name]][r]))
      if (any(lengths(per_unit) != 1)) {
        stop("'strata' must be constant within each subject.", call. = FALSE)
      }
      unit_strata <- unlist(per_unit)
    }
  }

  draw_units <- function() {
    if (is.null(unit_strata)) {
      return(sample.int(n_units, n_units, replace = TRUE))
    }
    unlist(lapply(sort(unique(unit_strata)), function(s) {
      idx <- which(unit_strata == s)
      idx[sample.int(length(idx), length(idx), replace = TRUE)]
    }), use.names = FALSE)
  }

  make_sample <- function(pick) {
    if (is.null(unit_rows)) {
      d <- data[pick, , drop = FALSE]
    } else {
      rows <- unit_rows[pick]
      d <- data[unlist(rows, use.names = FALSE), , drop = FALSE]
      d[[id_name]] <- rep(seq_along(pick), lengths(rows))
    }
    rownames(d) <- NULL
    d
  }

  # 3. Estimate on the original data ----------------------------------------
  t0 <- statistic(data, ...)
  if (!is.numeric(t0) || length(t0) == 0) {
    stop("'statistic' must return a non-empty numeric vector.", call. = FALSE)
  }
  t0 <- stats::setNames(as.numeric(t0), .ipwboot_names(t0))
  p <- length(t0)

  # 4. Bootstrap ------------------------------------------------------------
  t_mat <- matrix(NA_real_, nrow = R, ncol = p, dimnames = list(NULL, names(t0)))
  errors <- character()
  for (b in seq_len(R)) {
    d <- make_sample(draw_units())
    val <- tryCatch(statistic(d, ...), error = function(e) e)
    if (inherits(val, "error")) {
      errors <- c(errors, conditionMessage(val))
    } else if (!is.numeric(val) || length(val) != p) {
      errors <- c(errors, "'statistic' returned a result of a different length or type.")
    } else {
      t_mat[b, ] <- as.numeric(val)
    }
  }
  n_fail <- sum(!stats::complete.cases(t_mat))
  if (n_fail > 0) {
    warning(n_fail, " of ", R, " bootstrap samples failed and were left out. ",
            "First error: ", errors[1], call. = FALSE)
  }
  if (R - n_fail < 2) {
    stop("Fewer than 2 bootstrap samples succeeded. First error: ", errors[1], call. = FALSE)
  }

  res <- list(
    t0     = t0,
    t      = t_mat,
    se     = apply(t_mat, 2, stats::sd, na.rm = TRUE),
    ci     = .ipwboot_ci(t_mat, conf),
    R      = R,
    conf   = conf,
    n.fail = n_fail,
    errors = unique(errors),
    id     = id_name,
    strata = strata_name,
    call   = call_obj
  )
  class(res) <- c("ipwboot", "list")
  res
}

.ipwboot_names <- function(t0) {
  nm <- names(t0)
  if (is.null(nm)) nm <- if (length(t0) == 1) "estimate" else paste0("t", seq_along(t0))
  nm[is.na(nm) | nm == ""] <- paste0("t", which(is.na(nm) | nm == ""))
  nm
}

.ipwboot_ci <- function(t_mat, level) {
  a <- (1 - level) / 2
  ci <- t(apply(t_mat, 2, stats::quantile, probs = c(a, 1 - a), na.rm = TRUE, names = FALSE))
  dimnames(ci) <- list(colnames(t_mat), paste(format(100 * c(a, 1 - a), trim = TRUE), "%"))
  ci
}

#' @rdname ipwboot
#' @export
print.ipwboot <- function(x, digits = 4, ...) {
  unit <- if (is.null(x$id)) "rows" else paste0("subjects ('", x$id, "')")
  cat("Bootstrap of an inverse probability weighted analysis\n")
  cat(x$R, " bootstrap samples, resampling ", unit,
      if (!is.null(x$strata)) paste0(" within strata ('", x$strata, "')"), "\n", sep = "")
  if (x$n.fail > 0) cat(x$n.fail, " failed samples left out\n", sep = "")
  tab <- cbind(estimate = x$t0, se = x$se, x$ci)
  cat("\n")
  print(signif(tab, digits))
  invisible(x)
}

#' @rdname ipwboot
#' @importFrom stats confint
#' @export
confint.ipwboot <- function(object, parm, level = object$conf, ...) {
  ci <- .ipwboot_ci(object$t, level)
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}
