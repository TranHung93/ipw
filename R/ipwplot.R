#' Plot Inverse Probability Weights
#'
#' @description
#' Plots the distribution of inverse probability weights. Displays a density plot
#' for point treatment weights or boxplots across follow-up intervals for time-varying weights.
#'
#' @param weights Numerical vector of inverse probability weights.
#' @param timevar Optional numerical vector of follow-up time.
#' @param binwidth Interval width for grouping follow-up time when \code{timevar} is supplied.
#' @param logscale Logical; if \code{TRUE}, weights are plotted on a log scale. Default is \code{TRUE}.
#' @param xlab Horizontal axis label.
#' @param ylab Vertical axis label.
#' @param main Main plot title.
#' @param ref Logical; if \code{TRUE}, adds a reference line at 1 (or log(1) = 0).
#' @param ... Additional graphical parameters passed to \code{boxplot} or \code{plot}.
#' @author Willem M. van der Wal \email{willem@vanderwalresearch.com}, Ronald B. Geskus \email{rgeskus@oucru.org}
#' @export
#'
#' @examples
#' # Point treatment weights density plot
#' set.seed(123)
#' w_point <- runif(500, 0.5, 4)
#' ipwplot(weights = w_point, logscale = FALSE, main = "Density Plot of Point Treatment Weights")
#'
#' # Time-varying weights boxplot across follow-up strata
#' w_tv <- runif(1100, 0.2, 5)
#' time_tv <- rep(seq(0, 500, by = 50), each = 100)
#' ipwplot(weights = w_tv, timevar = time_tv, binwidth = 50, main = "Time-Varying Weights Boxplot")
ipwplot <- function(
    weights,
    timevar = NULL,
    binwidth = NULL,
    logscale = TRUE,
    xlab = NULL,
    ylab = NULL,
    main = "",
    ref = TRUE,
    ...)
{
  if (missing(weights) || !is.numeric(weights) || length(weights) == 0) {
    stop("'weights' must be a non-empty numeric vector.")
  }

  # Derive default axis labels if not provided
  weights_expr <- tryCatch(deparse(substitute(weights)), error = function(e) "weights")
  timevar_expr <- if (!is.null(timevar)) tryCatch(deparse(substitute(timevar)), error = function(e) "timevar") else NULL

  # Time-varying boxplot mode
  if (!is.null(timevar)) {
    if (!is.numeric(timevar) || length(timevar) != length(weights)) {
      stop("'timevar' must be a numeric vector of the same length as 'weights'.")
    }
    if (is.null(binwidth) || !is.numeric(binwidth) || binwidth <= 0) {
      stop("A positive numeric 'binwidth' must be specified when 'timevar' is provided.")
    }

    timevargrp <- round(timevar / binwidth) * binwidth
    if (is.null(xlab)) xlab <- timevar_expr
    if (is.null(ylab)) {
      ylab <- if (logscale) sprintf("log(%s)", weights_expr) else weights_expr
    }

    plot_data <- if (logscale) log(weights) else weights
    ref_val <- if (logscale) 0 else 1

    graphics::boxplot(plot_data ~ timevargrp, pch = 20, xlab = xlab, ylab = ylab, main = main, ...)
    if (ref) graphics::abline(h = ref_val, lty = 2)

    # Point treatment density plot mode
  } else {
    if (is.null(xlab)) {
      xlab <- if (logscale) sprintf("log(%s)", weights_expr) else weights_expr
    }
    if (is.null(ylab)) ylab <- "density"

    plot_data <- if (logscale) log(weights) else weights
    ref_val <- if (logscale) 0 else 1

    dens <- stats::density(plot_data)
    graphics::plot(dens, pch = 20, xlab = xlab, ylab = ylab, main = main, ...)
    if (ref) graphics::abline(v = ref_val, lty = 2)
  }
  invisible(NULL)
}
