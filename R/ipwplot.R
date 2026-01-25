#' Plot Inverse Probability Weights
#' @description
#' For time varying weights: display boxplots within strata of follow-up time. For point treatment weights: display density plot.
#'
#' @param weights numerical vector of inverse probability weights to plot.
#' @param timevar numerical vector representing follow-up time. When specified, boxplots within strata of follow-up time are displayed. When left unspecified, a density plot is displayed.
#' @param binwidth numerical value indicating the width of the intervals of follow-up time; for each interval a boxplot is made. Ignored when \code{timevar} is not specified.
#' @param logscale logical value. If \code{TRUE}, weights are plotted on a logarithmic scale.
#' @param xlab label for the horizontal axis.
#' @param ylab label for the vertical axis.
#' @param main main title for the plot.
#' @param ref logical value. If \code{TRUE}, a reference line is plotted at \code{y=1}.
#' @param ... additional arguments passed to \code{\link{boxplot}} (when \code{timevar} is specified) or \code{\link{plot}} (when \code{timevar} is not specified).
#'
#' @return A plot is displayed.
#' @author Willem M. van der Wal \email{willem@vanderwalresearch.com}, Ronald B. Geskus \email{rgeskus@oucru.org}
#' @references Van der Wal W.M. & Geskus R.B. (2011). ipw: An R Package for Inverse  Probability Weighting. \emph{Journal of Statistical Software}, \bold{43}(13), 1-23. \url{https://doi.org/10.18637/jss.v043.i13}
#' @seealso \code{\link{basdat}}, \code{\link{haartdat}}, \code{\link{ipwplot}}, \code{\link{ipwpoint}}, \code{\link{ipwtm}}, \code{\link{timedat}}, \code{\link{tstartfun}}.
#' @keywords hplot
#'
#' @export
#'
#' @examples #see ?ipwpoint and ?ipwtm for examples
ipwplot <- function(weights, timevar = NULL, binwidth = NULL, logscale = TRUE, xlab = NULL, ylab = NULL, main = "", ref = TRUE, ...){
  if (!is.null(timevar)){
    timevargrp <- round(timevar/binwidth)*binwidth
    if (is.null(xlab)) xlab <- deparse(match.call()$timevar, width.cutoff = 500)
    if (is.null(ylab) & logscale == FALSE) ylab <- deparse(match.call()$weights, width.cutoff = 500)
    if (is.null(ylab) & logscale == TRUE) ylab <- paste("log(", deparse(match.call()$weights, width.cutoff = 500), ")", sep = "")
    if(logscale == TRUE){
      boxplot(log(weights) ~ timevargrp, pch = 20, xlab = xlab, ylab = ylab, main = main, ...)
      if (ref == TRUE) abline(h = log(1), lty = 2)
    }
    if(logscale == FALSE){
      boxplot(weights ~ timevargrp, pch = 20, xlab = xlab, ylab = ylab, main = main, ...)
      if (ref == TRUE) abline(h = 1, lty = 2)
    }
  }
  if (is.null(timevar)){
    if (is.null(xlab) & logscale == FALSE) xlab <- deparse(match.call()$weights, width.cutoff = 500)
    if (is.null(xlab) & logscale == TRUE) xlab <- paste("log(", deparse(match.call()$weights, width.cutoff = 500), ")", sep = "")
    if (is.null(ylab)) ylab <- "density"
    if(logscale == TRUE){
      plot(density(log(weights)), pch = 20, xlab = xlab, ylab = ylab, main = main, ...)
      if (ref == TRUE) abline(v = log(1), lty = 2)
    }
    if(logscale == FALSE){
      plot(density(weights), pch = 20, xlab = xlab, ylab = ylab, main = main, ...)
      if (ref == TRUE) abline(v = 1, lty = 2)
    }
  }
}
