#' IQ, Income and Health
#'
#' @description
#' A simulated dataset containing measurements of IQ, monthly income, and
#' health scores for 1000 individuals.
#'
#' @docType data
#' @keywords datasets
#' @name healthdat
#' @usage data(healthdat)
#'
#' @format A data frame with 1000 rows, with each row corresponding to a
#'   separate individual. The following 4 variables are included:
#' \describe{
#'   \item{id}{Individual ID.}
#'   \item{iq}{IQ score.}
#'   \item{income}{Gross monthly income (EUR).}
#'   \item{health}{Health score (0-100).}
#' }
#'
#' @details
#' In these simulated data, IQ acts as a confounder for the causal effect
#' of income on health. This dataset is primarily used to illustrate
#' point-treatment inverse probability weighting.
#'
#'
#'
#' @author Willem M. van der Wal \email{willem@@vanderwalresearch.com},
#'   Ronald B. Geskus \email{rgeskus@@oucru.org}
#'
#' @references
#' Van der Wal W.M. & Geskus R.B. (2011). ipw: An R Package for Inverse Probability Weighting.
#' \emph{Journal of Statistical Software}, \bold{43}(13), 1-23. \doi{10.18637/jss.v043.i13}
#'
#' @seealso \code{\link{basdat}}, \code{\link{haartdat}}, \code{\link{ipwplot}},
#'   \code{\link{ipwpoint}}, \code{\link{ipwtm}}, \code{\link{timedat}},
#'   \code{\link{tstartfun}}
#'
#' @examples
#' # For an example of how to use this data with ipwpoint, see:
#' # ?ipwpoint
"healthdat"
