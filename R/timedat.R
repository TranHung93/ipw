#' HIV: TB and Survival (Longitudinal Measurements)
#'
#' @description
#' A simulated dataset containing time-varying CD4 measurements for 386
#' HIV-positive individuals. Corresponding baseline data, including timing
#' of tuberculosis and death, are available in \code{\link{basdat}}.
#'
#' @docType data
#' @keywords datasets
#' @name timedat
#' @usage data(timedat)
#'
#' @format A data frame with 6291 observations on the following 3 variables:
#' \describe{
#'   \item{id}{Patient ID.}
#'   \item{fuptime}{Follow-up time (days since HIV seroconversion).}
#'   \item{cd4count}{CD4 count measured at \code{fuptime}.}
#' }
#'
#' @details
#' These simulated data are used together with \code{\link{basdat}} in a
#' detailed causal modeling example using inverse probability weighting (IPW).
#' See \code{\link{ipwtm}} for the full example. Data were simulated using
#' the algorithm described in Van der Wal et al. (2009).
#'
#'
#'
#' @author Willem M. van der Wal \email{willem@@vanderwalresearch.com},
#'   Ronald B. Geskus \email{rgeskus@@oucru.org}
#'
#' @references
#' Cole, S.R. & Hernán, M.A. (2008). Constructing inverse probability weights
#' for marginal structural models. \emph{American Journal of Epidemiology},
#' \bold{168}(6), 656-664.
#'
#' Robins, J.M., Hernán, M.A. & Brumback, B.A. (2000). Marginal structural
#' models and causal inference in epidemiology. \emph{Epidemiology}, \bold{11}, 550-560.
#'
#' Van der Wal W.M. & Geskus R.B. (2011). ipw: An R Package for Inverse Probability Weighting.
#' \emph{Journal of Statistical Software}, \bold{43}(13), 1-23. \url{https://doi.org/10.18637/jss.v043.i13}
#'
#' Van der Wal W.M., Prins M., Lumbreras B. & Geskus R.B. (2009). A simple G-computation
#' algorithm to quantify the causal effect of a secondary illness on the progression
#' of a chronic disease. \emph{Statistics in Medicine}, \bold{28}(18), 2325-2337.
#'
#' @seealso \code{\link{basdat}}, \code{\link{haartdat}}, \code{\link{ipwplot}},
#'   \code{\link{ipwpoint}}, \code{\link{ipwtm}}, \code{\link{tstartfun}}
#'
#' @examples
#' # For an example of how to use these longitudinal measurements with basdat, see:
#' # ?ipwtm
"timedat"
