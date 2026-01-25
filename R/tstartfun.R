#' Compute Starting Time For Counting Process Notation
#' @description
#' Function to compute starting time for intervals of follow-up, when using the counting process notation. Within each unit under observation (usually individuals), computes starting time equal to:
#' \itemize{
#'   \item time of previous record when there is a previous record.
#'   \item -1 for first record.
#' }
#'
#' @param id numerical vector, uniquely identifying the units under observation, within which the longitudinal measurements are taken.
#' @param timevar numerical vector, representing follow-up time, starting at 0.
#' @param data dataframe containing \code{id} and \code{timevar}
#'
#' @return Numerical vector containing starting time for each record. In the same order as the records in \code{data}, to facilitate merging.
#' @section Missing values:
#' Currently, \code{id} and \code{timevar} should not contain missing values.
#'
#' @author Willem M. van der Wal \email{willem@vanderwalresearch.com}, Ronald B. Geskus \email{rgeskus@oucru.org}
#' @references
#' Van der Wal W.M. & Geskus R.B. (2011). ipw: An R Package for Inverse  Probability Weighting. \emph{Journal of Statistical Software}, \bold{43}(13), 1-23. \url{https://doi.org/10.18637/jss.v043.i13}.
#' @seealso \code{\link{basdat}}, \code{\link{haartdat}}, \code{\link{ipwplot}}, \code{\link{ipwpoint}}, \code{\link{ipwtm}}, \code{\link{timedat}}, \code{\link{tstartfun}}.
#' @export
#'
#' @examples
#' #data
#' mydata1 <- data.frame(
#'   patient = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2),
#'   time.days = c(14, 34, 41, 56, 72, 98, 0, 11, 28, 35))
#'
#' #compute starting time for each interval
#' mydata1$tstart <- tstartfun(patient, time.days, mydata1)
#'
#' #result
#' mydata1
#'
#' #see also ?ipwtm for example
#'
#' @keywords methods survival
tstartfun <- function(
    id,
    timevar,
    data)
{
  #save input
  tempcall <- match.call()
  #record original order of dataframe so that the output can be returned in the same order
  order.orig <- 1:nrow(data)
  order.orig <- order.orig[order(
    eval(parse(text = paste("data$", deparse(tempcall$id, width.cutoff = 500), sep = ""))),
    eval(parse(text = paste("data$", deparse(tempcall$timevar, width.cutoff = 500), sep = "")))
  )] #sort as below
  #sort dataframe on follow-up time within each individual, necessary for cumulative products below
  data <- data[order(
    eval(parse(text = paste("data$", deparse(tempcall$id, width.cutoff = 500), sep = ""))),
    eval(parse(text = paste("data$", deparse(tempcall$timevar, width.cutoff = 500), sep = "")))
  ),]
  #make new dataframe for newly computed variables, to prevent variable name conflicts
  tempdat <- data.frame(
    id = data[,as.character(tempcall$id)],
    timevar = data[,as.character(tempcall$timevar)]
  )
  #compute tstart
  #time of previous record when there is a previous record
  #-1 for first record
  tempfun <- function(x) {
    tempdif <- diff(c(min(x), x))
    tempdif[tempdif == 0] <- min(x) + 1
    return(x - tempdif)
  }
  tempdat$tstart <- unsplit(lapply(split(tempdat$timevar, tempdat$id), function(x)tempfun(x)), tempdat$id)
  #return results in the same order as the original input dataframe
  return(tstart = tempdat$tstart[order(order.orig)])
}
