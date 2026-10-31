#' Compute Starting Time For Counting Process Notation
#'
#' @description
#' Computes starting time for follow-up intervals in longitudinal data using
#' counting process notation. Within each subject, starting time equals -1 for the
#' initial record and the previous record's follow-up time for subsequent records.
#'
#' @param id Patient or unit identifier (unquoted symbol, string, or vector).
#' @param timevar Follow-up time variable (unquoted symbol, string, or vector).
#' @param data Data frame containing \code{id} and \code{timevar}.
#'
#' @return Numeric vector of starting times in the same row order as \code{data}.
#' @author Willem M. van der Wal \email{willem@vanderwalresearch.com}, Ronald B. Geskus \email{rgeskus@oucru.org}
#' @export
#'
#' @examples
#' mydata1 <- data.frame(
#'   patient = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2),
#'   time.days = c(14, 34, 41, 56, 72, 98, 0, 11, 28, 35)
#' )
#' mydata1$tstart <- tstartfun(patient, time.days, mydata1)
#' mydata1
tstartfun <- function(id, timevar, data) {
  call_obj <- match.call()

  if (missing(data) || !is.data.frame(data)) {
    stop("A valid 'data' frame must be provided.")
  }

  # Extract id and timevar vectors safely
  id_sub <- substitute(id)
  id_name <- if (is.name(id_sub) || is.character(id_sub)) as.character(id_sub) else deparse(id_sub)
  if (id_name %in% names(data)) {
    id_vec <- data[[id_name]]
  } else {
    id_vec <- eval(id_sub, data, parent.frame())
  }

  time_sub <- substitute(timevar)
  time_name <- if (is.name(time_sub) || is.character(time_sub)) as.character(time_sub) else deparse(time_sub)
  if (time_name %in% names(data)) {
    time_vec <- data[[time_name]]
  } else {
    time_vec <- eval(time_sub, data, parent.frame())
  }

  if (anyNA(id_vec) || anyNA(time_vec)) {
    stop("'id' and 'timevar' must not contain missing values.")
  }

  n <- nrow(data)
  orig_order <- seq_len(n)

  # Sort indices by id and timevar
  sort_idx <- order(id_vec, time_vec)
  sorted_id <- id_vec[sort_idx]
  sorted_time <- time_vec[sort_idx]

  # Compute tstart per subject using lag
  calc_subject_tstart <- function(t_seq) {
    if (length(t_seq) == 1) {
      return(-1)
    }
    c(-1, t_seq[-length(t_seq)])
  }

  sorted_tstart <- unsplit(
    lapply(split(sorted_time, sorted_id), calc_subject_tstart),
    sorted_id
  )

  # Restore original row ordering
  res_tstart <- numeric(n)
  res_tstart[sort_idx] <- sorted_tstart

  return(res_tstart)
}
