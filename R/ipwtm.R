#' Estimate Inverse Probability Weights (Time Varying)
#'
#' @description
#' Estimates inverse probability weights to fit marginal structural models with
#' time-varying exposures and time-varying confounders across follow-up intervals.
#' Supports binomial, survival, multinomial, ordinal, and gaussian families.
#'
#' @import MASS nnet survival geepack graphics methods stats
#' @param exposure Exposure variable (unquoted symbol, string, or vector).
#' @param family Model family: \code{"binomial"}, \code{"survival"}, \code{"multinomial"},
#'   \code{"ordinal"}, or \code{"gaussian"}.
#' @param link Link function for binomial or ordinal models.
#' @param numerator Formula for numerator model. If \code{NULL}, unstabilized weights are computed.
#' @param denominator Formula for denominator model.
#' @param id Patient/subject identifier variable.
#' @param tstart Start time variable for counting process format (required for \code{family = "survival"}).
#' @param timevar End time/follow-up time variable.
#' @param type Exposure switch type: \code{"first"}, \code{"cens"}, or \code{"all"}.
#' @param data Data frame containing all required variables.
#' @param corstr Correlation structure for GEE when \code{family = "gaussian"} (default: \code{"ar1"}).
#' @param trunc Optional truncation percentile between 0 and 0.5.
#' @param ... Additional arguments passed to model fitting engines.
#'
#' @return A list of class \code{"ipwtm"} containing:
#' \item{ipw.weights}{Numeric vector of inverse probability weights in original row order.}
#' \item{weights.trunc}{Numeric vector of truncated weights (if \code{trunc} is specified).}
#' \item{call}{Original function call.}
#' \item{selvar}{Selection indicator vector (1 = active observation, 0 = post-switch).}
#' \item{num.mod}{Numerator model object.}
#' \item{den.mod}{Denominator model object.}
#' @export
#'
#'
#' @examples
#' ########################################################################
#' #EXAMPLE 1
#'
#' #Load longitudinal data from HIV positive individuals.
#' data(haartdat)
#'
#' #CD4 is confounder for the effect of initiation of HAART therapy on mortality.
#' #Estimate inverse probability weights to correct for confounding.
#' #Exposure allocation model is Cox proportional hazards model.
#' temp <- ipwtm(
#'   exposure = haartind,
#'   family = "survival",
#'   numerator = ~ sex + age,
#'   denominator = ~ sex + age + cd4.sqrt,
#'   id = patient,
#'   tstart = tstart,
#'   timevar = fuptime,
#'   type = "first",
#'   data = haartdat)
#'
#' #plot inverse probability weights
#' graphics.off()
#' ipwplot(weights = temp$ipw.weights, timevar = haartdat$fuptime,
#'         binwidth = 100, ylim = c(-1.5, 1.5), main = "Stabilized inverse probability weights")
#'
#' #CD4 count has an effect both on dropout and mortality, which causes informative censoring.
#' #Use inverse probability of censoring weighting to correct for effect of CD4 on dropout.
#' #Use Cox proportional hazards model for dropout.
#' temp2 <- ipwtm(
#'   exposure = dropout,
#'   family = "survival",
#'   numerator = ~ sex + age,
#'   denominator = ~ sex + age + cd4.sqrt,
#'   id = patient,
#'   tstart = tstart,
#'   timevar = fuptime,
#'   type = "cens",
#'   data = haartdat)
#'
#' #plot inverse probability of censoring weights
#' graphics.off()
#' ipwplot(weights = temp2$ipw.weights, timevar = haartdat$fuptime,
#'         binwidth = 100, ylim = c(-1.5, 1.5),
#'         main = "Stabilized inverse probability of censoring weights")
#'
#' #MSM for the causal effect of initiation of HAART on mortality.
#' #Corrected both for confounding and informative censoring.
#' #With robust standard error obtained using cluster().
#' require(survival)
#' summary(coxph(Surv(tstart, fuptime, event) ~ haartind + cluster(patient),
#'               data = haartdat, weights = temp$ipw.weights*temp2$ipw.weights))
#'
#' #uncorrected model
#' summary(coxph(Surv(tstart, fuptime, event) ~ haartind, data = haartdat))
#'
#' ########################################################################
#' #EXAMPLE 2
#'
#' data(basdat)
#' data(timedat)
#'
#' #Aim: to model the causal effect of active tuberculosis (TB) on mortality.
#' #Longitudinal CD4 is a confounder as well as intermediate for the effect of TB.
#'
#' #process original measurements
#' #check for ties (not allowed)
#' table(duplicated(timedat[,c("id", "fuptime")]))
#' #take square root of CD4 because of skewness
#' timedat$cd4.sqrt <- sqrt(timedat$cd4count)
#' #add TB time to dataframe
#' timedat <- merge(timedat, basdat[,c("id", "Ttb")], by = "id", all.x = TRUE)
#' #compute TB status
#' timedat$tb.lag <- ifelse(with(timedat, !is.na(Ttb) & fuptime > Ttb), 1, 0)
#' #longitudinal CD4-model
#' require(nlme)
#' cd4.lme <- lme(cd4.sqrt ~ fuptime + tb.lag, random = ~ fuptime | id,
#'                data = timedat)
#'
#' #build new dataset:
#' #rows corresponding to TB-status switches, and individual end times
#' times <- sort(unique(c(basdat$Ttb, basdat$Tend)))
#' startstop <- data.frame(
#'   id = rep(basdat$id, each = length(times)),
#'   fuptime = rep(times, nrow(basdat)))
#' #add baseline data to dataframe
#' startstop <- merge(startstop, basdat, by = "id", all.x = TRUE)
#' #limit individual follow-up using Tend
#' startstop <- startstop[with(startstop, fuptime <= Tend),]
#' startstop$tstart <- tstartfun(id, fuptime, startstop) #compute tstart (?tstartfun)
#' #indicate TB status
#' startstop$tb <- ifelse(with(startstop, !is.na(Ttb) & fuptime >= Ttb), 1, 0)
#' #indicate TB status at previous time point
#' startstop$tb.lag <- ifelse(with(startstop, !is.na(Ttb) & fuptime > Ttb), 1, 0)
#' #indicate death
#' startstop$event <- ifelse(with(startstop, !is.na(Tdeath) & fuptime >= Tdeath),
#'                           1, 0)
#' #impute CD4, based on TB status at previous time point.
#' startstop$cd4.sqrt <- predict(cd4.lme,
#'                               newdata = data.frame(id = startstop$id,
#'                                                    fuptime = startstop$fuptime,
#'                                                    tb.lag = startstop$tb.lag))
#' #compute inverse probability weights
#' temp <- ipwtm(
#'   exposure = tb,
#'   family = "survival",
#'   numerator = ~ 1,
#'   denominator = ~ cd4.sqrt,
#'   id = id,
#'   tstart = tstart,
#'   timevar = fuptime,
#'   type = "first",
#'   data = startstop)
#' summary(temp$ipw.weights)
#' ipwplot(weights = temp$ipw.weights, timevar = startstop$fuptime, binwidth = 100)
#'
#' #models
#' #IPW-fitted MSM, using cluster() to obtain robust standard error estimate
#' require(survival)
#' summary(coxph(Surv(tstart, fuptime, event) ~ tb + cluster(id),
#'               data = startstop, weights = temp$ipw.weights))
#' #unadjusted
#' summary(coxph(Surv(tstart, fuptime, event) ~ tb, data = startstop))
#' #adjusted using conditioning: part of the effect of TB is adjusted away
#' summary(coxph(Surv(tstart, fuptime, event) ~ tb + cd4.sqrt, data = startstop))
#'
#' ## Not run:
#' #compute bootstrap CI for TB parameter (takes a few hours)
#' #taking into account the uncertainty introduced by modelling longitudinal CD4
#' #taking into account the uncertainty introduced by estimating the inverse probability weights
#' #robust with regard to weights unequal to 1
#' #  require(boot)
#' #  boot.fun <- function(data, index, data.tm){
#' #     data.samp <- data[index,]
#' #     data.samp$id.samp <- 1:nrow(data.samp)
#' #     data.tm.samp <- do.call("rbind", lapply(data.samp$id.samp, function(id.samp) {
#' #       cbind(data.tm[data.tm$id == data.samp$id[data.samp$id.samp == id.samp],],
#' #         id.samp = id.samp)
#' #       }
#' #     ))
#' #     cd4.lme <- lme(cd4.sqrt ~ fuptime + tb.lag, random = ~ fuptime | id.samp, data = data.tm.samp)
#' #     times <- sort(unique(c(data.samp$Ttb, data.samp$Tend)))
#' #     startstop.samp <- data.frame(id.samp = rep(data.samp$id.samp, each = length(times)),
#' #                                  fuptime = rep(times, nrow(data.samp)))
#' #     startstop.samp <- merge(startstop.samp, data.samp, by = "id.samp", all.x = TRUE)
#' #     startstop.samp <- startstop.samp[with(startstop.samp, fuptime <= Tend),]
#' #     startstop.samp$tstart <- tstartfun(id.samp, fuptime, startstop.samp)
#' #     startstop.samp$tb <- ifelse(with(startstop.samp, !is.na(Ttb) & fuptime >= Ttb), 1, 0)
#' #     startstop.samp$tb.lag <- ifelse(with(startstop.samp, !is.na(Ttb) & fuptime > Ttb), 1, 0)
#' #     startstop.samp$event <- ifelse(with(startstop.samp, !is.na(Tdeath) & fuptime >= Tdeath), 1, 0)
#' #     startstop.samp$cd4.sqrt <- predict(cd4.lme, newdata = data.frame(id.samp =
#' #       startstop.samp$id.samp, fuptime = startstop.samp$fuptime, tb.lag = startstop.samp$tb.lag))
#' #
#' #     return(coef(coxph(Surv(tstart, fuptime, event) ~ tb, data = startstop.samp,
#' #        weights = ipwtm(
#' #             exposure = tb,
#' #             family = "survival",
#' #             numerator = ~ 1,
#' #             denominator = ~ cd4.sqrt,
#' #             id = id.samp,
#' #             tstart = tstart,
#' #             timevar = fuptime,
#' #             type = "first",
#' #             data = startstop.samp)$ipw.weights))[1])
#' #     }
#' #  bootres <- boot(data = basdat, statistic = boot.fun, R = 999, data.tm = timedat)
#' #  bootres
#' #  boot.ci(bootres, type = "basic")
#' #
#' ## End(Not run)
#' @keywords htest models
ipwtm <- function(
    exposure,
    family = c("binomial", "survival", "multinomial", "ordinal", "gaussian"),
    link = NULL,
    numerator = NULL,
    denominator,
    id,
    tstart,
    timevar,
    type = c("first", "cens", "all"),
    data,
    corstr = "ar1",
    trunc = NULL,
    ...) {
  call_obj <- match.call()

  # 1. Capture Column Identifiers Without Caller Evaluation -----------------
  var_names <- list(
    exposure = if (!missing(exposure)) .get_col_name(substitute(exposure)) else NULL,
    id       = if (!missing(id)) .get_col_name(substitute(id)) else NULL,
    timevar  = if (!missing(timevar)) .get_col_name(substitute(timevar)) else NULL,
    tstart   = if (!missing(tstart)) .get_col_name(substitute(tstart)) else NULL
  )

  # 2. Input Validation -----------------------------------------------------
  family <- match.arg(family)
  type <- match.arg(type)

  .validate_inputs(
    family      = family,
    type        = type,
    link        = link,
    numerator   = numerator,
    denominator = denominator,
    var_names   = var_names,
    data        = data,
    trunc       = trunc
  )

  # 3. Sorting & Preserving Original Order ----------------------------------
  ord <- order(data[[var_names$id]], data[[var_names$timevar]])
  order_orig <- seq_len(nrow(data))[ord]
  sorted_data <- data[ord, , drop = FALSE]

  # 4. Selection Variable Generation (selvar) -------------------------------
  selvar <- .compute_selvar(
    exposure = sorted_data[[var_names$exposure]],
    id       = sorted_data[[var_names$id]],
    type     = type,
    family   = family
  )
  sorted_data$selvar <- selvar

  # 5. Model Fitting & Probability Estimation -------------------------------
  engine_args <- list(
    data        = sorted_data,
    var_names   = var_names,
    numerator   = numerator,
    denominator = denominator,
    type        = type,
    link        = link,
    corstr      = corstr,
    dots        = list(...)
  )

  res_engine <- switch(family,
                       "binomial"    = do.call(.fit_ipw_binomial, engine_args),
                       "survival"    = do.call(.fit_ipw_survival, engine_args),
                       "multinomial" = do.call(.fit_ipw_categorical, c(engine_args, list(is_ordinal = FALSE))),
                       "ordinal"     = do.call(.fit_ipw_categorical, c(engine_args, list(is_ordinal = TRUE))),
                       "gaussian"    = do.call(.fit_ipw_gaussian, engine_args)
  )

  # 6. Cumulative Weights & Safety Checks -----------------------------------
  weights_sorted <- res_engine$weights
  if (anyNA(weights_sorted) || any(is.infinite(weights_sorted))) {
    stop("NA's in weights!", call. = FALSE)
  }

  # 7. Weight Truncation ----------------------------------------------------
  trunc_sorted <- NULL
  if (!is.null(trunc)) {
    q_lo <- stats::quantile(weights_sorted, probs = 0 + trunc, na.rm = TRUE)
    q_hi <- stats::quantile(weights_sorted, probs = 1 - trunc, na.rm = TRUE)
    trunc_sorted <- pmin(pmax(weights_sorted, q_lo), q_hi)
  }

  # 8. Assemble Return Object in Original Row Order -------------------------
  orig_idx <- order(order_orig)
  res <- list(
    ipw.weights   = weights_sorted[orig_idx],
    weights.trunc = if (!is.null(trunc)) trunc_sorted[orig_idx] else NULL,
    call          = call_obj,
    selvar        = selvar[orig_idx],
    num.mod       = res_engine$num_mod,
    den.mod       = res_engine$den_mod
  )

  res <- res[!vapply(res, is.null, logical(1))]
  class(res) <- c("ipwtm", "list")
  res
}


# ==========================================================================
# Input Validation & Helper Functions
# ==========================================================================

.validate_inputs <- function(family, type, link, numerator, denominator, var_names, data, trunc) {
  if (missing(data) || !is.data.frame(data)) {
    stop("No data specified", call. = FALSE)
  }
  if (is.null(var_names$exposure) || !(var_names$exposure %in% names(data))) {
    stop("No exposure variable specified", call. = FALSE)
  }
  if (is.null(var_names$id) || !(var_names$id %in% names(data))) {
    stop("No patient id specified", call. = FALSE)
  }
  if (is.null(var_names$timevar) || !(var_names$timevar %in% names(data))) {
    stop("No timevar specified", call. = FALSE)
  }
  if (family == "survival") {
    if (is.null(var_names$tstart) || !(var_names$tstart %in% names(data))) {
      stop("No tstart specified, is necessary for family = \"survival\"", call. = FALSE)
    }
  }

  if (missing(denominator) || !inherits(denominator, "formula")) {
    stop("No denominator model specified", call. = FALSE)
  }
  if (!is.null(numerator) && !inherits(numerator, "formula")) {
    stop("Invalid numerator formula specified", call. = FALSE)
  }
  if (family == "binomial") {
    if (is.null(link) || !(link %in% c("logit", "probit", "cauchit", "log", "cloglog"))) {
      stop("No valid link function specified for family = binomial (\"logit\", \"probit\", \"cauchit\", \"log\", \"cloglog\")", call. = FALSE)
    }
  }
  if (family == "ordinal") {
    if (is.null(link) || !(link %in% c("logit", "probit", "cauchit", "cloglog"))) {
      stop("No valid link function specified for family = ordinal (\"logit\", \"probit\", \"cauchit\", \"cloglog\")", call. = FALSE)
    }
  }
  if (family %in% c("survival", "multinomial", "ordinal") && type == "all") {
    stop(paste0("Type \"all\" not yet implemented for family = \"", family, "\""), call. = FALSE)
  }
  if (family %in% c("multinomial", "ordinal", "gaussian") && type == "cens") {
    stop(paste0("Type \"cens\" not yet implemented for family = \"", family, "\""), call. = FALSE)
  }
  if (family == "gaussian") {
    if (type == "first") {
      stop("Type \"first\" not implemented for family = \"gaussian\"", call. = FALSE)
    }
    if (is.null(numerator)) {
      stop("Numerator necessary for family = \"gaussian\"", call. = FALSE)
    }
  }
  if (!is.null(trunc) && (!is.numeric(trunc) || length(trunc) != 1 || trunc < 0 || trunc > 0.5)) {
    stop("Invalid truncation percentage specified (0-0.5)", call. = FALSE)
  }
}

.get_col_name <- function(expr) {
  if (is.symbol(expr) || is.character(expr)) as.character(expr) else deparse(expr, width.cutoff = 500)
}

.build_model_formula <- function(lhs_str, rhs_formula) {
  rhs_terms <- attr(stats::terms(rhs_formula), "term.labels")
  intercept <- attr(stats::terms(rhs_formula), "intercept") == 1
  rhs_str <- if (length(rhs_terms) == 0) {
    if (intercept) "1" else "0"
  } else {
    paste(rhs_terms, collapse = " + ")
  }
  if (length(rhs_terms) > 0 && !intercept) rhs_str <- paste(rhs_str, "- 1")
  stats::as.formula(paste(lhs_str, "~", rhs_str), env = environment(rhs_formula))
}

.compute_selvar <- function(exposure, id, type, family) {
  if (type == "all") {
    return(rep(1, length(exposure)))
  }

  if (family %in% c("binomial", "survival")) {
    unsplit(lapply(split(exposure, id), function(x) {
      idx <- match(1, x)
      if (!is.na(idx)) c(rep(1, idx), rep(0, length(x) - idx)) else rep(1, length(x))
    }), id)
  } else {
    unsplit(lapply(split(exposure, id), function(x) {
      base_lvl <- sort(unique(x))[1]
      switches <- which(x != base_lvl)
      idx <- if (length(switches) > 0) switches[1] else NA_integer_
      if (!is.na(idx)) c(rep(1, idx), rep(0, length(x) - idx)) else rep(1, length(x))
    }), id)
  }
}


# ==========================================================================
# Modeling Engines
# ==========================================================================

# --- Binomial Engine -----------------------------------------------------
.fit_ipw_binomial <- function(data, var_names, numerator, denominator, type, link, dots, ...) {
  fam_link <- switch(link,
                     "logit"   = stats::binomial(link = "logit"),
                     "probit"  = stats::binomial(link = "probit"),
                     "cauchit" = stats::binomial(link = "cauchit"),
                     "log"     = stats::binomial(link = "log"),
                     "cloglog" = stats::binomial(link = "cloglog")
  )

  exp_col <- var_names$exposure
  id_col  <- var_names$id

  # Denominator model
  form_den <- .build_model_formula(exp_col, denominator)
  mod_den <- do.call(stats::glm, c(list(
    formula   = form_den,
    family    = fam_link,
    data      = data,
    subset    = if (type != "all") quote(selvar == 1) else NULL,
    na.action = stats::na.fail
  ), dots))

  pred_den <- stats::predict.glm(mod_den, type = "response")
  p_den <- rep(1, nrow(data))
  active <- (data$selvar == 1)

  if (type == "all") {
    p_den[data[[exp_col]] == 0] <- 1 - pred_den[data[[exp_col]] == 0]
    p_den[data[[exp_col]] == 1] <- pred_den[data[[exp_col]] == 1]
    w_den <- unlist(lapply(split(p_den, data[[id_col]]), cumprod))
  } else {
    p_den[active & data[[exp_col]] == 0] <- 1 - pred_den[data[[exp_col]][active] == 0]
    if (type == "first") {
      p_den[active & data[[exp_col]] == 1] <- pred_den[data[[exp_col]][active] == 1]
    } else if (type == "cens") {
      p_den[active & data[[exp_col]] == 1] <- 1 - pred_den[data[[exp_col]][active] == 1]
    }
    p_den[!active] <- 1
    w_den <- unlist(lapply(split(p_den, data[[id_col]]), cumprod))
  }

  # Numerator model
  mod_num <- NULL
  if (is.null(numerator)) {
    w_num <- rep(1, nrow(data))
  } else {
    form_num <- .build_model_formula(exp_col, numerator)
    mod_num <- do.call(stats::glm, c(list(
      formula   = form_num,
      family    = fam_link,
      data      = data,
      subset    = if (type != "all") quote(selvar == 1) else NULL,
      na.action = stats::na.fail
    ), dots))

    pred_num <- stats::predict.glm(mod_num, type = "response")
    p_num <- rep(1, nrow(data))

    if (type == "all") {
      p_num[data[[exp_col]] == 0] <- 1 - pred_num[data[[exp_col]] == 0]
      p_num[data[[exp_col]] == 1] <- pred_num[data[[exp_col]] == 1]
      w_num <- unlist(lapply(split(p_num, data[[id_col]]), cumprod))
    } else {
      p_num[active & data[[exp_col]] == 0] <- 1 - pred_num[data[[exp_col]][active] == 0]
      if (type == "first") {
        p_num[active & data[[exp_col]] == 1] <- pred_num[data[[exp_col]][active] == 1]
      } else if (type == "cens") {
        p_num[active & data[[exp_col]] == 1] <- 1 - pred_num[data[[exp_col]][active] == 1]
      }
      p_num[!active] <- 1
      w_num <- unlist(lapply(split(p_num, data[[id_col]]), cumprod))
    }
  }

  list(weights = w_num / w_den, num_mod = mod_num, den_mod = mod_den)
}

# --- Survival Engine -----------------------------------------------------
.fit_ipw_survival <- function(data, var_names, numerator, denominator, type, dots, ...) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required for family = 'survival'.", call. = FALSE)
  }

  exp_col   <- var_names$exposure
  id_col    <- var_names$id
  time_col  <- var_names$timevar
  surv_lhs  <- sprintf("survival::Surv(%s, %s, %s)", var_names$tstart, time_col, exp_col)
  u_times   <- sort(unique(data[[time_col]]))

  .compute_surv_p <- function(mod) {
    bh <- survival::basehaz(mod, centered = TRUE)
    df_bh <- data.frame(timevar = u_times)
    df_bh <- merge(df_bh, data.frame(timevar = bh$time, cum_haz = bh$hazard), by = "timevar", all.x = TRUE)

    if (is.na(df_bh$cum_haz[1])) df_bh$cum_haz[1] <- 0
    df_bh$cum_haz  <- stats::approx(df_bh$timevar, df_bh$cum_haz, xout = df_bh$timevar, method = "constant", rule = 2)$y
    df_bh$step_haz <- c(df_bh$cum_haz[1], diff(df_bh$cum_haz, 1))

    merged <- merge(
      data.frame(orig_seq = seq_len(nrow(data)), timevar = data[[time_col]]),
      df_bh[, c("timevar", "step_haz")],
      by = "timevar",
      all.x = TRUE
    )
    merged <- merged[order(merged$orig_seq), ]

    risk <- rep(NA_real_, nrow(data))
    risk[data$selvar == 1] <- stats::predict(mod, type = "risk", centered = TRUE)
    inst_haz <- merged$step_haz * risk

    p <- rep(1, nrow(data))
    s0 <- which(data$selvar == 1 & data[[exp_col]] == 0)
    p[s0] <- exp(-1 * inst_haz[s0])

    s1 <- which(data$selvar == 1 & data[[exp_col]] == 1)
    if (type == "first") {
      p[s1] <- 1 - exp(-1 * inst_haz[s1])
    } else if (type == "cens") {
      p[s1] <- exp(-1 * inst_haz[s1])
    }
    p[data$selvar == 0] <- 1
    p
  }

  # Denominator model
  form_den <- .build_model_formula(surv_lhs, denominator)
  mod_den <- do.call(survival::coxph, c(list(
    formula   = form_den,
    data      = data,
    subset    = quote(selvar == 1),
    na.action = stats::na.fail,
    method    = "efron",
    model     = TRUE
  ), dots))
  p_den <- .compute_surv_p(mod_den)
  w_den <- unsplit(lapply(split(p_den, data[[id_col]]), cumprod), data[[id_col]])

  # Numerator model
  mod_num <- NULL
  if (is.null(numerator)) {
    w_num <- rep(1, nrow(data))
  } else {
    form_num <- .build_model_formula(surv_lhs, numerator)
    mod_num <- do.call(survival::coxph, c(list(
      formula   = form_num,
      data      = data,
      subset    = quote(selvar == 1),
      na.action = stats::na.fail,
      method    = "efron",
      model     = TRUE
    ), dots))
    p_num <- .compute_surv_p(mod_num)
    w_num <- unsplit(lapply(split(p_num, data[[id_col]]), cumprod), data[[id_col]])
  }

  list(weights = w_num / w_den, num_mod = mod_num, den_mod = mod_den)
}

# --- Multinomial & Ordinal Engine ----------------------------------------
.fit_ipw_categorical <- function(data, var_names, numerator, denominator, is_ordinal, link, dots, ...) {
  exp_col <- var_names$exposure
  id_col  <- var_names$id
  u_exp   <- sort(unique(data[[exp_col]]))

  lhs_str <- if (is_ordinal) sprintf("as.factor(%s)", exp_col) else exp_col
  polr_method <- if (is_ordinal) {
    switch(link, logit = "logistic", probit = "probit", cloglog = "cloglog", cauchit = "cauchit")
  } else NULL

  .fit_and_extract_p <- function(rhs_form) {
    form <- .build_model_formula(lhs_str, rhs_form)
    mod <- if (is_ordinal) {
      if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' required.", call. = FALSE)
      do.call(MASS::polr, c(list(
        formula   = form,
        data      = data,
        method    = polr_method,
        subset    = quote(selvar == 1),
        na.action = stats::na.fail,
        model     = TRUE
      ), dots))
    } else {
      if (!requireNamespace("nnet", quietly = TRUE)) stop("Package 'nnet' required.", call. = FALSE)
      do.call(nnet::multinom, c(list(
        formula   = form,
        data      = data,
        subset    = quote(selvar == 1),
        na.action = stats::na.fail,
        trace     = FALSE
      ), dots))
    }

    probs <- as.data.frame(stats::predict(mod, type = "probs"))
    p_vec <- rep(1, nrow(data))
    for (i in seq_along(u_exp)) {
      match_rows <- which(data$selvar == 1 & data[[exp_col]] == u_exp[i])
      p_vec[match_rows] <- probs[data[[exp_col]][data$selvar == 1] == u_exp[i], i]
    }
    p_vec[data$selvar == 0] <- 1
    list(model = mod, p = p_vec)
  }

  den_res <- .fit_and_extract_p(denominator)
  if (is.null(numerator)) {
    p_ratio <- 1 / den_res$p
    num_mod <- NULL
  } else {
    num_res <- .fit_and_extract_p(numerator)
    p_ratio <- num_res$p / den_res$p
    num_mod <- num_res$model
  }

  weights <- unsplit(lapply(split(p_ratio, data[[id_col]]), cumprod), data[[id_col]])
  list(weights = weights, num_mod = num_mod, den_mod = den_res$model)
}

# --- Gaussian Engine -----------------------------------------------------
.fit_ipw_gaussian <- function(data, var_names, numerator, denominator, corstr, dots, ...) {
  if (!requireNamespace("geepack", quietly = TRUE)) {
    stop("Package 'geepack' is required for family = 'gaussian'.", call. = FALSE)
  }

  exp_col  <- var_names$exposure
  id_col   <- var_names$id
  time_col <- var_names$timevar

  .fit_gee <- function(rhs_form) {
    form <- .build_model_formula(exp_col, rhs_form)
    mod <- do.call(geepack::geeglm, c(list(
      formula = form,
      data    = data,
      id      = data[[id_col]],
      corstr  = corstr,
      waves   = data[[time_col]]
    ), dots))
    disp <- as.numeric(sqrt(summary(mod)$dispersion)[1])
    dens <- stats::dnorm(data[[exp_col]], mean = stats::predict(mod), sd = disp)
    list(model = mod, dens = dens)
  }

  num_res <- .fit_gee(numerator)
  den_res <- .fit_gee(denominator)

  weights <- unsplit(lapply(split(num_res$dens / den_res$dens, data[[id_col]]), cumprod), data[[id_col]])
  list(weights = weights, num_mod = num_res$model, den_mod = den_res$model)
}
