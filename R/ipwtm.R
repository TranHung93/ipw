#' Estimate Inverse Probability Weights (Time Varying)
#'
#' @description
#' Estimates inverse probability weights to fit marginal structural models with
#' time-varying exposures and time-varying confounders across follow-up intervals.
#' Supports binomial, survival, multinomial, ordinal, and gaussian families.
#'
#' @details
#' Standard errors: \code{ipwtm} only estimates the weights. Standard errors
#' for the marginal structural model have to be obtained when that model is fitted.
#' Robust standard errors, for example from \code{coxph()} with \code{cluster(id)}
#' as in Van der Wal and Geskus (2011), treat the weights as known and tend to be
#' conservative (too large). Austin (2016) compared naive, robust and bootstrap
#' standard errors in simulations and recommends the bootstrap, in which the weights
#' are re-estimated in each bootstrap sample. Austin studied a point treatment with
#' baseline confounders; for longitudinal data in long format the bootstrap must
#' resample subjects, not rows: keep all rows of a sampled subject, give each sampled
#' copy a new id, re-run \code{ipwtm} (for treatment and, if used, censoring weights)
#' on the bootstrap sample and refit the marginal structural model. See the cluster
#' bootstrap in Example 1 below; \code{\link{ipwboot}} does this resampling for you
#' (argument \code{id}). About 200 bootstrap samples are usually enough to estimate
#' a standard error.
#'
#' @param exposure Name of the exposure column in \code{data} (unquoted or as a string).
#' @param family Model family: \code{"binomial"}, \code{"survival"}, \code{"multinomial"},
#'   \code{"ordinal"}, or \code{"gaussian"}.
#' @param link Link function for binomial or ordinal models.
#' @param numerator Formula for numerator model. If \code{NULL}, unstabilized weights are computed.
#' @param denominator Formula for denominator model. For \code{family = "binomial"},
#'   smooth terms (\code{s()}, \code{te()}, \code{ti()}, \code{t2()}) are allowed;
#'   the model is then fitted with \code{mgcv::gam}.
#' @param id Name of the subject identifier column in \code{data}.
#' @param tstart Name of the interval start time column (required for \code{family = "survival"}).
#' @param timevar Name of the follow-up (interval end) time column.
#' @param type Exposure switch type: \code{"first"}, \code{"cens"}, or \code{"all"}.
#'   \code{"all"} is available for all families. With \code{family = "survival"}
#'   and \code{type = "all"}, \code{exposure} is an event indicator (1 when an
#'   exposure event, such as starting treatment, happens at the end of the
#'   interval) that may occur more than once per subject: an Andersen-Gill Cox
#'   model is fitted on all rows and every row contributes to the weights
#'   (Graffeo et al. 2018). For treatment that is started and stopped, compute
#'   weights separately for starting (on the rows off treatment) and for
#'   stopping (on the rows on treatment).
#' @param data Data frame containing all required variables.
#' @param corstr Correlation structure for GEE when \code{family = "gaussian"} (default: \code{"ar1"}).
#' @param trunc Optional truncation percentile between 0 and 0.5.
#' @param ... Additional arguments passed to model fitting engines.
#'
#' @references
#' Graffeo N., Latouche A., Geskus R.B. & Chevret S. (2018). Modeling time-varying
#' exposure using inverse probability of treatment weights. \emph{Biometrical
#' Journal}, \bold{60}(2), 323-332. \doi{10.1002/bimj.201600223}.
#'
#' Van der Wal W.M. & Geskus R.B. (2011). ipw: An R Package for Inverse Probability
#' Weighting. \emph{Journal of Statistical Software}, \bold{43}(13), 1-23.
#' \doi{10.18637/jss.v043.i13}.
#'
#' Austin, P.C. (2016). Variance estimation when using inverse probability of treatment weighting (IPTW) with survival analysis. \emph{Statistics in Medicine}, \bold{35}, 5642-5655. \doi{10.1002/sim.7084}.
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
#' #Bootstrap standard error for the MSM above (cluster bootstrap by patient).
#' #cluster() gives a robust standard error that treats the weights as known and
#' #tends to be conservative. Austin (2016) recommends the bootstrap, with the
#' #weights re-estimated in each bootstrap sample.
#' #With longitudinal data, resample patients, not rows: keep all rows of a
#' #sampled patient, give each sampled copy a new id, re-run both ipwtm() calls
#' #and refit the MSM.
#' #B = 20 keeps the example short; use about 200 in practice.
#' \donttest{
#' msm.fun <- function(d) {
#'   tw <- ipwtm(exposure = haartind, family = "survival",
#'               numerator = ~ sex + age, denominator = ~ sex + age + cd4.sqrt,
#'               id = patient, tstart = tstart, timevar = fuptime,
#'               type = "first", data = d)
#'   cw <- ipwtm(exposure = dropout, family = "survival",
#'               numerator = ~ sex + age, denominator = ~ sex + age + cd4.sqrt,
#'               id = patient, tstart = tstart, timevar = fuptime,
#'               type = "cens", data = d)
#'   d$w <- tw$ipw.weights * cw$ipw.weights
#'   coef(coxph(Surv(tstart, fuptime, event) ~ haartind, data = d,
#'              weights = w))["haartind"]
#' }
#' rows.pat <- split(seq_len(nrow(haartdat)), haartdat$patient)
#' set.seed(123)
#' B <- 20
#' est.boot <- numeric(B)
#' for (b in seq_len(B)) {
#'   pick <- sample(length(rows.pat), replace = TRUE)
#'   d.boot <- haartdat[unlist(rows.pat[pick]), ]
#'   d.boot$patient <- rep(seq_along(pick), lengths(rows.pat[pick]))
#'   est.boot[b] <- msm.fun(d.boot)
#' }
#' #log hazard ratio with bootstrap standard error
#' c(estimate = unname(msm.fun(haartdat)), boot.se = sd(est.boot))
#' #percentile interval (crude with B = 20)
#' quantile(est.boot, c(0.025, 0.975))
#' }
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
    stop("Invalid IPW weights calculated (NA or Inf detected). ",
         "Check for near-zero denominator probabilities.", call. = FALSE)
  }

  # 7. Weight Truncation ----------------------------------------------------
  trunc_sorted <- if (!is.null(trunc)) .truncate_weights(weights_sorted, trunc) else NULL

  # 8. Assemble Return Object in Original Row Order -------------------------
  call_repl <- list(
    data   = call_obj$data,
    id     = as.name(var_names$id),
    waves  = as.name(var_names$timevar),
    family = if (family == "binomial") call("binomial", link = link) else NULL
  )
  res_engine$num_mod <- .tidy_model_call(res_engine$num_mod, call_repl)
  res_engine$den_mod <- .tidy_model_call(res_engine$den_mod, call_repl)

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
  .check_link(family, link)
  if ((.uses_smooths(denominator) || .uses_smooths(numerator)) && family != "binomial") {
    stop("Smooth terms (s(), te(), ti(), t2()) are only supported for family = \"binomial\" ",
         "in ipwtm().", call. = FALSE)
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
  .check_trunc(trunc)

  # All model variables present and complete
  id_vars <- unlist(var_names[c("exposure", "id", "timevar",
                                 if (family == "survival") "tstart")])
  model_vars <- unique(c(id_vars, all.vars(denominator),
                         if (!is.null(numerator)) all.vars(numerator)))
  missing_vars <- setdiff(model_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  na_vars <- model_vars[vapply(model_vars, function(v) anyNA(data[[v]]), logical(1))]
  if (length(na_vars) > 0) {
    stop("Missing values detected in: ", paste(na_vars, collapse = ", "),
         ". 'ipwtm' requires complete cases.", call. = FALSE)
  }

  if (family %in% c("binomial", "survival") && !all(data[[var_names$exposure]] %in% c(0, 1))) {
    stop("Exposure must be coded as 0 or 1 for family = \"", family, "\".", call. = FALSE)
  }

  if (family == "survival") {
    for (f in list(numerator, denominator)) {
      if (!is.null(f) &&
          !is.null(attr(stats::terms(f, specials = "strata"), "specials")$strata)) {
        stop("strata() terms are not supported for family = \"survival\".", call. = FALSE)
      }
    }
  }
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
    # First switch away from the lowest exposure level in the whole data set
    # (the reference level is global, not per subject; as in ipw <= 1.2.1.1)
    ref_lvl <- sort(unique(exposure))[1]
    unsplit(lapply(split(exposure, id), function(x) {
      idx <- match(TRUE, x != ref_lvl)
      if (!is.na(idx)) c(rep(1, idx), rep(0, length(x) - idx)) else rep(1, length(x))
    }), id)
  }
}


# ==========================================================================
# Modeling Engines
# ==========================================================================

# --- Binomial Engine -----------------------------------------------------
.fit_ipw_binomial <- function(data, var_names, numerator, denominator, type, link, dots, ...) {
  fam_link <- stats::binomial(link = link)

  exp_col <- var_names$exposure
  id_col  <- var_names$id

  # Denominator model
  form_den <- .build_model_formula(exp_col, denominator)
  mod_den <- .fit_glm_or_gam(c(list(
    formula   = form_den,
    family    = fam_link,
    data      = data,
    subset    = if (type != "all") quote(selvar == 1) else NULL,
    na.action = stats::na.fail
  ), dots))

  pred_den <- as.numeric(stats::predict(mod_den, type = "response"))
  p_den <- rep(1, nrow(data))
  active <- (data$selvar == 1)

  if (type == "all") {
    p_den[data[[exp_col]] == 0] <- 1 - pred_den[data[[exp_col]] == 0]
    p_den[data[[exp_col]] == 1] <- pred_den[data[[exp_col]] == 1]
    w_den <- unsplit(lapply(split(p_den, data[[id_col]]), cumprod), data[[id_col]])
  } else {
    p_den[active & data[[exp_col]] == 0] <- 1 - pred_den[data[[exp_col]][active] == 0]
    if (type == "first") {
      p_den[active & data[[exp_col]] == 1] <- pred_den[data[[exp_col]][active] == 1]
    } else if (type == "cens") {
      p_den[active & data[[exp_col]] == 1] <- 1 - pred_den[data[[exp_col]][active] == 1]
    }
    p_den[!active] <- 1
    w_den <- unsplit(lapply(split(p_den, data[[id_col]]), cumprod), data[[id_col]])
  }

  # Numerator model
  mod_num <- NULL
  if (is.null(numerator)) {
    w_num <- rep(1, nrow(data))
  } else {
    form_num <- .build_model_formula(exp_col, numerator)
    mod_num <- .fit_glm_or_gam(c(list(
      formula   = form_num,
      family    = fam_link,
      data      = data,
      subset    = if (type != "all") quote(selvar == 1) else NULL,
      na.action = stats::na.fail
    ), dots))

    pred_num <- as.numeric(stats::predict(mod_num, type = "response"))
    p_num <- rep(1, nrow(data))

    if (type == "all") {
      p_num[data[[exp_col]] == 0] <- 1 - pred_num[data[[exp_col]] == 0]
      p_num[data[[exp_col]] == 1] <- pred_num[data[[exp_col]] == 1]
      w_num <- unsplit(lapply(split(p_num, data[[id_col]]), cumprod), data[[id_col]])
    } else {
      p_num[active & data[[exp_col]] == 0] <- 1 - pred_num[data[[exp_col]][active] == 0]
      if (type == "first") {
        p_num[active & data[[exp_col]] == 1] <- pred_num[data[[exp_col]][active] == 1]
      } else if (type == "cens") {
        p_num[active & data[[exp_col]] == 1] <- 1 - pred_num[data[[exp_col]][active] == 1]
      }
      p_num[!active] <- 1
      w_num <- unsplit(lapply(split(p_num, data[[id_col]]), cumprod), data[[id_col]])
    }
  }

  list(weights = w_num / w_den, num_mod = mod_num, den_mod = mod_den)
}

# --- Survival Engine -----------------------------------------------------
.fit_ipw_survival <- function(data, var_names, numerator, denominator, type, dots, ...) {
  exp_col   <- var_names$exposure
  id_col    <- var_names$id
  time_col  <- var_names$timevar
  start_col <- var_names$tstart
  surv_lhs  <- sprintf("survival::Surv(%s, %s, %s)", start_col, time_col, exp_col)
  u_times   <- sort(unique(data[[time_col]]))

  # type = "all" (Graffeo et al. 2018): Andersen-Gill model for recurrent
  # exposure events, fitted on all rows. For each row (tstart, timevar], the
  # probability of the observed status uses the cumulative baseline hazard over
  # that row's own interval: P(no event) = exp(-(L0(timevar) - L0(tstart)) * risk).
  .compute_surv_p_all <- function(mod) {
    bh <- survival::basehaz(mod, centered = TRUE)
    cum_haz_at <- function(t) {
      out <- stats::approx(bh$time, bh$hazard, xout = t, method = "constant", rule = 2)$y
      out[t < min(bh$time)] <- 0
      out
    }
    d_cum_haz <- cum_haz_at(data[[time_col]]) - cum_haz_at(data[[start_col]])
    h <- d_cum_haz * stats::predict(mod, type = "risk", centered = TRUE)
    ifelse(data[[exp_col]] == 1, 1 - exp(-h), exp(-h))
  }

  .compute_surv_p <- function(mod) {
    if (type == "all") return(.compute_surv_p_all(mod))
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
    subset    = if (type != "all") quote(selvar == 1) else NULL,
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
      subset    = if (type != "all") quote(selvar == 1) else NULL,
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

  lhs_str <- if (is_ordinal) sprintf("as.factor(%s)", exp_col) else exp_col
  polr_method <- if (is_ordinal) .polr_method(link) else NULL

  .fit_and_extract_p <- function(rhs_form) {
    form <- .build_model_formula(lhs_str, rhs_form)
    mod <- if (is_ordinal) {
      do.call(MASS::polr, .with_defaults(list(
        formula   = form,
        data      = data,
        method    = polr_method,
        subset    = quote(selvar == 1),
        na.action = stats::na.fail,
        model     = TRUE
      ), dots, list(Hess = TRUE)))
    } else {
      do.call(nnet::multinom, .with_defaults(list(
        formula   = form,
        data      = data,
        subset    = quote(selvar == 1),
        na.action = stats::na.fail
      ), dots, list(trace = FALSE, Hess = TRUE)))
    }

    active <- data$selvar == 1
    p_vec <- rep(1, nrow(data))
    p_vec[active] <- .prob_of_observed(stats::predict(mod, type = "probs"),
                                       data[[exp_col]][active], mod$lev)
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
