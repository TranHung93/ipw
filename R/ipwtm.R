#' Estimate Inverse Probability Weights (Time Varying)
#' @description
#' Estimate inverse probability weights to fit marginal structural models, with a time-varying exposure and time-varying confounders. Within each unit under observation this function computes inverse probability weights at each time point during follow-up. The exposure can be binomial, multinomial, ordinal or continuous. Both stabilized and unstabilized weights can be estimated.
#'
#' @param exposure vector, representing the exposure of interest. Both numerical and categorical variables can be used. A binomial exposure variable should be coded using values \code{0}/\code{1}.
#' @param family specifies a family of link functions, used to model the relationship between the variables in \code{numerator} or \code{denominator} and \code{exposure}, respectively. Alternatives are \code{"binomial"}, \code{"survival"}, \code{"multinomial"}, \code{"ordinal"} and \code{"gaussian"}. A specific link function is then chosen using the argument \code{link}, as explained below. Regression models are fitted using \code{\link{glm}}, \code{\link{coxph}}, \code{\link{multinom}}, \code{\link{polr}} or \code{\link{geeglm}}, respectively.
#' @param link specifies the specific link function between the variables in \code{numerator} or \code{denominator} and exposure, respectively. For \code{family="binomial"} (fitted using \code{\link{glm}}) alternatives are \code{"logit"}, \code{"probit"}, \code{"cauchit"}, \code{"log"} and \code{"cloglog"}. For \code{family="survival"} this argument is ignored, and Cox proportional hazards models are always used (fitted using \code{\link{coxph}}). For \code{family="multinomial"} this argument is ignored, and multinomial logistic regression models are always used (fitted using \code{\link{multinom}}). For \code{family=} \code{"ordinal"} (fitted using \code{\link{polr}}) alternatives are \code{"logit"}, \code{"probit"}, \code{"cauchit"}, and \code{"cloglog"}. For \code{family="gaussian"} this argument is ignored, and GEE models with an identity link are always used (fitted using \code{\link{geeglm}}.)
#' @param numerator is a formula, specifying the right-hand side of the model used to estimate the elements in the numerator of the inverse probability weights. When left unspecified, unstabilized weights with a numerator of 1 are estimated.
#' @param denominator is a formula, specifying the right-hand side of the model used to estimate the elements in the denominator of the inverse probability weights.
#' @param id vector, uniquely identifying the units under observation (typically patients) within which the longitudinal measurements are taken.
#' @param tstart numerical vector, representing the starting time of follow-up intervals, using the counting process notation. This argument is only needed when \code{family=} \code{"survival"}, otherwise it is ignored. The Cox proportional hazards models are fitted using counting process data. Since a switch in exposure level can occur at the start of follow-up, \code{tstart} should be negative for the first interval (with \code{timevar=0}) within each patient.
#' @param timevar numerical vector, representing follow-up time, starting at \code{0}. This variable is used as the end time of follow-up intervals, using the counting process notation, when \code{family="survival"}.
#' @param type specifies the type of exposure. Alternatives are \code{"first"}, \code{"cens"} and \code{"all"}. With \code{type="first"}, weights are estimated up to the first switch from the lowest exposure value (typically \code{0} or the first factor level) to any other value. After this switch, weights will then be constant. Such a weight is e.g. used when estimating the effect of ``initiation of HAART'' on mortality (see example 1 below). \code{type="first"} is currently only implemented for \code{"binomial"}, \code{"survival"}, \code{"multinomial"} and \code{"ordinal"} families. With \code{type="cens"} inverse probability of censoring weights (IPCW) are estimated as defined in appendix 1 in Cole & Hernán (2008). IPCW is illustrated in example 1 below. \code{type="cens"} is currently only implemented for \code{"binomial"} and \code{"survival"} families. With \code{type="all"}, all time points are used to estimate weights. \code{type="all"} is implemented only for the \code{"binomial"} and \code{"gaussian"} family.
#' @param data dataframe containing \code{exposure}, variables in \code{numerator} and \code{denominator}, \code{id}, \code{tstart} and \code{timevar}.
#' @param corstr correlation structure, only needed when using \code{family = "gaussian"}. Defaults to "ar1". See \code{\link{geeglm}} for details.
#' @param trunc optional truncation percentile (0-0.5). E.g. when \code{trunc = 0.01}, the left tail is truncated to the 1st percentile, and the right tail is truncated to the 99th percentile. When specified, both un-truncated and truncated weights are returned.
#' @param ... are further arguments passed to the function that is used to estimate the numerator and denominator models (the function is chosen using \code{family}).
#' @details
#' Within each unit under observation i (usually patients), this function computes inverse probability weights at each time point j during follow-up. These weights are the cumulative product over all previous time points up to j of the ratio of two probabilities:
#' \itemize{
#' \item the numerator contains at each time point the probability of the observed exposure level given observed values of stabilization factors and the observed exposure history up to the time point before j. These probabilities are estimated using the model regressing \code{exposure} on the terms in \code{numerator}, using the link function indicated by \code{family} and \code{link}.
#' \item the denominator contains at each time point the probability of the observed exposure level given the observed history of time varying confounders up to j, as well as the stabilization factors in the numerator and the observed exposure history up to the time point before j. These probabilities are estimated using the model regressing \code{exposure} on the terms in \code{denominator}, using the link function indicated by \code{family} and \code{link}.}
#'
#' When the models from which the elements in the numerator and denominator are predicted are correctly specified, and there is no unmeasured confounding, weighting observations ij by the inverse probability weights adjusts for confounding of the effect of the exposure of interest. On the weighted dataset a marginal structural model can then be fitted, quantifying the causal effect of the exposure on the outcome of interest.
#'
#' With \code{numerator} specified, stabilized weights are computed, otherwise unstabilized weights with a numerator of 1 are computed. With a continuous exposure, using \code{family = "gaussian"}, weights are computed using the ratio of predicted densities at each time point. Therefore, for \code{family = "gaussian"} only stabilized weights can be used, since unstabilized weights would have infinity variance.
#'
#'
#' @return A list containing the following elements:
#' \item{ipw.weights }{vector containing inverse probability weights for each observation. Returned in the same order as the observations in \code{data}, to facilitate merging.}
#' \item{weights.trunc }{vector containing truncated inverse probability weights, only returned when \code{trunc} is specified.}
#' \item{call }{the original function call.}
#' \item{selvar }{selection variable. With \code{type = "first"}, \code{selvar = 1} within each unit under observation, up to and including the first time point at which a switch from the lowest value of \code{exposure} to any other value is made, and \code{selvar = 0} after the first switch. For \code{type = "all"}, \code{selvar = 1} for all measurements. The numerator and denominator models are fitted only on observations with \code{selvar = 1}. Returned in the same order as observations in \code{data}, to facilitate merging.}
#' \item{num.mod }{the numerator model, only returned when \code{numerator} is specified.}
#' \item{den.mod }{the denominator model.}
#'
#' @section Missing values:
#' Currently, the \code{exposure} variable and the variables used in \code{numerator} and \code{denominator}, \code{id}, \code{tstart} and \code{timevar} should not contain missing values.
#'
#' @author Willem M. van der Wal \email{willem@vanderwalresearch.com}, Ronald B. Geskus \email{rgeskus@oucru.org}
#' @references Cole, S.R. & Hernán, M.A. (2008). Constructing inverse probability weights for marginal structural models. \emph{American Journal of Epidemiology}, \bold{168}(6), 656-664. \url{https://pubmed.ncbi.nlm.nih.gov:443/18682488/}.
#'
#' Robins, J.M., Hernán, M.A. & Brumback, B.A. (2000). Marginal structural models and causal inference in epidemiology. \emph{Epidemiology}, \bold{11}, 550-560. \url{https://pubmed.ncbi.nlm.nih.gov/10955408/}.
#'
#' Van der Wal W.M. & Geskus R.B. (2011). ipw: An R Package for Inverse  Probability Weighting. \emph{Journal of Statistical Software}, \bold{43}(13), 1-23. \url{https://doi.org/10.18637/jss.v043.i13}
#' @seealso \code{\link{basdat}}, \code{\link{haartdat}}, \code{\link{ipwplot}}, \code{\link{ipwpoint}}, \code{\link{ipwtm}}, \code{\link{timedat}}, \code{\link{tstartfun}}.
#' @export
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
#'         binwidth = 100, ylim = c(-1.5, 1.5), main = "Stabilized inverse probability of censoring weights")
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
#' startstop$cd4.sqrt <- predict(cd4.lme, newdata = data.frame(id = startstop$id,
#'                                                             fuptime = startstop$fuptime, tb.lag = startstop$tb.lag))
#'
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
    family,
    link,
    numerator = NULL,
    denominator,
    id,
    tstart,
    timevar,
    type,
    data,
    corstr = "ar1",
    trunc = NULL,
    ...)
{
  #save input
  tempcall <- match.call()
  #some basic input checks
  if (!("exposure" %in% names(tempcall))) stop("No exposure variable specified")
  if (!("family" %in% names(tempcall)) | ("family" %in% names(tempcall) & !(tempcall$family %in% c("binomial", "survival", "multinomial", "ordinal", "gaussian")))) stop("No valid family specified (\"binomial\", \"survival\", \"multinomial\", \"ordinal\", \"gaussian\")")
  if (tempcall$family == "binomial") {if(!(tempcall$link %in% c("logit", "probit", "cauchit", "log", "cloglog"))) stop("No valid link function specified for family = binomial (\"logit\", \"probit\", \"cauchit\", \"log\", \"cloglog\")")}
  if (tempcall$family == "ordinal" ) {if(!(tempcall$link %in% c("logit", "probit", "cauchit", "cloglog"))) stop("No valid link function specified for family = ordinal (\"logit\", \"probit\", \"cauchit\", \"cloglog\")")}
  if (!("denominator" %in% names(tempcall))) stop("No denominator model specified")
  if (!is.null(tempcall$numerator) & !is(eval(tempcall$numerator), "formula")) stop("Invalid numerator formula specified")
  if (!is.null(tempcall$denominator) & !is(eval(tempcall$denominator), "formula")) stop("Invalid denominator formula specified")
  if (!("id" %in% names(tempcall))) stop("No patient id specified")
  if (tempcall$family == "survival" & !("tstart" %in% names(tempcall))) stop("No tstart specified, is necessary for family = \"survival\"")
  if (!("timevar" %in% names(tempcall))) stop("No timevar specified")
  if (!("type" %in% names(tempcall))) stop("No type specified (\"first\" or \"all\")")
  if (!(tempcall$type %in% c("first", "all", "cens"))) stop("No type specified (\"first\", \"all\" or \"cens\")")
  if (tempcall$family %in% c("survival", "multinomial", "ordinal") & tempcall$type == "all") stop(paste("Type \"all\" not yet implemented for family = ", deparse(tempcall$family, width.cutoff = 500), sep = ""))
  if (tempcall$family %in% c("multinomial", "ordinal", "gaussian") & tempcall$type == "cens") stop(paste("Type \"cens\" not yet implemented for family = ", deparse(tempcall$family, width.cutoff = 500), sep = ""))
  if (tempcall$family %in% c("gaussian") & tempcall$type == "first") stop(paste("Type \"first\" not implemented for family = ", deparse(tempcall$family, width.cutoff = 500), sep = ""))
  if (tempcall$family %in% c("gaussian") & !("numerator" %in% names(tempcall))) stop("Numerator necessary for family = \"gaussian\"")
  if (!("data" %in% names(tempcall))) stop("No data specified")
  if (!is.null(tempcall$trunc)) {if(tempcall$trunc < 0 | tempcall$trunc > 0.5) stop("Invalid truncation percentage specified (0-0.5)")}
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
    timevar = data[,as.character(tempcall$timevar)],
    exposure = data[,as.character(tempcall$exposure)]
  )
  #make selection variable, time points up to first switch from lowest value, or all time points
  if (type %in% c("first", "cens") & (family == "binomial" | family == "survival"))
  {tempdat$selvar <- do.call("c", lapply(split(tempdat$exposure, tempdat$id),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))}
  if (type %in% c("first", "cens") & (family == "multinomial" | family == "ordinal")){
    z <- unique(tempdat$exposure)[unique(tempdat$exposure) != sort(unique(tempdat$exposure))[1]]
    min2 <- function(x)ifelse(min(is.na(unique(x))) == 1, NA, min(x, na.rm = TRUE))
    tempdat$selvar <- do.call("c", lapply(split(tempdat$exposure, tempdat$id),function(x)if (!is.na(min2(match(z, x)))) return(c(rep(1,min2(match(z, x))),rep(0,length(x)-min2(match(z, x))))) else return(rep(1,length(x)))))
  }
  if (type == "all")
  {tempdat$selvar <- rep(1, nrow(tempdat))}
  #weights binomial, type "first"
  if (tempcall$family == "binomial" & tempcall$type %in% c("first", "cens")) {
    if(tempcall$link == "logit") lf <- binomial(link = logit)
    if(tempcall$link == "probit") lf  <- binomial(link = probit)
    if(tempcall$link == "cauchit") lf  <- binomial(link = cauchit)
    if(tempcall$link == "log") lf  <- binomial(link = log)
    if(tempcall$link == "cloglog") lf  <- binomial(link = cloglog)
    if (is.null(tempcall$numerator)) tempdat$w.numerator <- 1
    else {
      mod1 <- glm(
        formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = ""))),
        family = lf,
        data = data,
        subset = tempdat$selvar == 1,
        na.action = na.fail,
        ...)
      tempdat$p.numerator <- vector("numeric", nrow(tempdat))
      tempdat$p.numerator[tempdat$exposure == 0 & tempdat$selvar == 1] <- 1 - predict.glm(mod1, type = "response")[tempdat$exposure[tempdat$selvar == 1] == 0]
      if(type == "first"){tempdat$p.numerator[tempdat$exposure == 1 & tempdat$selvar == 1] <- predict.glm(mod1, type = "response")[tempdat$exposure[tempdat$selvar == 1] == 1]}
      if(type == "cens"){tempdat$p.numerator[tempdat$exposure == 1 & tempdat$selvar == 1] <- 1 - predict.glm(mod1, type = "response")[tempdat$exposure[tempdat$selvar == 1] == 1]}
      tempdat$p.numerator[tempdat$selvar == 0] <- 1
      tempdat$w.numerator <- unlist(lapply(split(tempdat$p.numerator, tempdat$id), function(x)cumprod(x)))
      mod1$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = "")))
      mod1$call$family <- tempcall$link
      mod1$call$data <- tempcall$data
      mod1$call$subset <- paste("up to first instance of ", deparse(tempcall$exposure, width.cutoff = 500), " = 1 (selvar == 1)", sep = "")
    }
    mod2 <- glm(
      formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = ""))),
      family = lf,
      data = data,
      subset = tempdat$selvar == 1,
      na.action = na.fail,
      ...)
    tempdat$p.denominator <- vector("numeric", nrow(tempdat))
    tempdat$p.denominator[tempdat$exposure == 0 & tempdat$selvar == 1] <- 1 - predict.glm(mod2, type = "response")[tempdat$exposure[tempdat$selvar == 1] == 0]
    if(type == "first"){tempdat$p.denominator[tempdat$exposure == 1 & tempdat$selvar == 1] <- predict.glm(mod2, type = "response")[tempdat$exposure[tempdat$selvar == 1] == 1]}
    if(type == "cens"){tempdat$p.denominator[tempdat$exposure == 1 & tempdat$selvar == 1] <- 1 - predict.glm(mod2, type = "response")[tempdat$exposure[tempdat$selvar == 1] == 1]}
    tempdat$p.denominator[tempdat$selvar == 0] <- 1
    tempdat$w.denominator <- unlist(lapply(split(tempdat$p.denominator, tempdat$id), function(x)cumprod(x)))
    mod2$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = "")))
    mod2$call$family <- tempcall$link
    mod2$call$data <- tempcall$data
    mod2$call$subset <- paste("up to first instance of ", deparse(tempcall$exposure, width.cutoff = 500), " = 1 (selvar == 1)", sep = "")
    tempdat$ipw.weights <- tempdat$w.numerator/tempdat$w.denominator
  }
  #weights binomial, type "all"
  if (tempcall$family == "binomial" & tempcall$type == "all") {
    if(tempcall$link == "logit") lf <- binomial(link = logit)
    if(tempcall$link == "probit") lf  <- binomial(link = probit)
    if(tempcall$link == "cauchit") lf  <- binomial(link = cauchit)
    if(tempcall$link == "log") lf  <- binomial(link = log)
    if(tempcall$link == "cloglog") lf  <- binomial(link = cloglog)
    if (is.null(tempcall$numerator)) tempdat$w.numerator <- 1
    else {
      mod1 <- glm(
        formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = ""))),
        family = lf,
        data = data,
        na.action = na.fail,
        ...)
      tempdat$p.numerator <- vector("numeric", nrow(tempdat))
      tempdat$p.numerator[tempdat$exposure == 0] <- 1 - predict.glm(mod1, type = "response")[tempdat$exposure == 0]
      tempdat$p.numerator[tempdat$exposure == 1] <- predict.glm(mod1, type = "response")[tempdat$exposure == 1]
      tempdat$w.numerator <- unlist(lapply(split(tempdat$p.numerator, tempdat$id), function(x)cumprod(x)))
      mod1$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = "")))
      mod1$call$family <- tempcall$link
      mod1$call$data <- tempcall$data
    }
    mod2 <- glm(
      formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = ""))),
      family = lf,
      data = data,
      na.action = na.fail,
      ...)
    tempdat$p.denominator <- vector("numeric", nrow(tempdat))
    tempdat$p.denominator[tempdat$exposure == 0] <- 1 - predict.glm(mod2, type = "response")[tempdat$exposure == 0]
    tempdat$p.denominator[tempdat$exposure == 1] <- predict.glm(mod2, type = "response")[tempdat$exposure == 1]
    tempdat$w.denominator <- unlist(lapply(split(tempdat$p.denominator, tempdat$id), function(x)cumprod(x)))
    mod2$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = "")))
    mod2$call$family <- tempcall$link
    mod2$call$data <- tempcall$data
    tempdat$ipw.weights <- tempdat$w.numerator/tempdat$w.denominator
  }
  #weights Cox
  if (tempcall$family == "survival") {
    if (is.null(tempcall$numerator)) tempdat$w.numerator <- 1
    else {
      mod1 <- coxph(
        formula = eval(parse(text = paste("Surv(", deparse(tempcall$tstart), ", ", deparse(tempcall$timevar, width.cutoff = 500), ", ", deparse(tempcall$exposure, width.cutoff = 500), ") ", deparse(tempcall$numerator, width.cutoff = 500), sep = ""))),
        data = data,
        subset = tempdat$selvar == 1,
        na.action = na.fail,
        method = "efron",
        ...)
      bh <- basehaz(mod1, centered = TRUE)
      temp <- data.frame(timevar = sort(unique(tempdat$timevar)))
      temp <- merge(temp, data.frame(timevar = bh$time, bashaz.cum.numerator = bh$hazard), by = "timevar", all.x = TRUE);rm(bh)
      if (is.na(temp$bashaz.cum.numerator[temp$timevar == min(unique(tempdat$timevar))])) temp$bashaz.cum.numerator[temp$timevar == min(unique(tempdat$timevar))] <- 0
      temp$bashaz.cum.numerator <- approx(x = temp$timevar, y = temp$bashaz.cum.numerator, xout = temp$timevar, method = "constant", rule = 2)$y
      temp$bashaz.numerator[1] <- temp$bashaz.cum.numerator[1]
      temp$bashaz.numerator[2:nrow(temp)] <- diff(temp$bashaz.cum.numerator, 1)
      temp$bashaz.cum.numerator <- NULL
      tempdat <- merge(tempdat, temp, by = "timevar", all.x = TRUE);rm(temp)
      tempdat <- tempdat[order(tempdat$id, tempdat$timevar),]
      tempdat$risk.numerator[tempdat$selvar == 1] <-predict(mod1, type="risk", centered = TRUE)
      tempdat$hazard.numerator[tempdat$selvar == 1] <- with(tempdat[tempdat$selvar == 1,], bashaz.numerator*risk.numerator)
      tempdat$p.numerator[with(tempdat, selvar == 1 & exposure == 0)] <- with(tempdat[with(tempdat, selvar == 1 & exposure == 0),], exp(-1*bashaz.numerator*risk.numerator))
      if(type == "first"){tempdat$p.numerator[with(tempdat, selvar == 1 & exposure == 1)] <- 1 - with(tempdat[with(tempdat, selvar == 1 & exposure == 1),], exp(-1*bashaz.numerator*risk.numerator))}
      if(type == "cens"){tempdat$p.numerator[with(tempdat, selvar == 1 & exposure == 1)] <- with(tempdat[with(tempdat, selvar == 1 & exposure == 1),], exp(-1*bashaz.numerator*risk.numerator))}
      tempdat$p.numerator[tempdat$selvar == 0] <- 1
      tempdat$w.numerator <- unsplit(lapply(split(tempdat$p.numerator, tempdat$id), function(x)cumprod(x)), tempdat$id)
      mod1$call$formula <- eval(parse(text = paste("Surv(", deparse(tempcall$tstart), ", ", deparse(tempcall$timevar, width.cutoff = 500), ", ", deparse(tempcall$exposure, width.cutoff = 500), ") ", deparse(tempcall$numerator, width.cutoff = 500), sep = "")))
      mod1$call$data <- tempcall$data
    }
    mod2 <- coxph(
      formula = eval(parse(text = paste("Surv(", deparse(tempcall$tstart), ", ", deparse(tempcall$timevar, width.cutoff = 500), ", ", deparse(tempcall$exposure, width.cutoff = 500), ") ", deparse(tempcall$denominator, width.cutoff = 500), sep = ""))),
      data = data,
      subset = tempdat$selvar == 1,
      na.action = na.fail,
      method = "efron",
      ...)
    bh <- basehaz(mod2, centered = TRUE)
    temp <- data.frame(timevar = sort(unique(tempdat$timevar)))
    temp <- merge(temp, data.frame(timevar = bh$time, bashaz.cum.denominator = bh$hazard), by = "timevar", all.x = TRUE);rm(bh)
    if (is.na(temp$bashaz.cum.denominator[temp$timevar == min(unique(tempdat$timevar))])) temp$bashaz.cum.denominator[temp$timevar == min(unique(tempdat$timevar))] <- 0
    temp$bashaz.cum.denominator <- approx(x = temp$timevar, y = temp$bashaz.cum.denominator, xout = temp$timevar, method = "constant", rule = 2)$y
    temp$bashaz.denominator[1] <- temp$bashaz.cum.denominator[1]
    temp$bashaz.denominator[2:nrow(temp)] <- diff(temp$bashaz.cum.denominator, 1)
    temp$bashaz.cum.denominator <- NULL
    tempdat <- merge(tempdat, temp, by = "timevar", all.x = TRUE);rm(temp)
    tempdat <- tempdat[order(tempdat$id, tempdat$timevar),]
    tempdat$risk.denominator[tempdat$selvar == 1] <-predict(mod2, type="risk", centered = TRUE)
    tempdat$hazard.denominator[tempdat$selvar == 1] <- with(tempdat[tempdat$selvar == 1,], bashaz.denominator*risk.denominator)
    tempdat$p.denominator[with(tempdat, selvar == 1 & exposure == 0)] <- with(tempdat[with(tempdat, selvar == 1 & exposure == 0),], exp(-1*bashaz.denominator*risk.denominator))
    if(type == "first"){tempdat$p.denominator[with(tempdat, selvar == 1 & exposure == 1)] <- 1 - with(tempdat[with(tempdat, selvar == 1 & exposure == 1),], exp(-1*bashaz.denominator*risk.denominator))}
    if(type == "cens"){tempdat$p.denominator[with(tempdat, selvar == 1 & exposure == 1)] <- with(tempdat[with(tempdat, selvar == 1 & exposure == 1),], exp(-1*bashaz.denominator*risk.denominator))}
    tempdat$p.denominator[tempdat$selvar == 0] <- 1
    tempdat$w.denominator <- unsplit(lapply(split(tempdat$p.denominator, tempdat$id), function(x)cumprod(x)), tempdat$id)
    mod2$call$formula <- eval(parse(text = paste("Surv(", deparse(tempcall$tstart), ", ", deparse(tempcall$timevar, width.cutoff = 500), ", ", deparse(tempcall$exposure, width.cutoff = 500), ") ", deparse(tempcall$denominator, width.cutoff = 500), sep = "")))
    mod2$call$data <- tempcall$data
    mod2$call$subset <- paste("up to first instance of ", deparse(tempcall$exposure, width.cutoff = 500), " = 1 (selvar == 1)", sep = "")
    tempdat$ipw.weights <- tempdat$w.numerator/tempdat$w.denominator
  }
  #weights multinomial
  if (tempcall$family == "multinomial") {
    if (is.null(tempcall$numerator)) tempdat$p.numerator <- 1
    else {
      mod1 <- multinom(
        formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = ""))),
        data = data,
        subset = tempdat$selvar == 1,
        na.action = na.fail,
        ...)
      pred1 <- as.data.frame(predict(mod1, type = "probs"))
      tempdat$p.numerator[tempdat$selvar == 0] <- 1
      for (i in 1:length(unique(tempdat$exposure)))tempdat$p.numerator[with(tempdat, tempdat$selvar == 1 & exposure == sort(unique(tempdat$exposure))[i])] <- pred1[tempdat$exposure[tempdat$selvar == 1] == sort(unique(tempdat$exposure))[i],i]
      mod1$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = "")))
      mod1$call$data <- tempcall$data
      mod1$call$subset <- paste("up to first instance of ", deparse(tempcall$exposure, width.cutoff = 500), " = 1 (selvar == 1)", sep = "")
    }
    mod2 <- multinom(
      formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = ""))),
      data = data,
      subset = tempdat$selvar == 1,
      na.action = na.fail,
      ...)
    pred2 <- as.data.frame(predict(mod2, type = "probs"))
    tempdat$p.denominator[tempdat$selvar == 0] <- 1
    for (i in 1:length(unique(tempdat$exposure)))tempdat$p.denominator[with(tempdat, tempdat$selvar == 1 & exposure == sort(unique(tempdat$exposure))[i])] <- pred2[tempdat$exposure[tempdat$selvar == 1] == sort(unique(tempdat$exposure))[i],i]
    mod2$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = "")))
    mod2$call$data <- tempcall$data
    mod2$call$subset <- paste("up to first instance of ", deparse(tempcall$exposure, width.cutoff = 500), " = 1 (selvar == 1)", sep = "")
    tempdat$ipw.weights <- unsplit(lapply(split(with(tempdat, p.numerator/p.denominator), tempdat$id), function(x)cumprod(x)), tempdat$id)
  }
  #weights ordinal
  if (tempcall$family == "ordinal") {
    if(tempcall$link == "logit") m <- "logistic"
    if(tempcall$link == "probit") m  <- "probit"
    if(tempcall$link == "cloglog") m  <- "cloglog"
    if(tempcall$link == "cauchit") m  <- "cauchit"
    if (is.null(tempcall$numerator)) tempdat$p.numerator <- 1
    else {
      mod1 <- polr(
        formula = eval(parse(text = paste("as.factor(", deparse(tempcall$exposure, width.cutoff = 500), ")", deparse(tempcall$numerator, width.cutoff = 500), sep = ""))),
        data = data,
        method = m,
        subset = tempdat$selvar == 1,
        na.action = na.fail,
        ...)
      pred1 <- as.data.frame(predict(mod1, type = "probs"))
      tempdat$p.numerator[tempdat$selvar == 0] <- 1
      for (i in 1:length(unique(tempdat$exposure)))tempdat$p.numerator[with(tempdat, tempdat$selvar == 1 & exposure == sort(unique(tempdat$exposure))[i])] <- pred1[tempdat$exposure[tempdat$selvar == 1] == sort(unique(tempdat$exposure))[i],i]
      mod1$call$formula <- eval(parse(text = paste("as.factor(", deparse(tempcall$exposure, width.cutoff = 500), ")", deparse(tempcall$numerator, width.cutoff = 500), sep = "")))
      mod1$call$data <- tempcall$data
      mod1$call$method <- m
      mod1$call$subset <- paste("up to first instance of ", deparse(tempcall$exposure, width.cutoff = 500), " = 1 (selvar == 1)", sep = "")
    }
    mod2 <- polr(
      formula = eval(parse(text = paste("as.factor(", deparse(tempcall$exposure, width.cutoff = 500), ")", deparse(tempcall$denominator, width.cutoff = 500), sep = ""))),
      data = data,
      method = m,
      subset = tempdat$selvar == 1,
      na.action = na.fail,
      ...)
    pred2 <- as.data.frame(predict(mod2, type = "probs"))
    tempdat$p.denominator[tempdat$selvar == 0] <- 1
    for (i in 1:length(unique(tempdat$exposure)))tempdat$p.denominator[with(tempdat, tempdat$selvar == 1 & exposure == sort(unique(tempdat$exposure))[i])] <- pred2[tempdat$exposure[tempdat$selvar == 1] == sort(unique(tempdat$exposure))[i],i]
    mod2$call$formula <- eval(parse(text = paste("as.factor(", deparse(tempcall$exposure, width.cutoff = 500), ")", deparse(tempcall$denominator, width.cutoff = 500), sep = "")))
    mod2$call$data <- tempcall$data
    mod2$call$method <- m
    mod2$call$subset <- paste("up to first instance of ", deparse(tempcall$exposure, width.cutoff = 500), " = 1 (selvar == 1)", sep = "")
    tempdat$ipw.weights <- unsplit(lapply(split(with(tempdat, p.numerator/p.denominator), tempdat$id), function(x)cumprod(x)), tempdat$id)
  }
  #weights gaussian
  if (tempcall$family == "gaussian") {
    mod1 <- geeglm(
      formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = ""))),
      data = data,
      id = tempdat$id,
      corstr = tempcall$corstr,
      waves = tempdat$timevar,
      ...)
    mod1$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = "")))
    mod1$call$data <- tempcall$data
    mod1$call$id <- tempcall$id
    mod1$call$corstr <- tempcall$corstr
    mod1$call$waves <- tempcall$waves
    mod2 <- geeglm(
      formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = ""))),
      data = data,
      id = tempdat$id,
      corstr = tempcall$corstr,
      waves = tempdat$timevar,
      ...)
    mod2$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = "")))
    mod2$call$data <- tempcall$data
    mod2$call$id <- tempcall$id
    mod2$call$corstr <- tempcall$corstr
    mod2$call$waves <- tempcall$waves
    tempdat$kdens1 <- dnorm(tempdat$exposure, predict(mod1), as.numeric(sqrt(summary(mod1)$dispersion)[1]))
    tempdat$kdens2 <- dnorm(tempdat$exposure, predict(mod2), as.numeric(sqrt(summary(mod2)$dispersion)[1]))
    tempdat$ipw.weights <- unsplit(lapply(split(with(tempdat, kdens1/kdens2), tempdat$id), function(x)cumprod(x)), tempdat$id)
  }
  #check for NA's in weights
  if (sum(is.na(tempdat$ipw.weights)) > 0) stop ("NA's in weights!")
  #truncate weights, when trunc value is specified (0-0.5)
  if (!(is.null(tempcall$trunc))){
    tempdat$weights.trunc <- tempdat$ipw.weights
    tempdat$weights.trunc[tempdat$ipw.weights <= quantile(tempdat$ipw.weights, 0+trunc)] <- quantile(tempdat$ipw.weights, 0+trunc)
    tempdat$weights.trunc[tempdat$ipw.weights >  quantile(tempdat$ipw.weights, 1-trunc)] <- quantile(tempdat$ipw.weights, 1-trunc)
  }
  #return results in the same order as the original input dataframe
  if (is.null(tempcall$trunc)){
    if (is.null(tempcall$numerator)) return(list(ipw.weights = tempdat$ipw.weights[order(order.orig)], call = tempcall, selvar = tempdat$selvar[order(order.orig)], den.mod = mod2))
    else return(list(ipw.weights = tempdat$ipw.weights[order(order.orig)], call = tempcall, selvar = tempdat$selvar[order(order.orig)], num.mod = mod1, den.mod = mod2))
  }
  else{
    if (is.null(tempcall$numerator)) return(list(ipw.weights = tempdat$ipw.weights[order(order.orig)], weights.trunc = tempdat$weights.trunc[order(order.orig)], call = tempcall, selvar = tempdat$selvar[order(order.orig)], den.mod = mod2))
    else return(list(ipw.weights = tempdat$ipw.weights[order(order.orig)], weights.trunc = tempdat$weights.trunc[order(order.orig)], call = tempcall, selvar = tempdat$selvar[order(order.orig)], num.mod = mod1, den.mod = mod2))
  }
}
