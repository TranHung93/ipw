#' Estimate Inverse Probability Weights (Point Treatment)
#' @description
#' Estimate inverse probability weights to fit marginal structural models in a point treatment situation. The exposure for which we want to estimate the causal effect can be binomial, multinomial, ordinal or continuous. Both stabilized and unstabilized weights can be estimated.
#'
#'
#' @param exposure a vector, representing the exposure variable of interest. Both numerical and categorical variables can be used. A binomial exposure variable should be coded using values \code{0}/\code{1}.
#' @param family is used to specify a family of link functions, used to model the relationship between the variables in \code{numerator} or \code{denominator} and \code{exposure}, respectively. Alternatives are \code{"binomial"},\code{"multinomial"}, \code{"ordinal"} and \code{"gaussian"}. A specific link function is then chosen using the argument \code{link}, as explained below. Regression models are fitted using \code{\link{glm}}, \code{\link[nnet]{multinom}}, \code{\link[MASS]{polr}} or \code{\link{glm}}, respectively.
#' @param link specifies the link function between the variables in \code{numerator} or \code{denominator} and \code{exposure}, respectively. For \code{family = "binomial"} (fitted using \code{\link{glm}}) alternatives are \code{"logit"}, \code{"probit"}, \code{"cauchit"}, \code{"log"} and \code{"cloglog"}. For \code{family = "multinomial"} this argument is ignored, and multinomial logistic regression models are always used (fitted using \code{\link[nnet]{multinom}}). For \code{family = }\code{"ordinal"} (fitted using \code{\link[MASS]{polr}}) alternatives are \code{"logit"}, \code{"probit"}, \code{"cauchit"}, and \code{"cloglog"}. For \code{family = "gaussian"} this argument is ignored, and a linear regression model with identity link is always used (fitted using \code{\link{glm}}).
#' @param numerator is a formula, specifying the right-hand side of the model used to estimate the elements in the numerator of the inverse probability weights. When left unspecified, unstabilized weights with a numerator of 1 are estimated.
#' @param denominator is a formula, specifying the right-hand side of the model used to estimate the elements in the denominator of the inverse probability weights. This typically includes the variables specified in the numerator model, as well as confounders for which to correct.
#' @param data is a dataframe containing \code{exposure} and the variables used in \code{numerator} and \code{denominator}.
#' @param trunc optional truncation percentile (0-0.5). E.g. when \code{trunc = 0.01}, the left tail is truncated to the 1st percentile, and the right tail is truncated to the 99th percentile.When specified, both un-truncated and truncated weights are returned.
#' @param ... are further arguments passed to the function that is used to estimate the numerator and denominator models (the function is chosen using \code{family}).
#' @details
#' For each unit under observation, this function computes an inverse probability weight, which is the ratio of two probabilities:
#'   \itemize{
#'     \item the numerator contains the probability of the observed exposure level given observed values of stabilization factors (usually a set of baseline covariates). These probabilities are estimated using the model regressing \code{exposure} on the terms in \code{numerator}, using the link function indicated by \code{family} and \code{link}.
#'     \item the denominator contains the probability of the observed exposure level given the observed values of a set of confounders, as well as the stabilization factors in the numerator. These probabilities are estimated using the model regressing \code{exposure} on the terms in \code{denominator}, using the link function indicated by \code{family} and \code{link}.}
#'
#' When the models from which the elements in the numerator and denominator are predicted are correctly specified, and there is no unmeasured confounding, weighting the observations by the inverse probability weights adjusts for confounding of the effect of the exposure of interest. On the weighted dataset a marginal structural model can then be fitted, quantifying the causal effect of the exposure on the outcome of interest.
#'
#' With \code{numerator} specified, stabilized weights are computed, otherwise unstabilized weighs with a numerator of 1 are computed. With a continuous exposure, using \code{family = "gaussian"}, weights are computed using the ratio of predicted densities. Therefore, for \code{family = "gaussian"} only stabilized weights can be used, since unstabilized weights would have infinity variance.
#'
#' @return A list containing the following elements:
#' \item{ipw.weights }{is a vector containing inverse probability weights for each unit under observation. This vector is returned in the same order as the measurements contained in \code{data}, to facilitate merging.}
#' \item{weights.trunc }{is a vector containing truncated inverse probability weights for each unit under observation. This vector is only returned when \code{trunc} is specified.}
#' \item{call }{is the original function call to \code{ipwpoint}.}
#' \item{num.mod }{is the numerator model, only returned when \code{numerator} is specified.}
#' \item{den.mod }{is the denominator model.}
#' @section Missing values:
#' Currently, the \code{exposure} variable and the variables used in \code{numerator} and \code{denominator} should not contain missing values.
#' @author Willem M. van der Wal \email{willem@vanderwalresearch.com}, Ronald B. Geskus \email{rgeskus@oucru.org}
#' @references Cole, S.R. & Hernán, M.A. (2008). Constructing inverse probability weights for marginal structural models. \emph{American Journal of Epidemiology}, \bold{168}(6), 656-664.
#'
#' Robins, J.M., Hernán, M.A. & Brumback, B.A. (2000). Marginal structural models and causal inference in epidemiology. \emph{Epidemiology}, \bold{11}, 550-560.
#'
#' Van der Wal W.M. & Geskus R.B. (2011). ipw: An R Package for Inverse  Probability Weighting. \emph{Journal of Statistical Software}, \bold{43}(13), 1-23. \doi{10.18637/jss.v043.i13}.
#' @export
#' @seealso
#' \code{\link{basdat}}, \code{\link{haartdat}}, \code{\link{ipwplot}}, \code{\link{ipwpoint}}, \code{\link{ipwtm}}, \code{\link{timedat}}, \code{\link{tstartfun}}.
#'
#' @examples
#' # Simulate data with continuous confounder and outcome, binomial exposure.
#' # Marginal causal effect of exposure on outcome: 10.
#' n <- 1000
#' simdat <- data.frame(l = rnorm(n, 10, 5))
#' a.lin <- simdat$l - 10
#' pa <- exp(a.lin)/(1 + exp(a.lin))
#' simdat$a <- rbinom(n, 1, prob = pa)
#' simdat$y <- 10*simdat$a + 0.5*simdat$l + rnorm(n, -10, 5)
#' simdat[1:5,]
#'
#' # Estimate ipw weights.
#' temp <- ipwpoint(
#'   exposure = a,
#'   family = "binomial",
#'   link = "logit",
#'   numerator = ~ 1,
#'   denominator = ~ l,
#'   data = simdat)
#' summary(temp$ipw.weights)
#'
#' # Plot inverse probability weights
#' # ipwplot(weights = temp$ipw.weights, logscale = FALSE,
#' #         main = "Stabilized weights", xlim = c(0, 8))
#'
#' #Examine numerator and denominator models.
#' summary(temp$num.mod)
#' summary(temp$den.mod)
#'
#' #Paste inverse probability weights
#' simdat$sw <- temp$ipw.weights
#'
#' #Marginal structural model for the causal effect of a on y
#' #corrected for confounding by l using inverse probability weighting
#' #with robust standard error from the survey package.
#' if (requireNamespace("survey", quietly = TRUE)) {
#'   library(survey)
#'   msm <- svyglm(y ~ a,
#'                 design = svydesign(~1, weights = ~temp$ipw.weights,
#'                 data = simdat))
#'   summary(msm)
#' }
#' \dontrun{
#' # Compute basic bootstrap confidence interval
#' # require(boot)
#' # boot.fun <- function(dat, index){
#' #   coef(glm(
#' #       formula = y ~ a,
#' #       data = dat[index,],
#' #       weights = ipwpoint(
#' #           exposure = a,
#' #           family = "gaussian",
#' #           numerator = ~ 1,
#' #           denominator = ~ l,
#' #           data = dat[index,])$ipw.weights))[2]
#' #   }
#' # bootres <- boot(simdat, boot.fun, 499);bootres
#' # boot.ci(bootres, type = "basic")
#' }

ipwpoint <- function(
    exposure,
    family,
    link,
    numerator = NULL,
    denominator,
    data,
    trunc = NULL,
    ...)
{
  #save input
  tempcall <- match.call()
  #some basic input checks
  if (!("exposure" %in% names(tempcall))) stop("No exposure variable specified")
  if (!("family" %in% names(tempcall)) | ("family" %in% names(tempcall) & !(tempcall$family %in% c("binomial", "multinomial", "ordinal", "gaussian")))) stop("No valid family specified (\"binomial\", \"multinomial\", \"ordinal\", \"gaussian\")")
  if (tempcall$family == "binomial") {if(!(tempcall$link %in% c("logit", "probit", "cauchit", "log", "cloglog"))) stop("No valid link function specified for family = binomial (\"logit\", \"probit\", \"cauchit\", \"log\", \"cloglog\")")}
  if (tempcall$family == "ordinal" ) {if(!(tempcall$link %in% c("logit", "probit", "cauchit", "cloglog"))) stop("No valid link function specified for family = binomial (\"logit\", \"probit\", \"cauchit\", \"cloglog\")")}
  if (!("denominator" %in% names(tempcall))) stop("No denominator model specified")
  if (!is.null(tempcall$numerator) & !is(eval(tempcall$numerator), "formula")) stop("Invalid numerator formula specified")
  if (!is.null(tempcall$denominator) & !is(eval(tempcall$denominator), "formula")) stop("Invalid denominator formula specified")
  if (tempcall$family %in% c("gaussian") & !("numerator" %in% names(tempcall))) stop("Numerator necessary for family = \"gaussian\"")
  if (!("data" %in% names(tempcall))) stop("No data specified")
  if (!is.null(tempcall$trunc)) {if(tempcall$trunc < 0 | tempcall$trunc > 0.5) stop("Invalid truncation percentage specified (0-0.5)")}
  #make new dataframe for newly computed variables, to prevent variable name conflicts
  tempdat <- data.frame(
    exposure = data[,as.character(tempcall$exposure)]
  )
  #weights binomial
  if (tempcall$family == "binomial") {
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
      tempdat$w.numerator <- vector("numeric", nrow(tempdat))
      tempdat$w.numerator[tempdat$exposure == 0] <- 1 - predict.glm(mod1, type = "response")[tempdat$exposure == 0]
      tempdat$w.numerator[tempdat$exposure == 1] <- predict.glm(mod1, type = "response")[tempdat$exposure == 1]
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
    tempdat$w.denominator <- vector("numeric", nrow(tempdat))
    tempdat$w.denominator[tempdat$exposure == 0] <- 1 - predict.glm(mod2, type = "response")[tempdat$exposure == 0]
    tempdat$w.denominator[tempdat$exposure == 1] <- predict.glm(mod2, type = "response")[tempdat$exposure == 1]
    mod2$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = "")))
    mod2$call$family <- tempcall$link
    mod2$call$data <- tempcall$data
    tempdat$ipw.weights <- tempdat$w.numerator/tempdat$w.denominator
  }
  #weights multinomial
  if (tempcall$family == "multinomial") {
    if (is.null(tempcall$numerator)) tempdat$p.numerator <- 1
    else {
      mod1 <- multinom(
        formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = ""))),
        data = data,
        na.action = na.fail,
        ...)
      pred1 <- as.data.frame(predict(mod1, type = "probs"))
      tempdat$w.numerator <- vector("numeric", nrow(tempdat))
      for (i in 1:length(unique(tempdat$exposure)))tempdat$w.numerator[with(tempdat, exposure == sort(unique(tempdat$exposure))[i])] <- pred1[tempdat$exposure == sort(unique(tempdat$exposure))[i],i]
      mod1$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = "")))
      mod1$call$data <- tempcall$data
    }
    mod2 <- multinom(
      formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = ""))),
      data = data,
      na.action = na.fail,
      ...)
    pred2 <- as.data.frame(predict(mod2, type = "probs"))
    tempdat$w.denominator <- vector("numeric", nrow(tempdat))
    for (i in 1:length(unique(tempdat$exposure)))tempdat$w.denominator[with(tempdat, exposure == sort(unique(tempdat$exposure))[i])] <- pred2[tempdat$exposure == sort(unique(tempdat$exposure))[i],i]
    mod2$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = "")))
    mod2$call$data <- tempcall$data
    tempdat$ipw.weights <- tempdat$w.numerator/tempdat$w.denominator
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
        na.action = na.fail,
        ...)
      pred1 <- as.data.frame(predict(mod1, type = "probs"))
      tempdat$w.numerator <- vector("numeric", nrow(tempdat))
      for (i in 1:length(unique(tempdat$exposure)))tempdat$w.numerator[with(tempdat, exposure == sort(unique(tempdat$exposure))[i])] <- pred1[tempdat$exposure == sort(unique(tempdat$exposure))[i],i]
      mod1$call$formula <- eval(parse(text = paste("as.factor(", deparse(tempcall$exposure, width.cutoff = 500), ")", deparse(tempcall$numerator, width.cutoff = 500), sep = "")))
      mod1$call$data <- tempcall$data
      mod1$call$method <- m
    }
    mod2 <- polr(
      formula = eval(parse(text = paste("as.factor(", deparse(tempcall$exposure, width.cutoff = 500), ")", deparse(tempcall$denominator, width.cutoff = 500), sep = ""))),
      data = data,
      method = m,
      na.action = na.fail,
      ...)
    pred2 <- as.data.frame(predict(mod2, type = "probs"))
    tempdat$w.denominator <- vector("numeric", nrow(tempdat))
    for (i in 1:length(unique(tempdat$exposure)))tempdat$w.denominator[with(tempdat, exposure == sort(unique(tempdat$exposure))[i])] <- pred2[tempdat$exposure == sort(unique(tempdat$exposure))[i],i]
    mod2$call$formula <- eval(parse(text = paste("as.factor(", deparse(tempcall$exposure, width.cutoff = 500), ")", deparse(tempcall$denominator, width.cutoff = 500), sep = "")))
    mod2$call$data <- tempcall$data
    mod2$call$method <- m
    tempdat$ipw.weights <- tempdat$w.numerator/tempdat$w.denominator
  }
  #weights gaussian
  if (tempcall$family == "gaussian") {
    mod1 <- glm(
      formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = ""))),
      data = data,
      na.action = na.fail,
      ...)
    tempdat$w.numerator <- dnorm(tempdat$exposure, predict(mod1), sd(mod1$residuals))
    mod1$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = "")))
    mod1$call$data <- tempcall$data
    mod2 <- glm(
      formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = ""))),
      data = data,
      na.action = na.fail,
      ...)
    tempdat$w.denominator <- dnorm(tempdat$exposure, predict(mod2), sd(mod2$residuals))
    mod2$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = "")))
    mod2$call$data <- tempcall$data
    tempdat$ipw.weights <- tempdat$w.numerator/tempdat$w.denominator
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
    if (is.null(tempcall$numerator)) return(list(ipw.weights = tempdat$ipw.weights, call = tempcall, den.mod = mod2))
    else return(list(ipw.weights = tempdat$ipw.weights, call = tempcall, num.mod = mod1, den.mod = mod2))
  }
  else{
    if (is.null(tempcall$numerator)) return(list(ipw.weights = tempdat$ipw.weights, weights.trunc = tempdat$weights.trunc, call = tempcall, den.mod = mod2))
    else return(list(ipw.weights = tempdat$ipw.weights, weights.trunc = tempdat$weights.trunc, call = tempcall, num.mod = mod1, den.mod = mod2))
  }
}
