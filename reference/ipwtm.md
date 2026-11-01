# Estimate Inverse Probability Weights (Time Varying)

Estimates inverse probability weights to fit marginal structural models
with time-varying exposures and time-varying confounders across
follow-up intervals. Supports binomial, survival, multinomial, ordinal,
and gaussian families.

## Usage

``` r
ipwtm(
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
  ...
)
```

## Arguments

- exposure:

  Name of the exposure column in `data` (unquoted or as a string).

- family:

  Model family: `"binomial"`, `"survival"`, `"multinomial"`,
  `"ordinal"`, or `"gaussian"`.

- link:

  Link function for binomial or ordinal models.

- numerator:

  Formula for numerator model. If `NULL`, unstabilized weights are
  computed.

- denominator:

  Formula for denominator model. For `family = "binomial"`, smooth terms
  (`s()`, `te()`, `ti()`, `t2()`) are allowed; the model is then fitted
  with [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html).

- id:

  Name of the subject identifier column in `data`.

- tstart:

  Name of the interval start time column (required for
  `family = "survival"`).

- timevar:

  Name of the follow-up (interval end) time column.

- type:

  Exposure switch type: `"first"`, `"cens"`, or `"all"`. `"all"` is
  available for all families. With `family = "survival"` and
  `type = "all"`, `exposure` is an event indicator (1 when an exposure
  event, such as starting treatment, happens at the end of the interval)
  that may occur more than once per subject: an Andersen-Gill Cox model
  is fitted on all rows and every row contributes to the weights
  (Graffeo et al. 2018). For treatment that is started and stopped,
  compute weights separately for starting (on the rows off treatment)
  and for stopping (on the rows on treatment).

- data:

  Data frame containing all required variables.

- corstr:

  Correlation structure for GEE when `family = "gaussian"` (default:
  `"ar1"`).

- trunc:

  Optional truncation percentile between 0 and 0.5.

- ...:

  Additional arguments passed to model fitting engines.

## Value

A list of class `"ipwtm"` containing:

- ipw.weights:

  Numeric vector of inverse probability weights in original row order.

- weights.trunc:

  Numeric vector of truncated weights (if `trunc` is specified).

- call:

  Original function call.

- selvar:

  Selection indicator vector (1 = active observation, 0 = post-switch).

- num.mod:

  Numerator model object.

- den.mod:

  Denominator model object.

## References

Graffeo N., Latouche A., Geskus R.B. & Chevret S. (2018). Modeling
time-varying exposure using inverse probability of treatment weights.
*Biometrical Journal*, **60**(2), 323-332.
[doi:10.1002/bimj.201600223](https://doi.org/10.1002/bimj.201600223) .

Van der Wal W.M. & Geskus R.B. (2011). ipw: An R Package for Inverse
Probability Weighting. *Journal of Statistical Software*, **43**(13),
1-23. [doi:10.18637/jss.v043.i13](https://doi.org/10.18637/jss.v043.i13)
.

## Examples

``` r
########################################################################
#EXAMPLE 1

#Load longitudinal data from HIV positive individuals.
data(haartdat)

#CD4 is confounder for the effect of initiation of HAART therapy on mortality.
#Estimate inverse probability weights to correct for confounding.
#Exposure allocation model is Cox proportional hazards model.
temp <- ipwtm(
  exposure = haartind,
  family = "survival",
  numerator = ~ sex + age,
  denominator = ~ sex + age + cd4.sqrt,
  id = patient,
  tstart = tstart,
  timevar = fuptime,
  type = "first",
  data = haartdat)

#plot inverse probability weights
graphics.off()
ipwplot(weights = temp$ipw.weights, timevar = haartdat$fuptime,
        binwidth = 100, ylim = c(-1.5, 1.5), main = "Stabilized inverse probability weights")

#CD4 count has an effect both on dropout and mortality, which causes informative censoring.
#Use inverse probability of censoring weighting to correct for effect of CD4 on dropout.
#Use Cox proportional hazards model for dropout.
temp2 <- ipwtm(
  exposure = dropout,
  family = "survival",
  numerator = ~ sex + age,
  denominator = ~ sex + age + cd4.sqrt,
  id = patient,
  tstart = tstart,
  timevar = fuptime,
  type = "cens",
  data = haartdat)

#plot inverse probability of censoring weights
graphics.off()
ipwplot(weights = temp2$ipw.weights, timevar = haartdat$fuptime,
        binwidth = 100, ylim = c(-1.5, 1.5),
        main = "Stabilized inverse probability of censoring weights")

#MSM for the causal effect of initiation of HAART on mortality.
#Corrected both for confounding and informative censoring.
#With robust standard error obtained using cluster().
require(survival)
summary(coxph(Surv(tstart, fuptime, event) ~ haartind + cluster(patient),
              data = haartdat, weights = temp$ipw.weights*temp2$ipw.weights))
#> Call:
#> coxph(formula = Surv(tstart, fuptime, event) ~ haartind, data = haartdat, 
#>     weights = temp$ipw.weights * temp2$ipw.weights, cluster = patient)
#> 
#>   n= 19175, number of events= 31 
#> 
#>             coef exp(coef) se(coef) robust se      z Pr(>|z|)  
#> haartind -0.9363    0.3921   0.4299    0.4527 -2.068   0.0386 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>          exp(coef) exp(-coef) lower .95 upper .95
#> haartind    0.3921      2.551    0.1614    0.9521
#> 
#> Concordance= 0.597  (se = 0.029 )
#> Likelihood ratio test= 5.57  on 1 df,   p=0.02
#> Wald test            = 4.28  on 1 df,   p=0.04
#> Score (logrank) test = 5.07  on 1 df,   p=0.02,   Robust = 4.83  p=0.03
#> 
#>   (Note: the likelihood ratio and score tests assume independence of
#>      observations within a cluster, the Wald and robust score tests do not).

#uncorrected model
summary(coxph(Surv(tstart, fuptime, event) ~ haartind, data = haartdat))
#> Call:
#> coxph(formula = Surv(tstart, fuptime, event) ~ haartind, data = haartdat)
#> 
#>   n= 19175, number of events= 31 
#> 
#>             coef exp(coef) se(coef)      z Pr(>|z|)
#> haartind -0.5862    0.5565   0.4368 -1.342     0.18
#> 
#>          exp(coef) exp(-coef) lower .95 upper .95
#> haartind    0.5565      1.797    0.2364      1.31
#> 
#> Concordance= 0.569  (se = 0.032 )
#> Likelihood ratio test= 1.97  on 1 df,   p=0.2
#> Wald test            = 1.8  on 1 df,   p=0.2
#> Score (logrank) test = 1.85  on 1 df,   p=0.2
#> 

########################################################################
#EXAMPLE 2

data(basdat)
data(timedat)

#Aim: to model the causal effect of active tuberculosis (TB) on mortality.
#Longitudinal CD4 is a confounder as well as intermediate for the effect of TB.

#process original measurements
#check for ties (not allowed)
table(duplicated(timedat[,c("id", "fuptime")]))
#> 
#> FALSE 
#>  6291 
#take square root of CD4 because of skewness
timedat$cd4.sqrt <- sqrt(timedat$cd4count)
#add TB time to dataframe
timedat <- merge(timedat, basdat[,c("id", "Ttb")], by = "id", all.x = TRUE)
#compute TB status
timedat$tb.lag <- ifelse(with(timedat, !is.na(Ttb) & fuptime > Ttb), 1, 0)
#longitudinal CD4-model
require(nlme)
#> Loading required package: nlme
cd4.lme <- lme(cd4.sqrt ~ fuptime + tb.lag, random = ~ fuptime | id,
               data = timedat)

#build new dataset:
#rows corresponding to TB-status switches, and individual end times
times <- sort(unique(c(basdat$Ttb, basdat$Tend)))
startstop <- data.frame(
  id = rep(basdat$id, each = length(times)),
  fuptime = rep(times, nrow(basdat)))
#add baseline data to dataframe
startstop <- merge(startstop, basdat, by = "id", all.x = TRUE)
#limit individual follow-up using Tend
startstop <- startstop[with(startstop, fuptime <= Tend),]
startstop$tstart <- tstartfun(id, fuptime, startstop) #compute tstart (?tstartfun)
#indicate TB status
startstop$tb <- ifelse(with(startstop, !is.na(Ttb) & fuptime >= Ttb), 1, 0)
#indicate TB status at previous time point
startstop$tb.lag <- ifelse(with(startstop, !is.na(Ttb) & fuptime > Ttb), 1, 0)
#indicate death
startstop$event <- ifelse(with(startstop, !is.na(Tdeath) & fuptime >= Tdeath),
                          1, 0)
#impute CD4, based on TB status at previous time point.
startstop$cd4.sqrt <- predict(cd4.lme,
                              newdata = data.frame(id = startstop$id,
                                                   fuptime = startstop$fuptime,
                                                   tb.lag = startstop$tb.lag))
#compute inverse probability weights
temp <- ipwtm(
  exposure = tb,
  family = "survival",
  numerator = ~ 1,
  denominator = ~ cd4.sqrt,
  id = id,
  tstart = tstart,
  timevar = fuptime,
  type = "first",
  data = startstop)
summary(temp$ipw.weights)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.2117  0.9409  0.9908  1.0373  1.0125 13.6539 
ipwplot(weights = temp$ipw.weights, timevar = startstop$fuptime, binwidth = 100)

#models
#IPW-fitted MSM, using cluster() to obtain robust standard error estimate
require(survival)
summary(coxph(Surv(tstart, fuptime, event) ~ tb + cluster(id),
              data = startstop, weights = temp$ipw.weights))
#> Call:
#> coxph(formula = Surv(tstart, fuptime, event) ~ tb, data = startstop, 
#>     weights = temp$ipw.weights, cluster = id)
#> 
#>   n= 84494, number of events= 134 
#> 
#>      coef exp(coef) se(coef) robust se     z Pr(>|z|)   
#> tb 0.8127    2.2541   0.1901    0.2599 3.127  0.00177 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>    exp(coef) exp(-coef) lower .95 upper .95
#> tb     2.254     0.4436     1.354     3.751
#> 
#> Concordance= 0.57  (se = 0.024 )
#> Likelihood ratio test= 17.4  on 1 df,   p=3e-05
#> Wald test            = 9.78  on 1 df,   p=0.002
#> Score (logrank) test = 19.09  on 1 df,   p=1e-05,   Robust = 8.46  p=0.004
#> 
#>   (Note: the likelihood ratio and score tests assume independence of
#>      observations within a cluster, the Wald and robust score tests do not).
#unadjusted
summary(coxph(Surv(tstart, fuptime, event) ~ tb, data = startstop))
#> Call:
#> coxph(formula = Surv(tstart, fuptime, event) ~ tb, data = startstop)
#> 
#>   n= 84494, number of events= 134 
#> 
#>      coef exp(coef) se(coef)     z Pr(>|z|)    
#> tb 1.4954    4.4612   0.1811 8.257   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>    exp(coef) exp(-coef) lower .95 upper .95
#> tb     4.461     0.2242     3.128     6.362
#> 
#> Concordance= 0.636  (se = 0.024 )
#> Likelihood ratio test= 63.84  on 1 df,   p=1e-15
#> Wald test            = 68.17  on 1 df,   p=<2e-16
#> Score (logrank) test = 79.8  on 1 df,   p=<2e-16
#> 
#adjusted using conditioning: part of the effect of TB is adjusted away
summary(coxph(Surv(tstart, fuptime, event) ~ tb + cd4.sqrt, data = startstop))
#> Call:
#> coxph(formula = Surv(tstart, fuptime, event) ~ tb + cd4.sqrt, 
#>     data = startstop)
#> 
#>   n= 84494, number of events= 134 
#> 
#>              coef exp(coef) se(coef)      z Pr(>|z|)    
#> tb        0.24618   1.27913  0.24288  1.014    0.311    
#> cd4.sqrt -0.24444   0.78314  0.03313 -7.378 1.61e-13 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>          exp(coef) exp(-coef) lower .95 upper .95
#> tb          1.2791     0.7818    0.7947    2.0590
#> cd4.sqrt    0.7831     1.2769    0.7339    0.8357
#> 
#> Concordance= 0.738  (se = 0.025 )
#> Likelihood ratio test= 124  on 2 df,   p=<2e-16
#> Wald test            = 118.9  on 2 df,   p=<2e-16
#> Score (logrank) test = 136.6  on 2 df,   p=<2e-16
#> 

## Not run:
#compute bootstrap CI for TB parameter (takes a few hours)
#taking into account the uncertainty introduced by modelling longitudinal CD4
#taking into account the uncertainty introduced by estimating the inverse probability weights
#robust with regard to weights unequal to 1
#  require(boot)
#  boot.fun <- function(data, index, data.tm){
#     data.samp <- data[index,]
#     data.samp$id.samp <- 1:nrow(data.samp)
#     data.tm.samp <- do.call("rbind", lapply(data.samp$id.samp, function(id.samp) {
#       cbind(data.tm[data.tm$id == data.samp$id[data.samp$id.samp == id.samp],],
#         id.samp = id.samp)
#       }
#     ))
#     cd4.lme <- lme(cd4.sqrt ~ fuptime + tb.lag, random = ~ fuptime | id.samp, data = data.tm.samp)
#     times <- sort(unique(c(data.samp$Ttb, data.samp$Tend)))
#     startstop.samp <- data.frame(id.samp = rep(data.samp$id.samp, each = length(times)),
#                                  fuptime = rep(times, nrow(data.samp)))
#     startstop.samp <- merge(startstop.samp, data.samp, by = "id.samp", all.x = TRUE)
#     startstop.samp <- startstop.samp[with(startstop.samp, fuptime <= Tend),]
#     startstop.samp$tstart <- tstartfun(id.samp, fuptime, startstop.samp)
#     startstop.samp$tb <- ifelse(with(startstop.samp, !is.na(Ttb) & fuptime >= Ttb), 1, 0)
#     startstop.samp$tb.lag <- ifelse(with(startstop.samp, !is.na(Ttb) & fuptime > Ttb), 1, 0)
#     startstop.samp$event <- ifelse(with(startstop.samp, !is.na(Tdeath) & fuptime >= Tdeath), 1, 0)
#     startstop.samp$cd4.sqrt <- predict(cd4.lme, newdata = data.frame(id.samp =
#       startstop.samp$id.samp, fuptime = startstop.samp$fuptime, tb.lag = startstop.samp$tb.lag))
#
#     return(coef(coxph(Surv(tstart, fuptime, event) ~ tb, data = startstop.samp,
#        weights = ipwtm(
#             exposure = tb,
#             family = "survival",
#             numerator = ~ 1,
#             denominator = ~ cd4.sqrt,
#             id = id.samp,
#             tstart = tstart,
#             timevar = fuptime,
#             type = "first",
#             data = startstop.samp)$ipw.weights))[1])
#     }
#  bootres <- boot(data = basdat, statistic = boot.fun, R = 999, data.tm = timedat)
#  bootres
#  boot.ci(bootres, type = "basic")
#
## End(Not run)
```
