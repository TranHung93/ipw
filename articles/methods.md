# Methods behind the new features: formulas, code and rationale

Each new way of estimating weights added since ipw 1.3.2 comes down to
the same step: estimate the probability, or density, of the exposure
each person actually received, then take a ratio. This article gives,
for each new feature, the formula, the code in the package and why that
approach was chosen. Section 6 covers recurrent exposure events (Graffeo
et al. 2018), section 7
[`ipwprobtm()`](https://tranhung93.github.io/ipw/reference/ipwprobtm.md),
which turns probabilities you estimated yourself into time-varying
weights, and section 8 the standard errors of
[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md).
How to use the features is shown in the [main
vignette](https://tranhung93.github.io/ipw/articles/ipw.md); whether
they recover known effects is checked in the [simulation
study](https://tranhung93.github.io/ipw/articles/simulations.md).

## Notation and the shared formula

Every family only has to supply *f*, the probability or density of the
observed exposure; the weight formula itself is the same for all of
them.

- *A*: the exposure each person received; *L*: confounders in the
  denominator model; *V*: covariates in the numerator model (often
  none).
- *f*(*a* \| ·): the probability of *a* for a discrete exposure, the
  density of *a* for a continuous one.

Point treatment, stabilized and unstabilized:

``` math
SW_i = \frac{f(A_i \mid V_i)}{f(A_i \mid L_i)}, \qquad W_i = \frac{1}{f(A_i \mid L_i)}
```

Time-varying, for subject *i* at time *t*, with $`\bar{A}`$ and
$`\bar{L}`$ the history up to that interval:

``` math
SW_i(t) = \prod_{k \le t} \frac{f(A_{ik} \mid \bar{A}_{i,k-1}, V_i)}{f(A_{ik} \mid \bar{A}_{i,k-1}, \bar{L}_{ik})}
```

In
[`ipwpoint()`](https://tranhung93.github.io/ipw/reference/ipwpoint.md),
each family’s branch returns `dens` for the numerator and denominator
models, and one line forms the weight:

``` r

ipw_weights <- num$dens / den$dens   # num$dens = 1 when numerator = NULL
```

## 1. `estimand`: ATT, ATC and overlap weights for a binomial point treatment

The `estimand` argument changes which population the weights standardise
to, using the same fitted propensity score *e*(*L*) = P(*A* = 1 \| *L*)
from the denominator model.

**Formula.**

``` math
\text{ATE: } w = \frac{A}{e} + \frac{1-A}{1-e} \qquad \text{ATT: } w = A + (1-A)\frac{e}{1-e} \qquad \text{ATC: } w = A\frac{1-e}{e} + (1-A) \qquad \text{ATO: } w = A(1-e) + (1-A)\,e
```

**Code.**

``` r

ipw_weights <- switch(estimand,
  ATE = num$dens / den$dens,
  ATT = ifelse(exp_vec == 1, 1, den$p / (1 - den$p)),
  ATC = ifelse(exp_vec == 1, (1 - den$p) / den$p, 1),
  ATO = ifelse(exp_vec == 1, 1 - den$p, den$p)
)
```

`den$p` is the fitted P(*A* = 1 \| *L*) kept from the binomial
denominator model; `"ATT"`, `"ATC"` and `"ATO"` require
`numerator = NULL`.

**Why this approach.** Each estimand answers a different clinical
question; the formulas only change who gets how much weight.

- **ATE (default):** “What would happen if everyone were treated,
  compared with no one?” Each person is weighted to stand for the whole
  study population.
- **ATT:** “Did the treatment help the people who actually received it?”
  Treated people keep weight 1; untreated people are reweighted to look
  like the treated group. Useful when treatment is only realistic for
  some patients.
- **ATC:** “What would the treatment have done for the people who did
  not receive it?” The mirror image of the ATT: untreated people keep
  weight 1 and treated people are reweighted to look like the untreated
  group. Useful when the question is whether to extend a treatment to
  those not yet receiving it.
- **ATO (overlap):** “What is the effect in patients for whom doctors
  could have gone either way?” Every weight lies between 0 and 1, so no
  single person can dominate. This helps when some patients are almost
  always or almost never treated: there, ATE weights become huge and
  need truncation. With a logistic model, overlap weights also make the
  average of every covariate in the model exactly equal between the
  groups.
- **Why no numerator:** ATT, ATC and overlap weights are not a ratio of
  two probabilities, so stabilizing them has no meaning. ipw gives an
  error instead of silently ignoring a numerator.
- **Which to use:** choose by the clinical question first. If ATE
  weights are extreme, consider ATO, and state clearly which patients
  the result describes.

**Caveat.** The overlap population is defined by the propensity score
itself, so results describe that population, not a pre-specified one;
this should be stated when reporting ATO estimates.

### ATT and overlap weights for multinomial and ordinal exposures

For an exposure with levels $`k = 1, \dots, K`$, let $`P_i(k)`$ be the
fitted probability of level $`k`$ for person *i* (from
[`nnet::multinom`](https://rdrr.io/pkg/nnet/man/multinom.html) or
[`MASS::polr`](https://rdrr.io/pkg/MASS/man/polr.html)) and $`A_i`$ the
observed level.

**Formula.** With a focal level $`f`$ (argument `focal`), and for
generalized overlap weights (Li & Li 2019):

``` math
\text{ATT: } w_i = \frac{P_i(f)}{P_i(A_i)} \qquad
\text{ATO: } w_i = \frac{1 / P_i(A_i)}{\sum_{k=1}^{K} 1 / P_i(k)}
```

People at the focal level get ATT weight 1. With $`K = 2`$,
$`P_i(2) = e`$ and $`f = 2`$ these are the binomial formulas above: ATT
gives $`e/(1-e)`$ to the others, and ATO gives
$`(1/e) / (1/e + 1/(1-e)) = 1 - e`$ to level 2 and $`e`$ to level 1.

**Code.** The full matrix of fitted probabilities is built by the
internal helper `.prob_matrix()` (also used for the probability of the
observed level):

``` r

ipw_weights <- switch(estimand,
  ATT = ifelse(as.character(exp_vec) == focal, 1, den$pmat[, focal] / den$dens),
  ATO = (1 / den$dens) / rowSums(1 / den$pmat)
)
```

**Why this approach.** The ATT for a categorical exposure compares every
level with what would have happened to the people at the focal level, so
the weights reweight each other level to resemble the focal group. The
generalized overlap weights target the people with substantial
probability of every level; like binomial overlap weights they are
bounded, here by 1. `"ATC"` is not offered for categorical exposures:
with more than two levels there is no single control group, and the
effect in the controls is the ATT with `focal` set to the control level.

## 2. Skewed continuous exposures: `"gamma"` and `"lognormal"`

For a positive, right-skewed exposure (dose, duration, a biomarker), *f*
is a gamma or log-normal density instead of a normal one.

**Formula.** Gamma with log link, mean μ_(i) = exp(*L*_(i)ᵀβ) and shape
*k* = 1/φ, where φ is the dispersion estimated by the GLM:

``` math
f(a \mid L) = \frac{(k/\mu)^{k}}{\Gamma(k)}\, a^{k-1} e^{-k a / \mu}
```

Log-normal: a normal linear model for log *A* with mean μ_(i) =
*L*_(i)ᵀβ and standard deviation σ:

``` math
f(a \mid L) = \frac{1}{a \sigma \sqrt{2\pi}} \exp\!\left(-\frac{(\log a - \mu)^2}{2\sigma^2}\right)
```

The factor 1/*a* is the same in numerator and denominator, so the
stabilized log-normal weight equals the gaussian weight computed on log
*A*.

**Code.**

``` r

# gamma
mod   <- do.call(stats::glm, c(base, list(family = stats::Gamma(link = "log")), dots))
shape <- 1 / summary(mod)$dispersion
mu    <- as.numeric(stats::predict(mod, type = "response"))
dens  <- stats::dgamma(exp_vec, shape = shape, rate = shape / mu)

# lognormal: the model formula is log(exposure) ~ ...
mod  <- do.call(stats::glm, c(base, list(family = stats::gaussian()), dots))
dens <- stats::dlnorm(exp_vec, meanlog = as.numeric(stats::predict(mod)),
                      sdlog = stats::sd(stats::residuals(mod)))
```

The exposure must be positive, and a numerator model is required.

**Why this approach.**

- **The problem:** doses, durations and many lab values are always
  positive and have a long right tail: most values are small, a few are
  very large. The normal model is symmetric and allows negative values,
  so it fits the tails badly, and the tails are where the weights are
  largest.
- **Gamma:** suits exposures whose spread grows with their size, for
  example a duration that varies by about 20% whatever its typical
  length. Common for durations, doses and costs.
- **Log-normal:** suits exposures that look bell-shaped after taking the
  log, such as viral load and many biomarkers. In practice it is the
  gaussian model applied to log(exposure).
- **Which to use:** draw a histogram of log(exposure). If it looks
  bell-shaped, use log-normal; otherwise try gamma. You can also compare
  the AIC of the two denominator models.
- **Why the numerator is required:** for a continuous exposure, *f* is a
  density, not a probability. Its size depends on the units (mg or g),
  so 1/*f* has no fixed meaning and can be huge. Dividing by a numerator
  density, as already required for gaussian, cancels the units and keeps
  the weights around 1.

**Implementation notes.** The gamma shape uses the moment (Pearson)
estimate of φ from
[`summary.glm()`](https://rdrr.io/r/stats/summary.glm.html), not maximum
likelihood
([`MASS::gamma.shape()`](https://rdrr.io/pkg/MASS/man/gamma.shape.glm.html)
would be the alternative). The log-normal σ uses
[`sd()`](https://rdrr.io/r/stats/sd.html) of the residuals (divisor *n*
− 1, not *n* − *p*), to match the existing gaussian family.

**Caveat.** In the simulation study a small bias (4–6%) remained with
the correct family. Weights computed from the true densities had the
same bias, so it is a property of weighting with a skewed exposure in
finite samples, and it shrinks slowly as the sample grows.

## 3. Count exposures: `"poisson"` and `"negbin"`

For a count exposure (number of doses, visits), *f* is the probability
mass of the observed count, from a log-linear model for its mean.

**Formula.** Both models have a log link, μ_(i) = exp(*L*_(i)ᵀβ).
Poisson, then negative binomial with dispersion parameter θ (variance
μ + μ²/θ):

``` math
f(a \mid L) = \frac{e^{-\mu} \mu^{a}}{a!}
```

``` math
f(a \mid L) = \frac{\Gamma(a + \theta)}{\Gamma(\theta)\, a!} \left(\frac{\theta}{\theta + \mu}\right)^{\theta} \left(\frac{\mu}{\theta + \mu}\right)^{a}
```

**Code.**

``` r

# poisson
mod  <- .fit_glm_or_gam(c(base, list(family = stats::poisson()), dots))
dens <- stats::dpois(exp_vec, lambda = as.numeric(stats::predict(mod, type = "response")))

# negbin: theta estimated by glm.nb, separately for numerator and denominator
mod  <- do.call(MASS::glm.nb, c(base, dots))
dens <- stats::dnbinom(exp_vec, size = mod$theta,
                       mu = as.numeric(stats::predict(mod, type = "response")))
```

The exposure must be a non-negative integer.

**Why this approach.**

- **The problem:** a count (number of doses: 0, 1, 2, …) is neither
  yes/no nor a smooth measurement. A normal model allows impossible
  values such as −1.3 doses and gives poor probabilities for large
  counts, which is exactly where weights become large.
- **Poisson:** the standard model for counts. It assumes the spread of
  the counts equals their average.
- **Negative binomial:** real counts often vary more than that: many
  patients get no doses and a few get very many. A Poisson model then
  gives too small a chance to these unusual counts, and so too large a
  weight. The negative binomial adds one parameter for the extra spread.
- **Which to use:** start with Poisson; if the variance of the exposure
  is clearly larger than its mean, use the negative binomial. In the
  simulation study, the negative binomial also worked for Poisson counts
  (bias 3.0%), while Poisson weights for overdispersed counts were badly
  biased.
- **Why the numerator is optional:** for counts, *f* is a true
  probability between 0 and 1, as for a yes/no exposure, so unstabilized
  weights 1/*f* still make sense.

**Caveat.** Zero-inflated or truncated counts are not covered.

## 4. Smooth terms with mgcv

When a formula contains `s()`, `te()`, `ti()` or `t2()`, the exposure
model is fitted as a generalized additive model, so a confounder’s
effect can be nonlinear without hand-built splines.

**Formula.** With link *g* and smooth functions *f_(j)* represented by
penalized regression splines:

``` math
g\big(E[A \mid L]\big) = \beta_0 + \sum_j f_j(L_j) + L_{\text{lin}}^\top \gamma
```

Each *f_(j)* is fitted by minimising the penalized deviance below; λ_(j)
sets its smoothness and is estimated by mgcv (default
`method = "GCV.Cp"`; `method = "REML"` can be passed through `...`).

``` math
D(\beta) + \sum_j \lambda_j \int f_j''(x)^2 \, dx
```

The weights then follow the notation section, with *f* taken from the
GAM’s fitted values.

**Code.** One helper chooses [`glm()`](https://rdrr.io/r/stats/glm.html)
or [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html), and gives the
formula an environment in which `s()` and the others are mgcv’s, so
users need not attach mgcv:

``` r

.fit_glm_or_gam <- function(args) {
  if (.uses_smooths(args$formula)) {
    .require_mgcv()
    args$formula <- .gam_formula(args$formula)
    args$na.action <- NULL
    do.call(mgcv::gam, args)
  } else {
    do.call(stats::glm, args)
  }
}
```

**Why this approach.**

- **The problem:** a regression assumes each confounder has a
  straight-line effect (on the logit scale for a yes/no exposure). That
  is often wrong: the chance of starting treatment may rise steeply at
  low CD4 counts and then level off. A straight line then gives the
  wrong probabilities, and some confounding remains in the weighted
  analysis.
- **The old fix:** add splines or polynomials by hand. That means
  choosing knots or degrees, and it is easy to make the curve too stiff
  or too wiggly.
- **What a GAM does:** writing `s(cd4)` lets the data decide the shape
  of the curve, while a penalty keeps it smooth so it does not follow
  random noise.
- **Why mgcv:** it is the standard, well-tested R package for GAMs and
  is shipped with R as a recommended package, so most users already have
  it.
- **Why only some families:** for binomial, gaussian and Poisson, a GAM
  is a direct replacement for the GLM already used. Multinomial,
  ordinal, gamma and negative binomial would need different fitting
  code, so they are left for later.
- **Tip:** use `s()` for continuous confounders where a curved effect is
  plausible, and look at `plot(result$den.mod)` to see the fitted
  curves.

**Implementation note.** For gaussian exposures, the GAM uses σ =
`sqrt(sig2)` (residual sum of squares over *n* minus the effective
degrees of freedom), while the GLM path keeps `sd(residuals)` from
earlier versions, so the two are not on exactly the same footing.

## 5. `ipwtm()`: type “all” for multinomial and ordinal exposures

Every interval now contributes the probability of the observed category,
so the weight follows exposures that move between categories in both
directions.

**Formula.** With categories 1, …, *J*, the multinomial model is a
baseline-category logit; the ordinal model a cumulative link model with
link *g* (logit, probit, cauchit or cloglog):

``` math
P(A_{ik} = j \mid L_{ik}) = \frac{\exp(L_{ik}^\top \beta_j)}{\sum_{m=1}^{J} \exp(L_{ik}^\top \beta_m)}, \quad \beta_1 = 0
```

``` math
P(A_{ik} \le j \mid L_{ik}) = g^{-1}(\zeta_j - L_{ik}^\top \beta)
```

The weight is the time-varying product in the notation section, with
*f*(*A_(ik)* \| ·) the fitted probability of the category observed in
interval *k*.

**Code.** With `type = "all"` every row is selected (`selvar = 1`), so
the existing categorical engine fits
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) or
[`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html) on all rows:

``` r

p_vec[active] <- .prob_of_observed(stats::predict(mod, type = "probs"),
                                   data[[exp_col]][active], mod$lev)
weights <- unsplit(lapply(split(p_ratio, data[[id_col]]), cumprod), data[[id_col]])
```

**Why this approach.**

- **The problem:** some exposures change more than once. A patient’s
  adherence can be high, then low, then high again. With
  `type = "first"`, ipw stops following a person after their first
  change, so the later changes are ignored.
- **What it does:** every visit counts. At each visit the model asks:
  given this person’s history, how likely was the category they actually
  had? These chances are multiplied over time to give the weight.
- **Why this way:** it is the method ipw already uses for yes/no
  exposures with `type = "all"`, extended to more than two categories.
  The multinomial and ordinal models were already in the package, so
  little new code was needed.
- **Which model to use:** multinomial when the categories have no order
  (drug A, B or C); ordinal when they do (low, medium, high). The
  ordinal model uses one effect per covariate instead of one per
  category, so it needs less data.

**Caveat.** The models are pooled over time, so the user must put past
exposure (for example the previous category) and time in the formulas;
the function does not add them.

## 6. `ipwtm()`: type “all” for recurrent exposure events (survival)

Exposure events can now happen more than once per person, for example
repeated prescriptions or hospital admissions. This follows Graffeo,
Latouche, Geskus and Chevret (2018).

**Formula.** An Andersen–Gill Cox model for the exposure events is
fitted on all rows, with baseline cumulative hazard Λ₀(*t*). Row *k* of
person *i* covers the interval (*s_(ik)*, *t_(ik)*\] and has covariates
*L_(ik)*:

``` math
P(A_{ik} = 0 \mid L_{ik}) = \exp\big\{-\big[\Lambda_0(t_{ik}) - \Lambda_0(s_{ik})\big]\, e^{\beta^\top L_{ik}}\big\}
```

``` math
P(A_{ik} = 1 \mid L_{ik}) = 1 - P(A_{ik} = 0 \mid L_{ik})
```

The weight is the time-varying product in the notation section, with
*f*(*A_(ik)* \| ·) the probability of the status observed in row *k*
(event or no event).

**Code.** The model is `coxph(Surv(tstart, timevar, exposure) ~ ...)` on
all rows; each row’s probability uses the baseline hazard over that
row’s own interval:

``` r

bh <- survival::basehaz(mod, centered = TRUE)
d_cum_haz <- cum_haz_at(data[[time_col]]) - cum_haz_at(data[[start_col]])
h <- d_cum_haz * stats::predict(mod, type = "risk", centered = TRUE)
ifelse(data[[exp_col]] == 1, 1 - exp(-h), exp(-h))
```

**Why this approach.**

- **The problem:** with `type = "first"`, ipw stops following a person
  after their first exposure event, so later prescriptions or admissions
  are ignored.
- **What it does:** every row counts. For each interval the model asks:
  given this person’s covariates, how likely was an event (or no event)
  in this interval? These chances are multiplied over time to give the
  weight.
- **Why the Andersen–Gill model:** it is the standard Cox model for
  recurrent events and is already in the survival package, so ipw adds
  only the step from hazard to probability.
- **Why each row’s own interval:** in start-stop data cut at events,
  rows have different lengths for different people. Using Λ₀(*t_(ik)*) −
  Λ₀(*s_(ik)*) gives each row the hazard over exactly the time it
  covers, instead of assuming a shared grid of visits.
- **Simulation check:** the Cox coefficients were unbiased, the
  estimated weights were within 7–9% of the true weights, and the effect
  of the number of events was recovered with 0.8% bias, against 21% for
  regression adjustment and 79% unweighted ([simulation
  study](https://tranhung93.github.io/ipw/articles/simulations.md),
  scenarios 6a–6c).

**Caveat.** For an interval that ends in an event, the method uses the
probability of an event at some time during the interval, which is an
approximation (Graffeo et al. 2018). The model assumes events are
independent given the covariates, so the number of earlier events should
be in the formula when it affects the hazard.

## 7. `ipwprobtm()`: time-varying weights from user-supplied probabilities

[`ipwprobtm()`](https://tranhung93.github.io/ipw/reference/ipwprobtm.md)
separates the weight arithmetic from the modelling: the user brings *f*
from any model, and ipw multiplies the ratios over time within each
subject and truncates the weights. It is the last step of
[`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) without
the model fitting.

**Formula.** With `den` = *f*(*A_(ik)* \| ·) and `num` = *f*(*A_(ik)* \|
·) from the numerator model, per row (`num` = 1 when not given), the
weight is the cumulative product within each subject, ordered by time,
as in the notation section:

``` math
W_i(t) = \prod_{k \le t} \frac{num_{ik}}{den_{ik}}
```

**Code.**

``` r

ratio <- num / den
ord <- order(id, timevar)
cum <- unsplit(lapply(split(ratio[ord], id[ord]), cumprod), id[ord])
weights <- numeric(n)
weights[ord] <- cum        # back to the user's row order
```

It checks that `den` and `num` are positive and finite, that lengths
match, and that each id–time pair occurs once; `id` and `timevar` are
required.

**Why this approach.**

- **The problem:** ipw only fits regression models. Many analysts now
  estimate the chance of treatment with machine learning (random
  forests, boosting, SuperLearner), which can pick up complex patterns a
  regression misses. Until now they had to compute the weights by hand.
- **What it does:** you give, for each row, the probability of the
  treatment the person actually received. ipw does the rest: divides
  numerator by denominator, multiplies over time within each person, and
  truncates if asked.
- **Why this way:** building every learner into ipw would add many
  dependencies and break whenever those packages change. Accepting
  probabilities works with any method, now and in future, and keeps ipw
  small.
- **Why time-varying only:** for a point treatment the weight is simply
  `num / den`, one line of R, so a function adds nothing. For
  time-varying data the ratios must be multiplied within each person, in
  time order, and returned in the original row order; that is where
  mistakes happen.
- **Why it asks for `timevar`:** time-varying weights must be multiplied
  in time order. Asking for the time, instead of trusting the row order,
  avoids wrong weights when the data are not sorted.
- **When to use:** when you fitted the treatment model yourself, or want
  to compare several models on the same data.

**Caveats.** Flexible models can give probabilities near 0 or 1 and so
extreme weights; truncation (`trunc`) and cross-fitting help. The user
is responsible for setting probabilities to 1 after the first switch for
`type = "first"` or `"cens"`. Standard errors that ignore the estimation
of the probabilities are only approximate.

## 8. `ipwvcov()`: standard errors that account for the estimated weights

[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
(new, experimental) gives the variance of the coefficients of a marginal
structural model fitted with [`lm()`](https://rdrr.io/r/stats/lm.html),
[`glm()`](https://rdrr.io/r/stats/glm.html) or
[`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) and weights
from `ipwpoint(family = "binomial")`, with the estimation of the weights
taken into account (M-estimation: Lunceford & Davidian 2004; Stefanski &
Boos 2002).

**Formula.** Stack the parameters $`\theta = (\beta, \alpha, \gamma)`$
of the MSM, the denominator model and the (optional) numerator model.
For person *i*, with MSM design row $`X_i`$, mean $`\mu_i`$, variance
function $`V`$, and exposure-model design rows $`Z_i`$ (denominator,
fitted $`p_i`$) and $`Z^{num}_i`$ (numerator, fitted $`q_i`$):

``` math
\psi_i(\theta) = \begin{pmatrix}
U_i \\ S_i \\ R_i
\end{pmatrix} = \begin{pmatrix}
w_i X_i (y_i - \mu_i)\, \frac{d\mu_i/d\eta_i}{V(\mu_i)} \\
Z_i (a_i - p_i)\, \frac{dp_i/d\eta_i}{p_i(1-p_i)} \\
Z^{num}_i (a_i - q_i)\, \frac{dq_i/d\eta_i}{q_i(1-q_i)}
\end{pmatrix}
```

Because $`U_i`$ is proportional to $`w_i`$,
$`\partial U_i / \partial \alpha =
U_i d_i^T`$ with
$`d_i = \partial \log w_i / \partial \alpha = g(a_i, p_i)\,
(dp_i/d\eta_i)\, Z_i`$, where

| estimand | $`g(1, p)`$     | $`g(0, p)`$    |
|----------|-----------------|----------------|
| ATE      | $`-1/p`$        | $`1/(1-p)`$    |
| ATT      | $`0`$           | $`1/(p(1-p))`$ |
| ATC      | $`-1/(p(1-p))`$ | $`0`$          |
| ATO      | $`-1/(1-p)`$    | $`1/p`$        |

and, for stabilized ATE weights,
$`e_i = \partial \log w_i / \partial \gamma
= (a_i/q_i - (1-a_i)/(1-q_i))\, (dq_i/d\eta_i)\, Z^{num}_i`$. The bread
$`A = -\sum_i \partial \psi_i / \partial \theta`$ is block
upper-triangular:

``` math
A = \begin{pmatrix}
\sum_i w_i \frac{(d\mu_i/d\eta_i)^2}{V(\mu_i)} X_i X_i^T & -\sum_i U_i d_i^T & -\sum_i U_i e_i^T \\
0 & \sum_i \frac{(dp_i/d\eta_i)^2}{p_i(1-p_i)} Z_i Z_i^T & 0 \\
0 & 0 & \sum_i \frac{(dq_i/d\eta_i)^2}{q_i(1-q_i)} Z^{num}_i Z^{num\,T}_i
\end{pmatrix}
```

the meat is $`B = \sum_i \psi_i \psi_i^T`$, and the variance is
$`A^{-1} B A^{-T}`$;
[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
returns its $`\beta`$ block.

**Code.** Each block is a cross-product computed in base R; for example,
for the denominator model and the $`\beta`$-$`\alpha`$ block:

``` r

g <- switch(estimand,
  ATE = ifelse(a == 1, -1 / p, 1 / (1 - p)),
  ATT = ifelse(a == 1, 0, 1 / (p * (1 - p))),
  ATC = ifelse(a == 1, -1 / (p * (1 - p)), 0),
  ATO = ifelse(a == 1, -1 / (1 - p), 1 / p))
D <- den$Z * (g * den$dp)          # rows d_i
A[beta, alpha] <- -crossprod(U, D)
V <- solve(A) %*% crossprod(Psi) %*% t(solve(A))
```

**Why this approach.** Robust standard errors of the weighted MSM treat
the weights as known. Estimating the weights usually *reduces* the
variance of an ATE estimate, so these standard errors are conservative
for the ATE, but for the ATT, ATC and ATO they can be too small. The
stacked estimating equations account for this exactly (asymptotically)
and need no resampling. The bread uses expected (Fisher) information.
For canonical links (logit, gaussian identity) this equals the observed
information, so the sandwich is exact; for other links the two are only
asymptotically equivalent, and in the tests a probit or cloglog exposure
model gave standard errors up to 2% away from the exact
(numerical-Jacobian) sandwich.
[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
therefore requires logistic exposure models. For the marginal structural
model any [`glm()`](https://rdrr.io/r/stats/glm.html) link is accepted.
For a saturated model such as `y ~ a` with a binary exposure, the
expected and observed information are exactly equal in any sample,
whatever the link: the fitted mean in each exposure group is that
group’s weighted mean of `y`, so the extra term of the observed
information, $`\sum_i w_i X_i X_i^T (y_i - \mu_i)\, c(\eta_i)`$ with
$`c`$ the derivative of $`(d\mu/d\eta)/V(\mu)`$ with respect to
$`\eta`$, sums to zero within each group, where $`\eta_i`$ and hence
$`c(\eta_i)`$ are constant. For other MSMs with a non-canonical link
(for example a log link with covariates) the two are asymptotically
equivalent when the model is correctly specified. Dispersion parameters
cancel. No small-sample correction is applied. The package tests compare
the result with a numerical Jacobian of the stacked equations and with
`glm_weightit(vcov = "asympt")` from the WeightIt package.

**Cox marginal structural models.** For
[`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) with
right-censored data, one row per person, the estimating equation for
$`\beta`$ is the weighted partial likelihood score (Binder 1992; Lin &
Wei 1989), with the Breslow form for tied event times. With $`\delta_i`$
the event indicator, $`Y_j(t)`$ the at-risk indicator,
$`S_0(t) = \sum_j w_j Y_j(t)
e^{\beta^T X_j}`$, $`S_1(t) = \sum_j w_j Y_j(t) e^{\beta^T X_j} X_j`$
and $`\bar X(t) = S_1(t) / S_0(t)`$, person $`i`$’s contribution is

``` math
U_i = w_i \Big[ \delta_i \big(X_i - \bar X(t_i)\big) - \sum_k \frac{\delta_k w_k
Y_i(t_k)\, e^{\beta^T X_i} \big(X_i - \bar X(t_k)\big)}{S_0(t_k)} \Big]
= w_i\, r_i ,
```

where $`r_i`$, the score residual, is what
`residuals(msm, type = "score")` returns (without the weight). The bread
for $`\beta`$ is the weighted information
$`\sum_k \delta_k w_k \big[S_2(t_k)/S_0(t_k) - \bar X(t_k) \bar
X(t_k)^T\big]`$, the inverse of `msm$naive.var`.

The weights appear in each person’s own term *and* in the risk-set sums
$`S_0`$, $`S_1`$, so it is not obvious that the cross block keeps the
glm form. It does. The total score is
$`\sum_k \delta_k w_k (X_k - \bar X(t_k))`$, and its derivative with
respect to one weight $`w_j`$ is

``` math
\frac{\partial}{\partial w_j} \sum_i U_i = \delta_j \big(X_j - \bar X(t_j)\big) -
\sum_k \delta_k w_k \frac{Y_j(t_k)\, e^{\beta^T X_j} \big(X_j - \bar X(t_k)\big)}{S_0(t_k)} = r_j ,
```

because $`\partial \bar X(t_k) / \partial w_j = Y_j(t_k) e^{\beta^T X_j}
(X_j - \bar X(t_k)) / S_0(t_k)`$. Hence $`\partial \sum_i U_i / \partial
\alpha = \sum_j r_j\, \partial w_j / \partial \alpha^T = \sum_j U_j d_j^T`$,
the same expression as for [`glm()`](https://rdrr.io/r/stats/glm.html),
and [`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
computes it analytically. The tests compare the whole variance with
central finite differences of the stacked equations, with risk sets
computed by brute force; they agree to about $`10^{-7}`$. Without the
exposure-model blocks the sandwich is
$`I^{-1} \sum_i U_i U_i^T I^{-1}`$, which is exactly
[`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)’s robust
variance (`msm$var`): the standard error that treats the weights as
known.

With tied event times, Efron’s approximation changes the score and the
information in a way that ipw does not reproduce, so
[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
requires `ties = "breslow"` when event times are tied. Without tied
event times all tie methods give the same partial likelihood, and any
`ties` is accepted.

`WeightIt::coxph_weightit(vcov = "asympt")` gives slightly different
values. Its `vcov = "HC0"` equals
[`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)’s robust
variance exactly, and the weights agree to $`10^{-14}`$, so the
difference lies in how the dependence of the Cox score on the
propensity-model coefficients enters. In our checks the relative
difference was below 0.1% for ATE and overlap weights and shrank with
the sample size; for ATT and ATC weights WeightIt’s standard error was 1
to 4% larger and the difference did not shrink (n from 600 to 38 400).
[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
matches the numerical Jacobian of the exact stacked equations in all
cases.

**Time-varying weights from
[`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md).** For
`family = "binomial"` with logistic models, each model (denominator,
numerator; treatment, censoring) is a parameter block, and the person,
not the row, is the unit of independence. For row $`r`$ of person $`p`$,
[`ipwtm()`](https://tranhung93.github.io/ipw/reference/ipwtm.md) defines
the weight as a product over that person’s rows $`s \le r`$ in time
order,

``` math
w_r = \prod_{s \le r} \frac{f^{num}_s}{f^{den}_s}, \qquad
f_s = \begin{cases}
P(A_s = a_s) & \text{type "first" or "all", rows used by the model}\\
1 - p_s & \text{type "cens", rows used by the model (the censoring row included)}\\
1 & \text{rows after the first switch ("first", "cens")},
\end{cases}
```

where the rows used by the model are those with `selvar = 1` (all rows
for `type = "all"`). With logistic models,
$`\partial \log f_s / \partial \alpha =
(a_s - p_s) Z_s`$ for “first” and “all”, and $`-p_s Z_s`$ for “cens”, so

``` math
d_r = \frac{\partial \log w_r}{\partial \alpha} = -\sum_{s \le r,\ s \in p} \frac{\partial \log f^{den}_s}{\partial \alpha},
\qquad
e_r = \frac{\partial \log w_r}{\partial \gamma} = \sum_{s \le r,\ s \in p} \frac{\partial \log f^{num}_s}{\partial \gamma},
```

cumulative sums within person, which also carry the value forward after
the switch. The exposure-model score of person $`p`$ is
$`S_p = \sum_s (a_s - p_s)
Z_s`$ over the rows used by the model, and its information comes from
the fitted glm. For the MSM, $`U_r`$ is the row score (the glm score, or
for a Cox model on counting-process data `Surv(tstart, time, event)` the
weighted score residual with risk sets
$`\{j: t^{start}_j < t \le t_j\}`$), the cross blocks are
$`-\sum_r U_r d_r^T`$ and $`-\sum_r U_r e_r^T`$ (for the Cox model by
the same identity as above, now row by row), and the meat is
$`B = \sum_p \psi_p \psi_p^T`$ with $`\psi_p`$ the sum of all scores of
person $`p`$. A list of weight objects whose weights multiply, such as
treatment times censoring weights, adds one block per model.

The package tests compare the result with central finite differences of
the stacked equations, with the weights recomputed from the model
coefficients by independent code that follows the definition above
(agreement about $`10^{-7}`$, on shuffled rows), for a
[`glm()`](https://rdrr.io/r/stats/glm.html) MSM with “all” and “cens”
weights, a Cox MSM with “first” and “cens” weights and a subset of
`haartdat`. They also compare with
`WeightIt::glm_weightit(vcov = "asympt")` for a `weightitMSM()` object:
with visits at the same times for everyone and time-specific columns in
the pooled model, the weights are the same, and an MSM fitted on the
long data with one parameter per non-last row (so that only the last
rows carry information) has the same coefficients and variance as
WeightIt’s person-level model (agreement $`10^{-6}`$).

**Heavy-tailed weights.** Cumulative weights are products of many
factors. When some people have a very low probability of their observed
exposure history (a near-violation of positivity), the weights have a
heavy right tail: P(W \> x) falls off like $`x^{-\gamma}`$ for a tail
index $`\gamma`$. The weighted estimating equations then have a finite
variance only if $`\gamma > 2`$. The usual normal approximation also
relies on a finite third moment ($`\gamma > 3`$), and estimating the
variance well, from the sandwich or from the bootstrap, relies on a
finite fourth moment ($`\gamma > 4`$).

Simulations with individual visit times (see the comparison article)
estimated $`\gamma`$ with the Hill estimator on the top 5% of the final
cumulative weights. The more strongly treatment depended on the
confounder and on the previous treatment, the smaller the tail index:
about 3.9, 2.8 and 2.3 in three scenarios. Coverage of 95% intervals
from
[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md)
followed the tail index: about 95%, about 94% and 92 to 93%. In the
heaviest-tailed scenario coverage did not improve between 300 and 4800
people, while the largest weight kept growing with *n*. The cluster
bootstrap did no better (90 to 95%, on fewer data sets), and
[`ipwvcov()`](https://tranhung93.github.io/ipw/reference/ipwvcov.md) and
the bootstrap gave about the same standard errors throughout. The
formula is exact for the stacked estimating equations; the problem is
the weights, and it affects any normal-theory interval. Check the weight
distribution ([`summary()`](https://rdrr.io/r/base/summary.html),
[`ipwplot()`](https://tranhung93.github.io/ipw/reference/ipwplot.md),
the effective sample size) before trusting any standard error. If a few
people carry much of the weight, consider truncation, a less extreme
weight model or a different estimand.

**Caveats.** Only binary exposures with logistic
[`glm()`](https://rdrr.io/r/stats/glm.html) exposure models; truncated
weights (truncation at estimated percentiles is not a smooth function of
$`\alpha`$), smooth terms, other families (including Cox exposure models
in `ipwtm(family = "survival")`), Cox MSMs with strata or frailties, and
Efron ties with tied event times are not covered, and the bootstrap
remains the advice there.

## References

- Binder DA (1992). Fitting Cox’s proportional hazards models from
  survey data. *Biometrika* 79, 139–147.
- Cole SR, Hernán MA (2008). Constructing inverse probability weights
  for marginal structural models. *American Journal of Epidemiology*
  168, 656–664.
- Graffeo N, Latouche A, Geskus RB, Chevret S (2018). Modeling
  time-varying exposure using inverse probability of treatment weights.
  *Biometrical Journal* 60, 323–332.
- Hirano K, Imbens GW (2004). The propensity score with continuous
  treatments. In Gelman A, Meng XL (eds), *Applied Bayesian Modeling and
  Causal Inference from Incomplete-Data Perspectives*. Wiley.
- Imai K, van Dyk DA (2004). Causal inference with general treatment
  regimes: generalizing the propensity score. *Journal of the American
  Statistical Association* 99, 854–866.
- Li F, Li F (2019). Propensity score weighting for causal inference
  with multiple treatments. *Annals of Applied Statistics* 13,
  2389–2415.
- Li F, Morgan KL, Zaslavsky AM (2018). Balancing covariates via
  propensity score weighting. *Journal of the American Statistical
  Association* 113, 390–400.
- Lin DY, Wei LJ (1989). The robust inference for the Cox proportional
  hazards model. *Journal of the American Statistical Association* 84,
  1074–1078.
- Lunceford JK, Davidian M (2004). Stratification and weighting via the
  propensity score in estimation of causal treatment effects: a
  comparative study. *Statistics in Medicine* 23, 2937–2960.
- Robins JM, Hernán MA, Brumback B (2000). Marginal structural models
  and causal inference in epidemiology. *Epidemiology* 11, 550–560.
- Stefanski LA, Boos DD (2002). The calculus of M-estimation. *The
  American Statistician* 56, 29–38.
- van der Laan MJ, Polley EC, Hubbard AE (2007). Super learner.
  *Statistical Applications in Genetics and Molecular Biology* 6(1).
- van der Wal WM, Geskus RB (2011). ipw: An R package for inverse
  probability weighting. *Journal of Statistical Software* 43(13), 1–23.
- Venables WN, Ripley BD (2002). *Modern Applied Statistics with S*, 4th
  ed. Springer.
- Wood SN (2017). *Generalized Additive Models: An Introduction with R*,
  2nd ed. CRC Press.
