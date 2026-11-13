# Tests for ipwdiag(): weight diagnostics

sim_point <- function(n = 1000, seed = 16) {
  set.seed(seed)
  d <- data.frame(l = rnorm(n, 10, 5))
  d$a <- rbinom(n, 1, plogis(d$l - 10))
  d$y <- 10 * d$a + 0.5 * d$l + rnorm(n, -10, 5)
  d
}

# Evaluate expr, muffling its messages; return the value with the messages as attribute
with_messages <- function(expr) {
  msgs <- character()
  val <- withCallingHandlers(expr, message = function(m) {
    msgs <<- c(msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  structure(list(value = val, messages = msgs))
}

# A minimal ipwpoint-like object holding given weights
fake_point <- function(w) structure(list(ipw.weights = w), class = c("ipwpoint", "list"))

test_that("ipwdiag returns the documented structure for ipwpoint (binomial)", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  simdat <- sim_point()
  w <- ipwpoint(exposure = a, family = "binomial", link = "logit",
                numerator = ~ 1, denominator = ~ l, data = simdat)
  msm <- lm(y ~ a, data = simdat, weights = w$ipw.weights)
  d <- ipwdiag(w, msm = msm)
  expect_s3_class(d, "ipwdiag")
  expect_named(d, c("overlap", "share", "tail", "trunc", "summary"))
  expect_named(d$overlap, c("data", "density", "support", "exposure", "xlim", "extreme"))
  expect_equal(nrow(d$overlap$data), 1000)
  expect_equal(d$overlap$data$p, unname(fitted(w$den.mod)))
  expect_equal(d$overlap$exposure, "a")
  expect_named(d$share, c("data", "weights"))
  expect_named(d$tail, c("data", "alpha", "k", "threshold", "intercept", "slope"))
  expect_named(d$trunc, c("q", "estimate", "se", "lower", "upper"))
  expect_equal(d$trunc$q, c(0, 0.01, 0.025, 0.05, 0.1))
  expect_equal(attr(d$trunc, "coef"), "a")
  expect_named(d$summary, c("n", "unit", "overlap", "share", "tail", "trunc"))
  expect_named(d$summary$share, c("top1", "top5", "max", "mean", "ess", "ess.fraction"))
  p <- fitted(w$den.mod)
  expect_equal(d$summary$overlap[["proportion"]], mean(p < 0.05 | p > 0.95))
  # q = 0 reproduces the MSM
  expect_equal(d$trunc$estimate[1], unname(coef(msm)["a"]))
  expect_output(print(d), "Kish ESS")
  expect_output(print(d), "Truncation sensitivity for 'a'")
  expect_output(print(d), "Warning: tail index < 2")
  # invisible result
  expect_invisible(ipwdiag(w, which = "share"))
})

test_that("ipwdiag skips the overlap plot for a gaussian ipwpoint exposure", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  set.seed(3)
  d <- data.frame(l = rnorm(500))
  d$a <- d$l + rnorm(500)
  w <- ipwpoint(exposure = a, family = "gaussian", numerator = ~ 1, denominator = ~ l,
                data = d)
  out <- with_messages(ipwdiag(w))
  res <- out$value
  msgs <- out$messages
  expect_length(msgs, 2)
  expect_match(msgs[1], "overlap plot needs a binomial")
  expect_match(msgs[2], "no 'msm' given")
  expect_s3_class(res, "ipwdiag")
  expect_named(res, c("share", "tail", "summary"))
  expect_null(res$summary$overlap)
  expect_error(ipwdiag(w, which = "overlap"), "binomial exposure model")
})

test_that("ipwdiag works for ipwtm() binomial and survival weights on haartdat", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  data(haartdat, package = "ipw", envir = environment())
  wb <- ipwtm(exposure = haartind, family = "binomial", link = "logit",
              numerator = ~ sex + age, denominator = ~ sex + age + cd4.sqrt,
              id = patient, timevar = fuptime, type = "first", data = haartdat)
  db <- ipwdiag(wb, which = c("overlap", "share", "tail"))
  expect_s3_class(db, "ipwdiag")
  expect_named(db, c("overlap", "share", "tail", "summary"))
  # overlap uses the rows of the denominator model (selvar == 1)
  expect_equal(nrow(db$overlap$data), sum(wb$selvar == 1))
  expect_equal(db$overlap$data$exposure, haartdat$haartind[wb$selvar == 1])
  # share uses each person's final cumulative weight
  n_id <- length(unique(haartdat$patient))
  expect_equal(db$summary$n, n_id)
  last <- !duplicated(haartdat$patient, fromLast = TRUE)   # haartdat is sorted
  expect_equal(sort(db$share$weights), sort(wb$ipw.weights[last]))

  ws <- ipwtm(exposure = haartind, family = "survival",
              numerator = ~ sex + age, denominator = ~ sex + age + cd4.sqrt,
              id = patient, tstart = tstart, timevar = fuptime, type = "first",
              data = haartdat)
  out <- with_messages(ipwdiag(ws))
  ds <- out$value
  expect_match(out$messages[1], "overlap plot needs a binomial")
  expect_named(ds, c("share", "tail", "summary"))
  expect_equal(ds$summary$n, n_id)
  expect_equal(ds$summary$share[["max"]], max(ws$ipw.weights[last]))

  # final weights are found from the id and order attributes, whatever the row order
  shuf <- haartdat[sample(nrow(haartdat)), ]
  wb2 <- ipwtm(exposure = haartind, family = "binomial", link = "logit",
               numerator = ~ sex + age, denominator = ~ sex + age + cd4.sqrt,
               id = patient, timevar = fuptime, type = "first", data = shuf)
  db2 <- ipwdiag(wb2, which = "share")
  expect_equal(db2$share$weights, db$share$weights)
})

test_that("a list of weight objects uses the product of the weights", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  data(haartdat, package = "ipw", envir = environment())
  args <- list(family = "binomial", link = "logit", numerator = ~ sex + age,
               denominator = ~ sex + age + cd4.sqrt, id = quote(patient),
               timevar = quote(fuptime), data = haartdat)
  wt <- do.call(ipwtm, c(list(exposure = quote(haartind), type = "first"), args))
  wc <- do.call(ipwtm, c(list(exposure = quote(dropout), type = "cens"), args))
  d <- ipwdiag(list(wt, wc), which = c("overlap", "share"))
  last <- !duplicated(haartdat$patient, fromLast = TRUE)
  expect_equal(d$share$weights,
               sort((wt$ipw.weights * wc$ipw.weights)[last], decreasing = TRUE))
  # overlap from the first element
  expect_equal(d$overlap$exposure, "haartind")
  expect_error(ipwdiag(list(wt, fake_point(1:3))), "same number of weights")
})

test_that("share and Kish ESS match hand calculations", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  set.seed(1)
  w <- rexp(250) + 0.1
  d <- ipwdiag(fake_point(w), which = "share")
  ws <- sort(w, decreasing = TRUE)
  expect_equal(d$summary$share[["top1"]], sum(ws[1:3]) / sum(w))    # ceiling(2.5) = 3
  expect_equal(d$summary$share[["top5"]], sum(ws[1:13]) / sum(w))   # ceiling(12.5) = 13
  expect_equal(d$summary$share[["max"]], max(w))
  expect_equal(d$summary$share[["mean"]], mean(w))
  expect_equal(d$summary$share[["ess"]], sum(w)^2 / sum(w^2))
  expect_equal(d$summary$share[["ess.fraction"]], sum(w)^2 / sum(w^2) / 250)
  # the curve starts at (0, 0) and ends at (1, 1)
  expect_equal(unlist(d$share$data[1, ]), c(fraction = 0, share = 0))
  expect_equal(unlist(d$share$data[251, ]), c(fraction = 1, share = 1))
  # equal weights: on the diagonal, ESS = n
  de <- ipwdiag(fake_point(rep(2, 200)), which = "share")
  expect_equal(de$share$data$share, de$share$data$fraction)
  expect_equal(de$summary$share[["ess"]], 200)
})

test_that("the Hill estimate matches a hand calculation and recovers a Pareto tail", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  set.seed(2)
  w <- exp(rnorm(1000))
  d <- ipwdiag(fake_point(w), which = "tail")
  ws <- sort(w, decreasing = TRUE)
  k <- 50
  expect_equal(d$tail$k, k)
  expect_equal(d$tail$threshold, ws[k + 1])
  expect_equal(d$tail$alpha, 1 / mean(log(ws[1:k]) - log(ws[k + 1])))
  expect_equal(d$tail$data$log.weight, log(ws[1:k]))
  expect_equal(d$tail$data$log.rank, log((1:k) / 1000))
  expect_equal(d$tail$slope, -d$tail$alpha)
  # a different top fraction
  d10 <- ipwdiag(fake_point(w), which = "tail", top = 0.1)
  expect_equal(d10$tail$alpha, 1 / mean(log(ws[1:100]) - log(ws[101])))

  # Pareto weights with alpha = 3: P(W > x) = x^-3 for x >= 1
  set.seed(4)
  wp <- runif(100000)^(-1 / 3)
  dp <- ipwdiag(fake_point(wp), which = "tail")
  expect_equal(dp$tail$alpha, 3, tolerance = 0.05)   # SE about 3 / sqrt(5000) = 0.04
  expect_output(print(dp), "tail index < 4")
  # no warning for a light tail
  set.seed(5)
  dl <- ipwdiag(fake_point(runif(2000, 0.5, 1.5)), which = "tail")
  expect_true(dl$tail$alpha > 4)
  expect_false(any(grepl("Warning", capture.output(print(dl)))))
})

test_that("trunc: estimates equal refits with ipwpoint(trunc = q); SE equals svyglm's", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  simdat <- sim_point()
  w <- ipwpoint(exposure = a, family = "binomial", link = "logit",
                numerator = ~ 1, denominator = ~ l, data = simdat)
  msm <- glm(y ~ a + l, data = simdat, weights = w$ipw.weights)
  qs <- c(0, 0.02, 0.05)
  d <- ipwdiag(w, which = "trunc", msm = msm, trunc = qs)
  for (j in seq_along(qs)[-1]) {
    wq <- ipwpoint(exposure = a, family = "binomial", link = "logit", numerator = ~ 1,
                   denominator = ~ l, data = simdat, trunc = qs[j])$weights.trunc
    fit <- glm(y ~ a + l, data = simdat, weights = wq)
    expect_equal(d$trunc$estimate[j], unname(coef(fit)["a"]))
    if (requireNamespace("survey", quietly = TRUE)) {
      sd2 <- simdat
      sd2$wq <- wq
      svy <- survey::svyglm(y ~ a + l, design = survey::svydesign(~ 1, weights = ~ wq,
                                                                  data = sd2))
      expect_equal(d$trunc$se[j], unname(sqrt(diag(vcov(svy)))["a"]))
    }
  }
  expect_equal(d$trunc$lower, d$trunc$estimate - qnorm(0.975) * d$trunc$se)
  expect_equal(d$trunc$upper, d$trunc$estimate + qnorm(0.975) * d$trunc$se)
  # coef by name and position
  expect_equal(attr(ipwdiag(w, "trunc", msm = msm, coef = "l")$trunc, "coef"), "l")
  expect_equal(attr(ipwdiag(w, "trunc", msm = msm, coef = 1)$trunc, "coef"), "(Intercept)")
})

test_that("trunc: robust SE for a logistic glm equals svyglm's", {
  skip_if_not_installed("survey")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  simdat <- sim_point()
  simdat$yb <- as.numeric(simdat$y > 0)
  w <- ipwpoint(exposure = a, family = "binomial", link = "logit",
                numerator = ~ 1, denominator = ~ l, data = simdat)
  msm <- suppressWarnings(glm(yb ~ a, family = binomial, data = simdat,
                              weights = w$ipw.weights))
  d <- suppressWarnings(ipwdiag(w, "trunc", msm = msm, trunc = c(0, 0.05)))
  simdat$wq <- ipwpoint(exposure = a, family = "binomial", link = "logit", numerator = ~ 1,
                        denominator = ~ l, data = simdat, trunc = 0.05)$weights.trunc
  svy <- suppressWarnings(survey::svyglm(yb ~ a, family = quasibinomial,
                                         design = survey::svydesign(~ 1, weights = ~ wq,
                                                                    data = simdat)))
  expect_equal(d$trunc$estimate[2], unname(coef(svy)["a"]), tolerance = 1e-6)
  expect_equal(d$trunc$se[2], unname(sqrt(diag(vcov(svy)))["a"]), tolerance = 1e-6)
})

test_that("trunc: Cox MSM SE equals coxph(robust = TRUE)", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  simdat <- sim_point()
  set.seed(5)
  simdat$time <- rexp(nrow(simdat), 0.1 * exp(-0.5 * simdat$a + 0.05 * simdat$l))
  simdat$status <- as.numeric(simdat$time < 15)
  simdat$time <- pmin(simdat$time, 15)
  w <- ipwpoint(exposure = a, family = "binomial", link = "logit",
                numerator = ~ 1, denominator = ~ l, data = simdat)
  msm <- survival::coxph(survival::Surv(time, status) ~ a, data = simdat,
                         weights = w$ipw.weights)
  d <- ipwdiag(w, "trunc", msm = msm, trunc = c(0, 0.05))
  wq <- ipwpoint(exposure = a, family = "binomial", link = "logit", numerator = ~ 1,
                 denominator = ~ l, data = simdat, trunc = 0.05)$weights.trunc
  fit <- survival::coxph(survival::Surv(time, status) ~ a, data = simdat, weights = wq,
                         robust = TRUE)
  expect_equal(d$trunc$estimate[2], unname(coef(fit)))
  expect_equal(d$trunc$se[2], unname(sqrt(diag(fit$var))))
  expect_equal(attr(d$trunc, "coef"), "a")
})

test_that("trunc with ipwtm(): equals ipwtm(trunc = q); clustered by person", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  data(haartdat, package = "ipw", envir = environment())
  ws <- ipwtm(exposure = haartind, family = "survival",
              numerator = ~ sex + age, denominator = ~ sex + age + cd4.sqrt,
              id = patient, tstart = tstart, timevar = fuptime, type = "first",
              data = haartdat)
  wq <- ipwtm(exposure = haartind, family = "survival",
              numerator = ~ sex + age, denominator = ~ sex + age + cd4.sqrt,
              id = patient, tstart = tstart, timevar = fuptime, type = "first",
              data = haartdat, trunc = 0.05)$weights.trunc
  # Cox MSM without cluster(): the ipwtm() id is used
  msm <- survival::coxph(survival::Surv(tstart, fuptime, event) ~ haartind,
                         data = haartdat, weights = ws$ipw.weights)
  d <- ipwdiag(ws, "trunc", msm = msm, trunc = c(0, 0.05))
  fit <- survival::coxph(survival::Surv(tstart, fuptime, event) ~ haartind +
                           cluster(patient), data = haartdat, weights = wq)
  expect_equal(d$trunc$estimate[2], unname(coef(fit)))
  expect_equal(d$trunc$se[2], unname(sqrt(diag(fit$var))))
  # with cluster() in the model
  msm2 <- survival::coxph(survival::Surv(tstart, fuptime, event) ~ haartind +
                            cluster(patient), data = haartdat, weights = ws$ipw.weights)
  d2 <- ipwdiag(ws, "trunc", msm = msm2, trunc = c(0, 0.05))
  expect_equal(d2$trunc, d$trunc)

  # lm MSM on the long data: scores summed within person, as svyglm with ids = ~ patient
  skip_if_not_installed("survey")
  msm3 <- lm(event ~ haartind, data = haartdat, weights = ws$ipw.weights)
  d3 <- ipwdiag(ws, "trunc", msm = msm3, trunc = c(0, 0.05))
  hd <- haartdat
  hd$wq <- wq
  svy <- survey::svyglm(event ~ haartind,
                        design = survey::svydesign(ids = ~ patient, weights = ~ wq, data = hd))
  expect_equal(d3$trunc$estimate[2], unname(coef(svy)["haartind"]))
  expect_equal(d3$trunc$se[2], unname(sqrt(diag(vcov(svy)))["haartind"]))
})

test_that("trunc refits msm where it was created when the caller lacks the data", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  fit_here <- function() {
    dd <- sim_point(300, seed = 8)
    w <- ipwpoint(exposure = a, family = "binomial", link = "logit",
                  numerator = ~ 1, denominator = ~ l, data = dd)
    list(w = w, msm = lm(y ~ a, data = dd, weights = w$ipw.weights))
  }
  f <- fit_here()
  d <- ipwdiag(f$w, "trunc", msm = f$msm, trunc = c(0, 0.05))
  expect_equal(d$trunc$estimate[1], unname(coef(f$msm)["a"]))

  # data nowhere to be found: clear error
  dd2 <- sim_point(300, seed = 8)
  w2 <- ipwpoint(exposure = a, family = "binomial", link = "logit",
                 numerator = ~ 1, denominator = ~ l, data = dd2)
  msm2 <- lm(y ~ a, data = dd2, weights = w2$ipw.weights)
  rm(dd2)
  expect_error(ipwdiag(w2, "trunc", msm = msm2), "Could not refit 'msm'")
})

test_that("ipwdiag gives clear errors and messages for wrong input", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  simdat <- sim_point(300)
  w <- ipwpoint(exposure = a, family = "binomial", link = "logit",
                numerator = ~ 1, denominator = ~ l, data = simdat)
  expect_error(ipwdiag(simdat), "must be the result of ipwpoint\\(\\) or ipwtm\\(\\)")
  expect_error(ipwdiag(list()), "must be the result")
  expect_error(ipwdiag(fake_point(NULL)), "no valid weights")
  expect_error(ipwdiag(fake_point(c(1, NA))), "no valid weights")
  expect_error(ipwdiag(fake_point(c(1, 0, 2))), "must be positive")
  expect_error(ipwdiag(w, which = "boxplot"), "should be one of")
  expect_error(ipwdiag(w, top = 0.5), "'top' must be")
  expect_error(ipwdiag(w, top = 0), "'top' must be")
  expect_error(ipwdiag(w, trunc = 0.5), "'trunc' must be")
  expect_error(ipwdiag(w, trunc = -0.1), "'trunc' must be")
  expect_error(ipwdiag(w, trunc = NA_real_), "'trunc' must be")
  expect_error(ipwdiag(w, ask = NA), "'ask' must be")
  expect_error(ipwdiag(w, which = "trunc"), "needs 'msm'")
  msm <- lm(y ~ a, data = simdat, weights = w$ipw.weights)
  expect_error(ipwdiag(w, "trunc", msm = lm(y ~ a, data = simdat[1:200, ])),
               "fitted on 200 rows but there are 300 weights")
  expect_error(ipwdiag(w, "trunc", msm = msm, coef = "b"), "'coef' must be")
  expect_error(ipwdiag(w, "trunc", msm = msm, coef = 5), "'coef' must be")
  expect_error(ipwdiag(w, "trunc", msm = lm(y ~ 1, data = simdat, weights = w$ipw.weights)),
               "no coefficient other than the intercept")
  expect_error(ipwdiag(w, "trunc", msm = structure(list(), class = "svyglm")),
               "must be a model fitted with lm\\(\\), glm\\(\\)")
  expect_warning(ipwdiag(w, "trunc", msm = lm(y ~ a, data = simdat), trunc = 0),
                 "weights of 'msm' differ")
  # too few weights in the tail: skipped with a message
  expect_message(r <- ipwdiag(fake_point(runif(100) + 1), which = "tail"),
                 "at least 20 weights")
  expect_null(r$tail)
  expect_null(r$summary$tail)
  # old ipwtm() objects without attributes
  old <- structure(list(ipw.weights = c(1, 2)), class = c("ipwtm", "list"))
  expect_error(ipwdiag(old), "re-run ipwtm\\(\\)")
})

test_that("ipwdiag draws one plot per page with ask = TRUE and restores par()", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  simdat <- sim_point(300)
  w <- ipwpoint(exposure = a, family = "binomial", link = "logit",
                numerator = ~ 1, denominator = ~ l, data = simdat)
  op <- graphics::par(c("mfrow", "ask", "cex"))
  suppressMessages(ipwdiag(w, cex = 0.7))
  expect_equal(graphics::par(c("mfrow", "ask", "cex")), op)
  suppressMessages(ipwdiag(w, ask = TRUE))
  expect_equal(graphics::par(c("mfrow", "ask", "cex")), op)
  # a layout the user set up is kept
  graphics::par(mfrow = c(1, 3))
  suppressMessages(ipwdiag(w, which = c("overlap", "share", "tail")))
  expect_equal(graphics::par("mfrow"), c(1, 3))
})
