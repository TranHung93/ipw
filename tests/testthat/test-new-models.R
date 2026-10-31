# Tests for the models added after 1.3.2: each new option is checked against
# weights computed by hand from the underlying model.

# ---- ipwtm(): type = "all" for multinomial and ordinal ------------------------

make_long_categorical <- function(seed = 3, n_id = 100, n_t = 3) {
  set.seed(seed)
  dat <- data.frame(id = rep(seq_len(n_id), each = n_t), time = rep(seq_len(n_t), n_id))
  dat$l <- rnorm(nrow(dat))
  eta2 <- -0.3 + 0.8 * dat$l
  eta3 <- -0.3 + 1.6 * dat$l
  den <- 1 + exp(eta2) + exp(eta3)
  u <- runif(nrow(dat))
  dat$a <- factor(ifelse(u < 1 / den, "x", ifelse(u < (1 + exp(eta2)) / den, "y", "z")))
  dat$o <- factor(as.integer(dat$a), ordered = TRUE)
  dat
}

observed_prob <- function(probs, obs) {
  probs[cbind(seq_along(obs), match(as.character(obs), colnames(probs)))]
}

test_that("ipwtm type = 'all' works for a multinomial exposure", {
  dat <- make_long_categorical()
  res <- ipwtm(exposure = a, family = "multinomial", numerator = ~ 1, denominator = ~ l,
               id = id, timevar = time, type = "all", data = dat)

  expect_true(all(res$selvar == 1))
  den <- nnet::multinom(a ~ l, data = dat, trace = FALSE)
  num <- nnet::multinom(a ~ 1, data = dat, trace = FALSE)
  ratio <- observed_prob(predict(num, type = "probs"), dat$a) /
    observed_prob(predict(den, type = "probs"), dat$a)
  expect_equal(unname(res$ipw.weights), unname(ave(ratio, dat$id, FUN = cumprod)),
               tolerance = 1e-5)
})

test_that("ipwtm type = 'all' works for an ordinal exposure", {
  dat <- make_long_categorical()
  res <- ipwtm(exposure = o, family = "ordinal", link = "logit", numerator = ~ 1,
               denominator = ~ l, id = id, timevar = time, type = "all", data = dat)

  expect_true(all(res$selvar == 1))
  den <- MASS::polr(o ~ l, data = dat, method = "logistic", Hess = TRUE)
  num <- MASS::polr(o ~ 1, data = dat, method = "logistic", Hess = TRUE)
  ratio <- observed_prob(predict(num, type = "probs"), dat$o) /
    observed_prob(predict(den, type = "probs"), dat$o)
  expect_equal(unname(res$ipw.weights), unname(ave(ratio, dat$id, FUN = cumprod)),
               tolerance = 1e-5)
})

test_that("ipwtm still refuses type = 'all' for survival", {
  dat <- make_long_categorical()
  dat$event <- rbinom(nrow(dat), 1, 0.2)
  dat$tstart <- dat$time - 1
  expect_error(
    ipwtm(exposure = event, family = "survival", denominator = ~ l, id = id,
          tstart = tstart, timevar = time, type = "all", data = dat),
    "not yet implemented"
  )
})

# ---- ipwprob(): user-supplied probabilities ----------------------------------

test_that("ipwprob reproduces ipwpoint weights from the same probabilities", {
  set.seed(5)
  d <- data.frame(l = rnorm(200))
  d$a <- rbinom(200, 1, plogis(d$l))
  ref <- ipwpoint(exposure = a, family = "binomial", link = "logit",
                  numerator = ~ 1, denominator = ~ l, data = d)
  pd <- fitted(ref$den.mod)
  pn <- fitted(ref$num.mod)
  res <- ipwprob(den = ifelse(d$a == 1, pd, 1 - pd), num = ifelse(d$a == 1, pn, 1 - pn))

  expect_s3_class(res, "ipwprob")
  expect_equal(unname(res$ipw.weights), unname(ref$ipw.weights))
  expect_equal(ipwprob(den = c(0.5, 0.25))$ipw.weights, c(2, 4))
})

test_that("ipwprob cumulates within subjects over time, in input order", {
  long <- data.frame(id = c(1, 1, 1, 2, 2), time = c(1, 2, 3, 1, 2),
                     den = c(0.5, 0.8, 0.4, 0.9, 0.6), num = c(0.6, 0.6, 0.5, 0.7, 0.7))
  expected <- c(0.6 / 0.5, 0.6 / 0.5 * 0.6 / 0.8, 0.6 / 0.5 * 0.6 / 0.8 * 0.5 / 0.4,
                0.7 / 0.9, 0.7 / 0.9 * 0.7 / 0.6)
  res <- ipwprob(den = long$den, num = long$num, id = long$id, timevar = long$time)
  expect_equal(res$ipw.weights, expected)

  shuffled <- c(5, 2, 4, 1, 3)
  res_s <- ipwprob(den = long$den[shuffled], num = long$num[shuffled],
                   id = long$id[shuffled], timevar = long$time[shuffled])
  expect_equal(res_s$ipw.weights, expected[shuffled])
})

test_that("ipwprob truncates and rejects invalid input", {
  w <- ipwprob(den = c(0.1, 0.2, 0.5, 0.9, 0.95), trunc = 0.2)$weights.trunc
  raw <- 1 / c(0.1, 0.2, 0.5, 0.9, 0.95)
  q <- unname(quantile(raw, c(0.2, 0.8)))
  expect_equal(w, pmin(pmax(raw, q[1]), q[2]))

  expect_error(ipwprob(den = c(0.5, 0)), "positive")
  expect_error(ipwprob(den = c(0.5, NA)), "positive")
  expect_error(ipwprob(den = c(0.5, 0.4), num = 0.5), "same length")
  expect_error(ipwprob(den = c(0.5, 0.4), id = c(1, 1)), "timevar")
  expect_error(ipwprob(den = c(0.5, 0.4), timevar = c(1, 2)), "'id' is required")
  expect_error(ipwprob(den = c(0.5, 0.4), id = c(1, 1), timevar = c(1, 1)), "only once")
})

# ---- ipwpoint(): count exposures ---------------------------------------------

test_that("ipwpoint family = 'poisson' matches dpois weights", {
  set.seed(6)
  d <- data.frame(l = rnorm(300))
  d$a <- rpois(300, exp(0.3 + 0.5 * d$l))
  res <- ipwpoint(exposure = a, family = "poisson", numerator = ~ 1, denominator = ~ l, data = d)

  mn <- glm(a ~ 1, family = poisson, data = d)
  md <- glm(a ~ l, family = poisson, data = d)
  expected <- dpois(d$a, fitted(mn)) / dpois(d$a, fitted(md))
  expect_equal(unname(res$ipw.weights), unname(expected))
  expect_equal(mean(res$ipw.weights), 1, tolerance = 0.15)
  expect_identical(res$den.mod$call$family, quote(poisson()))

  unstab <- ipwpoint(exposure = a, family = "poisson", denominator = ~ l, data = d)
  expect_equal(unname(unstab$ipw.weights), unname(1 / dpois(d$a, fitted(md))))

  d$b <- d$a + 0.5
  expect_error(ipwpoint(exposure = b, family = "poisson", denominator = ~ l, data = d),
               "non-negative integers")
})

test_that("ipwpoint family = 'negbin' matches dnbinom weights", {
  set.seed(7)
  d <- data.frame(l = rnorm(400))
  d$a <- rnbinom(400, size = 2, mu = exp(0.5 + 0.5 * d$l))
  res <- ipwpoint(exposure = a, family = "negbin", numerator = ~ 1, denominator = ~ l, data = d)

  mn <- MASS::glm.nb(a ~ 1, data = d)
  md <- MASS::glm.nb(a ~ l, data = d)
  expected <- dnbinom(d$a, size = mn$theta, mu = fitted(mn)) /
    dnbinom(d$a, size = md$theta, mu = fitted(md))
  expect_s3_class(res$den.mod, "negbin")
  expect_equal(unname(res$ipw.weights), unname(expected), tolerance = 1e-6)
})

# ---- ipwpoint(): skewed continuous exposures ---------------------------------

test_that("ipwpoint family = 'gamma' matches dgamma weights", {
  set.seed(8)
  d <- data.frame(l = rnorm(300))
  d$a <- rgamma(300, shape = 3, rate = 3 / exp(0.5 + 0.3 * d$l))
  res <- ipwpoint(exposure = a, family = "gamma", numerator = ~ 1, denominator = ~ l, data = d)

  dens <- function(m) {
    shape <- 1 / summary(m)$dispersion
    dgamma(d$a, shape = shape, rate = shape / fitted(m))
  }
  mn <- glm(a ~ 1, family = Gamma(link = "log"), data = d)
  md <- glm(a ~ l, family = Gamma(link = "log"), data = d)
  expect_equal(unname(res$ipw.weights), unname(dens(mn) / dens(md)))

  expect_error(ipwpoint(exposure = a, family = "gamma", denominator = ~ l, data = d),
               "numerator")
  d$z <- d$a - 1
  expect_error(ipwpoint(exposure = z, family = "gamma", numerator = ~ 1,
                        denominator = ~ l, data = d), "positive")
})

test_that("ipwpoint family = 'lognormal' matches dlnorm weights", {
  set.seed(9)
  d <- data.frame(l = rnorm(300))
  d$a <- exp(rnorm(300, 0.5 * d$l, 0.5))
  res <- ipwpoint(exposure = a, family = "lognormal", numerator = ~ 1, denominator = ~ l, data = d)

  dens <- function(m) dlnorm(d$a, meanlog = fitted(m), sdlog = sd(residuals(m)))
  mn <- glm(log(a) ~ 1, data = d)
  md <- glm(log(a) ~ l, data = d)
  expect_equal(unname(res$ipw.weights), unname(dens(mn) / dens(md)))
  expect_equal(mean(res$ipw.weights), 1, tolerance = 0.2)
})

# ---- Smooth terms (mgcv) -----------------------------------------------------

test_that("ipwpoint fits smooth terms with mgcv for a binomial exposure", {
  skip_if_not_installed("mgcv")
  s <- mgcv::s
  set.seed(10)
  d <- data.frame(l = rnorm(300))
  d$a <- rbinom(300, 1, plogis(sin(2 * d$l)))
  res <- ipwpoint(exposure = a, family = "binomial", link = "logit",
                  numerator = ~ 1, denominator = ~ s(l), data = d)

  expect_s3_class(res$den.mod, "gam")
  g <- mgcv::gam(a ~ s(l), family = binomial(), data = d)
  p <- as.numeric(fitted(g))
  pm <- mean(d$a)
  expected <- ifelse(d$a == 1, pm, 1 - pm) / ifelse(d$a == 1, p, 1 - p)
  expect_equal(unname(res$ipw.weights), expected, tolerance = 1e-6)
})

test_that("ipwpoint fits smooth terms with mgcv for a gaussian exposure", {
  skip_if_not_installed("mgcv")
  s <- mgcv::s
  set.seed(11)
  d <- data.frame(l = rnorm(300))
  d$a <- sin(2 * d$l) + rnorm(300, 0, 0.5)
  res <- ipwpoint(exposure = a, family = "gaussian", numerator = ~ 1,
                  denominator = ~ s(l), data = d)

  g <- mgcv::gam(a ~ s(l), data = d)
  m0 <- glm(a ~ 1, data = d)
  expected <- dnorm(d$a, fitted(m0), sd(residuals(m0))) /
    dnorm(d$a, as.numeric(fitted(g)), sqrt(g$sig2))
  expect_equal(unname(res$ipw.weights), unname(expected), tolerance = 1e-6)
})

test_that("ipwtm fits smooth terms with mgcv for a binomial exposure", {
  skip_if_not_installed("mgcv")
  s <- mgcv::s
  set.seed(12)
  dat <- data.frame(id = rep(1:100, each = 3), time = rep(1:3, 100))
  dat$l <- rnorm(300)
  dat$a <- rbinom(300, 1, plogis(sin(2 * dat$l)))
  res <- ipwtm(exposure = a, family = "binomial", link = "logit", numerator = ~ 1,
               denominator = ~ s(l), id = id, timevar = time, type = "all", data = dat)

  expect_s3_class(res$den.mod, "gam")
  g <- mgcv::gam(a ~ s(l), family = binomial(), data = dat)
  p <- as.numeric(fitted(g))
  pm <- mean(dat$a)
  ratio <- ifelse(dat$a == 1, pm, 1 - pm) / ifelse(dat$a == 1, p, 1 - p)
  expect_equal(unname(res$ipw.weights), unname(ave(ratio, dat$id, FUN = cumprod)),
               tolerance = 1e-6)
})

test_that("smooth terms are refused where they are not supported", {
  dat <- make_long_categorical()
  expect_error(
    ipwpoint(exposure = a, family = "multinomial", denominator = ~ s(l), data = dat),
    "Smooth terms"
  )
  expect_error(
    ipwtm(exposure = a, family = "multinomial", denominator = ~ s(l), id = id,
          timevar = time, type = "first", data = dat),
    "Smooth terms"
  )
})

# ---- Estimands for a binomial point treatment --------------------------------

test_that("ipwpoint estimand = 'ATT' and 'ATO' give the standard weights", {
  set.seed(13)
  d <- data.frame(l = rnorm(300))
  d$a <- rbinom(300, 1, plogis(0.8 * d$l))
  p <- fitted(glm(a ~ l, family = binomial, data = d))

  att <- ipwpoint(exposure = a, family = "binomial", link = "logit",
                  denominator = ~ l, data = d, estimand = "ATT")
  expect_equal(unname(att$ipw.weights), unname(ifelse(d$a == 1, 1, p / (1 - p))))

  ato <- ipwpoint(exposure = a, family = "binomial", link = "logit",
                  denominator = ~ l, data = d, estimand = "ATO")
  expect_equal(unname(ato$ipw.weights), unname(ifelse(d$a == 1, 1 - p, p)))

  ate <- ipwpoint(exposure = a, family = "binomial", link = "logit",
                  denominator = ~ l, data = d)
  expect_equal(unname(ate$ipw.weights), unname(1 / ifelse(d$a == 1, p, 1 - p)))
})

test_that("ipwpoint refuses ATT/ATO for other families or with a numerator", {
  set.seed(14)
  d <- data.frame(l = rnorm(100))
  d$a <- rbinom(100, 1, 0.5)
  d$c <- rpois(100, 2)
  expect_error(ipwpoint(exposure = c, family = "poisson", denominator = ~ l, data = d,
                        estimand = "ATT"), "only available")
  expect_error(ipwpoint(exposure = a, family = "binomial", link = "logit", numerator = ~ 1,
                        denominator = ~ l, data = d, estimand = "ATO"), "must be NULL")
  expect_error(ipwpoint(exposure = a, family = "binomial", link = "logit",
                        denominator = ~ l, data = d, estimand = "XYZ"), "should be one of")
})
