test_that("ipwboot resampling rows reproduces a hand-written bootstrap", {
  set.seed(31)
  d <- data.frame(y = rnorm(40))
  stat <- function(dd) mean(dd$y)

  set.seed(5)
  res <- ipwboot(d, stat, R = 30)

  set.seed(5)
  manual <- vapply(1:30, function(b) mean(d$y[sample.int(40, 40, replace = TRUE)]), numeric(1))

  expect_s3_class(res, "ipwboot")
  expect_equal(unname(res$t0), mean(d$y))
  expect_equal(unname(res$t[, 1]), manual)
  expect_equal(unname(res$se), sd(manual))
  expect_equal(unname(res$ci[1, ]), unname(quantile(manual, c(0.025, 0.975))))
  expect_identical(colnames(res$t), "estimate")
  expect_identical(res$n.fail, 0L)
})

test_that("ipwboot with id resamples whole subjects and gives copies new ids", {
  set.seed(32)
  d <- data.frame(id = rep(1:10, times = 1:10), v = 1)   # subject i has i rows
  stat <- function(dd) {
    c(n.ids = length(unique(dd$id)),
      rows.ok = all(table(dd$id) %in% 1:10),
      ids.consecutive = identical(sort(unique(dd$id)), 1:10))
  }
  res <- ipwboot(d, stat, id = id, R = 25)

  expect_true(all(res$t[, "n.ids"] == 10))          # duplicates count as separate subjects
  expect_true(all(res$t[, "rows.ok"] == 1))          # each copy keeps all its rows
  expect_true(all(res$t[, "ids.consecutive"] == 1))  # new ids 1..10
  expect_identical(res$id, "id")

  # the same with the id column given as a string
  set.seed(7); a <- ipwboot(d, function(dd) nrow(dd), id = id, R = 10)
  set.seed(7); b <- ipwboot(d, function(dd) nrow(dd), id = "id", R = 10)
  expect_equal(a$t, b$t)
})

test_that("ipwboot keeps the number of subjects per stratum fixed", {
  d <- data.frame(id = rep(1:12, each = 2), g = rep(c("a", "b"), times = c(8, 16)))
  stat <- function(dd) c(a = sum(dd$g == "a"), b = sum(dd$g == "b"))
  set.seed(33)
  res <- ipwboot(d, stat, id = id, strata = g, R = 20)
  expect_true(all(res$t[, "a"] == 8))
  expect_true(all(res$t[, "b"] == 16))

  d$bad <- rep(c("x", "y"), 12)
  expect_error(ipwboot(d, stat, id = id, strata = bad, R = 5), "constant within each subject")
})

test_that("ipwboot records failed bootstrap samples", {
  set.seed(34)
  d <- data.frame(y = rnorm(30), g = rbinom(30, 1, 0.5))
  calls <- 0
  stat <- function(dd) {
    calls <<- calls + 1
    if (calls %% 4 == 0) stop("model did not converge")
    mean(dd$y)
  }
  expect_warning(res <- ipwboot(d, stat, R = 20), "bootstrap samples failed")
  expect_identical(res$n.fail, 5L)                 # calls 4, 8, ..., 20 (call 1 is t0)
  expect_identical(res$errors, "model did not converge")
  expect_equal(unname(res$se), sd(res$t[, 1], na.rm = TRUE))

  expect_warning(
    res2 <- ipwboot(d, function(dd) if (nrow(dd) > 0 && runif(1) < 0.3) 1:2 else 1, R = 30),
    "different length"
  )
  expect_true(res2$n.fail > 0)
})

test_that("ipwboot works with ipwpoint and gives a sensible standard error", {
  set.seed(35)
  n <- 400
  d <- data.frame(l = rnorm(n))
  d$a <- rbinom(n, 1, plogis(d$l))
  d$y <- 2 * d$a + d$l + rnorm(n)
  stat <- function(dd) {
    w <- ipwpoint(exposure = a, family = "binomial", link = "logit",
                  numerator = ~ 1, denominator = ~ l, data = dd)$ipw.weights
    coef(lm(y ~ a, data = dd, weights = w))
  }
  set.seed(36)
  res <- ipwboot(d, stat, R = 60)
  expect_identical(colnames(res$t), c("(Intercept)", "a"))
  expect_true(all(res$se > 0))
  expect_true(res$ci["a", 1] < 2 && res$ci["a", 2] > 2)
  ci90 <- confint(res, "a", level = 0.9)
  expect_identical(dimnames(ci90), list("a", c("5 %", "95 %")))
  expect_equal(unname(ci90[1, ]), unname(quantile(res$t[, "a"], c(0.05, 0.95))))
  expect_output(print(res), "60 bootstrap samples, resampling rows")
})

test_that("ipwboot works with ipwtm and subject resampling", {
  data(haartdat, package = "ipw")
  stat <- function(dd) {
    w <- ipwtm(exposure = haartind, family = "survival", numerator = ~ sex + age,
               denominator = ~ sex + age + cd4.sqrt, id = patient, tstart = tstart,
               timevar = fuptime, type = "first", data = dd)$ipw.weights
    coef(survival::coxph(survival::Surv(tstart, fuptime, event) ~ haartind,
                         data = dd, weights = w))
  }
  set.seed(37)
  res <- ipwboot(haartdat, stat, id = patient, R = 3)
  expect_identical(dim(res$t), c(3L, 1L))
  expect_true(is.finite(res$se) && res$se > 0)
  expect_output(print(res), "resampling subjects \\('patient'\\)")
})

test_that("ipwboot rejects invalid input", {
  d <- data.frame(id = 1:5, y = 1:5)
  f <- function(dd) mean(dd$y)
  expect_error(ipwboot(list(y = 1), f), "data")
  expect_error(ipwboot(d, "mean"), "function")
  expect_error(ipwboot(d, f, R = 1), "'R'")
  expect_error(ipwboot(d, f, R = 10.5), "'R'")
  expect_error(ipwboot(d, f, conf = 1), "'conf'")
  expect_error(ipwboot(d, f, id = nope), "not found")
  expect_error(ipwboot(d, function(dd) "a"), "numeric")
  d$id[2] <- NA
  expect_error(ipwboot(d, f, id = id), "missing values")
})
