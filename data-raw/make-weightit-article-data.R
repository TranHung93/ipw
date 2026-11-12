# Data for vignettes/articles/ipw-and-weightit.Rmd.
#
# Bootstrap standard errors of the MSM coefficient for the binary point
# treatment example, for each estimand: the weights are re-estimated with
# ipwpoint() in each bootstrap sample. Everything else in the article is
# computed when the article is built.
# Rerun when the example data or the number of bootstrap samples change.
#
# Needs: the installed development ipw.

library(ipw)
out_file <- file.path("vignettes", "articles", "ipw-and-weightit.rds")

# Same data as in the article (section 1)
make_binary <- function() {
  set.seed(2026)
  n <- 1000
  d <- data.frame(l = stats::rnorm(n), x = stats::rbinom(n, 1, 0.4))
  d$a <- stats::rbinom(n, 1, stats::plogis(-0.3 + 0.8 * d$l + 0.5 * d$x))
  d$y <- 2 * d$a + d$l + d$x + stats::rnorm(n)
  d
}
d <- make_binary()

estimands <- c("ATE", "ATT", "ATC", "ATO")
R <- 2000
effect <- function(dd, est) {
  w <- ipwpoint(exposure = a, family = "binomial", link = "logit",
                denominator = ~ l + x, data = dd, estimand = est)$ipw.weights
  unname(stats::coef(stats::lm(y ~ a, data = dd, weights = w))["a"])
}
set.seed(1)
idx <- replicate(R, sample(nrow(d), replace = TRUE), simplify = FALSE)
boot <- sapply(estimands, function(est) vapply(idx, function(i) effect(d[i, ], est), numeric(1)))

saveRDS(list(R = R, boot_se = apply(boot, 2, stats::sd),
             data_check = sum(d$y), ipw_version = as.character(utils::packageVersion("ipw")),
             date = Sys.Date()),
        out_file)
cat("Saved", out_file, "\n")
