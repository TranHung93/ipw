# Results of the old CRAN versions for vignettes/articles/bug-fixes.Rmd.
#
# Runs each bug case in ipw 1.2.1.1 and ipw 1.3.2 (installed from CRAN into
# separate libraries, each in its own R process via callr) and saves the
# data, the calls and the old results. The article recomputes the correct
# weights and the current version's weights itself, so only the old results
# are stored. Rerun this script only when a case is added or changed.
#
# Needs: callr, remotes, survival, MASS, and the installed development ipw.

stopifnot(requireNamespace("callr"), requireNamespace("remotes"))
out_file <- file.path("vignettes", "articles", "bug-fixes-old-versions.rds")

# ---- data --------------------------------------------------------------------
set.seed(20261111)
n <- 500

# Binomial exposure whose probability depends on follow-up length through an offset
d_offset <- data.frame(l = rnorm(n), fup = runif(n, 0.5, 4))
d_offset$off <- log(d_offset$fup)
d_offset$a <- rbinom(n, 1, plogis(-1.5 + 0.8 * d_offset$l + d_offset$off))

# Binomial exposure, only covariates besides the exposure (for ~ .)
d_dot <- data.frame(l1 = rnorm(n), l2 = rbinom(n, 1, 0.4))
d_dot$a <- rbinom(n, 1, plogis(-0.3 + 0.8 * d_dot$l1 - 0.6 * d_dot$l2))

# 0/1 exposure stored as a factor with levels in the order "1", "0"
d_factor <- d_dot[, c("l1", "l2")]
a_num <- rbinom(n, 1, plogis(-0.3 + 0.8 * d_factor$l1))
d_factor$a <- factor(a_num, levels = c(1, 0))
d_factor$a_num <- a_num

# Ordinal exposure given as text
d_ordinal <- data.frame(l = rnorm(n))
lat <- 0.9 * d_ordinal$l + rlogis(n)
d_ordinal$dose <- c("low", "mid", "high")[cut(lat, c(-Inf, -0.5, 1, Inf), labels = FALSE)]

# Survival exposure, subjects with their own visit times (type = "first")
rows <- lapply(seq_len(300), function(i) {
  tt <- sort(unique(round(cumsum(runif(6, 5, 15)))))
  data.frame(id = i, time = tt, l = rnorm(length(tt)))
})
d_surv <- do.call(rbind, rows)
d_surv$tstart <- ave(d_surv$time, d_surv$id, FUN = function(t) c(0, head(t, -1)))
d_surv$a <- rbinom(nrow(d_surv), 1,
                   1 - exp(-0.02 * (d_surv$time - d_surv$tstart) * exp(0.7 * d_surv$l)))
d_surv <- d_surv[ave(d_surv$a, d_surv$id, FUN = function(x) cumsum(cumsum(x)) <= 1) == 1, ]
rownames(d_surv) <- NULL

# Time-varying binomial exposure with a covariate that happens to be called selvar
d_long <- data.frame(id = rep(1:200, each = 5), time = rep(1:5, 200))
d_long$l <- rnorm(nrow(d_long))
d_long$selvar <- rnorm(nrow(d_long))
d_long$a <- rbinom(nrow(d_long), 1, plogis(-2 + d_long$l + 0.8 * d_long$selvar))

# The same data with three rows entered twice
d_dup <- rbind(d_long, d_long[c(1, 2, 3), ])

# haartdat for the ties case (common time grid, many tied event times)
data("haartdat", package = "ipw", envir = environment())

# ---- cases -------------------------------------------------------------------
cases <- list(
  offset = list(
    data = d_offset,
    call = quote(ipwpoint(exposure = a, family = "binomial", link = "logit",
                          denominator = ~ l + offset(off), data = data))),
  dot = list(
    data = d_dot,
    call = quote(ipwpoint(exposure = a, family = "binomial", link = "logit",
                          denominator = ~ ., data = data))),
  factor = list(
    data = d_factor,
    call = quote(ipwpoint(exposure = a, family = "binomial", link = "logit",
                          denominator = ~ l1 + l2, data = data))),
  ordinal = list(
    data = d_ordinal,
    call = quote(ipwpoint(exposure = dose, family = "ordinal", link = "logit",
                          denominator = ~ l, data = data))),
  survival = list(
    data = d_surv,
    call = quote(ipwtm(exposure = a, family = "survival", denominator = ~ l,
                       id = id, tstart = tstart, timevar = time,
                       type = "first", data = data))),
  selvar = list(
    data = d_long,
    call = quote(ipwtm(exposure = a, family = "binomial", link = "logit",
                       denominator = ~ l + selvar, id = id, timevar = time,
                       type = "first", data = data))),
  duplicates = list(
    data = d_dup,
    call = quote(ipwtm(exposure = a, family = "binomial", link = "logit",
                       denominator = ~ l, id = id, timevar = time,
                       type = "first", data = data))),
  ties = list(
    data = haartdat,
    call = quote(ipwtm(exposure = haartind, family = "survival",
                       numerator = ~ sex + age, denominator = ~ sex + age + cd4.sqrt,
                       id = patient, tstart = tstart, timevar = fuptime,
                       type = "first", data = data, ties = "breslow")))
)

# ---- run in a given library -------------------------------------------------
run_cases <- function(cases, lib) {
  callr::r(function(cases) {
    suppressPackageStartupMessages({ library(ipw); library(survival) })
    res <- lapply(cases, function(cs) {
      warns <- character()
      fit <- withCallingHandlers(
        tryCatch(eval(cs$call, list(data = cs$data)), error = function(e) e),
        warning = function(w) {
          warns <<- c(warns, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
      err <- inherits(fit, "error")
      list(weights  = if (err) NULL else as.numeric(fit$ipw.weights),
           error    = if (err) conditionMessage(fit) else NA_character_,
           warnings = unique(warns))
    })
    list(version = as.character(utils::packageVersion("ipw")), results = res)
  }, args = list(cases), libpath = c(lib, .libPaths()))
}

libs <- c("1.2.1.1" = file.path(tempdir(), "ipw-1.2.1.1"),
          "1.3.2"   = file.path(tempdir(), "ipw-1.3.2"))
for (v in names(libs)) {
  dir.create(libs[[v]], showWarnings = FALSE)
  remotes::install_version("ipw", version = v, lib = libs[[v]], upgrade = "never",
                           dependencies = FALSE, repos = "https://cloud.r-project.org",
                           quiet = TRUE)
}

old <- lapply(libs, function(lib) run_cases(cases, lib))
got <- unname(vapply(old, `[[`, "", "version"))
cat("Versions run:", got, "\n")
stopifnot(identical(got, names(libs)))

# Strip the large haartdat copy: the article loads it from the package
cases$ties$data <- NULL
saveRDS(list(cases = cases, old = old,
             made = format(Sys.Date()), R = R.version.string),
        out_file, compress = "xz")

# ---- print a summary (results of the current version too, for checking) ------
dev <- run_cases(c(cases[names(cases) != "ties"],
                   list(ties = list(data = haartdat, call = cases$ties$call))), .libPaths()[1])
show <- function(r) if (!is.na(r$error)) paste("ERROR:", r$error) else
  sprintf("ok, n = %d, range %.4g-%.4g%s", length(r$weights), min(r$weights), max(r$weights),
          if (length(r$warnings)) paste(" | WARNING:", paste(r$warnings, collapse = " / ")) else "")
for (nm in names(cases)) {
  cat("\n==", nm, "==\n")
  for (v in names(old)) cat(sprintf("  %-8s %s\n", v, show(old[[v]]$results[[nm]])))
  cat(sprintf("  %-8s %s\n", dev$version, show(dev$results[[nm]])))
}
cat("\nSaved", out_file, "(", file.size(out_file), "bytes )\n")
