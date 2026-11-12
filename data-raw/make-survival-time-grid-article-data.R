# Data for vignettes/articles/survival-weights-time-grid.Rmd.
#
# Installs every version of ipw from the CRAN archive (plus the current CRAN
# release), each into its own library, runs the same survival cases in each
# (in a separate R process with callr), and keeps the source of ipwtm() from
# three old versions. Saves data, calls, results and code.
# Rerun only when a case is added or changed.
#
# Needs: callr, survival, and the installed development ipw.

stopifnot(requireNamespace("callr"))
suppressPackageStartupMessages(library(survival))
out_file <- file.path("vignettes", "articles", "survival-weights-time-grid.rds")
cran <- "https://cloud.r-project.org"

# ---- versions on CRAN -------------------------------------------------------
idx <- readLines(url(paste0(cran, "/src/contrib/Archive/ipw/")), warn = FALSE)
files <- unique(regmatches(idx, regexpr("ipw_[0-9][0-9.-]*\\.tar\\.gz", idx)))
current <- available.packages(repos = cran)["ipw", "Version"]
tarballs <- c(paste0(cran, "/src/contrib/Archive/ipw/", files),
              paste0(cran, "/src/contrib/ipw_", current, ".tar.gz"))
versions <- sub("^ipw_(.*)\\.tar\\.gz$", "\\1", basename(tarballs))
ord <- order(numeric_version(gsub("-", ".", versions)))
tarballs <- tarballs[ord]; versions <- versions[ord]
cat("CRAN versions found:", versions, "\n")

install_one <- function(tb, v) {
  lib <- file.path(tempdir(), paste0("ipw-", v))
  dir.create(lib, showWarnings = FALSE)
  msgs <- character()
  log <- utils::capture.output(res <- withCallingHandlers(
    tryCatch(utils::install.packages(tb, repos = NULL, type = "source", lib = lib, quiet = TRUE),
             error = function(e) e),
    warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") }),
    type = "message")
  ok <- file.exists(file.path(lib, "ipw", "DESCRIPTION"))
  list(lib = lib, ok = ok,
       message = if (ok) NA_character_ else
         paste(c(if (inherits(res, "error")) conditionMessage(res), msgs, tail(log, 3)), collapse = " | "))
}
installed <- Map(install_one, tarballs, versions)
names(installed) <- versions

# ---- data -------------------------------------------------------------------
data("haartdat", package = "ipw", envir = environment())

# Subjects with their own visit times. Rows run up to and including the first
# exposure (type = "first" data); x is a baseline covariate, l time-varying.
make_own_times <- function(n, seed) {
  set.seed(seed)
  rows <- lapply(seq_len(n), function(i) {
    tt <- sort(unique(round(cumsum(runif(6, 5, 15)))))
    data.frame(id = i, time = tt, x = rbinom(1, 1, 0.5), l = rnorm(length(tt)))
  })
  d <- do.call(rbind, rows)
  d$tstart <- ave(d$time, d$id, FUN = function(t) c(0, head(t, -1)))
  d$a <- rbinom(nrow(d), 1,
                1 - exp(-0.02 * (d$time - d$tstart) * exp(0.7 * d$l + 0.3 * d$x)))
  d <- d[ave(d$a, d$id, FUN = function(z) cumsum(cumsum(z)) <= 1) == 1, ]
  # censoring (drop-out) at the end of a row, for type = "cens"
  d$dropout <- rbinom(nrow(d), 1,
                      1 - exp(-0.01 * (d$time - d$tstart) * exp(0.5 * d$l)))
  rownames(d) <- NULL
  d
}
d_own <- make_own_times(300, 20261111)
# for the censoring weights: follow-up ends at the first drop-out
d_cens <- d_own[ave(d_own$dropout, d_own$id, FUN = function(z) cumsum(cumsum(z)) <= 1) == 1, ]
rownames(d_cens) <- NULL

# The same data split at every time that occurs in the data (one common grid,
# as in Example 2 of ?ipwtm). Covariates are constant within an original row;
# an exposure at the end of an original row is put on its last piece.
split_common_grid <- function(d) {
  grid <- sort(unique(c(d$tstart, d$time)))
  out <- lapply(seq_len(nrow(d)), function(r) {
    cuts <- grid[grid > d$tstart[r] & grid < d$time[r]]
    s <- c(d$tstart[r], cuts); e <- c(cuts, d$time[r]); k <- length(s)
    data.frame(id = d$id[r], tstart = s, time = e, x = d$x[r], l = d$l[r],
               a = c(rep(0, k - 1), d$a[r]), orig_row = r, last_piece = seq_len(k) == k)
  })
  do.call(rbind, out)
}
d_grid <- split_common_grid(d_own)

# Small example with 3 patients for the worked example in the article
d_tiny <- data.frame(
  id     = c(1, 1, 2, 3, 3),
  tstart = c(0, 10, 0, 0, 5),
  time   = c(10, 20, 15, 5, 20),
  l      = c(0.5, 1.2, -0.3, 0.1, 0.8),
  a      = c(0, 1, 1, 0, 0))

# ---- cases ------------------------------------------------------------------
cases <- list(
  haartdat_first = list(
    data = NULL,  # haartdat, loaded from the package
    call = quote(ipwtm(exposure = haartind, family = "survival",
                       numerator = ~ sex + age, denominator = ~ sex + age + cd4.sqrt,
                       id = patient, tstart = tstart, timevar = fuptime,
                       type = "first", data = data))),
  haartdat_cens = list(
    data = NULL,
    call = quote(ipwtm(exposure = dropout, family = "survival",
                       numerator = ~ sex + age, denominator = ~ sex + age + cd4.sqrt,
                       id = patient, tstart = tstart, timevar = fuptime,
                       type = "cens", data = data))),
  own_first = list(
    data = d_own,
    call = quote(ipwtm(exposure = a, family = "survival", denominator = ~ x + l,
                       id = id, tstart = tstart, timevar = time,
                       type = "first", data = data))),
  own_first_stab = list(
    data = d_own,
    call = quote(ipwtm(exposure = a, family = "survival", numerator = ~ x,
                       denominator = ~ x + l, id = id, tstart = tstart,
                       timevar = time, type = "first", data = data))),
  own_cens = list(
    data = d_cens,
    call = quote(ipwtm(exposure = dropout, family = "survival", numerator = ~ x,
                       denominator = ~ x + l, id = id, tstart = tstart,
                       timevar = time, type = "cens", data = data))),
  common_grid_first = list(
    data = d_grid,
    call = quote(ipwtm(exposure = a, family = "survival", denominator = ~ x + l,
                       id = id, tstart = tstart, timevar = time,
                       type = "first", data = data))),
  tiny = list(
    data = d_tiny,
    call = quote(ipwtm(exposure = a, family = "survival", denominator = ~ l,
                       id = id, tstart = tstart, timevar = time,
                       type = "first", data = data)))
)

# ---- run --------------------------------------------------------------------
run_cases <- function(cases, lib) {
  callr::r(function(cases) {
    suppressPackageStartupMessages({ library(ipw); library(survival) })
    res <- lapply(cases, function(cs) {
      warns <- character()
      fit <- withCallingHandlers(
        tryCatch(eval(cs$call, list(data = cs$data)), error = function(e) e),
        warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
      err <- inherits(fit, "error")
      list(weights = if (err) NULL else as.numeric(fit$ipw.weights),
           error = if (err) conditionMessage(fit) else NA_character_,
           warnings = unique(warns))
    })
    list(version = as.character(utils::packageVersion("ipw")), results = res)
  }, args = list(cases), libpath = c(lib, .libPaths()))
}
with_data <- function(cases) {
  for (nm in c("haartdat_first", "haartdat_cens")) cases[[nm]]$data <- haartdat
  cases
}

by_version <- lapply(names(installed), function(v) {
  inst <- installed[[v]]
  if (!inst$ok) return(list(version = v, installed = FALSE, message = inst$message, results = NULL))
  r <- tryCatch(run_cases(with_data(cases), inst$lib), error = function(e)
    list(version = v, results = NULL, message = conditionMessage(e)))
  c(r, list(installed = TRUE))
})
names(by_version) <- versions

dev <- run_cases(with_data(cases), .libPaths()[1])

# ---- source code of ipwtm() in old versions (for the article) -------------
code_versions <- c("1.0", "1.2.1.1", "1.3.2")
old_code <- lapply(code_versions, function(v) {
  tb <- tarballs[versions == v]
  f <- file.path(tempdir(), basename(tb))
  utils::download.file(tb, f, quiet = TRUE, mode = "wb")
  ex <- file.path(tempdir(), paste0("src-", v))
  utils::untar(f, files = "ipw/R/ipwtm.R", exdir = ex)
  readLines(file.path(ex, "ipw", "R", "ipwtm.R"), warn = FALSE)
})
names(old_code) <- code_versions
for (v in code_versions) {
  x <- old_code[[v]]
  hit <- grep("basehaz|bashaz|diff\\(|step_haz|cum_haz|merge\\(", x)
  cat("\n==== ipwtm.R in", v, "(", length(x), "lines ); lines near the baseline hazard ====\n")
  if (length(hit)) {
    show <- sort(unique(unlist(lapply(hit, function(h) max(1, h - 4):min(length(x), h + 4)))))
    cat(sprintf("%4d  %s", show, x[show]), sep = "\n")
  }
}

saveRDS(list(cases = cases, by_version = by_version, old_code = old_code,
             made = format(Sys.Date()), R = R.version.string),
        out_file, compress = "xz")

# ---- print a summary --------------------------------------------------------
ref_w <- function(d, exposure, num, den, type) {
  ord <- order(d$id, d$tstart)
  stopifnot(identical(ord, seq_len(nrow(d))))
  sel <- ave(d[[exposure]], d$id, FUN = function(z) cumsum(cumsum(z)) <= 1) == 1
  p <- function(rhs) {
    if (is.null(rhs)) return(rep(1, nrow(d)))
    f <- update(rhs, as.formula(paste("Surv(tstart, time,", exposure, ") ~ .")))
    dsel <- d[sel, ]
    fit <- coxph(f, data = dsel, model = TRUE)
    bh <- basehaz(fit, centered = TRUE)
    H <- function(t) { h <- approx(bh$time, bh$hazard, t, method = "constant", rule = 2)$y; h[t < min(bh$time)] <- 0; h }
    out <- rep(1, nrow(d))
    hz <- (H(d$time[sel]) - H(d$tstart[sel])) * predict(fit, type = "risk")
    # censoring weights use P(not censored) on every row, also the censored one
    out[sel] <- if (type == "cens") exp(-hz) else ifelse(d[[exposure]][sel] == 1, 1 - exp(-hz), exp(-hz))
    out
  }
  ave(p(num) / p(den), d$id, FUN = cumprod)
}
refs <- list(
  own_first = ref_w(d_own, "a", NULL, ~ x + l, "first"),
  own_first_stab = ref_w(d_own, "a", ~ x, ~ x + l, "first"),
  own_cens = ref_w(d_cens, "dropout", ~ x, ~ x + l, "cens"),
  common_grid_first = ref_w(d_grid, "a", NULL, ~ x + l, "first"))
for (v in c(names(by_version), "dev")) {
  r <- if (v == "dev") dev else by_version[[v]]
  cat(sprintf("\n== %s ==\n", if (v == "dev") dev$version else v))
  if (is.null(r$results)) { cat("  not run:", r$message, "\n"); next }
  for (nm in names(cases)) {
    x <- r$results[[nm]]
    s <- if (!is.na(x$error)) paste("ERROR:", x$error) else if (nm %in% names(refs))
      sprintf("max rel diff from correct %.3g", max(abs(x$weights / refs[[nm]] - 1))) else
      sprintf("ok, range %.4g-%.4g", min(x$weights), max(x$weights))
    cat(sprintf("  %-18s %s\n", nm, s))
  }
}
cat("\nSaved", out_file, "(", file.size(out_file), "bytes )\n")
