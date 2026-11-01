# Compare ipwtm(family = "survival", type = "all") with SVall(), the function
# published with Graffeo, Latouche, Geskus & Chevret (2018, Biometrical
# Journal 60, 323-332), on the authors' own data.
#
# Their code and data are not part of ipw. Unzip the paper's supporting
# information and set `graffeo_dir` to its Code_and_Data folder. Then, from the
# ipw package folder:
#   Rscript dev/compare-graffeo.R "path/to/Code_and_Data"
# or set graffeo_dir below and source() it after devtools::load_all().

args <- commandArgs(trailingOnly = TRUE)
graffeo_dir <- if (length(args)) args[1] else "path/to/Code_and_Data"
stopifnot(dir.exists(graffeo_dir))

library(survival)
if (!"ipw" %in% loadedNamespaces() || is.null(getNamespaceInfo("ipw", "path"))) {
  devtools::load_all(".")
}

their <- new.env()
sys.source(file.path(graffeo_dir, "code", "functions", "processing.data.R"), envir = their)
sys.source(file.path(graffeo_dir, "code", "functions", "SVall.R"), envir = their)
environment(their$process.data) <- their
their$tstartfun <- ipw::tstartfun
load(file.path(graffeo_dir, "data", "dataExample.RData"), envir = their)

compare <- function(label, ours, theirs, id, time) {
  o <- order(id, time)
  d <- abs(ours[o] - theirs[o])
  cat(sprintf("%-38s max |difference| = %.3g   (rows: %d)\n", label, max(d), length(d)))
  invisible(max(d))
}
pct <- function(w) round(c(min = min(w), p01 = unname(quantile(w, 0.01)), mean = mean(w),
                           p99 = unname(quantile(w, 0.99)), max = max(w)), 4)

# ---- Section 3.3: dataExample, one exposure process --------------------------
startstop <- their$process.data(data = their$dataExample, exposure = expo, id = id,
                                tstart = tstart, fuptime = fuptime,
                                timefixedcov = x, timevarconf = xt)

svall <- their$SVall(exposure = expo, numerator = ~ x, denominator = ~ x + xt,
                     id = id, tstart = tstart, timevar = fuptime, data = startstop, trunc = 0.01)
ours <- ipwtm(exposure = expo, family = "survival", numerator = ~ x, denominator = ~ x + xt,
              id = id, tstart = tstart, timevar = fuptime, type = "all",
              data = startstop, trunc = 0.01)

cat("\nSection 3.3, dataExample (stabilized weights)\n")
compare("ipw.weights", ours$ipw.weights, svall$ipw.weights, startstop$id, startstop$fuptime)
compare("weights.trunc", ours$weights.trunc, svall$weights.trunc, startstop$id, startstop$fuptime)
cat("Paper:  min 0.1349, 1st pct 0.4553, mean 1.0013, 99th pct 2.4511, max 11.0541\n")
cat("SVall: ", format(pct(svall$ipw.weights)), "\n")
cat("ipwtm: ", format(pct(ours$ipw.weights)), "\n")

svall_u <- their$SVall(exposure = expo, denominator = ~ x + xt, id = id, tstart = tstart,
                       timevar = fuptime, data = startstop)
ours_u <- ipwtm(exposure = expo, family = "survival", denominator = ~ x + xt, id = id,
                tstart = tstart, timevar = fuptime, type = "all", data = startstop)
compare("unstabilized ipw.weights", ours_u$ipw.weights, svall_u$ipw.weights,
        startstop$id, startstop$fuptime)

# ---- Section 3.4: one simulated data set, initiation and discontinuation -----
sim_files <- c("mySimrecTDC.R", "CumHazRange.R", "invH.R", "invHrecstep.R", "generEvent.R")
for (f in sim_files) sys.source(file.path(graffeo_dir, "code", "functions", "DataGeneration", f),
                                envir = their)
sys.source(file.path(graffeo_dir, "code", "functions", "processing.data.duration.2.R"), envir = their)
for (fn in ls(their)) if (is.function(their[[fn]])) environment(their[[fn]]) <- their

set.seed(103 + 2015)
sim <- their$mySimrecTDC(N = 200, beta.x = 0.5, beta.e = 1.1, beta.x.TD = -1.2,
                         beta.x.TD2 = -0.8, par.rec = c(1, 2), dfree = 1 / 12)$ipwtab
sim <- their$generEvent(sim, HR_e = 1.3, HR_xt = 1.2, 5, 5)
ss <- their$process.data.dur.2(data = sim, exposureI = Init, exposureT = Term, id = id,
                               tstart = tstart, fuptime = fuptime, timefixedcov1 = x,
                               timevarconf1 = xt, expTrtInterest = expo, event = event)

cat("\nSection 3.4, one simulated data set\n")
for (proc in c("Init", "Term")) {
  sub <- ss[!is.na(ss[[proc]]), ]
  # both functions take the exposure as an unquoted column name
  sv <- eval(bquote(their$SVall(exposure = .(as.name(proc)), numerator = ~ x,
                                denominator = ~ x + xt, id = id, tstart = tstart,
                                timevar = fuptime, data = sub)))
  ot <- eval(bquote(ipwtm(exposure = .(as.name(proc)), family = "survival", numerator = ~ x,
                          denominator = ~ x + xt, id = id, tstart = tstart,
                          timevar = fuptime, type = "all", data = sub)))
  compare(paste(proc, "ipw.weights"), ot$ipw.weights, sv$ipw.weights, sub$id, sub$fuptime)
}
