# Compare ipwtm() results across all CRAN versions of ipw
#
# Example: unstabilized weights for HAART initiation (JSS paper, p. 14;
# GitHub issue #4).
#
# Each archived version is installed from the CRAN archive into its own
# temporary library and run in a separate R process (via callr), so your
# installed ipw is not touched. All versions use the same haartdat (taken from
# your installed ipw) and the same survival package, so any difference
# comes from ipw itself.
#
# Run from a terminal:
#   Rscript compare-ipw-versions.R
# or source() it in a fresh R session. If you run it in an RStudio session
# where ipw is already loaded, restart R first (Session > Restart R).
# Needs: callr, survival (and an internet connection).

if (!requireNamespace("callr", quietly = TRUE)) install.packages("callr")

cran <- "https://cloud.r-project.org"

# ---- Versions to test -------------------------------------------------------
archive_url <- paste0(cran, "/src/contrib/Archive/ipw/")
archive_html <- tryCatch(readLines(archive_url, warn = FALSE), error = function(e) character())
tarballs <- unique(regmatches(archive_html, regexpr("ipw_[0-9][0-9.-]*\\.tar\\.gz", archive_html)))
if (length(tarballs) == 0) {
  # Fallback if the archive index can't be read
  tarballs <- paste0("ipw_", c("1.0", paste0("1.0-", 1:11), "1.2", "1.2.1", "1.2.1.1", "1.3.0"),
                     ".tar.gz")
}
versions <- sub("^ipw_(.*)\\.tar\\.gz$", "\\1", tarballs)
# Sort versions and tarballs together so urls[i] always belongs to versions[i]
ord <- order(package_version(gsub("-", ".", versions)))
versions <- versions[ord]
tarballs <- tarballs[ord]
urls <- paste0(archive_url, tarballs)
message("Archived versions found: ", paste(versions, collapse = ", "))

# ---- The example, run in a fresh R process ---------------------------------
# haartdat is read in a subprocess so ipw is never loaded in this session
# (on Windows, R refuses to install a package that is loaded).
haartdat <- callr::r(function() {
  e <- new.env()
  utils::data("haartdat", package = "ipw", envir = e)
  e$haartdat
}, show = FALSE)

run_example <- function(haartdat, lib) {
  suppressPackageStartupMessages({
    library(survival)
    # Load ipw from the version's own library when given, never by search path
    if (length(lib)) library(ipw, lib.loc = lib) else library(ipw)
  })
  res <- ipwtm(exposure = haartind, family = "survival", denominator = ~ cd4.sqrt,
               id = patient, tstart = tstart, timevar = fuptime, type = "first",
               data = haartdat)
  ns_path <- getNamespaceInfo("ipw", "path")
  list(version = read.dcf(file.path(ns_path, "DESCRIPTION"), fields = "Version")[1, 1],
       path    = ns_path,
       weights = as.numeric(res$ipw.weights))
}

run_in_lib <- function(lib) {
  callr::r(run_example, args = list(haartdat = haartdat, lib = lib),
           libpath = c(lib, .libPaths()), show = FALSE)
}

# ---- Current installed version (reference) ---------------------------------
current <- run_in_lib(character())
message("Reference: installed ipw ", current$version)

# ---- Old versions ----------------------------------------------------------
results <- list()
results[[current$version]] <- current$weights

for (i in seq_along(versions)) {
  v <- versions[i]
  if (v == current$version) next
  message("ipw ", v, " ...", appendLF = FALSE)
  lib <- file.path(tempdir(), paste0("ipw-", v))
  dir.create(lib, showWarnings = FALSE)
  ok <- tryCatch({
    # Install in a fresh R process where ipw is not loaded
    callr::r(function(url, lib) {
      utils::install.packages(url, repos = NULL, type = "source", lib = lib, quiet = TRUE)
    }, args = list(url = urls[i], lib = lib), show = FALSE)
    if (!dir.exists(file.path(lib, "ipw"))) stop("installation failed")
    out <- run_in_lib(lib)
    # Compare as versions: "1.0-1" and "1.0.1" are the same version
    if (package_version(out$version) != package_version(v)) {
      stop("loaded ipw ", out$version, " from ", out$path, " instead of ", v)
    }
    results[[v]] <- out$weights
    TRUE
  }, error = function(e) {
    message(" failed: ", conditionMessage(e))
    results[[v]] <<- conditionMessage(e)
    FALSE
  })
  if (ok) message(" done")
}

# ---- Summary table ---------------------------------------------------------
ref <- current$weights
rows <- lapply(names(results), function(v) {
  w <- results[[v]]
  if (!is.numeric(w)) {
    return(data.frame(version = v, Min = NA, Q1 = NA, Median = NA, Mean = NA, Q3 = NA,
                      Max = NA, max_abs_diff = NA, identical = NA,
                      note = substr(w, 1, 60)))
  }
  s <- as.numeric(summary(w))
  same_len <- length(w) == length(ref)
  data.frame(version = v, Min = s[1], Q1 = s[2], Median = s[3], Mean = s[4], Q3 = s[5],
             Max = s[6],
             max_abs_diff = if (same_len) max(abs(w - ref)) else NA,
             identical = if (same_len) isTRUE(all.equal(w, ref, tolerance = 1e-10)) else FALSE,
             note = if (same_len) "" else "different length")
})
tab <- do.call(rbind, rows)
tab <- tab[order(package_version(gsub("-", ".", tab$version))), ]
rownames(tab) <- NULL

cat("\nUnstabilized weights, haartdat, ipwtm(family = \"survival\", type = \"first\")\n")
cat("Reference for max_abs_diff / identical: ipw", current$version, "\n")
cat("survival version:", as.character(utils::packageVersion("survival")),
    "| R:", R.version.string, "\n\n")
print(tab, digits = 7, row.names = FALSE)

write.csv(tab, "ipw-version-comparison.csv", row.names = FALSE)
message("\nSaved: ipw-version-comparison.csv")
