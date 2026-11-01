# ==========================================================================
# Internal helpers shared by ipwpoint() and ipwtm()
# ==========================================================================

.valid_links <- list(
  binomial = c("logit", "probit", "cauchit", "log", "cloglog"),
  ordinal  = c("logit", "probit", "cauchit", "cloglog")
)

.check_link <- function(family, link) {
  valid <- .valid_links[[family]]
  if (is.null(valid)) return(invisible(TRUE))
  if (is.null(link) || length(link) != 1 || !(link %in% valid)) {
    stop("Invalid 'link' specified for family = \"", family, "\". Must be one of: ",
         paste(sprintf('"%s"', valid), collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.polr_method <- function(link) {
  switch(link, logit = "logistic", probit = "probit", cloglog = "cloglog", cauchit = "cauchit")
}

.check_trunc <- function(trunc) {
  if (!is.null(trunc) &&
      (!is.numeric(trunc) || length(trunc) != 1 || is.na(trunc) || trunc < 0 || trunc > 0.5)) {
    stop("'trunc' must be a single numeric percentile between 0 and 0.5.", call. = FALSE)
  }
  invisible(TRUE)
}

# Column name from an unquoted symbol or a string
.get_col_name <- function(expr) {
  if (is.symbol(expr) || is.character(expr)) as.character(expr) else deparse(expr, width.cutoff = 500)
}

# Build "lhs ~ rhs" from a one-sided formula, keeping its environment
.build_model_formula <- function(lhs_str, rhs_formula) {
  tt <- stats::terms(rhs_formula)
  rhs_terms <- attr(tt, "term.labels")
  intercept <- attr(tt, "intercept") == 1
  rhs_str <- if (length(rhs_terms) == 0) {
    if (intercept) "1" else "0"
  } else {
    paste(rhs_terms, collapse = " + ")
  }
  if (length(rhs_terms) > 0 && !intercept) rhs_str <- paste(rhs_str, "- 1")
  stats::as.formula(paste(lhs_str, "~", rhs_str), env = environment(rhs_formula))
}

# Truncate weights at the trunc and 1 - trunc quantiles
.truncate_weights <- function(w, trunc) {
  q <- stats::quantile(w, probs = c(trunc, 1 - trunc), na.rm = TRUE, names = FALSE)
  pmin(pmax(w, q[1]), q[2])
}

# Probability of the observed category from predict(<multinom/polr>, type = "probs").
# Columns are matched by level name, never by position. Handles the vector that
# multinom returns for a two-level exposure (P of the second level).
.prob_of_observed <- function(probs, obs, levels_all) {
  if (is.null(dim(probs))) {
    if (length(levels_all) == 2L && length(probs) == length(obs)) {
      probs <- cbind(1 - probs, probs)
    } else {
      probs <- matrix(probs, nrow = length(obs), byrow = TRUE)
    }
    colnames(probs) <- levels_all
  }
  probs <- as.matrix(probs)
  idx <- match(as.character(obs), colnames(probs))
  if (anyNA(idx)) {
    stop("Observed exposure levels do not match fitted probability levels: ",
         paste(unique(as.character(obs)[is.na(idx)]), collapse = ", "), call. = FALSE)
  }
  probs[cbind(seq_len(nrow(probs)), idx)]
}

# Add default arguments unless the user supplied them through `...`
.with_defaults <- function(args, dots, defaults) {
  defaults <- defaults[setdiff(names(defaults), names(dots))]
  c(args, defaults, dots)
}

# Replace data copies that do.call() inlines into a model's call with the
# user's own expressions, so print()/summary() of the model stay readable.
.tidy_model_call <- function(mod, replacements) {
  if (is.null(mod) || is.null(mod$call)) return(mod)
  for (nm in names(replacements)) {
    if (!is.null(mod$call[[nm]]) && !is.null(replacements[[nm]])) {
      mod$call[[nm]] <- replacements[[nm]]
    }
  }
  mod
}

# ---- Generalized additive models (mgcv) --------------------------------------

.smooth_fns <- c("s", "te", "ti", "t2")

# TRUE if a formula calls an mgcv smooth such as s(x), te(x, z) or mgcv::s(x).
# Only function calls count: a covariate that happens to be named s, te, ti or
# t2 is an ordinary variable.
.uses_smooths <- function(f) {
  if (is.null(f)) return(FALSE)
  is_smooth_call <- function(e) {
    if (!is.call(e)) return(FALSE)
    fn <- e[[1]]
    if (is.name(fn) && as.character(fn) %in% .smooth_fns) return(TRUE)
    if (is.call(fn) && identical(fn[[1]], as.name("::")) &&
        as.character(fn[[3]]) %in% .smooth_fns) return(TRUE)
    any(vapply(as.list(e)[-1], is_smooth_call, logical(1)))
  }
  is_smooth_call(f)
}

.require_mgcv <- function() {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is needed for smooth terms (s(), te(), ti(), t2()). ",
         "Please install it.", call. = FALSE)
  }
}

# Give a formula an environment in which s(), te(), ti() and t2() are mgcv's,
# so smooth terms work without the user attaching mgcv.
.gam_formula <- function(form) {
  env <- new.env(parent = environment(form))
  for (fn in .smooth_fns) assign(fn, getExportedValue("mgcv", fn), envir = env)
  environment(form) <- env
  form
}

# Fit a GLM, or a GAM when the formula has smooth terms. `args` holds formula,
# data, family and any further arguments for glm()/gam().
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
