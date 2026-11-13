#' Diagnostic Plots for Inverse Probability Weights
#'
#' @description
#' Plots that check whether inverse probability weights can be trusted:
#' overlap of the estimated exposure probabilities (positivity), the share of
#' the total weight held by the people with the largest weights, the tail of
#' the weight distribution, and the sensitivity of a marginal structural
#' model's estimate to truncating the weights. The numbers behind every plot
#' are returned (invisibly), so the plots can be redrawn with any graphics
#' package.
#'
#' @details
#' \describe{
#' \item{\code{"overlap"}}{Densities of the estimated probability of exposure
#'   level 1 from the denominator model, by exposure group (binomial exposure
#'   models only; for \code{\link{ipwtm}} the rows the denominator model was
#'   fitted on). The shaded ranges hold values from one exposure group only:
#'   people there have (almost) no counterparts with the other exposure, a
#'   sign of (near) positivity violations. The summary gives the proportion
#'   of people (rows) with a probability below 0.05 or above 0.95. For
#'   \code{ipwtm()} these are probabilities per interval, which are small
#'   when the intervals are short, so that proportion is not by itself a
#'   warning sign there; look at whether the two groups overlap. For other
#'   exposure families this plot is skipped with a message.}
#' \item{\code{"share"}}{The weights sorted from largest to smallest, with the
#'   cumulative share of the total weight plotted against the fraction of
#'   people. With equal weights the curve lies on the diagonal; a curve that
#'   rises steeply at the left means that a few people carry much of the
#'   weight. For \code{ipwtm()} each person's final cumulative weight (the
#'   weight of the last row in time order) is used. The summary gives the
#'   share of the total weight held by the top 1\% and top 5\% (the largest
#'   \code{ceiling(0.01 * n)} and \code{ceiling(0.05 * n)} weights), the
#'   largest weight, the mean weight and Kish's effective sample size
#'   \code{sum(w)^2 / sum(w^2)}, as a number and as a fraction of \code{n}.
#'   Stabilized (ATE) weights should have a mean close to 1; a mean far from
#'   1 suggests a misspecified exposure model or positivity problems (ATT,
#'   ATC and overlap weights have no such target).}
#' \item{\code{"tail"}}{For the largest \code{top} fraction of the weights
#'   (the same weights as for \code{"share"}), the log of the rank divided by
#'   \code{n} (the empirical probability of a weight at least this large)
#'   against the log of the weight. A Pareto-type tail,
#'   \eqn{P(W > x) \propto x^{-\alpha}}, gives a straight line with slope
#'   \eqn{-\alpha}. The Hill estimate of the tail index uses the \eqn{k}
#'   largest weights \eqn{w_{(1)} \ge \dots \ge w_{(k)}}:
#'   \eqn{\hat\alpha = 1 / \{k^{-1} \sum_{i=1}^{k} \log(w_{(i)} / w_{(k+1)})\}},
#'   with \eqn{k} = \code{floor(top * n)}; its line is drawn. At least 20
#'   weights are needed in the tail, otherwise the plot is skipped with a
#'   message. \eqn{\alpha < 2} means that the weights have an infinite
#'   variance. With \eqn{\alpha} between 2 and 4 the variance is finite but
#'   the third or fourth moment is not, and normal-theory confidence
#'   intervals, with standard errors from \code{\link{ipwvcov}} and from the
#'   bootstrap alike, can be too narrow (see "Heavy-tailed weights" in
#'   \code{\link{ipwvcov}} and in the article "ipw and WeightIt" on the
#'   package website). The summary warns when \eqn{\hat\alpha < 4}.}
#' \item{\code{"trunc"}}{Needs \code{msm}. For each level \eqn{q} in
#'   \code{trunc} the weights are truncated at their \eqn{q} and \eqn{1 - q}
#'   quantiles, exactly as \code{ipwpoint()} and \code{ipwtm()} do with their
#'   \code{trunc} argument (\eqn{q = 0} leaves them unchanged), \code{msm} is
#'   refitted with these weights with \code{stats::update()}, and the
#'   estimate of coefficient \code{coef} is plotted with a 95\% confidence
#'   interval against \eqn{q}. An estimate that hardly changes means that the
#'   extreme weights do not drive the result; a large change means that they
#'   do (truncation trades bias for variance, and changes slightly what is
#'   estimated). The intervals are robust (sandwich) intervals that treat the
#'   weights as known: for \code{lm()} and \code{glm()} the HC0 sandwich
#'   with the small-sample factor \eqn{g/(g - 1)} that \code{survey::svyglm()}
#'   uses (\eqn{g} the number of people; for long data from \code{ipwtm()}
#'   the scores are summed within person), and for \code{coxph()} its robust
#'   variance (\code{robust = TRUE}, or clustered by person for long data:
#'   an existing \code{cluster()} term is kept, otherwise the \code{ipwtm()}
#'   id is used). Truncation at estimated quantiles is not accounted for;
#'   for proper intervals with truncated weights use the bootstrap,
#'   re-estimating and truncating the weights in each sample (section 5 of
#'   \code{vignette("ipw")}). With a list of weight objects the product of
#'   the weights is truncated, which differs from truncating each set of
#'   weights with its own \code{trunc} argument.}
#' }
#'
#' Covariate balance is checked with the cobalt package (see "Covariate
#' balance" in \code{vignette("ipw")}), and \code{\link{ipwplot}} shows the
#' distribution of time-varying weights over follow-up. The diagnostics use
#' the untruncated weights (\code{ipw.weights}).
#'
#' @param object Result of \code{\link{ipwpoint}} or \code{\link{ipwtm}}, or a
#'   list of such results whose weights multiply (e.g. treatment times
#'   censoring weights from \code{ipwtm()} on the same long data, as in
#'   \code{\link{ipwvcov}}). For a list, \code{"overlap"} uses the first
#'   element (put the treatment weights first), and the other plots the
#'   product of the weights.
#' @param which Plots to draw: any of \code{"overlap"}, \code{"share"},
#'   \code{"tail"} and \code{"trunc"}. Plots that cannot be drawn (overlap for
#'   a non-binomial exposure model, truncation without \code{msm}, too few
#'   weights in the tail) are skipped with a message; asking for
#'   \code{"overlap"} or \code{"trunc"} explicitly in these cases is an error.
#' @param msm Optional marginal structural model, fitted with \code{lm()},
#'   \code{glm()} or \code{survival::coxph()} on the same rows, in the same
#'   order, as the weights (one row per person for \code{ipwpoint()}, all rows
#'   of the long data for \code{ipwtm()}). Needed for \code{"trunc"}. Its data
#'   must be available where \code{ipwdiag()} is called, or where the model
#'   was fitted, because the model is refitted.
#' @param trunc Truncation levels for \code{"trunc"}: numbers in
#'   \eqn{[0, 0.5)}.
#' @param coef The coefficient of \code{msm} to show in \code{"trunc"}, by
#'   name or position. Default: the first coefficient other than the
#'   intercept for \code{lm()} and \code{glm()}, the first for \code{coxph()}.
#' @param top Fraction of the largest weights used for \code{"tail"}, in
#'   \eqn{(0, 0.5)}.
#' @param ask Logical. If \code{TRUE}, each plot is drawn on its own page and
#'   R asks before showing the next one. If \code{FALSE} (default) and more
#'   than one plot is drawn, they share one page (side by side for two, a 2 x 2
#'   layout for three or four), unless a multi-panel layout was already set
#'   with \code{par(mfrow = )}. The graphical parameters are restored on exit.
#' @param ... Graphical parameters set with \code{\link[graphics]{par}} while
#'   the plots are drawn, e.g. \code{cex = 0.8}.
#'
#' @return Invisibly, a list of class \code{"ipwdiag"} with one element per
#'   plot that was drawn, holding its data, and an element \code{summary}.
#'   \code{print()} shows the summary.
#' \item{overlap}{List: \code{data} (data frame with the estimated
#'   probability \code{p} and the \code{exposure}, one row per row of the
#'   denominator model), \code{density} (list of two data frames \code{x},
#'   \code{y}, for exposure 0 and 1), \code{support} (the range where both
#'   groups have values), \code{exposure} (its name), \code{xlim} (the
#'   plotted range) and \code{extreme} (as in the summary).}
#' \item{share}{List: \code{data} (data frame with \code{fraction}, the
#'   fraction of people, and \code{share}, the cumulative share of the total
#'   weight, starting at 0, 0) and \code{weights} (the weights used, sorted
#'   from largest to smallest).}
#' \item{tail}{List: \code{data} (data frame with \code{log.weight} and
#'   \code{log.rank}, the log of rank / n, for the \code{k} largest weights),
#'   \code{alpha} (the Hill estimate), \code{k}, \code{threshold}
#'   (\eqn{w_{(k+1)}}), and the \code{intercept} and \code{slope}
#'   (\eqn{-\alpha}) of the drawn line.}
#' \item{trunc}{Data frame with columns \code{q}, \code{estimate}, \code{se},
#'   \code{lower} and \code{upper}; attribute \code{"coef"}.}
#' \item{summary}{List with \code{n} (number of people), \code{unit},
#'   \code{overlap} (proportion with a probability below 0.05 or above 0.95,
#'   and the number of rows), \code{share} (named vector \code{top1},
#'   \code{top5}, \code{max}, \code{mean}, \code{ess}, \code{ess.fraction}),
#'   \code{tail} (\code{alpha}, \code{k}) and \code{trunc} (the table above),
#'   for the plots that were computed.}
#'
#' @author Hung Thai Tran \email{hungtt@oucru.org}
#' @references
#' Cole, S.R. & Hernán, M.A. (2008). Constructing inverse probability weights
#' for marginal structural models. \emph{American Journal of Epidemiology},
#' \bold{168}(6), 656-664. \doi{10.1093/aje/kwn164}.
#'
#' Austin, P.C. & Stuart, E.A. (2015). Moving towards best practice when using
#' inverse probability of treatment weighting (IPTW) using the propensity
#' score to estimate causal treatment effects in observational studies.
#' \emph{Statistics in Medicine}, \bold{34}(28), 3661-3679.
#' \doi{10.1002/sim.6607}.
#'
#' Hill, B.M. (1975). A simple general approach to inference about the tail of
#' a distribution. \emph{Annals of Statistics}, \bold{3}(5), 1163-1174.
#' \doi{10.1214/aos/1176343247}.
#' @seealso \code{\link{ipwplot}}, \code{\link{ipwpoint}}, \code{\link{ipwtm}},
#'   \code{\link{ipwvcov}}.
#' @export
#'
#' @examples
#' # Point treatment: binary exposure a, confounder l, outcome y
#' set.seed(16)
#' n <- 1000
#' simdat <- data.frame(l = rnorm(n, 10, 5))
#' simdat$a <- rbinom(n, 1, plogis(simdat$l - 10))
#' simdat$y <- 10 * simdat$a + 0.5 * simdat$l + rnorm(n, -10, 5)
#' w <- ipwpoint(exposure = a, family = "binomial", link = "logit",
#'               numerator = ~ 1, denominator = ~ l, data = simdat)
#' msm <- lm(y ~ a, data = simdat, weights = w$ipw.weights)
#' d <- ipwdiag(w, msm = msm)
#' d
#' # the numbers behind the plots
#' head(d$share$data)
#' d$trunc
#'
#' # Time-varying weights (haartdat): treatment times censoring weights
#' data(haartdat)
#' w.haart <- ipwtm(exposure = haartind, family = "binomial", link = "logit",
#'                  numerator = ~ sex + age, denominator = ~ sex + age + cd4.sqrt,
#'                  id = patient, timevar = fuptime, type = "first",
#'                  data = haartdat)
#' w.cens <- ipwtm(exposure = dropout, family = "binomial", link = "logit",
#'                 numerator = ~ sex + age, denominator = ~ sex + age + cd4.sqrt,
#'                 id = patient, timevar = fuptime, type = "cens",
#'                 data = haartdat)
#' ipwdiag(list(w.haart, w.cens), which = c("overlap", "share", "tail"))
ipwdiag <- function(
    object,
    which = c("overlap", "share", "tail", "trunc"),
    msm = NULL,
    trunc = c(0, 0.01, 0.025, 0.05, 0.1),
    coef = NULL,
    top = 0.05,
    ask = FALSE,
    ...)
{
  # 1. Checks ------------------------------------------------------------------
  explicit <- !missing(which)
  if (!is.character(which) || length(which) == 0) {
    stop("'which' must be one or more of \"overlap\", \"share\", \"tail\", \"trunc\".",
         call. = FALSE)
  }
  which <- match.arg(which, several.ok = TRUE)
  if (!is.numeric(top) || length(top) != 1 || is.na(top) || top <= 0 || top >= 0.5) {
    stop("'top' must be a single number between 0 and 0.5.", call. = FALSE)
  }
  if (!is.numeric(trunc) || length(trunc) == 0 || anyNA(trunc) ||
      any(trunc < 0) || any(trunc >= 0.5)) {
    stop("'trunc' must be a numeric vector with values in [0, 0.5).", call. = FALSE)
  }
  trunc <- sort(unique(trunc))
  if (!is.logical(ask) || length(ask) != 1 || is.na(ask)) {
    stop("'ask' must be TRUE or FALSE.", call. = FALSE)
  }
  wt <- .diag_weights(object)

  # 2. Computations ---------------------------------------------------------------
  res <- list()
  smry <- list(n = length(wt$final), unit = wt$unit)

  if ("overlap" %in% which) {
    ov <- .diag_overlap(wt$first, explicit)
    if (!is.null(ov)) {
      res$overlap <- ov
      smry$overlap <- ov$extreme
    }
  }
  if ("share" %in% which) {
    res$share <- .diag_share(wt$final)
    smry$share <- res$share$summary
    res$share$summary <- NULL
  }
  if ("tail" %in% which) {
    tl <- .diag_tail(wt$final, top)
    if (!is.null(tl)) {
      res$tail <- tl
      smry$tail <- c(alpha = tl$alpha, k = tl$k)
    }
  }
  if ("trunc" %in% which) {
    if (is.null(msm)) {
      if (explicit) {
        stop("which = \"trunc\" needs 'msm', the marginal structural model.", call. = FALSE)
      }
      message("ipwdiag(): no 'msm' given, so the truncation plot is skipped.")
    } else {
      res$trunc <- .diag_trunc(msm, wt, trunc, coef, parent.frame())
      smry$trunc <- res$trunc
    }
  }
  res$summary <- smry
  class(res) <- "ipwdiag"

  # 3. Plots ----------------------------------------------------------------------
  drawn <- intersect(c("overlap", "share", "tail", "trunc"), names(res))
  if (length(drawn) > 0) {
    new_par <- list(...)
    if (ask) {
      new_par$ask <- TRUE
    } else if (length(drawn) > 1 && prod(graphics::par("mfrow")) == 1) {
      new_par$mfrow <- if (length(drawn) == 2) c(1, 2) else c(2, 2)
      # tighter margins for small panels, unless given in ...
      if (is.null(new_par$mar)) new_par$mar <- c(3.6, 3.6, 2.6, 1) + 0.1
      if (is.null(new_par$mgp)) new_par$mgp <- c(2.3, 0.7, 0)
    }
    if (length(new_par) > 0) {
      old_par <- graphics::par(new_par)
      on.exit(graphics::par(old_par), add = TRUE)
    }
    for (p in drawn) {
      switch(p,
             overlap = .plot_overlap(res$overlap),
             share   = .plot_share(res$share, smry$share, wt$unit),
             tail    = .plot_tail(res$tail, top),
             trunc   = .plot_trunc(res$trunc))
    }
  }
  invisible(res)
}

#' @export
print.ipwdiag <- function(x, digits = 3, ...) {
  s <- x$summary
  cat("Weight diagnostics for", s$n, s$unit, "\n\n")
  if (!is.null(s$overlap)) {
    cat(sprintf("Overlap:    %.1f%% of %d %s have an estimated probability < 0.05 or > 0.95\n",
                100 * s$overlap[["proportion"]], s$overlap[["rows"]],
                if (s$overlap[["rows"]] == s$n) "people" else "rows"))
  }
  if (!is.null(s$share)) {
    sh <- s$share
    cat(sprintf("Weights:    mean %s, largest %s\n",
                format(sh[["mean"]], digits = digits), format(sh[["max"]], digits = digits)))
    cat(sprintf("Share:      top 1%% hold %.1f%%, top 5%% hold %.1f%% of the total weight\n",
                100 * sh[["top1"]], 100 * sh[["top5"]]))
    cat(sprintf("Kish ESS:   %.1f (%.1f%% of n)\n", sh[["ess"]], 100 * sh[["ess.fraction"]]))
  }
  if (!is.null(s$tail)) {
    cat(sprintf("Tail index: %s (Hill estimate from the largest %d weights)\n",
                format(s$tail[["alpha"]], digits = digits), as.integer(s$tail[["k"]])))
    if (s$tail[["alpha"]] < 2) {
      cat("            Warning: tail index < 2: the weights have an infinite variance;\n",
          "            no normal-theory interval can be trusted (see ?ipwvcov,\n",
          "            \"Heavy-tailed weights\").\n", sep = "")
    } else if (s$tail[["alpha"]] < 4) {
      cat("            Warning: tail index < 4: normal-theory intervals (ipwvcov() and\n",
          "            the bootstrap alike) can be too narrow (see ?ipwvcov,\n",
          "            \"Heavy-tailed weights\").\n", sep = "")
    }
  }
  if (!is.null(s$trunc)) {
    cat("\nTruncation sensitivity for '", attr(s$trunc, "coef"),
        "' (95% CI: robust, weights treated as known):\n", sep = "")
    tab <- s$trunc
    attr(tab, "coef") <- NULL
    print(tab, digits = digits, row.names = FALSE)
  }
  invisible(x)
}

# ---- Weights -----------------------------------------------------------------

# The row weights (product for a list), each person's final weight, the
# object used for the overlap plot, and for long data the id of each row.
.diag_weights <- function(object) {
  if (inherits(object, c("ipwpoint", "ipwtm"))) {
    objects <- list(object)
  } else if (.is_weight_list(object)) {
    objects <- object
  } else {
    stop("'object' must be the result of ipwpoint() or ipwtm(), or a list of such ",
         "results whose weights multiply.", call. = FALSE)
  }
  for (o in objects) {
    w <- o$ipw.weights
    if (!is.numeric(w) || length(w) == 0 || anyNA(w) || any(!is.finite(w))) {
      stop("'object' has no valid weights (ipw.weights must be finite numbers).",
           call. = FALSE)
    }
    if (any(w <= 0)) stop("All weights must be positive.", call. = FALSE)
  }
  N <- length(objects[[1]]$ipw.weights)
  if (any(vapply(objects, function(o) length(o$ipw.weights), integer(1)) != N)) {
    stop("All weight objects in 'object' must have the same number of weights ",
         "(one per row of the same data).", call. = FALSE)
  }
  w <- Reduce(`*`, lapply(objects, function(o) o$ipw.weights))

  is_tm <- vapply(objects, inherits, logical(1), what = "ipwtm")
  id <- NULL
  if (any(is_tm)) {
    tm <- objects[[which(is_tm)[1]]]
    id <- attr(tm, "id")
    ord <- attr(tm, "sort.order")
    if (is.null(id) || is.null(ord)) {
      stop("An ipwtm() result has no 'id' or 'sort.order' attribute; re-run ipwtm() ",
           "with this version of ipw.", call. = FALSE)
    }
    for (o in objects[is_tm]) {
      if (!identical(attr(o, "id"), id)) {
        stop("The ipwtm() results in 'object' must come from the same long data ",
             "(same rows, same order, same id).", call. = FALSE)
      }
    }
    # last row of each person in time order
    last <- ord[!duplicated(id[ord], fromLast = TRUE)]
    final <- w[last]
    unit <- "people (final cumulative weight of each person)"
  } else {
    final <- w
    unit <- "people"
  }
  list(rows = w, final = final, id = id, first = objects[[1]], unit = unit)
}

# ---- Overlap -----------------------------------------------------------------

.diag_overlap <- function(obj, explicit) {
  mod <- obj$den.mod
  fam <- if (is.list(mod) && !is.null(mod$family)) mod$family$family else NULL
  if (!identical(fam, "binomial")) {
    msg <- "The overlap plot needs a binomial exposure model (family = \"binomial\")."
    if (explicit) stop(msg, call. = FALSE)
    message("ipwdiag(): ", msg, " It is skipped.")
    return(NULL)
  }
  p <- as.numeric(mod$fitted.values)
  a <- as.numeric(mod$y)
  exp_name <- attr(obj, "exposure")
  if (is.null(exp_name)) exp_name <- deparse(stats::formula(mod)[[2]])
  # plotting range: the range of the probabilities, a little wider, within [0, 1]
  pad <- 0.02 * max(diff(range(p)), 0.05)
  xlim <- c(max(0, min(p) - pad), min(1, max(p) + pad))
  dens <- lapply(c(0, 1), function(g) {
    pg <- p[a == g]
    if (length(pg) < 2) return(data.frame(x = numeric(0), y = numeric(0)))
    d <- stats::density(pg, from = xlim[1], to = xlim[2])
    data.frame(x = d$x, y = d$y)
  })
  names(dens) <- c("0", "1")
  support <- c(lower = max(min(p[a == 0]), min(p[a == 1])),
               upper = min(max(p[a == 0]), max(p[a == 1])))
  list(data = data.frame(p = p, exposure = a), density = dens, support = support,
       exposure = exp_name, xlim = xlim,
       extreme = c(proportion = mean(p < 0.05 | p > 0.95), rows = length(p)))
}

.plot_overlap <- function(ov) {
  dens <- ov$density
  ymax <- max(c(dens[["0"]]$y, dens[["1"]]$y, 0.1))
  xl <- ov$xlim
  graphics::plot(NA, xlim = xl, ylim = c(0, ymax * 1.08), xaxs = "i",
                 xlab = sprintf("Estimated P(%s = 1)", ov$exposure), ylab = "Density",
                 main = "Overlap")
  usr <- graphics::par("usr")
  sup <- ov$support
  shade <- "#00000022"
  if (sup[["lower"]] < sup[["upper"]]) {
    graphics::rect(xl[1], usr[3], sup[["lower"]], usr[4], col = shade, border = NA)
    graphics::rect(sup[["upper"]], usr[3], xl[2], usr[4], col = shade, border = NA)
  } else {
    graphics::rect(xl[1], usr[3], xl[2], usr[4], col = shade, border = NA)
  }
  cols <- c("#1F6FB4", "#C8501E")
  graphics::lines(dens[["0"]]$x, dens[["0"]]$y, col = cols[1], lwd = 2)
  graphics::lines(dens[["1"]]$x, dens[["1"]]$y, col = cols[2], lwd = 2, lty = 2)
  graphics::legend("top", legend = paste(ov$exposure, "=", 0:1), col = cols, lwd = 2,
                   lty = 1:2, bty = "n", cex = 0.85, horiz = TRUE)
  graphics::mtext("shaded: values in one group only", side = 3, line = 0.2,
                  cex = 0.8 * graphics::par("cex"))
}

# ---- Share of the total weight ----------------------------------------------

.diag_share <- function(w) {
  n <- length(w)
  ws <- sort(w, decreasing = TRUE)
  cum <- cumsum(ws) / sum(ws)
  top_share <- function(f) cum[ceiling(f * n)]
  list(data = data.frame(fraction = c(0, seq_len(n) / n), share = c(0, cum)),
       weights = ws,
       summary = c(top1 = top_share(0.01), top5 = top_share(0.05), max = ws[1],
                   mean = mean(w), ess = sum(w)^2 / sum(w^2),
                   ess.fraction = sum(w)^2 / sum(w^2) / n))
}

.plot_share <- function(sh, s, unit) {
  d <- sh$data
  graphics::plot(d$fraction, d$share, type = "l", lwd = 2, col = "#1F6FB4",
                 xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i",
                 xlab = "Fraction of people (largest weights first)",
                 ylab = "Share of the total weight", main = "Weight share")
  graphics::abline(0, 1, lty = 3, col = "grey40")
  graphics::abline(v = c(0.01, 0.05), lty = 3, col = "grey70")
  .inner_text(c(sprintf("top 1%%: %.1f%%, top 5%%: %.1f%%",
                         100 * s[["top1"]], 100 * s[["top5"]]),
                 sprintf("Kish ESS %.0f (%.0f%% of n)", s[["ess"]],
                         100 * s[["ess.fraction"]]),
                 "dotted line: equal weights"), side = 1)
}

# Lines of text inside the plot region, right-aligned at the top (side = 3) or
# bottom (side = 1). mtext() positions them in lines of text, so they keep
# their spacing when a recorded plot is replayed on a device of another size.
.inner_text <- function(txt, side, cex = 0.8) {
  k <- length(txt)
  pos <- if (side == 3) seq_len(k) else rev(seq_len(k))
  for (i in seq_len(k)) {
    graphics::mtext(txt[i], side = side, line = -1.1 * pos[i] - 0.2, adj = 0.97,
                    cex = cex * graphics::par("cex"))
  }
}

# ---- Tail index -----------------------------------------------------------------

.hill <- function(w, k) {
  ws <- sort(w, decreasing = TRUE)
  1 / mean(log(ws[seq_len(k)] / ws[k + 1]))
}

.diag_tail <- function(w, top) {
  n <- length(w)
  k <- floor(top * n)
  if (k < 20) {
    message("ipwdiag(): the tail plot needs at least 20 weights in the top ",
            100 * top, "% (there are ", k, "); it is skipped. Increase 'top' or use ",
            "more data.")
    return(NULL)
  }
  ws <- sort(w, decreasing = TRUE)
  thr <- ws[k + 1]
  if (thr == ws[1]) {
    message("ipwdiag(): the largest weights are all equal; the tail plot is skipped.")
    return(NULL)
  }
  alpha <- .hill(w, k)
  # Pareto tail above the threshold: P(W >= x) = (k / n) (x / thr)^(-alpha)
  list(data = data.frame(log.weight = log(ws[seq_len(k)]),
                         log.rank = log(seq_len(k) / n)),
       alpha = alpha, k = k, threshold = thr,
       intercept = log(k / n) + alpha * log(thr), slope = -alpha)
}

.plot_tail <- function(tl, top) {
  d <- tl$data
  graphics::plot(d$log.weight, d$log.rank, pch = 20, col = "#1F6FB4",
                 xlab = "log(weight)", ylab = "log(rank / n)",
                 main = sprintf("Tail of the weights (top %g%%)", 100 * top))
  graphics::abline(tl$intercept, tl$slope, lty = 2, lwd = 1.5)
  .inner_text(c(sprintf("Hill tail index %.2f (dashed line)", tl$alpha),
                 if (tl$alpha < 4) "< 4: intervals may be too narrow"), side = 3)
}

# ---- Truncation sensitivity ---------------------------------------------------

.diag_trunc <- function(msm, wt, trunc, coef, caller) {
  if (!(identical(class(msm), "lm") || identical(class(msm), c("glm", "lm")) ||
        identical(class(msm), "coxph"))) {
    stop("'msm' must be a model fitted with lm(), glm() or survival::coxph() ",
         "(svyglm(), glm.nb(), geeglm() and penalized Cox models are not supported; ",
         "fit the model with weights = <the weights>).", call. = FALSE)
  }
  is_cox <- inherits(msm, "coxph")
  N <- length(wt$rows)
  n_msm <- if (is_cox) msm$n else length(msm$residuals)
  if (n_msm != N) {
    stop("'msm' was fitted on ", n_msm, " rows but there are ", N, " weights. Fit it ",
         "on the same rows, in the same order, as the weights.", call. = FALSE)
  }
  cf <- stats::coef(msm)
  if (is.null(coef)) {
    nm <- if (is_cox) names(cf) else setdiff(names(cf), "(Intercept)")
    if (length(nm) == 0) stop("'msm' has no coefficient other than the intercept.", call. = FALSE)
    coef <- nm[1]
  } else if (is.numeric(coef) && length(coef) == 1 && coef >= 1 && coef <= length(cf)) {
    coef <- names(cf)[coef]
  } else if (!(is.character(coef) && length(coef) == 1 && coef %in% names(cf))) {
    stop("'coef' must be the name or position of one coefficient of 'msm': ",
         paste(names(cf), collapse = ", "), ".", call. = FALSE)
  }
  msm_w <- if (inherits(msm, "glm")) msm$prior.weights else msm$weights
  if (is.null(msm_w) || !isTRUE(all.equal(as.numeric(msm_w), as.numeric(wt$rows)))) {
    warning("The weights of 'msm' differ from the weights in 'object' (or 'msm' has ",
            "no weights); the truncation plot refits 'msm' with (truncated) weights ",
            "from 'object'.", call. = FALSE)
  }

  # How the robust variance of a Cox model is obtained
  extra <- list()
  if (is_cox) {
    has_cluster <- !is.null(attr(stats::terms(msm), "specials")$cluster) ||
      !is.null(msm$call$cluster)
    if (!has_cluster) {
      if (!is.null(wt$id)) extra$cluster <- wt$id else extra$robust <- TRUE
    }
  }

  rows <- lapply(trunc, function(q) {
    wq <- .truncate_weights(wt$rows, q)
    fit <- .refit_msm(msm, c(list(weights = wq), extra), caller)
    if (is_cox) {
      v <- fit$var
      dimnames(v) <- list(names(stats::coef(fit)), names(stats::coef(fit)))
    } else {
      v <- .sandwich_glm(fit, wq, wt$id)
    }
    est <- unname(stats::coef(fit)[coef])
    se <- sqrt(v[coef, coef])
    data.frame(q = q, estimate = est, se = se,
               lower = est - stats::qnorm(0.975) * se, upper = est + stats::qnorm(0.975) * se)
  })
  tab <- do.call(rbind, rows)
  attr(tab, "coef") <- coef
  tab
}

# Refit 'msm' with new arguments (inlined as values), evaluating the call where
# ipwdiag() was called or, failing that, where the model formula was created.
.refit_msm <- function(msm, args, caller) {
  cl <- msm$call
  if (is.null(cl)) stop("'msm' has no call to refit.", call. = FALSE)
  for (nm in names(args)) cl[[nm]] <- args[[nm]]
  envs <- list(caller, environment(stats::formula(msm)))
  err <- NULL
  for (e in envs) {
    if (is.null(e)) next
    fit <- tryCatch(eval(cl, e), error = function(x) { err <<- conditionMessage(x); NULL })
    if (!is.null(fit)) break
  }
  if (is.null(fit)) {
    stop("Could not refit 'msm' with truncated weights: ", err, "\nipwdiag() refits ",
         "the model with stats::update(); its data must be available where ipwdiag() ",
         "is called (fit 'msm' with data = <a data frame>).", call. = FALSE)
  }
  n_old <- if (inherits(msm, "coxph")) msm$n else length(msm$residuals)
  n_new <- if (inherits(fit, "coxph")) fit$n else length(fit$residuals)
  if (n_new != n_old) {
    stop("Refitting 'msm' with truncated weights used ", n_new, " rows instead of ",
         n_old, ".", call. = FALSE)
  }
  fit
}

# Sandwich variance (HC0) of an lm() or glm() fitted with weights w, with the
# factor g / (g - 1) that survey::svyglm() uses; for long data the scores are
# summed within person (g persons), otherwise g is the number of rows.
.sandwich_glm <- function(fit, w, id = NULL) {
  X <- stats::model.matrix(fit)
  if (inherits(fit, "glm")) {
    U <- X * (fit$weights * fit$residuals)   # working weights x working residuals
    A <- crossprod(X * fit$weights, X)
  } else {
    U <- X * (w * fit$residuals)
    A <- crossprod(X * w, X)
  }
  if (!is.null(id)) U <- rowsum(U, id, reorder = FALSE)
  g <- nrow(U)
  A_inv <- solve(A)
  V <- A_inv %*% crossprod(U) %*% A_inv * g / (g - 1)
  dimnames(V) <- list(colnames(X), colnames(X))
  V
}

.plot_trunc <- function(tab) {
  x <- seq_len(nrow(tab))
  ylim <- range(c(tab$lower, tab$upper))
  graphics::plot(x, tab$estimate, pch = 19, col = "#1F6FB4", ylim = ylim,
                 xlim = c(0.6, length(x) + 0.4), xaxt = "n",
                 xlab = "Truncation level q",
                 ylab = paste("Estimate of", attr(tab, "coef")), main = "Truncation sensitivity")
  graphics::axis(1, at = x, labels = format(tab$q, drop0trailing = TRUE))
  graphics::arrows(x, tab$lower, x, tab$upper, angle = 90, code = 3,
                   length = 0.04, col = "#1F6FB4")
  graphics::abline(h = tab$estimate[1], lty = 3, col = "grey40")
  graphics::mtext("95% CI: robust, weights treated as known", side = 3, line = 0.2,
                  cex = 0.8 * graphics::par("cex"))
}
