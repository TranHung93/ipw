# ==========================================================================
# ipwvcov() for time-varying weights from ipwtm() (stage B)
#
# Called by ipwvcov() when 'object' is an ipwtm() result or a list of
# ipwtm()/ipwpoint() results whose weights multiply to the MSM's weights.
# The unit of independence is the person (the ipwtm() id): all estimating
# equations are summed within person before the meat is formed.
# ==========================================================================

# TRUE for a list of ipwtm()/ipwpoint() results (not itself one of them)
.is_weight_list <- function(object) {
  is.list(object) && !inherits(object, c("ipwtm", "ipwpoint")) && length(object) > 0 &&
    all(vapply(object, function(o) inherits(o, c("ipwtm", "ipwpoint")), logical(1)))
}

.ipwvcov_tm <- function(objects, msm) {
  if (inherits(objects, "ipwtm")) objects <- list(objects)
  is_tm <- vapply(objects, inherits, logical(1), what = "ipwtm")
  if (!any(is_tm)) {
    stop("A list for 'object' needs at least one ipwtm() result; for ipwpoint() ",
         "weights alone, pass the ipwpoint() result itself.", call. = FALSE)
  }
  for (o in objects[is_tm]) .check_tm_object(o)
  for (o in objects[!is_tm]) {
    if (is.null(attr(o, "estimand"))) {
      stop("An ipwpoint() result has no 'estimand' attribute; re-run ipwpoint() ",
           "with this version of ipw.", call. = FALSE)
    }
    .check_vcov_model(o$den.mod, "denominator")
    if (!is.null(o$num.mod)) .check_vcov_model(o$num.mod, "numerator")
  }

  # Persons: the id of the ipwtm() data, the same for every ipwtm() object
  id <- attr(objects[[which(is_tm)[1]]], "id")
  for (o in objects[is_tm]) {
    if (!identical(attr(o, "id"), id)) {
      stop("The ipwtm() results in 'object' must come from the same long data ",
           "(same rows, same order, same id).", call. = FALSE)
    }
  }
  N <- length(id)
  for (o in objects) {
    if (length(o$ipw.weights) != N) {
      stop("All weight objects in 'object' must have one weight per row of the ",
           "long data (", N, " rows).", call. = FALSE)
    }
  }

  # MSM checks
  is_cox <- !missing(msm) && inherits(msm, "coxph")
  if (missing(msm) || !(is_cox || inherits(msm, "lm")) || inherits(msm, c("mlm", "gam"))) {
    stop("'msm' must be a model fitted with lm() or glm(), or with coxph() for ",
         "counting-process data Surv(tstart, time, event).", call. = FALSE)
  }
  if (anyNA(stats::coef(msm))) {
    stop("'msm' has aliased (NA) coefficients; remove them from the model.", call. = FALSE)
  }
  if (is_cox) .check_coxph_msm_tm(msm)
  X <- stats::model.matrix(msm)
  if (nrow(X) != N) {
    stop("'msm' was fitted on ", nrow(X), " rows but the long data have ", N,
         " rows. Fit it on all rows of the long data, in the same order.", call. = FALSE)
  }
  w <- Reduce(`*`, lapply(objects, function(o) o$ipw.weights))
  msm_w <- if (inherits(msm, "glm")) msm$prior.weights else msm$weights
  if (is.null(msm_w)) msm_w <- rep(1, N)
  msm_w <- as.numeric(msm_w)
  if (!isTRUE(all.equal(msm_w, as.numeric(w)))) {
    has_trunc <- vapply(objects, function(o) !is.null(o$weights.trunc), logical(1))
    if (any(has_trunc)) {
      w_tr <- Reduce(`*`, lapply(objects, function(o)
        if (!is.null(o$weights.trunc)) o$weights.trunc else o$ipw.weights))
      if (isTRUE(all.equal(msm_w, as.numeric(w_tr)))) {
        stop("Truncated weights are not supported: truncation at estimated ",
             "percentiles is not a smooth function of the model coefficients. Use ",
             "the cluster bootstrap (vignette(\"ipw\"), section 5).", call. = FALSE)
      }
    }
    stop("The weights of 'msm' must equal the product of the ipw.weights of the ",
         "objects in 'object', on the same rows in the same order.", call. = FALSE)
  }

  # MSM scores (rows) and information
  if (is_cox) {
    cox <- .coxph_parts_counting(msm, w, X)
    U <- cox$U
    A_bb <- cox$A
  } else {
    if (inherits(msm, "glm")) {
      mu <- msm$fitted.values
      resid <- msm$y - mu
      mu_eta <- msm$family$mu.eta(msm$linear.predictors)
      var_mu <- msm$family$variance(mu)
    } else {
      resid <- msm$residuals
      mu_eta <- rep(1, N)
      var_mu <- rep(1, N)
    }
    U <- X * (w * resid * mu_eta / var_mu)
    A_bb <- crossprod(X * (w * mu_eta^2 / var_mu), X)
  }

  # Exposure (and censoring) models: one parameter block per model
  blocks <- list()
  for (o in objects) {
    blocks <- c(blocks, if (inherits(o, "ipwtm")) .tm_blocks(o) else .point_blocks(o))
  }

  # Stacked sandwich, scores summed within person
  k <- ncol(X)
  Psi <- do.call(cbind, c(list(U), lapply(blocks, `[[`, "S")))
  Psi_person <- rowsum(Psi, id, reorder = FALSE)
  sizes <- vapply(blocks, function(b) ncol(b$A), integer(1))
  A <- matrix(0, ncol(Psi), ncol(Psi))
  A[seq_len(k), seq_len(k)] <- A_bb
  start <- k
  for (j in seq_along(blocks)) {
    idx <- start + seq_len(sizes[j])
    A[idx, idx] <- blocks[[j]]$A
    A[seq_len(k), idx] <- -crossprod(U, blocks[[j]]$D)
    start <- start + sizes[j]
  }
  A_inv <- solve(A)
  V <- A_inv %*% crossprod(Psi_person) %*% t(A_inv)
  V <- V[seq_len(k), seq_len(k), drop = FALSE]
  dimnames(V) <- list(names(stats::coef(msm)), names(stats::coef(msm)))
  V
}

# An ipwtm() result that ipwvcov() can handle: binomial family, logistic
# glm models (no GAM), any type, made with this version of ipw.
.check_tm_object <- function(o) {
  if (is.null(attr(o, "sort.order")) || is.null(attr(o, "type"))) {
    stop("An ipwtm() result has no 'sort.order' attribute; re-run ipwtm() with ",
         "this version of ipw.", call. = FALSE)
  }
  if (!identical(attr(o, "family"), "binomial")) {
    stop("ipwvcov() supports ipwtm() weights only for family = 'binomial' (it has '",
         attr(o, "family"), "'). For Cox (survival), multinomial, ordinal and ",
         "gaussian exposure models, use the cluster bootstrap that re-estimates the ",
         "weights (vignette(\"ipw\"), section 5).", call. = FALSE)
  }
  for (which in c("den.mod", "num.mod")) {
    mod <- o[[which]]
    if (is.null(mod)) next
    lab <- if (which == "den.mod") "denominator" else "numerator"
    if (inherits(mod, "gam")) {
      stop("ipwvcov() does not support smooth terms (GAM) in the ", lab, " model ",
           "of ipwtm(). Use the cluster bootstrap (vignette(\"ipw\"), section 5).",
           call. = FALSE)
    }
    .check_vcov_model(mod, lab)
  }
  invisible(TRUE)
}

# Parameter blocks of one ipwtm() result, in the row order of the data.
# Row factor f_s of the weight (sorted rows, selvar = 1 for type "first" and
# "cens", all rows for "all"): P(observed exposure) for "first" and "all",
# 1 - p for "cens" (the censoring row included), and 1 on other rows. The
# weight is prod(num factors) / prod(den factors) up to each row, so
#   d log w_r / d alpha = - cumsum_{s <= r, same id} d log f_s(den) / d alpha,
# with d log f_s / d alpha = (a_s - p_s) Z_s ("first", "all") or
# -p_s Z_s ("cens") for a logistic model (likewise + for the numerator).
.tm_blocks <- function(o) {
  ord <- attr(o, "sort.order")
  type <- attr(o, "type")
  N <- length(ord)
  id_sorted <- attr(o, "id")[ord]
  sel_sorted <- o$selvar[ord]
  fit_rows <- if (type == "all") seq_len(N) else which(sel_sorted == 1)

  one <- function(mod, sign) {
    Z <- stats::model.matrix(mod)
    if (nrow(Z) != length(fit_rows)) {
      stop("Internal error: the ipwtm() model rows do not match its selection ",
           "indicator.", call. = FALSE)
    }
    p <- as.numeric(mod$fitted.values)
    a <- as.numeric(mod$y)
    kz <- ncol(Z)
    S_sorted <- matrix(0, N, kz)
    S_sorted[fit_rows, ] <- Z * (a - p)
    G <- matrix(0, N, kz)
    G[fit_rows, ] <- Z * (if (type == "cens") -p else a - p)
    D_sorted <- sign * apply(G, 2, function(g) stats::ave(g, id_sorted, FUN = cumsum))
    D_sorted <- matrix(D_sorted, nrow = N)
    S <- D <- matrix(0, N, kz)
    S[ord, ] <- S_sorted
    D[ord, ] <- D_sorted
    list(S = S, D = D, A = crossprod(Z * (p * (1 - p)), Z))
  }
  out <- list(one(o$den.mod, -1))
  if (!is.null(o$num.mod)) out <- c(out, list(one(o$num.mod, +1)))
  out
}

# Parameter blocks of an ipwpoint() result used as one factor of the weights
# (one weight per row of the long data), as in ipwvcov() for ipwpoint().
.point_blocks <- function(o) {
  estimand <- attr(o, "estimand")
  den <- .binomial_parts(o$den.mod)
  a <- den$a
  p <- den$p
  g <- switch(estimand,
    ATE = ifelse(a == 1, -1 / p, 1 / (1 - p)),
    ATT = ifelse(a == 1, 0, 1 / (p * (1 - p))),
    ATC = ifelse(a == 1, -1 / (p * (1 - p)), 0),
    ATO = ifelse(a == 1, -1 / (1 - p), 1 / p)
  )
  out <- list(list(S = den$S, D = den$Z * (g * den$dp), A = den$A))
  if (!is.null(o$num.mod)) {
    num <- .binomial_parts(o$num.mod)
    q <- num$p
    out <- c(out, list(list(S = num$S, A = num$A,
                            D = num$Z * (ifelse(a == 1, 1 / q, -1 / (1 - q)) * num$dp))))
  }
  out
}

# A coxph() MSM for long data: counting-process Surv(tstart, time, event); no
# strata, frailties or time-transforms; cluster() is allowed (it does not
# change the coefficients; ipwvcov() always clusters on the ipwtm() id).
.check_coxph_msm_tm <- function(msm) {
  if (!identical(class(msm), "coxph")) {
    stop("ipwvcov() supports coxph() models without penalized terms (frailty(), ",
         "pspline()) and with at least one covariate.", call. = FALSE)
  }
  if (!identical(attr(msm$y, "type"), "counting")) {
    stop("With ipwtm() weights, the coxph() model must use counting-process data, ",
         "Surv(tstart, time, event), on the rows of the long data.", call. = FALSE)
  }
  specials <- attr(stats::terms(msm), "specials")
  used <- names(Filter(Negate(is.null), specials[c("strata", "frailty", "tt")]))
  if (length(used) > 0 || !is.null(msm$strata) || !is.null(msm$call$tt)) {
    stop("ipwvcov() does not support strata(), frailty() or tt() terms in a ",
         "coxph() model.", call. = FALSE)
  }
  event_times <- msm$y[msm$y[, 3] == 1, 2]
  if (anyDuplicated(event_times) && !identical(msm$method, "breslow")) {
    stop("The coxph() model has tied event times and uses ties = '", msm$method,
         "'. ipwvcov() needs ties = \"breslow\" when event times are tied.",
         call. = FALSE)
  }
  invisible(TRUE)
}

# Weighted Cox partial likelihood (Breslow) for counting-process rows. Row j
# is at risk at event time u when tstart_j < u <= time_j. With
# r_j = w_j exp(eta_j), S0(u) = sum_{at risk} r_j, S1(u) = sum r_j X_j,
# Xbar = S1 / S0, H0(t) = sum_{event times u <= t} dW(u) / S0(u) and H1 the
# same with dW(u) Xbar(u) / S0(u), row j's score is w_j times
#   res_j = status_j (X_j - Xbar(time_j))
#           - exp(eta_j) [X_j (H0(time_j) - H0(tstart_j)) - (H1(time_j) - H1(tstart_j))].
.coxph_parts_counting <- function(msm, w, X) {
  start <- msm$y[, 1]
  stop <- msm$y[, 2]
  status <- msm$y[, 3]
  X <- matrix(X, nrow = nrow(X))
  k <- ncol(X)
  eta <- msm$linear.predictors
  r <- w * exp(eta)
  ut <- sort(unique(stop[status == 1]))
  m <- length(ut)
  a_idx <- rep(seq_len(k), times = k)
  b_idx <- rep(seq_len(k), each = k)

  # sum of the rows of M over observations with times >= u, for each u in ut
  sum_ge <- function(times, M) {
    M <- as.matrix(M)
    o <- order(times)
    suffix <- rbind(apply(M[o, , drop = FALSE][rev(seq_along(o)), , drop = FALSE], 2, cumsum)[
      rev(seq_along(o)), , drop = FALSE], 0)
    pos <- findInterval(ut, times[o], left.open = TRUE)   # number of times < u
    suffix[pos + 1, , drop = FALSE]
  }
  at_risk_sum <- function(M) sum_ge(stop, M) - sum_ge(start, M)

  S0 <- at_risk_sum(r)[, 1]
  S1 <- matrix(at_risk_sum(X * r), nrow = m)
  S2 <- matrix(at_risk_sum(X[, a_idx, drop = FALSE] * X[, b_idx, drop = FALSE] * r), nrow = m)
  dW <- as.numeric(tapply(w * status, factor(stop, levels = ut), sum, default = 0))
  dW[is.na(dW)] <- 0
  Xbar <- S1 / S0

  info <- matrix(colSums(dW * (S2 / S0 - Xbar[, a_idx, drop = FALSE] * Xbar[, b_idx, drop = FALSE])), k, k)
  H0 <- c(0, cumsum(dW / S0))
  H1 <- rbind(0, matrix(apply(dW * Xbar / S0, 2, cumsum), nrow = m))
  i_stop <- findInterval(stop, ut) + 1      # event times <= time_j, +1 for the leading 0
  i_start <- findInterval(start, ut) + 1    # event times <= tstart_j
  ev <- match(stop, ut)
  Xbar_own <- matrix(0, length(stop), k)
  Xbar_own[status == 1, ] <- Xbar[ev[status == 1], , drop = FALSE]
  res <- status * (X - Xbar_own) -
    exp(eta) * (X * (H0[i_stop] - H0[i_start]) - (H1[i_stop, , drop = FALSE] - H1[i_start, , drop = FALSE]))
  list(U = res * w, A = info, res = res)
}
