#' GMM estimation of treatment effects under staggered adoption
#'
#' Estimates cohort-by-time average treatment effects (CATTs) under staggered
#' treatment adoption by the generalized method of moments, following Arora
#' and Bijani (2026). Only clean two-by-two difference-in-differences
#' comparisons -- against never-treated or not-yet-treated control cohorts --
#' enter the moment system. Comparisons against already-treated controls are
#' exact linear combinations of the clean ones and are therefore never
#' formed.
#'
#' @param data A data frame with one row per unit-period.
#' @param yname Character. Name of the outcome column. Must be numeric;
#'   factor, character and logical outcomes are rejected rather than coerced.
#' @param tname Character. Name of the calendar-time column. Must be
#'   integer-valued, and the set of periods present in `data` (pooled across
#'   units) must be consecutive.
#' @param idname Character. Name of the unit identifier column.
#' @param gname Character. Name of the cohort column, giving each unit's
#'   first treated period, with `0` for units that are never treated. Must be
#'   constant within unit.
#' @param weighting Character. The covariance model used to form the optimal
#'   GMM weight. One of `"pooled_toeplitz"` (the default: a single stationary
#'   autocovariance sequence shared by all groups), `"cohort_toeplitz"` (a
#'   separate stationary sequence per cohort), or `"unrestricted"` (a full,
#'   unrestricted within-cohort covariance). These correspond to GMM-T,
#'   GMM-HT and GMM-U in the paper.
#' @param never_treated Logical or `NULL`. Whether to use never-treated units
#'   as controls in addition to not-yet-treated cohorts. `NULL` (the default)
#'   detects them automatically. `TRUE` requires that at least one exists and
#'   errors otherwise. `FALSE` excludes them, and warns if any are present.
#' @param covar Optional character vector of baseline (pre-treatment,
#'   time-invariant) covariate columns. When supplied, the outcome-regression
#'   adjustment of Section 4.5 of the paper is applied, so that parallel
#'   trends need hold only conditional on `covar`.
#' @param max_iter Integer. Maximum number of iterated-GMM reweighting steps.
#'   The default is generous because the iteration can converge linearly at a
#'   slow rate: on `sim_panel`, `weighting = "cohort_toeplitz"` needs 86
#'   steps. Iterating to convergence is not required for asymptotic
#'   efficiency, since the two-step and iterated estimators are
#'   asymptotically equivalent, so a small value is a legitimate choice when
#'   speed matters; the simulations in Arora and Bijani (2026) used 10.
#' @param tol Numeric. Convergence tolerance on the largest absolute change
#'   in any CATT estimate between successive iterations. Note that the
#'   aggregate ATT typically stabilises several decimal places earlier than
#'   the individual effects do.
#'
#' @section Data requirements:
#' The following conditions are checked before any estimation, and each
#' raises an error naming the offending rows, units, periods or cohorts:
#'
#' * `data` is a data frame and every named column exists.
#' * The outcome is numeric; time and cohort are integer-valued.
#' * No missing unit id, period or cohort. Missing outcomes are permitted.
#' * Each unit-period combination appears exactly once.
#' * The cohort is constant within unit.
#' * No period is absent from the whole panel. Individual units may be
#'   missing individual periods.
#' * No unit is treated at or before the first observed period. Such units
#'   have no pre-treatment period, so none of their effects is identified,
#'   and silently including them contaminates the aggregate.
#' * At least one cohort adopts strictly inside the observed window.
#'
#' Units whose cohort falls after the last observed period are never treated
#' within the window and are therefore pooled with the never-treated group as
#' clean controls, rather than discarded.
#'
#' @section Unbalanced panels:
#' Missing outcomes are supported. Group-period means are taken over the
#' units observed at that period; autocovariances are divided by the number
#' of residual pairs actually observed at each lag; and the variance of a
#' group-period mean accounts for the overlap between the unit sets
#' contributing at each pair of periods. When the panel is balanced all three
#' reduce exactly to the complete-data formulas. Validity requires that
#' missingness be independent of the outcome given group and period;
#' selection on outcomes is not addressed by any weighting scheme.
#'
#' @section Partial identification:
#' A cohort-by-time effect with no clean comparison is not estimable. Such
#' cells are reported with `estimate = NA` and `identified = FALSE`, never as
#' a precisely-estimated zero. Where an unidentified cell carries positive
#' aggregation weight the corresponding aggregate is `NA`, and an
#' identified-subset aggregate renormalised over the estimable cells is
#' reported alongside it, with the share of weight it covers.
#'
#' @section Convergence:
#' Three states are reported separately. `converged` is `TRUE` only when the
#' largest change in any CATT fell below `tol`. `solve_ok` records the weaker
#' fact that at least one weighted solve succeeded. `termination` gives the
#' reason the iteration stopped.
#'
#' When `solve_ok` is `FALSE`, no weighted step was ever completed and the
#' returned estimates are the identity-weighted seed rather than a GMM
#' estimate under the requested weighting. That always raises a warning.
#'
#' The requirement is one of rank. Each cohort's estimated covariance is
#' built from that cohort's residual vectors, so its rank cannot exceed
#' \eqn{\min(N_g, T-1)}. Under `weighting = "unrestricted"` the weight's
#' total rank is therefore bounded by \eqn{\sum_g \min(N_g, T-1)}, and the
#' normal equations are singular unless that comfortably exceeds the number
#' of identified effects. The bound is necessary rather than sufficient:
#' equality leaves no margin and fails in practice. Where it binds, the
#' warning reports the number of moment directions the weight actually
#' supported.
#'
#' When `solve_ok` is `TRUE` but `converged` is `FALSE`, the iteration ran
#' and did not settle; the estimates are the last valid iterate. That warns
#' for the two Toeplitz weightings but not for `"unrestricted"`, where slow
#' convergence is routine.
#'
#' @section Covariates and the weighting matrix:
#' The covariate adjustment enters the moment vector only. The residuals from
#' which the covariance model is estimated are taken after sweeping out unit
#' and period fixed effects, which absorb any time-invariant covariate
#' exactly. Where a covariate has a time-varying loading -- the case in which
#' the outcome-regression adjustment does anything at all -- that component
#' remains in the residuals and inflates the estimated persistence. Point
#' estimates are unaffected, since the weight governs efficiency only, but
#' the cost is visible in the iteration: on `sim_panel`, whose covariates
#' carry time-varying loadings, `weighting = "cohort_toeplitz"` needs 86
#' steps to converge, against 11 on the same design with those terms removed.
#'
#' @return An object of class `staggered_gmm`: a list with components
#'   `call`, `weighting`, `coefficients`, `vcov`, `catt`, `aggregate`,
#'   `weights`, `convergence`, `design`, `controls` and `validation`. See
#'   [summary.staggered_gmm()] and `vignette("staggeredGMM")`.
#'
#' @references
#' Arora, P. and Bijani, R. (2026). "Estimating Treatment Effects under
#' Staggered Timing and Non-Spherical Errors." \doi{10.2139/ssrn.6558759}
#'
#' @seealso [gmm_j_test()] for the over-identification test of parallel
#'   trends and no anticipation.
#'
#' @examples
#' fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
#'                      idname = "unit_id", gname = "cohort")
#' fit
#'
#' @export
gmm_staggered <- function(data, yname, tname, idname, gname,
                          weighting = c("pooled_toeplitz",
                                        "cohort_toeplitz",
                                        "unrestricted"),
                          never_treated = NULL,
                          covar    = NULL,
                          max_iter = 100L,
                          tol      = 1e-6) {

  cl        <- match.call()
  weighting <- match.arg(weighting)

  if (!is.numeric(max_iter) || length(max_iter) != 1L || is.na(max_iter) ||
      max_iter < 1 || max_iter != trunc(max_iter)) {
    stop("`max_iter` must be a single whole number of at least 1.",
         call. = FALSE)
  }
  max_iter <- as.integer(max_iter)
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`tol` must be a single positive finite number.", call. = FALSE)
  }
  eig_tol <- 1e-8

  ## ---- design -------------------------------------------------------------
  v <- validate_staggered_panel(data, yname, tname, idname, gname,
                                covar = covar, never_treated = never_treated)

  cells  <- build_cells(v$treated_cohorts, v$T_max)
  N_beta <- nrow(cells)

  contrasts_df <- enumerate_contrasts(cells, v$treated_cohorts, v$T_min,
                                      v$has_never)
  if (nrow(contrasts_df) == 0L) {
    stop("No clean comparison exists for this design: every treated cohort ",
         "lacks both a never-treated control and a later not-yet-treated ",
         "cohort. No treatment effect is identified.", call. = FALSE)
  }

  coh_ymeans <- group_period_means(v$Y_mat, v$grp_of, v$in_design,
                                   v$all_groups)

  g_pos  <- match(contrasts_df$g, v$all_groups)
  c_pos  <- match(contrasts_df$c, v$all_groups)
  tp_idx <- contrasts_df$t_post - v$T_min + 1L
  tr_idx <- contrasts_df$t_pre  - v$T_min + 1L

  ## ---- drop comparisons touching an unobserved group-period ---------------
  usable <- is.finite(coh_ymeans[cbind(g_pos, tp_idx)]) &
            is.finite(coh_ymeans[cbind(g_pos, tr_idx)]) &
            is.finite(coh_ymeans[cbind(c_pos, tp_idx)]) &
            is.finite(coh_ymeans[cbind(c_pos, tr_idx)])
  n_dropped <- sum(!usable)
  if (n_dropped > 0L) {
    contrasts_df <- contrasts_df[usable, , drop = FALSE]
    g_pos  <- g_pos[usable]
    c_pos  <- c_pos[usable]
    tp_idx <- tp_idx[usable]
    tr_idx <- tr_idx[usable]
    warning(sprintf(
      "%d of %d clean comparison(s) reference a group-period with no observed ",
      n_dropped, n_dropped + nrow(contrasts_df)),
      "unit and have been dropped from the moment system.", call. = FALSE)
  }
  if (nrow(contrasts_df) == 0L) {
    stop("Every clean comparison references a group-period with no observed ",
         "unit; no treatment effect is identified.", call. = FALSE)
  }
  N_2x2 <- nrow(contrasts_df)

  ## ---- moment vector ------------------------------------------------------
  treated_change <- coh_ymeans[cbind(g_pos, tp_idx)] -
                    coh_ymeans[cbind(g_pos, tr_idx)]

  cov_n_fallback <- 0L
  if (is.null(v$covar)) {
    control_term <- coh_ymeans[cbind(c_pos, tp_idx)] -
                    coh_ymeans[cbind(c_pos, tr_idx)]
  } else {
    coh_covmeans <- group_covariate_means(v$X_unit, v$grp_of, v$in_design,
                                          v$all_groups)
    dcov <- build_delta_cov(contrasts_df, v$grp_of, v$in_design, v$Y_mat,
                            v$X_unit, v$T_min, coh_ymeans, coh_covmeans,
                            g_pos, c_pos, tp_idx, tr_idx)
    control_term   <- dcov$control_term
    cov_n_fallback <- dcov$n_fallback
    if (cov_n_fallback > 0L) {
      warning(sprintf(
        "%d of %d comparison(s) could not identify the covariate adjustment ",
        cov_n_fallback, N_2x2),
        "(too few complete-case control units, or a rank-deficient fit) and ",
        "reverted to the unconditional comparison for those contrasts only.",
        call. = FALSE)
    }
  }
  Delta <- treated_change - control_term

  ## ---- identification -----------------------------------------------------
  Q_H <- build_Q_H(contrasts_df$focal_col, N_beta)
  n_comparisons <- as.integer(colSums(Q_H))
  id_idx <- which(n_comparisons > 0L)
  if (length(id_idx) == 0L) {
    stop("No cohort-by-time effect is anchored by a clean comparison; ",
         "nothing is identified.", call. = FALSE)
  }
  identified <- logical(N_beta)
  identified[id_idx] <- TRUE

  if (length(id_idx) < N_beta) {
    bad <- cells[!identified, c("g", "t"), drop = FALSE]
    warning(sprintf(
      "%d of %d cohort-by-time effect(s) have no clean comparison and are ",
      N_beta - length(id_idx), N_beta),
      "reported as NA, not as an estimated zero. Affected: ",
      .preview(sprintf("(g=%d, t=%d)", bad$g, bad$t), 8L), ".",
      if (!v$has_never)
        " A never-treated control group, if available, would identify more."
      else "",
      call. = FALSE)
  }

  Q_S     <- Q_H[, id_idx, drop = FALSE]
  QtQ     <- crossprod(Q_S)
  QtQ_inv <- diag(1 / diag(QtQ), nrow = ncol(Q_S))   # Q_S is a selection matrix
  beta_id <- as.numeric(QtQ_inv %*% crossprod(Q_S, Delta))

  U_fac   <- build_U_factor(contrasts_df, v$all_groups, v$T_min, v$TT)
  red_fac <- gmm_reduced_factor(U_fac, Q_S, Delta)

  ## ---- treatment-effect map for the residualisation -----------------------
  obs_cohort <- v$cohort_of[v$unit_row]
  obs_time   <- v$time_idx + v$T_min - 1L
  cell_id <- match(paste(obs_cohort, obs_time, sep = "\r"),
                   paste(cells$g, cells$t, sep = "\r"))
  has_cell <- !is.na(cell_id)

  ## ---- iterated GMM -------------------------------------------------------
  converged   <- FALSE
  solve_ok    <- FALSE
  termination <- "max_iter_reached"
  n_iter      <- 0L
  max_change  <- NA_real_
  QtAQ_out    <- NULL
  W_last      <- NULL
  n_kept      <- NA_integer_
  n_zero_lag  <- 0L

  for (iter in seq_len(max_iter)) {
    n_iter   <- iter
    beta_old <- beta_id

    beta_full <- rep(NA_real_, N_beta)
    beta_full[id_idx] <- beta_id
    tau <- numeric(length(v$y))
    tau[has_cell] <- beta_full[cell_id[has_cell]]
    y_adj <- v$y - tau        # NA at observations in unidentified cells

    R  <- fe_residual_matrix(y_adj, v$unit_row, v$time_idx, v$N_units, v$TT)
    wt <- build_weight_matrix(R, v$obs_mat, v$grp_of, v$in_design,
                              v$all_groups, weighting)
    W_last     <- wt$W
    n_zero_lag <- wt$n_zero_lag

    red <- gmm_reduced_solve(red_fac, W_last, tol_rel = eig_tol)
    if (is.null(red)) {
      termination <- "singular_weight_matrix"
      break
    }
    # Recorded before the solve: when the solve fails, this is the diagnostic
    # that explains why -- the weight supports fewer independent moment
    # directions than there are effects to identify.
    n_kept <- red$n_kept

    beta_new <- tryCatch(solve(red$QtAQ, red$QtAD), error = function(e) NULL)
    if (is.null(beta_new)) {
      termination <- "singular_normal_equations"
      break
    }
    if (anyNA(beta_new) || !all(is.finite(beta_new))) {
      termination <- "non_finite_update"
      break
    }

    solve_ok   <- TRUE
    QtAQ_out   <- red$QtAQ
    max_change <- max(abs(as.numeric(beta_new) - beta_old))
    beta_id    <- as.numeric(beta_new)

    if (is.finite(max_change) && max_change < tol) {
      converged   <- TRUE
      termination <- "converged"
      break
    }
  }

  ## ---- variance -----------------------------------------------------------
  if (solve_ok) {
    V_id <- tryCatch(solve(QtAQ_out),
                     error = function(e) MASS::ginv(as.matrix(QtAQ_out)))
  } else {
    V_id <- gmm_seed_sandwich(red_fac, W_last, QtQ_inv)
  }
  V_id <- (V_id + t(V_id)) / 2

  if (!solve_ok) {
    # No weighted step was ever completed, so what is being returned is the
    # identity-weighted seed rather than a GMM estimate under the requested
    # weighting. This is a different situation from an efficient iteration
    # that ran but did not settle, and it is always worth surfacing --
    # including under "unrestricted", where a merely unconverged run is not.
    warning(sprintf(
      "No weighted GMM step could be completed (termination: %s). ",
      termination),
      sprintf("The reported estimates are the identity-weighted seed, not a %s-weighted ",
              weighting),
      "estimator, and standard errors use the sandwich fallback. ",
      if (!is.na(n_kept))
        sprintf(paste0("The estimated weighting matrix supports %d independent moment ",
                       "direction(s), but %d effect(s) must be identified. "),
                n_kept, length(id_idx))
      else "",
      "A cohort-specific covariance carries more parameters than a short ",
      "panel with few units per cohort can support; a more parsimonious ",
      "`weighting`, or more units, is needed for this design.",
      call. = FALSE)
  } else if (!converged && !identical(weighting, "unrestricted")) {
    warning(sprintf(
      "The efficient GMM reweighting did not converge (termination: %s) ",
      termination),
      sprintf("after %d of at most %d iteration(s), with a final change of %.2e. ",
              n_iter, max_iter, max_change),
      "Estimates are the last valid iterate. Consider a larger `max_iter`, ",
      "a looser `tol`, or a more parsimonious `weighting`.", call. = FALSE)
  }

  ## ---- assemble -----------------------------------------------------------
  beta_full <- rep(NA_real_, N_beta)
  beta_full[id_idx] <- beta_id
  se_full <- rep(NA_real_, N_beta)
  se_full[id_idx] <- sqrt(pmax(0, diag(V_id)))

  cell_names <- sprintf("g%d:t%d", cells$g, cells$t)
  names(beta_full) <- cell_names
  dimnames(V_id)   <- list(cell_names[id_idx], cell_names[id_idx])

  w <- aggregation_weights(cells, v$N_g, v$treated_cohorts)
  agg <- list(
    CW = aggregate_att(w$CW, beta_id, V_id, id_idx, N_beta),
    EW = aggregate_att(w$EW, beta_id, V_id, id_idx, N_beta)
  )

  catt <- data.frame(
    g             = cells$g,
    t             = cells$t,
    event_time    = cells$t - cells$g,
    estimate      = as.numeric(beta_full),
    std_error     = se_full,
    identified    = identified,
    n_comparisons = n_comparisons,
    row.names     = NULL
  )

  out <- list(
    call         = cl,
    weighting    = weighting,
    coefficients = beta_full,
    vcov         = V_id,
    catt         = catt,
    aggregate    = agg,
    weights      = list(CW = w$CW, EW = w$EW),
    convergence  = list(converged   = converged,
                        solve_ok    = solve_ok,
                        termination = termination,
                        n_iter      = n_iter,
                        max_change  = max_change,
                        tol         = tol,
                        max_iter    = max_iter),
    design       = list(N_units          = v$N_units,
                        T_periods        = v$TT,
                        treated_cohorts  = v$treated_cohorts,
                        n_cohorts        = length(v$treated_cohorts),
                        N_g              = v$N_g,
                        N_beta           = N_beta,
                        N_identified     = length(id_idx),
                        N_2x2            = N_2x2,
                        n_dropped_2x2    = n_dropped,
                        moment_rank      = red_fac$rank,
                        n_kept           = n_kept,
                        n_zero_lag       = n_zero_lag),
    controls     = list(has_never        = v$has_never,
                        n_never_zero     = v$n_never_zero,
                        n_never_late     = v$n_never_late,
                        covar            = v$covar,
                        cov_n_fallback   = cov_n_fallback),
    validation   = list(balanced         = v$balanced,
                        n_missing_cells  = v$n_missing_cells,
                        n_obs            = length(v$y)),
    .internal    = list(coh_ymeans      = coh_ymeans,
                        W               = W_last,
                        all_groups      = v$all_groups,
                        treated_cohorts = v$treated_cohorts,
                        has_never       = v$has_never,
                        T_min           = v$T_min,
                        T_max           = v$T_max,
                        TT              = v$TT,
                        N_gt            = v$N_gt)
  )
  class(out) <- "staggered_gmm"
  out
}
