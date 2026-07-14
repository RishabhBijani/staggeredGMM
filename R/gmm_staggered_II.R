#' @title GMM II: Staggered Treatment Effect Estimator (cohort-specific Toeplitz weighting)
#'
#' @description Estimates cohort x time-specific average treatment effects
#'   (CATTs) under staggered adoption by generalized method of moments (GMM),
#'   using only clean 2x2 DiD comparisons as moment conditions. Each cohort
#'   has its own Toeplitz autocovariance sequence, estimated from that
#'   cohort's own residuals. See [gmm_staggered_I()] (pooled Toeplitz) and
#'   [gmm_staggered_III()] (fully unrestricted within-cohort covariance) for
#'   the other two weightings.
#'
#'   Arora, P. and Bijani, R. (2026). "Efficient Estimation of Treatment
#'   Effects under Staggered Adoption: GMM Approach." Ashoka University
#'   Economics Discussion Paper 163.
#'   <https://dp.ashoka.edu.in/ash/wpaper/paper163_0.pdf>
#'
#'   Arora, P. and Bijani, R. (2026). "Estimating Treatment Effects under
#'   Staggered Timing and Non-Spherical Errors." Available at SSRN:
#'   <https://ssrn.com/abstract=6558759>.
#'
#' @inheritParams gmm_staggered_I
#'
#' @return A list with the same components as [gmm_staggered_I()]. The
#'   cohort-specific Toeplitz sequences underlying this estimator's weighting
#'   are used internally but are not part of the returned value.
#'
#' @section Warnings: `Q_H` (the incidence matrix mapping clean comparisons
#'   to CATTs) can be rank deficient under some designs, most commonly when
#'   a cohort has no later cohort available as a not-yet-treated control and
#'   `has_nt = FALSE`. When this happens, the affected CATTs have no clean
#'   comparison at all, so the fallback reports `beta_hat = 0, SE = 0` for
#'   exactly these cells (not a precisely-estimated null), and a warning
#'   names the exact affected cells. `has_nt = TRUE` with no `gname == 0`
#'   units present in `data` is treated as a likely data-coding error and
#'   raises an error rather than silently proceeding as `has_nt = FALSE`. If
#'   the efficient GMM reweighting does not converge (`eff_ok = FALSE`), a
#'   warning is raised, since non-convergence is unusual for this estimator;
#'   `SE_catt` and `SE_CW` still use a valid sandwich fallback. When `covar`
#'   is supplied, any (control cohort, pre-period, post-period) triple whose
#'   covariate adjustment cannot be identified reverts to the unconditional
#'   comparison for that triple only, and a warning reports how many
#'   contrasts were affected; a separate warning names any units whose
#'   covariate value is not actually constant over time.
#'
#' @examples
#' \donttest{
#' res <- gmm_staggered_II(sim_panel, y_var = "y", t_var = "year",
#'                          id_var = "unit_id", gname = "cohort",
#'                          has_nt = TRUE)
#' res$ATT_CW
#' res$eff_ok
#' }
#'
#' @export

gmm_staggered_II <- function(data, y_var, t_var, id_var, gname,
                             has_nt   = FALSE,
                             max_iter = 10L,
                             covar    = NULL) {

  if (max_iter < 1L) stop("max_iter must be at least 1.", call. = FALSE)

  tol     <- 1e-6
  eig_tol <- 1e-8

  Y_raw  <- as.numeric(data[[y_var]])
  t_raw  <- as.integer(data[[t_var]])
  id_raw <- data[[id_var]]
  g_raw  <- as.integer(data[[gname]])

  unit_ids  <- sort(unique(id_raw))
  N_units   <- length(unit_ids)
  T_min     <- min(t_raw)
  T_max     <- max(t_raw)
  TT        <- T_max - T_min + 1L
  unit_row  <- match(id_raw, unit_ids)
  time_idx  <- t_raw - T_min + 1L

  cohort_of   <- as.integer(tapply(g_raw, unit_row, function(x) x[1L]))
  start_yrs   <- sort(unique(cohort_of[cohort_of > 0L & cohort_of <= T_max]))
  K           <- length(start_yrs)

  if (has_nt && sum(cohort_of == 0L, na.rm = TRUE) == 0L) {
    stop("has_nt = TRUE requested, but no units have gname == 0 in data. ",
         "Set has_nt = FALSE, or check that never-treated units are coded as 0.",
         call. = FALSE)
  }

  all_cohorts <- if (has_nt) c(0L, start_yrs) else start_yrs

  N_g_lookup <- setNames(
    vapply(start_yrs, function(g) sum(cohort_of == g, na.rm = TRUE), integer(1L)),
    as.character(start_yrs)
  )
  if (has_nt)
    N_g_lookup["0"] <- sum(cohort_of == 0L, na.rm = TRUE)

  Y_mat <- matrix(NA_real_, nrow = N_units, ncol = TT)
  Y_mat[cbind(unit_row, time_idx)] <- Y_raw

  dt_long <- data.table(
    unit     = id_raw,
    unit_row = unit_row,
    time     = t_raw,
    time_idx = time_idx,
    Y        = Y_raw,
    g        = g_raw
  )
  setorder(dt_long, unit, time)

  cells        <- build_cells(start_yrs, T_max)
  N_beta       <- nrow(cells)
  contrasts_df <- enumerate_contrasts_clean(cells, all_cohorts, start_yrs,
                                            T_min, has_nt)
  N_2x2      <- nrow(contrasts_df)
  meta_focal <- contrasts_df$g
  meta_ctrl  <- contrasts_df$c
  meta_tp    <- contrasts_df$t_post
  meta_tr    <- contrasts_df$t_pre

  coh_row_map    <- setNames(seq_along(all_cohorts), as.character(all_cohorts))
  coh_ymeans_mat <- matrix(NA_real_, nrow = length(all_cohorts), ncol = TT)
  for (g in all_cohorts) {
    idx <- which(cohort_of == g)
    if (length(idx) > 0L)
      coh_ymeans_mat[coh_row_map[as.character(g)], ] <-
        colMeans(Y_mat[idx, , drop = FALSE], na.rm = TRUE)
  }

  tp_idx_v <- meta_tp - T_min + 1L
  tr_idx_v <- meta_tr - T_min + 1L
  g_row_v  <- as.integer(coh_row_map[as.character(meta_focal)])
  c_row_v  <- as.integer(coh_row_map[as.character(meta_ctrl)])

  no_covar <- is.null(covar) || length(covar) == 0L || identical(covar, "")

  if (no_covar) {
    Delta <- (coh_ymeans_mat[cbind(g_row_v, tp_idx_v)] -
                coh_ymeans_mat[cbind(g_row_v, tr_idx_v)]) -
             (coh_ymeans_mat[cbind(c_row_v, tp_idx_v)] -
                coh_ymeans_mat[cbind(c_row_v, tr_idx_v)])
  } else {
    xformla <- stats::reformulate(covar)
    X_unit  <- build_baseline_X(data, xformla, unit_ids, unit_row)

    coh_covmeans_mat <- matrix(NA_real_, nrow = length(all_cohorts),
                               ncol = ncol(X_unit))
    for (g in all_cohorts) {
      idx <- which(cohort_of == g)
      if (length(idx) > 0L)
        coh_covmeans_mat[coh_row_map[as.character(g)], ] <-
          colMeans(X_unit[idx, , drop = FALSE], na.rm = TRUE)
    }

    dcov <- build_delta_cov(contrasts_df, cohort_of, Y_mat, X_unit, T_min,
                             coh_ymeans_mat, coh_covmeans_mat,
                             g_row_v, c_row_v, tp_idx_v, tr_idx_v)
    if (dcov$n_fallback > 0L) {
      warning(sprintf(
        "%d of %d contrasts could not identify the covariate adjustment (too few complete-case units, or a rank-deficient fit) and reverted to the unconditional comparison for that contrast only.",
        dcov$n_fallback, N_2x2), call. = FALSE)
    }

    Delta <- (coh_ymeans_mat[cbind(g_row_v, tp_idx_v)] -
                coh_ymeans_mat[cbind(g_row_v, tr_idx_v)]) -
             dcov$control_term
  }

  Q_H      <- build_Q_H(contrasts_df, N_2x2, N_beta)
  QtQ      <- crossprod(Q_H)
  QtQ_inv  <- tryCatch(solve(QtQ), error = function(e) ginv(QtQ))
  Q_H_rank <- qr(Q_H)$rank
  beta_hat <- as.numeric(QtQ_inv %*% crossprod(Q_H, Delta))

  if (Q_H_rank < N_beta) {
    zero_cols <- which(colSums(abs(Q_H)) == 0)
    bad_cells <- cells[zero_cols, c("g", "t")]
    warning(sprintf(
      "Q_H is rank deficient (rank %d of %d): %d CATT(s) have no clean comparison at all under this design, so the fallback reports beta_hat = 0, SE = 0 for these cells (not a precisely-estimated null). Affected: %s.%s",
      Q_H_rank, N_beta, length(zero_cols),
      paste(sprintf("(g=%d, t=%d)", bad_cells$g, bad_cells$t), collapse = ", "),
      if (!has_nt) " Consider has_nt = TRUE if a never-treated group is available." else ""
    ), call. = FALSE)
  }

  n_grp         <- length(all_cohorts)
  units_by_grp  <- lapply(all_cohorts, function(g) which(cohort_of == g))
  active_by_grp <- lapply(all_cohorts,
                          function(g) which(meta_focal == g | meta_ctrl == g))
  sign_by_grp   <- lapply(seq_len(n_grp), function(gi) {
    idx <- active_by_grp[[gi]]
    ifelse(meta_focal[idx] == all_cohorts[gi], 1L, -1L)
  })

  T_g_count <- setNames(
    vapply(start_yrs, function(g) sum(cells$g == g), integer(1L)),
    as.character(start_yrs)
  )
  N_total_post <- sum(vapply(start_yrs, function(g)
    N_g_lookup[as.character(g)] * T_g_count[as.character(g)], numeric(1L)))

  w_CW <- numeric(N_beta)
  for (ci in seq_len(N_beta)) {
    g_c      <- as.character(cells$g[ci])
    w_CW[ci] <- N_g_lookup[g_c] / N_total_post
  }
  w_CW <- w_CW / sum(w_CW)

  eff_ok   <- FALSE
  QtAQ_out <- NULL
  W_out    <- NULL
  n_iter   <- 0L

  if (N_2x2 > 0L) {
  m_cols <- n_grp * TT
  U_fac  <- matrix(0.0, nrow = N_2x2, ncol = m_cols)
  for (gi in seq_len(n_grp)) {
    active <- active_by_grp[[gi]]
    if (length(active) == 0L) next
    sgn  <- sign_by_grp[[gi]]
    cols <- ((gi - 1L) * TT + 1L):(gi * TT)
    tp_i <- meta_tp[active] - T_min + 1L
    tr_i <- meta_tr[active] - T_min + 1L
    U_fac[cbind(active, cols[tp_i])] <- U_fac[cbind(active, cols[tp_i])] + sgn
    U_fac[cbind(active, cols[tr_i])] <- U_fac[cbind(active, cols[tr_i])] - sgn
  }
  red_fac <- gmm_reduced_factor(U_fac, Q_H, Delta)

  for (iter in seq_len(max_iter)) {
    n_iter   <- iter
    beta_old <- beta_hat

    idx <- match(paste(dt_long$g, dt_long$time),
                 paste(cells$g,   cells$t))
    dt_long[, tau_hat := ifelse(is.na(idx), 0, beta_hat[idx])]
    dt_long[, Y_adj   := Y - tau_hat]

    fe_mod <- feols(Y_adj ~ 1 | unit + time, data = dt_long)
    R_mat  <- matrix(NA_real_, nrow = N_units, ncol = TT)
    R_mat[cbind(dt_long$unit_row, dt_long$time_idx)] <- residuals(fe_mod)

    W <- matrix(0.0, nrow = m_cols, ncol = m_cols)

    for (gi in seq_len(n_grp)) {
      g_val   <- all_cohorts[gi]
      N_g     <- N_g_lookup[as.character(g_val)]
      units_g <- units_by_grp[[gi]]
      active  <- active_by_grp[[gi]]
      cols    <- ((gi - 1L) * TT + 1L):(gi * TT)

      if (length(active) == 0L || length(units_g) == 0L) next

      rm_g      <- R_mat[units_g, , drop = FALSE]
      sigma_d_g <- numeric(TT)
      for (d in 0L:(TT - 1L)) {
        r1 <- seq_len(TT - d)
        r2 <- (1L + d):TT
        sigma_d_g[d + 1L] <- sum(rm_g[, r1] * rm_g[, r2], na.rm = TRUE) /
          (N_g * (TT - d))
      }
      W[cols, cols] <- toeplitz(sigma_d_g) / N_g
    }
    W_out <- W

    red <- gmm_reduced_solve(red_fac, W, tol_rel = eig_tol)
    if (is.null(red)) break

    QtAQ     <- red$QtAQ
    QtAD     <- red$QtAD
    beta_new <- tryCatch(solve(QtAQ, QtAD), error = function(e) NULL)
    if (is.null(beta_new)) break

    QtAQ_out <- QtAQ
    eff_ok   <- TRUE
    beta_hat <- beta_new

    if (max(abs(beta_hat - beta_old)) < tol) break
  }

  if (eff_ok) {
    V_beta <- tryCatch(
      solve(QtAQ_out),
      error = function(e) ginv(as.matrix(QtAQ_out))
    )
  } else {
    M_out    <- outer(red_fac$d, red_fac$d) * crossprod(red_fac$v, W_out %*% red_fac$v)
    M_out    <- (M_out + t(M_out)) / 2
    QtOmegaQ <- crossprod(red_fac$uQ, M_out %*% red_fac$uQ)
    V_beta   <- QtQ_inv %*% QtOmegaQ %*% QtQ_inv
  }
  } else {
    V_beta <- matrix(0.0, N_beta, N_beta)
  }

  SE_catt <- sqrt(pmax(0, diag(V_beta)))
  ATT_CW  <- as.numeric(crossprod(w_CW, beta_hat))
  SE_CW   <- sqrt(max(0, as.numeric(crossprod(w_CW, V_beta %*% w_CW))))

  if (!eff_ok) {
    warning(sprintf(
      "eff_ok = FALSE: the efficient GMM reweighting did not converge within max_iter = %d iterations; SE_catt and SE_CW use the sandwich fallback from the identity-weighted seed rather than the efficient weighting. Consider increasing max_iter, or checking for near-collinear or very small cohorts.",
      max_iter
    ), call. = FALSE)
  }

  list(
    beta_hat = beta_hat,
    SE_catt  = SE_catt,
    catt_out = data.frame(g = cells$g, t = cells$t,
                          beta_hat = beta_hat, SE = SE_catt),
    ATT_CW   = ATT_CW,
    SE_CW    = SE_CW,
    w_CW     = w_CW,
    cells    = cells,
    eff_ok   = eff_ok,
    n_iter   = n_iter,
    N_beta   = N_beta,
    N_2x2    = N_2x2,
    Q_H_rank = Q_H_rank
  )
}