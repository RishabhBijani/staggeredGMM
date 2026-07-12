
#' @noRd
build_baseline_X <- function(data, xformla, unit_ids, unit_row) {
  mf      <- stats::model.frame(xformla, data = data, na.action = stats::na.pass)
  X_full  <- model.matrix(xformla, data = mf)
  n_units <- length(unit_ids)

  first_row_of_unit <- match(seq_len(n_units), unit_row)
  X_unit <- X_full[first_row_of_unit, , drop = FALSE]

  bad_units <- integer(0)
  for (j in seq_len(ncol(X_full))) {
    by_unit <- split(X_full[, j], unit_row)
    bad     <- vapply(by_unit, function(v) length(unique(v)) > 1L, logical(1L))
    if (any(bad)) bad_units <- union(bad_units, as.integer(names(by_unit)[bad]))
  }
  if (length(bad_units) > 0L) {
    bad_units <- sort(bad_units)
    warning(sprintf(
      "build_baseline_X: %d unit(s) have a covariate that is not constant across their observed rows; using each affected unit's first-observed value. Affected units: %s.",
      length(bad_units), paste(unit_ids[bad_units], collapse = ", ")
    ), call. = FALSE)
  }

  rownames(X_unit) <- as.character(unit_ids)
  X_unit
}
#' @noRd
build_delta_cov <- function(contrasts_df, cohort_of, Y_mat, X_unit, T_min,
                             coh_ymeans_mat, coh_covmeans_mat,
                             g_row_v, c_row_v, tp_idx_v, tr_idx_v) {

  meta_ctrl <- contrasts_df$c
  meta_tp   <- contrasts_df$t_post
  meta_tr   <- contrasts_df$t_pre
  N_2x2     <- nrow(contrasts_df)
  n_cov     <- ncol(X_unit)

  triple_key <- paste(meta_ctrl, meta_tr, meta_tp, sep = "_")
  uniq_keys  <- unique(triple_key)

  control_term <- numeric(N_2x2)
  n_fallback   <- 0L

  for (key in uniq_keys) {
    rows   <- which(triple_key == key)
    c_val  <- meta_ctrl[rows[1L]]
    tr_val <- meta_tr[rows[1L]]
    tp_val <- meta_tp[rows[1L]]

    units_c <- which(cohort_of == c_val)
    tr_i <- tr_val - T_min + 1L
    tp_i <- tp_val - T_min + 1L
    dY_c <- Y_mat[units_c, tp_i] - Y_mat[units_c, tr_i]
    X_c  <- X_unit[units_c, , drop = FALSE]

    complete <- complete.cases(X_c) & !is.na(dY_c)
    X_c  <- X_c[complete, , drop = FALSE]
    dY_c <- dY_c[complete]

    fittable  <- nrow(X_c) > n_cov
    gamma_hat <- NULL
    if (fittable) {
      fit <- tryCatch(lm.fit(x = X_c, y = dY_c), error = function(e) NULL)
      fittable <- !is.null(fit) && fit$rank == n_cov && !anyNA(fit$coefficients)
      if (fittable) gamma_hat <- fit$coefficients
    }

    if (fittable) {
      control_term[rows] <- as.numeric(
        coh_covmeans_mat[g_row_v[rows], , drop = FALSE] %*% gamma_hat)
    } else {
      n_fallback <- n_fallback + length(rows)
      control_term[rows] <- coh_ymeans_mat[cbind(c_row_v[rows], tp_idx_v[rows])] -
                             coh_ymeans_mat[cbind(c_row_v[rows], tr_idx_v[rows])]
    }
  }

  list(control_term = control_term, n_fallback = n_fallback)
}