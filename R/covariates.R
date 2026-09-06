# Baseline-covariate (outcome-regression) adjustment, Section 4.5 / Eq. (35)
# of the paper.
#
# Eq. (35) is
#
#   beta_tilde^{cm}_{g,g+k} = (1/N_g) sum_{i in g} [ (Y_{i,g+k} - Y_{i,g-m})
#                              - (Mhat_{c,g+k}(X_i) - Mhat_{c,g-m}(X_i)) ]
#
# with Mhat_{c,t}(x) = x' beta_{c,t} fitted on control cohort c alone. Writing
# gamma_c = beta_{c,g+k} - beta_{c,g-m}, that difference is obtained in one
# regression of the change (Y_{i,g+k} - Y_{i,g-m}) on X_i within cohort c --
# the form used by Callaway and Sant'Anna (2021). Because Mhat is linear in
# X_i and X_i is time-invariant, averaging over cohort g's units commutes
# with the linear map:
#
#   (1/N_g) sum_{i in g} X_i' gamma_c = Xbar_g' gamma_c
#
# exactly. So Eq. (35) collapses, with no approximation, to
#
#   Delta_tilde_s = (Ybar_{g,g+k} - Ybar_{g,g-m}) - Xbar_g' gamma_c
#
# requiring one regression per distinct (control cohort, pre-period,
# post-period) triple, reused across every focal cohort pairing with it.


# Extract the baseline covariate design, one row per unit.
#
# Baseline-ness is asserted, not assumed. Where a covariate is not constant
# within a unit, the unit's value at its EARLIEST observed period is used --
# not the first row in the supplied order, which would depend on how the
# caller happened to sort the data. `time_idx` is required for exactly this
# reason.
#
# `data` must already be in the canonical (unit, time) order, and `unit_row`
# and `time_idx` must be aligned to its rows.
# @noRd
build_baseline_X <- function(data, xformla, unit_ids, unit_row, time_idx) {
  mf     <- stats::model.frame(xformla, data = data,
                               na.action = stats::na.pass)
  X_full <- stats::model.matrix(xformla, data = mf)
  if (nrow(X_full) != length(unit_row)) {
    stop("Internal inconsistency: the covariate design has ", nrow(X_full),
         " rows but the panel has ", length(unit_row),
         ". Please report this with a reproducible example.", call. = FALSE)
  }

  n_units <- length(unit_ids)

  # Earliest observed period per unit, resolved on time rather than row order.
  ord_first <- order(unit_row, time_idx)
  first_row <- ord_first[!duplicated(unit_row[ord_first])]
  first_row <- first_row[order(unit_row[first_row])]
  X_unit    <- X_full[first_row, , drop = FALSE]

  bad_units <- integer(0)
  for (j in seq_len(ncol(X_full))) {
    by_unit <- split(X_full[, j], unit_row)
    bad <- vapply(by_unit, function(v) {
      v <- v[!is.na(v)]
      length(v) > 0L && length(unique(v)) > 1L
    }, logical(1L))
    if (any(bad)) bad_units <- union(bad_units, as.integer(names(by_unit)[bad]))
  }
  if (length(bad_units) > 0L) {
    bad_units <- sort(bad_units)
    warning(sprintf(
      "%d unit(s) have a `covar` value that is not constant over time. ",
      length(bad_units)),
      "`covar` is documented as a baseline (pre-treatment, time-invariant) ",
      "covariate; each affected unit's value at its earliest observed ",
      sprintf("period has been used. Affected units: %s.",
              .preview(unit_ids[bad_units])),
      call. = FALSE)
  }

  rownames(X_unit) <- as.character(unit_ids)
  X_unit
}


# Group-level mean covariate vectors, aligned to `all_groups`.
# @noRd
group_covariate_means <- function(X_unit, grp_of, in_design, all_groups) {
  n_grp <- length(all_groups)
  out <- matrix(NA_real_, nrow = n_grp, ncol = ncol(X_unit),
                dimnames = list(as.character(all_groups), colnames(X_unit)))
  for (gi in seq_len(n_grp)) {
    rows <- which(grp_of == all_groups[gi] & in_design)
    if (length(rows) > 0L)
      out[gi, ] <- colMeans(X_unit[rows, , drop = FALSE], na.rm = TRUE)
  }
  out[!is.finite(out)] <- NA_real_
  out
}


# Covariate-adjusted control-side term of Delta.
#
# gamma_c is fitted once per distinct (control group, pre-period,
# post-period) triple, on that control group's complete-case units only, and
# reused across every contrast row sharing the triple.
#
# Where gamma_c cannot be identified -- fewer complete-case units than
# covariate columns, or a rank-deficient fit -- THAT TRIPLE alone reverts to
# the unconditional control term, i.e. exactly what Delta would have been
# without covariates for contrasts using that triple. The count of affected
# contrast rows is returned so the caller can report it.
# @noRd
build_delta_cov <- function(contrasts_df, grp_of, in_design, Y_mat, X_unit,
                            T_min, coh_ymeans, coh_covmeans,
                            g_pos, c_pos, tp_idx, tr_idx) {

  n_2x2 <- nrow(contrasts_df)
  n_cov <- ncol(X_unit)
  control_term <- numeric(n_2x2)
  n_fallback   <- 0L
  if (n_2x2 == 0L) return(list(control_term = control_term, n_fallback = 0L))

  triple_key <- paste(contrasts_df$c, contrasts_df$t_pre,
                      contrasts_df$t_post, sep = "\r")
  by_triple  <- split(seq_len(n_2x2), triple_key)

  for (rows in by_triple) {
    r1     <- rows[1L]
    c_val  <- contrasts_df$c[r1]
    tr_i   <- tr_idx[r1]
    tp_i   <- tp_idx[r1]

    units_c <- which(grp_of == c_val & in_design)
    dY_c <- Y_mat[units_c, tp_i] - Y_mat[units_c, tr_i]
    X_c  <- X_unit[units_c, , drop = FALSE]

    keep <- stats::complete.cases(X_c) & !is.na(dY_c)
    X_c  <- X_c[keep, , drop = FALSE]
    dY_c <- dY_c[keep]

    gamma_hat <- NULL
    if (nrow(X_c) > n_cov) {
      fit <- tryCatch(stats::lm.fit(x = X_c, y = dY_c),
                      error = function(e) NULL)
      if (!is.null(fit) && isTRUE(fit$rank == n_cov) &&
          !anyNA(fit$coefficients)) {
        gamma_hat <- fit$coefficients
      }
    }

    if (!is.null(gamma_hat)) {
      control_term[rows] <- as.numeric(
        coh_covmeans[g_pos[rows], , drop = FALSE] %*% gamma_hat)
    } else {
      n_fallback <- n_fallback + length(rows)
      control_term[rows] <-
        coh_ymeans[cbind(c_pos[rows], tp_idx[rows])] -
        coh_ymeans[cbind(c_pos[rows], tr_idx[rows])]
    }
  }

  list(control_term = control_term, n_fallback = n_fallback)
}
