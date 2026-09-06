# Design objects. Everything here is a function of the treatment timing and
# the observation window only: none of it depends on the outcome, so all of
# it is computed once and reused across GMM iterations.


# The (cohort, post-period) cells whose CATTs are the parameters.
#
# @return data frame with columns g, t, col_id.
# @noRd
build_cells <- function(treated_cohorts, T_max) {
  g <- rep(treated_cohorts, times = T_max - treated_cohorts + 1L)
  t <- unlist(lapply(treated_cohorts, function(gg) seq.int(gg, T_max)),
              use.names = FALSE)
  data.frame(g = as.integer(g), t = as.integer(t),
             col_id = seq_along(g))
}


# Enumerate clean 2x2 comparisons: never-treated (NT) and not-yet-treated
# (NYT) controls only. Already-treated controls are exact linear combinations
# of these (Lemma 3.2 of the paper) and are therefore never formed.
#
# For each cell (g, t) and each pre-period s in [T_min, g-1]:
#   - one NT comparison against the never-treated pool, if one exists;
#   - one NYT comparison against every treated cohort h with h > t, which is
#     still untreated at the post-period t (and hence also at s < g <= t).
#
# Two passes: count exactly, then fill. Nothing is over-allocated.
#
# @return data frame with columns g, t_post, c, t_pre, type, focal_col.
#   `c = 0` denotes the never-treated pool.
# @noRd
enumerate_contrasts <- function(cells, treated_cohorts, T_min, has_never) {
  n_cells <- nrow(cells)
  if (n_cells == 0L) return(.empty_contrasts())

  n_pre   <- cells$g - T_min                       # pre-periods available
  n_later <- vapply(cells$t, function(tt) sum(treated_cohorts > tt), integer(1L))
  n_ctrl  <- n_later + as.integer(has_never)
  per_cell <- n_pre * n_ctrl
  total    <- sum(per_cell)
  if (total == 0L) return(.empty_contrasts())

  g_out    <- integer(total)
  tp_out   <- integer(total)
  c_out    <- integer(total)
  tr_out   <- integer(total)
  type_out <- character(total)
  fc_out   <- integer(total)

  pos <- 0L
  for (i in seq_len(n_cells)) {
    if (per_cell[i] == 0L) next
    g  <- cells$g[i]
    tt <- cells$t[i]
    fc <- cells$col_id[i]
    pre <- seq.int(T_min, g - 1L)
    ctrl <- c(if (has_never) 0L else integer(0),
              treated_cohorts[treated_cohorts > tt])
    n_here <- length(pre) * length(ctrl)
    idx <- (pos + 1L):(pos + n_here)

    # control varies slowest, pre-period fastest
    c_rep  <- rep(ctrl, each = length(pre))
    pre_rep <- rep(pre, times = length(ctrl))

    g_out[idx]    <- g
    tp_out[idx]   <- tt
    c_out[idx]    <- c_rep
    tr_out[idx]   <- pre_rep
    type_out[idx] <- ifelse(c_rep == 0L, "NT", "NYT")
    fc_out[idx]   <- fc
    pos <- pos + n_here
  }

  data.frame(g = g_out, t_post = tp_out, c = c_out, t_pre = tr_out,
             type = type_out, focal_col = fc_out,
             stringsAsFactors = FALSE)
}


.empty_contrasts <- function() {
  data.frame(g = integer(0), t_post = integer(0), c = integer(0),
             t_pre = integer(0), type = character(0), focal_col = integer(0),
             stringsAsFactors = FALSE)
}


# Incidence matrix mapping comparisons to CATTs. With clean comparisons only,
# every row carries a single +1 in its focal cell's column.
# @noRd
build_Q_H <- function(focal_col, N_beta) {
  n <- length(focal_col)
  Q <- matrix(0, nrow = n, ncol = N_beta)
  if (n > 0L) Q[cbind(seq_len(n), focal_col)] <- 1
  Q
}


# Factor U of the moment covariance Omega = U W U'. Row s of U expresses
# comparison s as a signed combination of the stacked group-by-period means:
# +1 on the focal group at the post period, -1 at the pre period, and the
# mirror image with opposite sign on the control group.
#
# U depends only on the design, so it is built once and its SVD reused across
# all GMM iterations.
# @noRd
build_U_factor <- function(contrasts_df, all_groups, T_min, TT) {
  n_2x2 <- nrow(contrasts_df)
  n_grp <- length(all_groups)
  U <- matrix(0, nrow = n_2x2, ncol = n_grp * TT)
  if (n_2x2 == 0L) return(U)

  grp_pos <- match(contrasts_df$g, all_groups)
  ctl_pos <- match(contrasts_df$c, all_groups)
  tp_i    <- contrasts_df$t_post - T_min + 1L
  tr_i    <- contrasts_df$t_pre  - T_min + 1L
  rows    <- seq_len(n_2x2)

  U[cbind(rows, (grp_pos - 1L) * TT + tp_i)] <-
    U[cbind(rows, (grp_pos - 1L) * TT + tp_i)] + 1
  U[cbind(rows, (grp_pos - 1L) * TT + tr_i)] <-
    U[cbind(rows, (grp_pos - 1L) * TT + tr_i)] - 1
  U[cbind(rows, (ctl_pos - 1L) * TT + tp_i)] <-
    U[cbind(rows, (ctl_pos - 1L) * TT + tp_i)] - 1
  U[cbind(rows, (ctl_pos - 1L) * TT + tr_i)] <-
    U[cbind(rows, (ctl_pos - 1L) * TT + tr_i)] + 1
  U
}


# Group-by-period outcome means, one row per group in `all_groups`.
# Cells with no observed unit are NaN and are handled by the caller.
# @noRd
group_period_means <- function(Y_mat, grp_of, in_design, all_groups) {
  n_grp <- length(all_groups)
  out <- matrix(NA_real_, nrow = n_grp, ncol = ncol(Y_mat),
                dimnames = list(as.character(all_groups), NULL))
  for (gi in seq_len(n_grp)) {
    rows <- which(grp_of == all_groups[gi] & in_design)
    if (length(rows) > 0L)
      out[gi, ] <- colMeans(Y_mat[rows, , drop = FALSE], na.rm = TRUE)
  }
  out[!is.finite(out)] <- NA_real_   # colMeans of an all-NA column gives NaN
  out
}
