# Aggregation of the cohort-by-time effects into scalar ATTs.
#
# Two aggregation schemes are reported. Both are weighted averages of the
# same CATT vector; they differ in what population they target.
#
#   CW ("treated-observation weighting")
#       w_{g,t} = N_g / sum_h N_h T_h
#     Every treated unit-period receives equal weight, so a cohort observed
#     for more post-treatment periods contributes more total weight. This is
#     the aggregation used throughout Arora and Bijani (2026), and it
#     coincides with the "simple" aggregation of Callaway and Sant'Anna
#     (2021).
#
#   EW ("cohort-equal weighting")
#       w_{g,t} = (1/K) / T_g
#     Each of the K treated cohorts receives total weight 1/K, spread equally
#     over its own post-treatment window, so exposure length does not affect
#     a cohort's influence.
#
# Neither is more correct than the other; they answer different questions.
# Users wanting a third scheme can weight the returned CATT table directly.
#
# PARTIAL IDENTIFICATION. A CATT with no clean comparison is not estimable,
# and is reported as NA rather than 0. Where such a cell carries positive
# weight, the corresponding aggregate is NA as well: the requested target is
# not identified by the design. An identified-subset aggregate, renormalised
# over the estimable cells only, is reported alongside it together with the
# share of weight it covers.


# @return list(CW = numeric(N_beta), EW = numeric(N_beta)); each sums to 1.
# @noRd
aggregation_weights <- function(cells, N_g, treated_cohorts) {
  N_beta <- nrow(cells)
  g_chr  <- as.character(cells$g)

  T_g <- table(cells$g)                          # post-periods per cohort
  T_g_of <- as.numeric(T_g[g_chr])
  N_g_of <- as.numeric(N_g[g_chr])

  denom_cw <- sum(as.numeric(N_g[names(T_g)]) * as.numeric(T_g))
  w_CW <- N_g_of / denom_cw

  K <- length(treated_cohorts)
  w_EW <- (1 / K) / T_g_of

  # Both are exact by construction; renormalise only against accumulated
  # floating-point error.
  list(CW = w_CW / sum(w_CW), EW = w_EW / sum(w_EW))
}


# Aggregate a weight vector against the estimated CATTs.
#
# @param w full-length weight vector, one entry per cell.
# @param beta_id estimates for the identified cells only.
# @param V_id covariance of `beta_id`.
# @param id_idx indices (into the full cell set) of the identified cells.
# @param N_beta total number of cells.
# @noRd
aggregate_att <- function(w, beta_id, V_id, id_idx, N_beta) {
  identified <- logical(N_beta)
  identified[id_idx] <- TRUE

  unidentified_weight <- sum(w[!identified])
  fully_identified <- unidentified_weight <= 0

  w_id <- w[id_idx]
  share <- sum(w_id)

  if (share > 0) {
    w_sub <- w_id / share
    est_sub <- sum(w_sub * beta_id)
    var_sub <- as.numeric(crossprod(w_sub, V_id %*% w_sub))
    se_sub  <- sqrt(max(0, var_sub))
  } else {
    est_sub <- NA_real_
    se_sub  <- NA_real_
  }

  if (fully_identified) {
    est <- est_sub
    se  <- se_sub
  } else {
    est <- NA_real_
    se  <- NA_real_
  }

  list(estimate               = est,
       std_error              = se,
       estimate_identified    = est_sub,
       std_error_identified   = se_sub,
       identified_weight_share = share,
       weights                = w)
}
