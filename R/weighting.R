# Estimation of the moment covariance model, and construction of the
# block-diagonal weight W = Cov(stacked group-by-period means).
#
# UNBALANCED PANELS. Under a complete panel, Cov(Ybar_{g,a}, Ybar_{g,b}) =
# sigma_{|a-b|} / N_g. Under missingness the cohort means at a and b are
# taken over different unit sets, and the correct expression is
#
#   Cov(Ybar_{g,a}, Ybar_{g,b}) = sigma_{|a-b|} * N_{g,ab} / (N_{g,a} N_{g,b})
#
# where N_{g,a} counts units of group g observed at a, and N_{g,ab} counts
# those observed at BOTH a and b. Writing K_g for that matrix of ratios, the
# weight block is Sigma_block * K_g (elementwise). When the panel is balanced
# N_{g,ab} = N_{g,a} = N_{g,b} = N_g and K_g collapses to the constant 1/N_g,
# reproducing the balanced formula exactly. There is therefore one code path,
# not two, and every balanced test also exercises the general machinery.
#
# The same logic applies to the autocovariances themselves: each is divided
# by the number of residual pairs actually observed at that lag, not by the
# nominal count the panel would have if complete.
#
# All of this is valid under missingness that is independent of the outcome
# given (group, period). Selection on outcomes is not addressed by any
# weighting scheme.


# Residual matrix from the two-way fixed-effect regression of the
# treatment-adjusted outcome.
#
# `residuals(fit, na.rm = FALSE)` is essential, not cosmetic. The default
# `na.rm = TRUE` returns only the observations used in the fit, which is
# shorter than the panel whenever any outcome is missing; assigning that
# short vector into the full index recycles it and displaces every residual
# after the first missing row.
# @noRd
fe_residual_matrix <- function(y_adj, unit_row, time_idx, N_units, TT) {
  df <- data.frame(y_adj = y_adj,
                   unit  = factor(unit_row),
                   time  = factor(time_idx))
  fit <- fixest::feols(y_adj ~ 1 | unit + time, data = df, notes = FALSE)
  res <- stats::residuals(fit, na.rm = FALSE)
  if (length(res) != length(y_adj)) {
    stop("Internal inconsistency: the fixed-effect fit returned ",
         length(res), " residuals for ", length(y_adj),
         " observations. Please report this with a reproducible example.",
         call. = FALSE)
  }
  R <- matrix(NA_real_, nrow = N_units, ncol = TT)
  R[cbind(unit_row, time_idx)] <- res
  R
}


# Pairwise-complete autocovariance sequence from a residual matrix.
# Lags with no observed pair are returned as 0 and counted.
# @noRd
autocov_sequence <- function(R) {
  TT  <- ncol(R)
  obs <- !is.na(R)
  Rz  <- R
  Rz[!obs] <- 0

  sigma  <- numeric(TT)
  n_pair <- integer(TT)
  for (d in seq.int(0L, TT - 1L)) {
    i1 <- seq_len(TT - d)
    i2 <- (1L + d):TT
    both <- obs[, i1, drop = FALSE] & obs[, i2, drop = FALSE]
    np <- sum(both)
    n_pair[d + 1L] <- np
    if (np > 0L)
      sigma[d + 1L] <- sum(Rz[, i1, drop = FALSE] * Rz[, i2, drop = FALSE]) / np
  }
  list(sigma = sigma, n_pair = n_pair)
}


# Unrestricted within-group covariance, pairwise-complete.
# @noRd
unrestricted_cov <- function(R) {
  TT  <- ncol(R)
  obs <- !is.na(R)
  Rz  <- R
  Rz[!obs] <- 0
  num <- crossprod(Rz)                 # TT x TT
  den <- crossprod(obs + 0)            # units observed at both periods
  S <- matrix(0, TT, TT)
  ok <- den > 0
  S[ok] <- num[ok] / den[ok]
  S
}


# Overlap matrix K_g[a, b] = N_{g,ab} / (N_{g,a} N_{g,b}).
# Entries involving a period with no observed unit are 0: no moment can
# reference such a period, because its group-period mean does not exist.
# @noRd
overlap_matrix <- function(obs_g) {
  n_a  <- colSums(obs_g)
  n_ab <- crossprod(obs_g + 0)
  den  <- outer(n_a, n_a)
  K <- matrix(0, ncol(obs_g), ncol(obs_g))
  ok <- den > 0
  K[ok] <- n_ab[ok] / den[ok]
  K
}


# Block-diagonal weight matrix over the stacked group-by-period means.
#
# @param R residual matrix, N_units x TT, NA where unobserved.
# @param obs_mat observation indicator, same shape.
# @param weighting one of the three covariance models.
# @return list(W, sigma, n_zero_lag) where `sigma` is the pooled
#   autocovariance sequence for "pooled_toeplitz" and NULL otherwise.
# @noRd
build_weight_matrix <- function(R, obs_mat, grp_of, in_design, all_groups,
                                weighting) {
  TT     <- ncol(R)
  n_grp  <- length(all_groups)
  W      <- matrix(0, nrow = n_grp * TT, ncol = n_grp * TT)
  sigma_pooled <- NULL
  n_zero_lag   <- 0L

  if (identical(weighting, "pooled_toeplitz")) {
    ac <- autocov_sequence(R[in_design, , drop = FALSE])
    sigma_pooled <- ac$sigma
    n_zero_lag   <- sum(ac$n_pair == 0L)
  }

  for (gi in seq_len(n_grp)) {
    rows <- which(grp_of == all_groups[gi] & in_design)
    if (length(rows) == 0L) next
    cols  <- ((gi - 1L) * TT + 1L):(gi * TT)
    R_g   <- R[rows, , drop = FALSE]
    obs_g <- obs_mat[rows, , drop = FALSE]

    Sigma_g <- switch(
      weighting,
      pooled_toeplitz = stats::toeplitz(sigma_pooled),
      cohort_toeplitz = {
        ac <- autocov_sequence(R_g)
        n_zero_lag <- n_zero_lag + sum(ac$n_pair == 0L)
        stats::toeplitz(ac$sigma)
      },
      unrestricted    = unrestricted_cov(R_g),
      stop("Unknown weighting scheme: ", weighting, call. = FALSE)
    )

    W[cols, cols] <- Sigma_g * overlap_matrix(obs_g)
  }

  list(W = W, sigma = sigma_pooled, n_zero_lag = n_zero_lag)
}
