test_that("autocovariances match a hand computation on a complete matrix", {
  R <- matrix(c(1, 2, 3,
                4, 5, 6), nrow = 2, byrow = TRUE)
  ac <- staggeredGMM:::autocov_sequence(R)

  expect_equal(ac$n_pair, c(6L, 4L, 2L))
  expect_equal(ac$sigma[1], sum(R^2) / 6)
  expect_equal(ac$sigma[2], (1*2 + 2*3 + 4*5 + 5*6) / 4)
  expect_equal(ac$sigma[3], (1*3 + 4*6) / 2)
})


test_that("autocovariance denominators count observed pairs only", {
  R <- matrix(c(1, 2, 3,
                4, NA, 6), nrow = 2, byrow = TRUE)
  ac <- staggeredGMM:::autocov_sequence(R)

  # lag 0: five observed cells; lag 1: only (1,2) and (2,3) from row 1;
  # lag 2: both rows.
  expect_equal(ac$n_pair, c(5L, 2L, 2L))
  expect_equal(ac$sigma[1], (1 + 4 + 9 + 16 + 36) / 5)
  expect_equal(ac$sigma[2], (1*2 + 2*3) / 2)
  expect_equal(ac$sigma[3], (1*3 + 4*6) / 2)
})


test_that("a lag with no observed pair yields zero, not NaN", {
  R <- matrix(c(1, NA,
                2, NA), nrow = 2, byrow = TRUE)
  ac <- staggeredGMM:::autocov_sequence(R)
  expect_equal(ac$n_pair[2], 0L)
  expect_equal(ac$sigma[2], 0)
  expect_false(anyNA(ac$sigma))
})


test_that("the overlap matrix collapses to 1/N_g when balanced", {
  obs <- matrix(TRUE, nrow = 5, ncol = 4)
  K <- staggeredGMM:::overlap_matrix(obs)
  expect_equal(K, matrix(1 / 5, 4, 4))
})


test_that("the overlap matrix reflects the observation pattern", {
  obs <- matrix(TRUE, nrow = 4, ncol = 3)
  obs[1, 2] <- FALSE          # unit 1 missing at period 2
  K <- staggeredGMM:::overlap_matrix(obs)

  # n_a = (4, 3, 4); n_ab[1,2] = 3 units observed at both 1 and 2.
  expect_equal(K[1, 1], 4 / (4 * 4))
  expect_equal(K[2, 2], 3 / (3 * 3))
  expect_equal(K[1, 2], 3 / (4 * 3))
  expect_equal(K[1, 3], 4 / (4 * 4))
})


test_that("a period with no observed unit contributes zero, not Inf", {
  obs <- matrix(TRUE, nrow = 3, ncol = 3)
  obs[, 2] <- FALSE
  K <- staggeredGMM:::overlap_matrix(obs)
  expect_true(all(is.finite(K)))
  expect_true(all(K[, 2] == 0))
  expect_true(all(K[2, ] == 0))
})


test_that("unrestricted covariance is pairwise complete and symmetric", {
  R <- matrix(c(1, 2,
                3, NA,
                5, 6), nrow = 3, byrow = TRUE)
  S <- staggeredGMM:::unrestricted_cov(R)
  expect_equal(S, t(S))
  expect_equal(S[1, 1], (1 + 9 + 25) / 3)
  expect_equal(S[2, 2], (4 + 36) / 2)
  expect_equal(S[1, 2], (1*2 + 5*6) / 2)   # unit 2 excluded from the pair
})


test_that("the balanced weight block reproduces the complete-data formula", {
  # W block should be toeplitz(sigma) / N_g exactly when nothing is missing.
  set.seed(3)
  R <- matrix(rnorm(20), nrow = 5, ncol = 4)
  obs <- matrix(TRUE, 5, 4)
  grp <- rep(1L, 5)

  wt <- staggeredGMM:::build_weight_matrix(
    R, obs, grp_of = grp, in_design = rep(TRUE, 5),
    all_groups = 1L, weighting = "pooled_toeplitz")

  ac <- staggeredGMM:::autocov_sequence(R)
  expect_equal(wt$W, stats::toeplitz(ac$sigma) / 5)
})


test_that("the residual matrix is aligned when outcomes are missing", {
  # The regression is exactly the one the estimator runs internally. With
  # na.rm = FALSE the residual vector is panel-length and NA exactly where
  # the outcome is; the default would be shorter and recycle.
  skip_if_not_installed("fixest")
  set.seed(5)
  N <- 6L; TT <- 5L
  unit_row <- rep(seq_len(N), each = TT)
  time_idx <- rep(seq_len(TT), times = N)
  y <- rnorm(N * TT)
  y[c(3L, 11L, 24L)] <- NA_real_

  R <- staggeredGMM:::fe_residual_matrix(y, unit_row, time_idx, N, TT)
  expect_equal(dim(R), c(N, TT))
  na_cells <- which(is.na(R[cbind(unit_row, time_idx)]))
  expect_equal(na_cells, which(is.na(y)))
})
