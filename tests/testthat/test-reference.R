# Cross-implementation check of Proposition 1 / Eq. (29): optimal GMM on the
# full stack of clean comparisons, using the generalized inverse of the
# singular N_2x2 x N_2x2 moment covariance, must reproduce the reduced-space
# solve exactly. The naive route is written here from the definition and is
# deliberately slow, so it is run only on a small design.

reconstruct_system <- function(dat) {
  v  <- staggeredGMM:::validate_staggered_panel(dat, "y", "year",
                                                "unit_id", "cohort")
  cells <- staggeredGMM:::build_cells(v$treated_cohorts, v$T_max)
  cd <- staggeredGMM:::enumerate_contrasts(cells, v$treated_cohorts,
                                           v$T_min, v$has_never)
  ym <- staggeredGMM:::group_period_means(v$Y_mat, v$grp_of, v$in_design,
                                          v$all_groups)

  g_pos <- match(cd$g, v$all_groups)
  c_pos <- match(cd$c, v$all_groups)
  tp    <- cd$t_post - v$T_min + 1L
  tr    <- cd$t_pre  - v$T_min + 1L

  Delta <- (ym[cbind(g_pos, tp)] - ym[cbind(g_pos, tr)]) -
           (ym[cbind(c_pos, tp)] - ym[cbind(c_pos, tr)])

  list(Delta = Delta,
       Q     = staggeredGMM:::build_Q_H(cd$focal_col, nrow(cells)),
       U     = staggeredGMM:::build_U_factor(cd, v$all_groups, v$T_min, v$TT),
       cells = cells)
}


small_design <- function(seed = 91L) {
  make_test_panel(cohorts = c(3L, 5L), n_per = 5L, n_never = 5L,
                  T_total = 7L, sd = 0.4, seed = seed)
}


test_that("the reduced-space solve matches naive full-stack optimal GMM", {
  dat <- small_design()

  for (w in c("pooled_toeplitz", "cohort_toeplitz", "unrestricted")) {
    fit <- suppressWarnings(
      gmm_staggered(dat, "y", "year", "unit_id", "cohort", weighting = w)
    )
    skip_if_not(fit$convergence$solve_ok,
                paste("weighted solve did not succeed for", w))

    r <- reconstruct_system(dat)
    Omega <- r$U %*% fit$.internal$W %*% t(r$U)
    Omega <- (Omega + t(Omega)) / 2

    Op <- MASS::ginv(Omega, tol = 1e-8)
    A  <- crossprod(r$Q, Op %*% r$Q)
    b  <- crossprod(r$Q, Op %*% r$Delta)
    beta_naive <- as.numeric(solve(A, b))

    expect_equal(beta_naive, unname(coef(fit)), tolerance = 1e-6, info = w)
  }
})


test_that("the naive and reduced-space variances agree", {
  dat <- small_design(seed = 92L)
  fit <- suppressWarnings(
    gmm_staggered(dat, "y", "year", "unit_id", "cohort")
  )
  skip_if_not(fit$convergence$solve_ok)

  r <- reconstruct_system(dat)
  Omega <- r$U %*% fit$.internal$W %*% t(r$U)
  Omega <- (Omega + t(Omega)) / 2
  Op <- MASS::ginv(Omega, tol = 1e-8)
  V_naive <- solve(crossprod(r$Q, Op %*% r$Q))

  expect_equal(unname(V_naive), unname(vcov(fit)), tolerance = 1e-6)
})


test_that("the moment covariance has the rank the factorisation implies", {
  dat <- small_design(seed = 93L)
  fit <- suppressWarnings(
    gmm_staggered(dat, "y", "year", "unit_id", "cohort")
  )
  r <- reconstruct_system(dat)
  Omega <- r$U %*% fit$.internal$W %*% t(r$U)
  Omega <- (Omega + t(Omega)) / 2

  ev <- eigen(Omega, symmetric = TRUE, only.values = TRUE)$values
  numerical_rank <- sum(ev > 1e-8 * max(ev))
  expect_lte(numerical_rank, fit$design$moment_rank)
  expect_equal(fit$design$moment_rank, 2L * (7L - 1L))
})


test_that("the identity-weighted seed is the per-cell mean of the moments", {
  # With clean comparisons only, Q is a selection matrix, so (Q'Q)^{-1}Q'Delta
  # is just the average moment for each cell. This pins down the starting
  # value the iteration departs from.
  dat <- small_design(seed = 94L)
  r <- reconstruct_system(dat)
  seed <- as.numeric(solve(crossprod(r$Q), crossprod(r$Q, r$Delta)))
  focal <- max.col(r$Q, ties.method = "first")
  manual <- as.numeric(tapply(r$Delta, focal, mean))
  expect_equal(seed, manual, tolerance = 1e-12)
})


test_that("a failed weighted step returns exactly the identity-weighted seed", {
  # When no weighted solve completes, the reported estimate must be the
  # seed -- the plain average of the clean comparisons for each cell --
  # and not something in between. sim_panel under the unrestricted
  # weighting is that case.
  fit <- suppressWarnings(
    gmm_staggered(sim_panel, yname = "y", tname = "year",
                  idname = "unit_id", gname = "cohort",
                  weighting = "unrestricted"))
  expect_false(fit$convergence$solve_ok)

  v  <- staggeredGMM:::validate_staggered_panel(sim_panel, "y", "year",
                                                "unit_id", "cohort")
  cells <- staggeredGMM:::build_cells(v$treated_cohorts, v$T_max)
  cd <- staggeredGMM:::enumerate_contrasts(cells, v$treated_cohorts,
                                           v$T_min, v$has_never)
  ym <- staggeredGMM:::group_period_means(v$Y_mat, v$grp_of, v$in_design,
                                          v$all_groups)
  gp <- match(cd$g, v$all_groups); cp <- match(cd$c, v$all_groups)
  tp <- cd$t_post - v$T_min + 1L;  tr <- cd$t_pre - v$T_min + 1L
  Delta <- (ym[cbind(gp, tp)] - ym[cbind(gp, tr)]) -
           (ym[cbind(cp, tp)] - ym[cbind(cp, tr)])

  seed <- as.numeric(tapply(Delta, cd$focal_col, mean))
  expect_equal(unname(coef(fit)), seed, tolerance = 1e-10)

  w <- staggeredGMM:::aggregation_weights(cells, v$N_g, v$treated_cohorts)$CW
  expect_equal(fit$aggregate$CW$estimate, sum(w * seed), tolerance = 1e-10)
})
