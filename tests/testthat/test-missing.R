# Missingness tests.
#
# Severity depends on the pattern. A whole unit dropping out shifts the
# residual vector by exactly one unit's worth of rows, which leaves pooled
# statistics almost unchanged and can hide a misalignment entirely.
# Scattered missingness does not. Both are tested.

test_that("scattered missingness leaves estimates near the truth", {
  truth <- sim_panel_truth()
  dat <- inject_mcar(sim_panel, frac = 0.10, seed = 7L)

  fit <- suppressWarnings(
    gmm_staggered(dat, "y", "year", "unit_id", "cohort")
  )
  expect_false(fit$validation$balanced)
  expect_gt(fit$validation$n_missing_cells, 0L)

  cw <- fit$aggregate$CW
  expect_true(is.finite(cw$estimate))
  expect_lt(abs(cw$estimate - truth$ATT_CW), 6 * cw$std_error)
})


test_that("a whole unit dropping out is handled", {
  truth <- sim_panel_truth()
  dat <- sim_panel
  dat$y[dat$unit_id == 1L] <- NA_real_

  fit <- suppressWarnings(
    gmm_staggered(dat, "y", "year", "unit_id", "cohort")
  )
  cw <- fit$aggregate$CW
  expect_true(is.finite(cw$estimate))
  expect_lt(abs(cw$estimate - truth$ATT_CW), 6 * cw$std_error)
})


test_that("estimates are insensitive to modest scattered missingness", {
  base <- gmm_staggered(sim_panel, "y", "year", "unit_id", "cohort")
  miss <- suppressWarnings(
    gmm_staggered(inject_mcar(sim_panel, frac = 0.05, seed = 13L),
                  "y", "year", "unit_id", "cohort")
  )
  # Not identical -- fewer observations -- but the same order of magnitude,
  # which a misaligned residual matrix would not deliver.
  expect_lt(abs(base$aggregate$CW$estimate - miss$aggregate$CW$estimate),
            3 * base$aggregate$CW$std_error)
})


test_that("comparisons touching an unobserved group-period are dropped", {
  dat <- sim_panel
  dat$y[dat$cohort == 13L & dat$year == 15L] <- NA_real_

  # Two warnings fire here: dropped comparisons, then the CATT that loses
  # all of its own. Assert both rather than letting one leak as noise.
  expect_warning(
    expect_warning(fit <- gmm_staggered(dat, "y", "year", "unit_id", "cohort"),
                   "no observed"),
    "no clean comparison")
  expect_gt(fit$design$n_dropped_2x2, 0L)
  # The affected CATT loses its own comparisons but the rest survive.
  expect_true(all(is.finite(fit$catt$estimate[fit$catt$identified])))
})


test_that("the balanced path is a special case of the general one", {
  # The pairwise-complete formulas must reproduce the complete-data ones
  # exactly, so that a balanced panel exercises the same code.
  set.seed(17)
  R <- matrix(rnorm(40), nrow = 8, ncol = 5)
  ac <- staggeredGMM:::autocov_sequence(R)

  manual <- vapply(0:4, function(d) {
    i1 <- seq_len(5 - d); i2 <- (1 + d):5
    sum(R[, i1] * R[, i2]) / (8 * (5 - d))
  }, numeric(1L))

  expect_equal(ac$sigma, manual)
  expect_equal(ac$n_pair, as.integer(8 * (5 - 0:4)))
})


test_that("an entirely unobserved outcome column is caught upstream", {
  # Every unit missing at one period is a whole-panel gap in the outcome,
  # not ordinary unbalancedness: the period's group means do not exist.
  dat <- sim_panel
  dat$y[dat$year == 5L] <- NA_real_
  expect_warning(
    fit <- gmm_staggered(dat, "y", "year", "unit_id", "cohort"),
    "no observed"
  )
  expect_gt(fit$design$n_dropped_2x2, 0L)
})
