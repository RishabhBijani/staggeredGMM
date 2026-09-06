test_that("a noiseless panel recovers its CATTs exactly", {
  # With no idiosyncratic error every clean 2x2 equals its target CATT, so
  # Delta = Q beta exactly and any weighting matrix returns beta. This tests
  # the moment construction and the solve independently of the covariance
  # model.
  dat <- make_test_panel(effects = tiny_effects(), sd = 0)

  for (w in c("pooled_toeplitz", "cohort_toeplitz", "unrestricted")) {
    fit <- suppressWarnings(
      gmm_staggered(dat, "y", "year", "unit_id", "cohort", weighting = w)
    )
    expect_equal(unname(coef(fit)), c(1, 2, 3), tolerance = 1e-8,
                 info = w)
    expect_equal(fit$catt$g, c(3L, 3L, 4L), info = w)
    expect_equal(fit$catt$t, c(3L, 4L, 4L), info = w)
    expect_equal(fit$catt$event_time, c(0L, 1L, 0L), info = w)
    expect_true(all(fit$catt$identified), info = w)
  }
})


test_that("the comparison count per cell is reported correctly", {
  dat <- make_test_panel(effects = tiny_effects(), sd = 0.05)
  fit <- suppressWarnings(gmm_staggered(dat, "y", "year", "unit_id", "cohort"))
  # (3,3): 2 pre-periods x {never, cohort 4}; (3,4): 2 x {never};
  # (4,4): 3 pre-periods x {never}.
  expect_equal(fit$catt$n_comparisons, c(4L, 2L, 3L))
  expect_equal(fit$design$N_2x2, 9L)
})


test_that("sim_panel estimates land near their analytic targets", {
  truth <- sim_panel_truth()
  fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                       idname = "unit_id", gname = "cohort")

  expect_true(fit$convergence$converged)
  expect_equal(fit$design$N_identified, 90L)

  cw <- fit$aggregate$CW
  ew <- fit$aggregate$EW
  expect_lt(abs(cw$estimate - truth$ATT_CW), 5 * cw$std_error)
  expect_lt(abs(ew$estimate - truth$ATT_EW), 5 * ew$std_error)

  # The two targets genuinely differ under heterogeneous exposure.
  expect_gt(abs(truth$ATT_CW - truth$ATT_EW), 0.9)
})


test_that("the three weightings agree on the target but not the numbers", {
  fit1 <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                        idname = "unit_id", gname = "cohort",
                        weighting = "pooled_toeplitz")
  fit2 <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                        idname = "unit_id", gname = "cohort",
                        weighting = "cohort_toeplitz")

  expect_equal(fit1$catt[, c("g", "t")], fit2$catt[, c("g", "t")])
  # Over-identified system: a different weight gives a different estimate.
  expect_false(isTRUE(all.equal(unname(coef(fit1)), unname(coef(fit2)))))
})


test_that("convergence and solve success are reported separately", {
  fit <- suppressWarnings(
    gmm_staggered(sim_panel, yname = "y", tname = "year",
                  idname = "unit_id", gname = "cohort", max_iter = 1L)
  )
  expect_true(fit$convergence$solve_ok)
  expect_false(fit$convergence$converged)
  expect_equal(fit$convergence$termination, "max_iter_reached")
  expect_equal(fit$convergence$n_iter, 1L)
  expect_true(is.finite(fit$convergence$max_change))
  expect_gt(fit$convergence$max_change, fit$convergence$tol)
})


test_that("non-convergence warns except under the unrestricted weighting", {
  expect_warning(
    gmm_staggered(sim_panel, yname = "y", tname = "year",
                  idname = "unit_id", gname = "cohort", max_iter = 1L),
    "did not converge"
  )

  # For the unrestricted case the weighted step must actually complete,
  # which needs enough units per cohort to give the weight sufficient rank.
  dat <- make_test_panel(cohorts = c(3L, 4L), n_per = 5L, n_never = 5L,
                         T_total = 4L, sd = 0.3, seed = 55L)
  probe <- suppressWarnings(
    gmm_staggered(dat, "y", "year", "unit_id", "cohort",
                  weighting = "unrestricted", max_iter = 1L))
  skip_if_not(probe$convergence$solve_ok,
              "unrestricted completed no weighted step on the probe design")

  expect_no_warning(
    gmm_staggered(dat, "y", "year", "unit_id", "cohort",
                  weighting = "unrestricted", max_iter = 1L)
  )
})


test_that("failing to complete any weighted step always warns", {
  # sim_panel has 60 units and 90 effects, so the unrestricted weight can
  # support at most sum_g min(N_g, T-1) = 60 moment directions and the
  # normal equations are singular from the first iteration.
  expect_warning(
    fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                         idname = "unit_id", gname = "cohort",
                         weighting = "unrestricted"),
    "identity-weighted seed"
  )
  expect_false(fit$convergence$solve_ok)
  expect_equal(fit$convergence$n_iter, 1L)
  # The rank bound here is sum_g min(N_g, T-1) = 6 * 10 = 60.
  expect_true(fit$design$n_kept < fit$design$N_identified)
})


test_that("unidentified cells are NA, not zero, and warn", {
  # No never-treated group and no later cohort: only (3,3) is identified.
  dat <- make_test_panel(cohorts = c(3L, 4L), n_per = 3L, n_never = 0L,
                         T_total = 4L, effects = tiny_effects(), sd = 0.05,
                         seed = 51L)
  expect_warning(fit <- gmm_staggered(dat, "y", "year", "unit_id", "cohort"),
                 "no clean comparison")

  expect_equal(fit$catt$identified, c(TRUE, FALSE, FALSE))
  expect_false(is.na(fit$catt$estimate[1]))
  expect_true(all(is.na(fit$catt$estimate[2:3])))
  expect_true(all(is.na(fit$catt$std_error[2:3])))
  expect_equal(fit$design$N_identified, 1L)
})


test_that("aggregates are NA when an unidentified cell carries weight", {
  dat <- make_test_panel(cohorts = c(3L, 4L), n_per = 3L, n_never = 0L,
                         T_total = 4L, effects = tiny_effects(), sd = 0.05,
                         seed = 52L)
  fit <- suppressWarnings(
    gmm_staggered(dat, "y", "year", "unit_id", "cohort")
  )

  cw <- fit$aggregate$CW
  expect_true(is.na(cw$estimate))
  expect_true(is.na(cw$std_error))
  expect_false(is.na(cw$estimate_identified))
  expect_gt(cw$identified_weight_share, 0)
  expect_lt(cw$identified_weight_share, 1)
})


test_that("partial identification does not disable the efficient weighting", {
  # Dropping the unidentified columns keeps the normal equations nonsingular,
  # so one unestimable cell must not force every other cell back to the
  # identity-weighted seed.
  dat <- make_test_panel(cohorts = c(3L, 4L), n_per = 6L, n_never = 0L,
                         T_total = 6L, sd = 0.3, seed = 53L)
  fit <- suppressWarnings(
    gmm_staggered(dat, "y", "year", "unit_id", "cohort")
  )
  expect_true(fit$convergence$solve_ok)
  expect_lt(fit$design$N_identified, fit$design$N_beta)
})


test_that("coef, vcov and confint are consistent", {
  fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                       idname = "unit_id", gname = "cohort")
  cf <- coef(fit)
  vc <- vcov(fit)
  ci <- confint(fit)

  expect_equal(length(cf), fit$design$N_beta)
  expect_equal(dim(vc), c(fit$design$N_identified, fit$design$N_identified))
  expect_equal(rownames(vc), colnames(vc))
  expect_equal(dim(ci), c(fit$design$N_beta, 2L))
  expect_true(all(ci[, 1] < cf, na.rm = TRUE))
  expect_true(all(ci[, 2] > cf, na.rm = TRUE))
  expect_error(confint(fit, level = 1.5), "between 0 and 1")
})


test_that("print and summary run without error", {
  fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                       idname = "unit_id", gname = "cohort")
  expect_output(print(fit), "Staggered-adoption GMM estimator")
  expect_output(print(fit), "pooled Toeplitz")
  expect_output(summary(fit), "Cohort-by-time effects")
})
