test_that("covariate adjustment runs and changes the estimate", {
  fit0 <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                        idname = "unit_id", gname = "cohort")
  fit1 <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                        idname = "unit_id", gname = "cohort",
                        covar = c("x1", "x2"))

  expect_equal(fit1$controls$covar, c("x1", "x2"))
  expect_equal(fit1$design$N_beta, fit0$design$N_beta)
  expect_false(isTRUE(all.equal(unname(coef(fit0)), unname(coef(fit1)))))
})


test_that("an empty covariate specification is the same as none", {
  fit0 <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                        idname = "unit_id", gname = "cohort")
  for (empty in list(NULL, character(0), "")) {
    fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                         idname = "unit_id", gname = "cohort", covar = empty)
    expect_equal(unname(coef(fit)), unname(coef(fit0)), tolerance = 1e-12)
    expect_null(fit$controls$covar)
  }
})


test_that("a missing covariate column is reported", {
  expect_error(
    gmm_staggered(sim_panel, yname = "y", tname = "year",
                  idname = "unit_id", gname = "cohort", covar = "nope"),
    "not found in `data`"
  )
})


test_that("a covariate may not double as a key column", {
  expect_error(
    gmm_staggered(sim_panel, yname = "y", tname = "year",
                  idname = "unit_id", gname = "cohort", covar = "year"),
    "must not name the outcome, time, id or cohort"
  )
})


test_that("the baseline is resolved on time, not on row order", {
  # A covariate that varies within unit must resolve to the EARLIEST period's
  # value regardless of how the caller sorted the rows.
  dat <- make_test_panel(cohorts = c(3L, 4L), n_per = 3L, n_never = 3L,
                         T_total = 5L, sd = 0.2, seed = 71L)
  dat$x <- ifelse(dat$year == 1L, 10, 99)      # baseline value is 10

  asc  <- dat[order(dat$unit_id, dat$year), ]
  desc <- dat[order(dat$unit_id, -dat$year), ]

  v_asc <- suppressWarnings(
    staggeredGMM:::validate_staggered_panel(asc, "y", "year", "unit_id",
                                            "cohort", covar = "x"))
  v_desc <- suppressWarnings(
    staggeredGMM:::validate_staggered_panel(desc, "y", "year", "unit_id",
                                            "cohort", covar = "x"))

  expect_true(all(v_asc$X_unit[, "x"] == 10))
  expect_equal(v_asc$X_unit, v_desc$X_unit)
})


test_that("a time-varying covariate warns", {
  dat <- make_test_panel(cohorts = c(3L, 4L), n_per = 3L, n_never = 3L,
                         T_total = 5L, sd = 0.2, seed = 72L)
  # Must vary across units as well as within them: a covariate whose
  # baseline is identical for every unit is collinear with the intercept
  # and would trigger the rank-deficiency fallback instead.
  set.seed(721)
  xu <- stats::rnorm(length(unique(dat$unit_id)))
  dat$x <- xu[dat$unit_id] + 0.1 * dat$year
  expect_warning(
    gmm_staggered(dat, "y", "year", "unit_id", "cohort", covar = "x"),
    "not constant over time"
  )
})


test_that("a genuinely time-invariant covariate does not warn", {
  dat <- make_test_panel(cohorts = c(3L, 4L), n_per = 3L, n_never = 3L,
                         T_total = 5L, sd = 0.2, seed = 73L)
  set.seed(74)
  xu <- stats::rnorm(length(unique(dat$unit_id)))
  dat$x <- xu[dat$unit_id]
  expect_no_warning(
    gmm_staggered(dat, "y", "year", "unit_id", "cohort", covar = "x")
  )
})


test_that("the covariate design keeps an intercept", {
  dat <- make_test_panel(cohorts = c(3L, 4L), n_per = 3L, n_never = 3L,
                         T_total = 5L, sd = 0.2, seed = 75L)
  set.seed(76)
  xu <- stats::rnorm(length(unique(dat$unit_id)))
  dat$x <- xu[dat$unit_id]

  v <- staggeredGMM:::validate_staggered_panel(dat, "y", "year", "unit_id",
                                               "cohort", covar = "x")
  expect_true("(Intercept)" %in% colnames(v$X_unit))
  expect_equal(nrow(v$X_unit), v$N_units)
})


test_that("adjustment falls back per triple when unidentifiable", {
  # Two covariates but only two control units: gamma cannot be identified,
  # so every affected comparison reverts to the unconditional form.
  dat <- make_test_panel(cohorts = c(3L, 4L), n_per = 2L, n_never = 2L,
                         T_total = 5L, sd = 0.2, seed = 77L)
  set.seed(78)
  nu <- length(unique(dat$unit_id))
  dat$xa <- stats::rnorm(nu)[dat$unit_id]
  dat$xb <- stats::rnorm(nu)[dat$unit_id]

  fit <- suppressWarnings(
    gmm_staggered(dat, "y", "year", "unit_id", "cohort",
                  covar = c("xa", "xb"))
  )
  expect_gt(fit$controls$cov_n_fallback, 0L)
  expect_true(all(is.finite(fit$catt$estimate[fit$catt$identified])))
})
