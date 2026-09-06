# The specification test is reverse-engineered from Section 4.4 and Appendix
# E of the paper, so the restriction counts are the primary check: on
# sim_panel's design (T = 33, cohorts at 10, 13, 16, 19, 22, never-treated
# pool) Appendix E reports 70 full-set and 14 local-window restrictions.
# Both follow from a base period t0 = g_min - 1 = 9, which is excluded from
# the moment set by construction and so removes one moment from the earliest
# cohort's window.

test_that("the base period is the last period before any adoption", {
  fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                       idname = "unit_id", gname = "cohort")
  j <- gmm_j_test(fit)
  expect_equal(j$base_period, 9L)
})


test_that("the full restriction count matches Appendix E", {
  fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                       idname = "unit_id", gname = "cohort")
  j <- gmm_j_test(fit, type = "full")
  expect_equal(j$n_moments, 70L)

  # Independent recount from the definition.
  cohorts <- c(10L, 13L, 16L, 19L, 22L)
  manual <- sum(vapply(cohorts, function(g) {
    tt <- seq.int(1L, g - 1L); sum(tt != 9L)
  }, integer(1L)))
  expect_equal(j$n_moments, manual)
})


test_that("the local-window restriction count matches Appendix E", {
  fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                       idname = "unit_id", gname = "cohort")
  j <- gmm_j_test(fit, type = "local", window = 3L)
  expect_equal(j$n_moments, 14L)

  cohorts <- c(10L, 13L, 16L, 19L, 22L)
  manual <- sum(vapply(cohorts, function(g) {
    tt <- seq.int(max(1L, g - 3L), g - 1L); sum(tt != 9L)
  }, integer(1L)))
  expect_equal(j$n_moments, manual)
})


test_that("full and local moment sets nest", {
  fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                       idname = "unit_id", gname = "cohort")
  jf <- gmm_j_test(fit, type = "full")
  jl <- gmm_j_test(fit, type = "local", window = 3L)

  key_f <- paste(jf$moments$g, jf$moments$t)
  key_l <- paste(jl$moments$g, jl$moments$t)
  expect_true(all(key_l %in% key_f))
  expect_gt(jf$n_moments, jl$n_moments)
})


test_that("every placebo moment is pre-treatment and excludes the base", {
  fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                       idname = "unit_id", gname = "cohort")
  j <- gmm_j_test(fit, type = "full")
  expect_true(all(j$moments$t < j$moments$g))
  expect_true(all(j$moments$t != j$base_period))
})


test_that("the statistic and p-value are well formed", {
  fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                       idname = "unit_id", gname = "cohort")
  j <- suppressWarnings(gmm_j_test(fit))

  expect_true(is.finite(j$statistic))
  expect_gte(j$statistic, 0)
  expect_gte(j$p_value, 0)
  expect_lte(j$p_value, 1)
  # An estimated Toeplitz sequence need not be positive definite, so the
  # degrees of freedom may be truncated below the moment count.
  expect_gt(j$df, 0L)
  expect_lte(j$df, j$n_moments)
  expect_equal(j$p_value,
               stats::pchisq(j$statistic, df = j$df, lower.tail = FALSE))
})


test_that("the test does not reject under the null", {
  # sim_panel satisfies parallel trends and no anticipation by construction.
  fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                       idname = "unit_id", gname = "cohort")
  j <- suppressWarnings(gmm_j_test(fit, type = "local", window = 3L))
  expect_gt(j$p_value, 0.01)
})


test_that("the test rejects a differential pre-trend", {
  dat <- sim_panel
  coh_index <- match(dat$cohort, c(0L, 10L, 13L, 16L, 19L, 22L)) - 1L
  dat$y <- dat$y + 0.5 * coh_index * dat$year   # cohort-specific linear trend

  fit <- suppressWarnings(
    gmm_staggered(dat, "y", "year", "unit_id", "cohort")
  )
  j <- suppressWarnings(gmm_j_test(fit, type = "local", window = 3L))
  expect_lt(j$p_value, 0.01)
})


test_that("the test requires a never-treated group", {
  dat <- make_test_panel(cohorts = c(3L, 4L, 5L), n_per = 4L, n_never = 0L,
                         T_total = 7L, sd = 0.3, seed = 61L)
  fit <- suppressWarnings(
    gmm_staggered(dat, "y", "year", "unit_id", "cohort")
  )
  expect_error(gmm_j_test(fit), "requires a never-treated control group")
})


test_that("argument checking is strict", {
  fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                       idname = "unit_id", gname = "cohort")
  expect_error(gmm_j_test(fit, window = 0), "at least 1")
  expect_error(gmm_j_test(fit, window = 2.5), "whole number")
  expect_error(gmm_j_test(sim_panel), "must be a fitted")
})


test_that("the result prints", {
  fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                       idname = "unit_id", gname = "cohort")
  j <- suppressWarnings(gmm_j_test(fit))
  expect_output(print(j), "Specification test")
  expect_output(print(j), "df =")
})
