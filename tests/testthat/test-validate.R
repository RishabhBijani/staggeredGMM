test_that("a clean panel estimates without condition", {
  expect_silent(
    gmm_staggered(sim_panel, "y", "year", "unit_id", "cohort")
  )
})


test_that("non-data-frame input is rejected", {
  expect_error(
    gmm_staggered(as.matrix(1:10), "y", "year", "unit_id", "cohort"),
    "must be a data frame"
  )
})


test_that("missing and malformed column-name arguments are rejected", {
  dat <- make_test_panel(effects = tiny_effects())
  expect_error(gmm_staggered(dat, "nope", "year", "unit_id", "cohort"),
               "not found in `data`")
  expect_error(gmm_staggered(dat, c("y", "y"), "year", "unit_id", "cohort"),
               "single non-missing column name")
  expect_error(gmm_staggered(dat, "y", "y", "unit_id", "cohort"),
               "four distinct")
})


test_that("outcome type is checked rather than coerced", {
  dat <- make_test_panel(effects = tiny_effects())

  fac <- dat; fac$y <- factor(round(fac$y))
  expect_error(gmm_staggered(fac, "y", "year", "unit_id", "cohort"),
               "must be numeric")

  chr <- dat; chr$y <- as.character(chr$y)
  expect_error(gmm_staggered(chr, "y", "year", "unit_id", "cohort"),
               "must be numeric")

  lgl <- dat; lgl$y <- dat$y > 0
  expect_error(gmm_staggered(lgl, "y", "year", "unit_id", "cohort"),
               "must be numeric")
})


test_that("non-integer time and cohort are rejected", {
  dat <- make_test_panel(effects = tiny_effects())

  frac <- dat; frac$year <- frac$year + 0.5
  expect_error(gmm_staggered(frac, "y", "year", "unit_id", "cohort"),
               "integer-valued")

  fcoh <- dat; fcoh$cohort <- as.character(fcoh$cohort)
  expect_error(gmm_staggered(fcoh, "y", "year", "unit_id", "cohort"),
               "integer-valued")
})


test_that("missing keys are rejected but missing outcomes are not", {
  dat <- make_test_panel(effects = tiny_effects(), sd = 0.05)

  bad <- dat; bad$year[3] <- NA_integer_
  expect_error(gmm_staggered(bad, "y", "year", "unit_id", "cohort"),
               "missing value")

  bad2 <- dat; bad2$cohort[5] <- NA_integer_
  expect_error(gmm_staggered(bad2, "y", "year", "unit_id", "cohort"),
               "missing value")

  ok <- dat; ok$y[4] <- NA_real_
  expect_no_error(gmm_staggered(ok, "y", "year", "unit_id", "cohort"))
})


test_that("duplicate unit-period rows are rejected", {
  dat <- make_test_panel(effects = tiny_effects())
  dup <- rbind(dat, dat[1, , drop = FALSE])
  expect_error(gmm_staggered(dup, "y", "year", "unit_id", "cohort"),
               "exactly once")
})


test_that("a cohort varying within unit is rejected", {
  dat <- make_test_panel(effects = tiny_effects())
  dat$cohort[dat$unit_id == 1 & dat$year == 2] <- 4L
  expect_error(gmm_staggered(dat, "y", "year", "unit_id", "cohort"),
               "constant within unit")
})


test_that("a period absent from the whole panel is rejected", {
  dat <- make_test_panel(cohorts = c(3L, 5L), n_per = 3L, n_never = 3L,
                         T_total = 6L, sd = 0.05)
  gap <- dat[dat$year != 2L, ]
  expect_error(gmm_staggered(gap, "y", "year", "unit_id", "cohort"),
               "must be consecutive")
})


test_that("a period missing for some units only is accepted", {
  dat <- make_test_panel(cohorts = c(3L, 5L), n_per = 4L, n_never = 4L,
                         T_total = 6L, sd = 0.05, seed = 21L)
  dat <- dat[!(dat$unit_id == 1L & dat$year == 2L), ]
  expect_no_error(
    suppressWarnings(gmm_staggered(dat, "y", "year", "unit_id", "cohort"))
  )
})


test_that("always-treated units raise an error naming the cohort", {
  dat <- make_test_panel(cohorts = c(1L, 3L), n_per = 3L, n_never = 3L,
                         T_total = 5L, sd = 0.05)
  expect_error(gmm_staggered(dat, "y", "year", "unit_id", "cohort"),
               "treated at or before the first observed period")
})


test_that("cohorts adopting after the window become clean controls", {
  base <- make_test_panel(cohorts = c(3L, 4L), n_per = 3L, n_never = 3L,
                          T_total = 5L, effects = tiny_effects(), sd = 0.05,
                          seed = 31L)
  late <- base
  late$cohort[late$cohort == 0L] <- 99L    # never treated, coded as late

  f0 <- gmm_staggered(base, "y", "year", "unit_id", "cohort")
  f1 <- gmm_staggered(late, "y", "year", "unit_id", "cohort")

  expect_true(f1$controls$has_never)
  expect_equal(f1$controls$n_never_late, 3L)
  expect_equal(f1$controls$n_never_zero, 0L)
  expect_equal(unname(coef(f0)), unname(coef(f1)), tolerance = 1e-10)
})


test_that("never_treated = TRUE errors when none exists", {
  dat <- make_test_panel(cohorts = c(3L, 4L, 5L), n_per = 3L, n_never = 0L,
                         T_total = 6L, sd = 0.05)
  expect_error(
    gmm_staggered(dat, "y", "year", "unit_id", "cohort",
                  never_treated = TRUE),
    "no unit is never treated"
  )
})


test_that("never_treated = FALSE warns when never-treated units exist", {
  dat <- make_test_panel(cohorts = c(3L, 4L, 5L), n_per = 3L, n_never = 3L,
                         T_total = 6L, sd = 0.05)
  # Dropping the controls also leaves cells unidentified, which warns too.
  suppressWarnings(expect_warning(
    gmm_staggered(dat, "y", "year", "unit_id", "cohort",
                  never_treated = FALSE),
    "will not be used as controls"
  ))
})


test_that("max_iter and tol are validated", {
  dat <- make_test_panel(effects = tiny_effects(), sd = 0.05)
  expect_error(gmm_staggered(dat, "y", "year", "unit_id", "cohort",
                             max_iter = 0), "at least 1")
  expect_error(gmm_staggered(dat, "y", "year", "unit_id", "cohort",
                             max_iter = 2.5), "whole number")
  expect_error(gmm_staggered(dat, "y", "year", "unit_id", "cohort",
                             tol = -1), "positive finite")
})


test_that("a single period is rejected", {
  dat <- make_test_panel(effects = tiny_effects())
  one <- dat[dat$year == 1L, ]
  expect_error(gmm_staggered(one, "y", "year", "unit_id", "cohort"),
               "two distinct time periods")
})
