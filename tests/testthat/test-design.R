test_that("build_cells enumerates every (cohort, post-period) pair", {
  cells <- staggeredGMM:::build_cells(c(3L, 4L), 4L)
  expect_equal(nrow(cells), 3L)
  expect_equal(cells$g, c(3L, 3L, 4L))
  expect_equal(cells$t, c(3L, 4L, 4L))
  expect_equal(cells$col_id, 1:3)
})


test_that("contrast enumeration matches an independent count", {
  # Independent recount, written from the definition rather than the code
  # under test: for each cell one comparison per (pre-period, control).
  count_contrasts <- function(cohorts, T_min, T_max, has_never) {
    total <- 0L
    for (g in cohorts) {
      n_pre <- g - T_min
      for (tt in seq.int(g, T_max)) {
        n_ctrl <- sum(cohorts > tt) + as.integer(has_never)
        total <- total + n_pre * n_ctrl
      }
    }
    total
  }

  cells <- staggeredGMM:::build_cells(c(3L, 4L), 4L)
  cd <- staggeredGMM:::enumerate_contrasts(cells, c(3L, 4L), 1L, TRUE)
  expect_equal(nrow(cd), count_contrasts(c(3L, 4L), 1L, 4L, TRUE))
  expect_equal(nrow(cd), 9L)

  cd_no_nt <- staggeredGMM:::enumerate_contrasts(cells, c(3L, 4L), 1L, FALSE)
  expect_equal(nrow(cd_no_nt), count_contrasts(c(3L, 4L), 1L, 4L, FALSE))
  expect_equal(nrow(cd_no_nt), 2L)
})


test_that("only clean controls are enumerated", {
  cells <- staggeredGMM:::build_cells(c(3L, 4L), 4L)
  cd <- staggeredGMM:::enumerate_contrasts(cells, c(3L, 4L), 1L, TRUE)

  # Control is never-treated, or a cohort not yet treated at the post period.
  expect_true(all(cd$c == 0L | cd$c > cd$t_post))
  # Pre-period always strictly before the focal cohort's adoption.
  expect_true(all(cd$t_pre < cd$g))
  expect_setequal(unique(cd$type), c("NT", "NYT"))
})


test_that("Q_H is a selection matrix with one entry per row", {
  cells <- staggeredGMM:::build_cells(c(3L, 4L), 4L)
  cd <- staggeredGMM:::enumerate_contrasts(cells, c(3L, 4L), 1L, TRUE)
  Q <- staggeredGMM:::build_Q_H(cd$focal_col, nrow(cells))

  expect_equal(dim(Q), c(nrow(cd), nrow(cells)))
  expect_true(all(rowSums(Q) == 1))
  expect_true(all(Q %in% c(0, 1)))
})


test_that("the moment space has the dimension the theory predicts", {
  # With a never-treated group and a balanced panel the clean comparisons
  # span the interaction subspace, of dimension G(T-1).
  dat <- make_test_panel(effects = tiny_effects(), sd = 0.05, seed = 41L)
  fit <- suppressWarnings(
    gmm_staggered(dat, "y", "year", "unit_id", "cohort"))
  expect_equal(fit$design$moment_rank, 2L * (4L - 1L))

  fit2 <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                        idname = "unit_id", gname = "cohort")
  expect_equal(fit2$design$N_beta, 90L)
  expect_equal(fit2$design$moment_rank, 5L * (33L - 1L))
})


test_that("sim_panel's comparison count matches an independent recount", {
  count_contrasts <- function(cohorts, T_min, T_max, has_never) {
    total <- 0L
    for (g in cohorts) {
      n_pre <- g - T_min
      for (tt in seq.int(g, T_max)) {
        total <- total + n_pre * (sum(cohorts > tt) + as.integer(has_never))
      }
    }
    total
  }
  fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                       idname = "unit_id", gname = "cohort")
  expect_equal(fit$design$N_2x2,
               count_contrasts(c(10L, 13L, 16L, 19L, 22L), 1L, 33L, TRUE))
})


test_that("group-period means ignore never-treated units of the wrong group", {
  dat <- make_test_panel(effects = tiny_effects(), sd = 0)
  v <- staggeredGMM:::validate_staggered_panel(dat, "y", "year", "unit_id",
                                               "cohort")
  ym <- staggeredGMM:::group_period_means(v$Y_mat, v$grp_of, v$in_design,
                                          v$all_groups)
  expect_equal(dim(ym), c(3L, 4L))
  expect_equal(rownames(ym), c("0", "3", "4"))

  manual <- mean(dat$y[dat$cohort == 3L & dat$year == 3L])
  # ym["3", 3] keeps the row dimname: R drops to a length-1 vector but
  # retains the one non-NULL dimnames component. Compare values only.
  expect_equal(unname(ym["3", 3]), manual)
})
