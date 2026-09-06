test_that("both weight vectors sum to one", {
  fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                       idname = "unit_id", gname = "cohort")
  expect_equal(sum(fit$weights$CW), 1, tolerance = 1e-12)
  expect_equal(sum(fit$weights$EW), 1, tolerance = 1e-12)
})


test_that("CW weights each treated observation equally", {
  cells <- staggeredGMM:::build_cells(c(10L, 13L), 15L)
  N_g <- c("10" = 20L, "13" = 5L)
  w <- staggeredGMM:::aggregation_weights(cells, N_g, c(10L, 13L))

  T_g <- c("10" = 6L, "13" = 3L)     # 10..15 and 13..15
  denom <- 20 * 6 + 5 * 3
  expect_equal(w$CW[cells$g == 10L], rep(20 / denom, 6))
  expect_equal(w$CW[cells$g == 13L], rep(5 / denom, 3))

  # Cohort totals are proportional to N_g * T_g, so exposure length matters.
  expect_equal(sum(w$CW[cells$g == 10L]), 120 / denom)
  expect_equal(sum(w$CW[cells$g == 13L]), 15 / denom)
})


test_that("EW gives every cohort the same total weight", {
  cells <- staggeredGMM:::build_cells(c(10L, 13L), 15L)
  N_g <- c("10" = 20L, "13" = 5L)
  w <- staggeredGMM:::aggregation_weights(cells, N_g, c(10L, 13L))

  expect_equal(sum(w$EW[cells$g == 10L]), 0.5)
  expect_equal(sum(w$EW[cells$g == 13L]), 0.5)
  # Within a cohort, spread equally over its own post-treatment window.
  expect_equal(w$EW[cells$g == 10L], rep(0.5 / 6, 6))
  expect_equal(w$EW[cells$g == 13L], rep(0.5 / 3, 3))
})


test_that("CW and EW coincide only when exposure is equal", {
  # Equal exposure windows: the two schemes must agree.
  cells <- data.frame(g = c(5L, 5L, 6L, 6L), t = c(5L, 6L, 6L, 7L),
                      col_id = 1:4)
  N_g <- c("5" = 4L, "6" = 4L)
  w <- staggeredGMM:::aggregation_weights(cells, N_g, c(5L, 6L))
  expect_equal(w$CW, w$EW)
})


test_that("aggregates match a direct weighted sum", {
  fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                       idname = "unit_id", gname = "cohort")
  manual_cw <- sum(fit$weights$CW * fit$catt$estimate)
  manual_ew <- sum(fit$weights$EW * fit$catt$estimate)
  expect_equal(fit$aggregate$CW$estimate, manual_cw)
  expect_equal(fit$aggregate$EW$estimate, manual_ew)
})


test_that("the aggregate standard error is the quadratic form in vcov", {
  fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                       idname = "unit_id", gname = "cohort")
  w <- fit$weights$CW
  manual <- sqrt(as.numeric(crossprod(w, vcov(fit) %*% w)))
  expect_equal(fit$aggregate$CW$std_error, manual, tolerance = 1e-10)
})


test_that("the identified-subset aggregate renormalises correctly", {
  w <- c(0.5, 0.3, 0.2)
  beta_id <- c(2, 4)
  V_id <- diag(c(1, 1))
  a <- staggeredGMM:::aggregate_att(w, beta_id, V_id, id_idx = c(1L, 3L),
                                    N_beta = 3L)

  expect_true(is.na(a$estimate))
  expect_equal(a$identified_weight_share, 0.7)
  expect_equal(a$estimate_identified, (0.5 * 2 + 0.2 * 4) / 0.7)
})


test_that("a fully identified design reports the aggregate itself", {
  w <- c(0.5, 0.5)
  a <- staggeredGMM:::aggregate_att(w, c(2, 4), diag(2), id_idx = c(1L, 2L),
                                    N_beta = 2L)
  expect_equal(a$estimate, 3)
  expect_equal(a$estimate_identified, 3)
  expect_equal(a$identified_weight_share, 1)
})
