# Invariances that would expose indexing or ordering errors.

test_that("results do not depend on row order", {
  set.seed(101)
  shuffled <- sim_panel[sample(nrow(sim_panel)), ]

  a <- gmm_staggered(sim_panel, "y", "year", "unit_id", "cohort")
  b <- gmm_staggered(shuffled, "y", "year", "unit_id", "cohort")

  expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-10)
  expect_equal(a$aggregate$CW$estimate, b$aggregate$CW$estimate,
               tolerance = 1e-10)
  expect_equal(unname(vcov(a)), unname(vcov(b)), tolerance = 1e-10)
})


test_that("results do not depend on the unit id labels", {
  set.seed(102)
  relabel <- sample(sort(unique(sim_panel$unit_id)))
  shifted <- sim_panel
  shifted$unit_id <- relabel[match(sim_panel$unit_id,
                                   sort(unique(sim_panel$unit_id)))]

  a <- gmm_staggered(sim_panel, "y", "year", "unit_id", "cohort")
  b <- gmm_staggered(shifted, "y", "year", "unit_id", "cohort")
  expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-10)
})


test_that("character unit ids give the same answer as numeric ones", {
  chr <- sim_panel
  chr$unit_id <- sprintf("unit%03d", chr$unit_id)

  a <- gmm_staggered(sim_panel, "y", "year", "unit_id", "cohort")
  b <- gmm_staggered(chr, "y", "year", "unit_id", "cohort")
  expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-10)
})


test_that("column order in `data` is irrelevant", {
  rev_cols <- sim_panel[, rev(names(sim_panel))]
  a <- gmm_staggered(sim_panel, "y", "year", "unit_id", "cohort")
  b <- gmm_staggered(rev_cols, "y", "year", "unit_id", "cohort")
  expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-10)
})


test_that("shifting calendar time by a constant is irrelevant", {
  shifted <- sim_panel
  shifted$year   <- shifted$year + 1000L
  shifted$cohort <- ifelse(shifted$cohort == 0L, 0L, shifted$cohort + 1000L)

  a <- gmm_staggered(sim_panel, "y", "year", "unit_id", "cohort")
  b <- gmm_staggered(shifted, "y", "year", "unit_id", "cohort")
  expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-10)
  expect_equal(a$design$moment_rank, b$design$moment_rank)
})


test_that("extra unused columns are ignored", {
  extra <- sim_panel
  extra$junk_chr <- "x"
  extra$junk_na  <- NA
  a <- gmm_staggered(sim_panel, "y", "year", "unit_id", "cohort")
  b <- gmm_staggered(extra, "y", "year", "unit_id", "cohort")
  expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-10)
})


test_that("a data.table is accepted", {
  skip_if_not_installed("data.table")
  dt <- data.table::as.data.table(sim_panel)
  a <- gmm_staggered(sim_panel, "y", "year", "unit_id", "cohort")
  b <- gmm_staggered(dt, "y", "year", "unit_id", "cohort")
  expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-10)
})


test_that("any data.frame subclass is accepted", {
  # Same mechanism as the data.table test above -- the validator calls
  # is.data.frame() then as.data.frame() -- but with no dependency, so this
  # runs everywhere.
  sub <- sim_panel
  class(sub) <- c("my_frame", "data.frame")
  a <- gmm_staggered(sim_panel, "y", "year", "unit_id", "cohort")
  b <- gmm_staggered(sub, "y", "year", "unit_id", "cohort")
  expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-10)
})
