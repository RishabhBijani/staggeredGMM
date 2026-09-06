# Deterministic panel builders shared across the test files.
#
# make_test_panel() produces a balanced staggered panel whose treatment
# effects are set exactly, so that with sd = 0 every clean 2x2 comparison
# recovers its target CATT to machine precision and the estimator's output
# can be checked against hand-computed values rather than a snapshot.

make_test_panel <- function(cohorts   = c(3L, 4L),
                            n_per     = 2L,
                            n_never   = 2L,
                            T_total   = 4L,
                            effects   = NULL,
                            sd        = 0,
                            seed      = 1L) {
  set.seed(seed)
  g_vec <- c(rep(as.integer(cohorts), each = n_per), rep(0L, n_never))
  N     <- length(g_vec)

  unit <- rep(seq_len(N), each = T_total)
  time <- rep(seq_len(T_total), times = N)
  g    <- g_vec[unit]

  alpha  <- stats::rnorm(N)
  lambda <- stats::rnorm(T_total)

  tau <- numeric(length(unit))
  if (!is.null(effects)) {
    m <- match(paste(g, time, sep = "\r"),
               paste(effects$g, effects$t, sep = "\r"))
    tau[!is.na(m)] <- effects$value[m[!is.na(m)]]
  }

  eps <- if (sd > 0) stats::rnorm(length(unit), sd = sd) else 0

  data.frame(unit_id = unit,
             year    = time,
             cohort  = g,
             y       = alpha[unit] + lambda[time] + tau + eps)
}


# The three-cell design used for the exact-recovery tests:
#   cohort 3 treated at t = 3 and 4, cohort 4 treated at t = 4.
tiny_effects <- function() {
  data.frame(g = c(3L, 3L, 4L), t = c(3L, 4L, 4L), value = c(1, 2, 3))
}


# Scattered (MCAR) missingness in the outcome. Distinct from block or
# whole-unit dropout, which can leave pooled statistics almost unchanged.
inject_mcar <- function(dat, frac = 0.10, seed = 7L) {
  set.seed(seed)
  n <- as.integer(nrow(dat) * frac)
  dat$y[sample(nrow(dat), n)] <- NA_real_
  dat
}


# Analytic CATTs and aggregate targets for sim_panel's data-generating
# process, recomputed here rather than taken from the package documentation.
sim_panel_truth <- function() {
  cohorts <- c(10L, 13L, 16L, 19L, 22L)
  beta    <- c(-16, -12, -10, -9, -2)
  rate    <- c(0.01, 0.04, 0.08, 0.10, 0.07)
  T_total <- 33L

  cells <- do.call(rbind, lapply(seq_along(cohorts), function(i) {
    g <- cohorts[i]
    data.frame(g = g, t = seq.int(g, T_total),
               catt = beta[i] * (1 + rate[i])^(seq.int(g, T_total) - g))
  }))

  T_g <- table(cells$g)
  # Equal cohort sizes, so treated-observation weighting is a plain mean.
  att_cw <- mean(cells$catt)
  att_ew <- mean(vapply(cohorts,
                        function(g) mean(cells$catt[cells$g == g]),
                        numeric(1L)))
  list(cells = cells, ATT_CW = att_cw, ATT_EW = att_ew, T_g = T_g)
}
