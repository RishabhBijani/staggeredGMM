
library(data.table)
library(MASS)

set.seed(312844)

cohort_size     <- 10L
n_never         <- 10L
treatment_times <- c(10, 13, 16, 19, 22)
n_cohorts       <- length(treatment_times)
N_total         <- n_cohorts * cohort_size + n_never
T_total         <- 33L

rho       <- 0.5
Sigma_ar1 <- toeplitz(rho^(0:(T_total - 1L)))

beta_het <- c(-16, -12, -10, -9, -2)
r_het    <- c(0.01, 0.04, 0.08, 0.10, 0.07)

unit_cohort <- c(rep(treatment_times, each = cohort_size), rep(0L, n_never))
unit_id     <- rep(seq_len(N_total), each = T_total)
time_id     <- rep(seq_len(T_total), times = N_total)
g_vec       <- unit_cohort[unit_id]

tau_vec <- numeric(N_total * T_total)
for (ci in seq_len(n_cohorts)) {
  g_c  <- treatment_times[ci]
  mask <- (g_vec == g_c) & (time_id >= g_c)
  tau_vec[mask] <- beta_het[ci] * (1 + r_het[ci])^(time_id[mask] - g_c)
}

alpha  <- rnorm(N_total)
lambda <- rnorm(T_total)
eps    <- as.vector(t(mvrnorm(N_total, rep(0, T_total), Sigma_ar1)))

cohort_idx_unit <- match(unit_cohort, treatment_times)
cohort_idx_unit[is.na(cohort_idx_unit)] <- 0L

x1_unit <- rnorm(N_total, mean = 0.25 * cohort_idx_unit, sd = 1)
x2_unit <- rbinom(N_total, size = 1, prob = plogis(-0.3 + 0.15 * cohort_idx_unit))

x1_vec <- x1_unit[unit_id]
x2_vec <- x2_unit[unit_id]

t_centered <- time_id - (T_total + 1) / 2
theta1_t   <- 0.06 * t_centered
theta2_t   <- 0.5  * t_centered / T_total

y_vec <- alpha[unit_id] + lambda[time_id] + tau_vec + eps +
  x1_vec * theta1_t + x2_vec * theta2_t

sim_panel <- data.frame(
  unit_id = unit_id,
  year    = time_id,
  cohort  = as.integer(g_vec),
  y       = y_vec,
  x1      = x1_vec,
  x2      = x2_vec
)

compute_true_effects <- function() {
  cells <- do.call(rbind, lapply(seq_len(n_cohorts), function(ci) {
    g_c <- treatment_times[ci]
    data.frame(g = g_c, t = g_c:T_total)
  }))
  cells$true_catt <- with(cells, {
    ci <- match(g, treatment_times)
    beta_het[ci] * (1 + r_het[ci])^(t - g)
  })

  T_g_count    <- table(cells$g)
  N_total_post <- sum(cohort_size * T_g_count)
  w_CW <- cohort_size / N_total_post
  w_EW <- (1 / n_cohorts) / T_g_count[as.character(cells$g)]

  list(
    cells    = cells,
    ATT_CW   = sum(w_CW * cells$true_catt),
    ATT_EW   = sum(as.numeric(w_EW) * cells$true_catt)
  )
}

true_effects <- compute_true_effects()
cat(sprintf("True ATT_CW = %.6f, ATT_EW = %.6f\n",
            true_effects$ATT_CW, true_effects$ATT_EW))

usethis::use_data(sim_panel, overwrite = TRUE)
