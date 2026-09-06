#' Specification test for parallel trends and no anticipation
#'
#' A serial-correlation robust joint test of the identifying assumptions,
#' following Section 4.4 of Arora and Bijani (2026).
#'
#' @details
#' The test is built on the never-treated event-study basis of Eq. (27) of
#' the paper. Fixing a base pre-period \eqn{t_0}, each treated cohort
#' \eqn{g} and period \eqn{t \ne t_0} defines
#'
#' \deqn{\Delta_{g,t} = (\bar{Y}_{g,t} - \bar{Y}_{g,t_0}) -
#'       (\bar{Y}_{\infty,t} - \bar{Y}_{\infty,t_0}).}
#'
#' Those with \eqn{t \ge g} identify the cohort-by-time effects. Those with
#' \eqn{t < g} are placebo restrictions whose expectation is zero only under
#' parallel trends and no anticipation. Collecting a set of \eqn{J} placebo
#' moments into \eqn{\hat{p}} with covariance \eqn{\hat\Omega_p} implied by
#' the fitted serial-correlation model,
#'
#' \deqn{J = \hat{p}' \hat\Omega_p^{-1} \hat{p} \;\to\; \chi^2_J}
#'
#' under the null. Because the placebo moments do not involve the treatment
#' effects, no degrees of freedom are subtracted for estimated parameters.
#'
#' The base period is \eqn{t_0 = g_{\min} - 1}, the last period before any
#' cohort adopts, so that every moment has an untreated reference.
#'
#' A large number of restrictions degrades the chi-squared approximation in
#' small panels. `type = "local"` therefore uses only the placebo moments
#' within `window` periods of each cohort's adoption, which holds its
#' nominal size across sample sizes at essentially no cost in power, since a
#' differential trend or an anticipation effect is already visible just
#' before treatment. `type = "full"` uses every placebo moment and is
#' well-sized only when the number of units is comfortably large.
#'
#' Estimation and testing are deliberately decoupled: the estimate uses the
#' full moment system while the test uses the local window. Using a
#' pre-trends test to gate an estimate distorts the estimator's conditional
#' distribution (Roth, 2022), the more so the stronger their correlation.
#' The test also speaks only to pre-treatment periods. Post-treatment
#' parallel trends, the assumption that actually identifies the effects, is
#' untestable, so a clean result is necessary but not sufficient.
#'
#' @param object A fitted [gmm_staggered()] object.
#' @param type Character. `"local"` (default) restricts the test to placebo
#'   moments within `window` periods before each cohort's adoption; `"full"`
#'   uses all of them.
#' @param window Integer. Width of the local pre-window, in periods. Ignored
#'   when `type = "full"`.
#'
#' @return An object of class `staggered_gmm_jtest`, a list with components
#'   `statistic`, `df`, `p_value`, `type`, `window`, `base_period`,
#'   `n_moments` and `moments` (the placebo moment table).
#'
#' @references
#' Arora, P. and Bijani, R. (2026). "Estimating Treatment Effects under
#' Staggered Timing and Non-Spherical Errors." \doi{10.2139/ssrn.6558759}
#'
#' Roth, J. (2022). "Pretest with Caution: Event-Study Estimates after
#' Testing for Parallel Trends." \emph{American Economic Review: Insights}
#' 4(3), 305-322.
#'
#' @examples
#' fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
#'                      idname = "unit_id", gname = "cohort")
#' gmm_j_test(fit)
#'
#' @export
gmm_j_test <- function(object, type = c("local", "full"), window = 3L) {

  if (!inherits(object, "staggered_gmm")) {
    stop("`object` must be a fitted `staggered_gmm` object, as returned by ",
         "gmm_staggered().", call. = FALSE)
  }
  type <- match.arg(type)
  if (!is.numeric(window) || length(window) != 1L || is.na(window) ||
      window < 1 || window != trunc(window)) {
    stop("`window` must be a single whole number of at least 1.",
         call. = FALSE)
  }
  window <- as.integer(window)

  it <- object$.internal
  if (!isTRUE(it$has_never)) {
    stop("The specification test requires a never-treated control group: the ",
         "event-study basis it is built on differences every cohort against ",
         "the same never-treated reference. This design has none.",
         call. = FALSE)
  }
  if (is.null(it$W)) {
    stop("The fitted object carries no weighting matrix, so the covariance ",
         "of the placebo moments cannot be formed.", call. = FALSE)
  }

  T_min <- it$T_min
  T_max <- it$T_max
  TT    <- it$TT
  cohorts    <- it$treated_cohorts
  all_groups <- it$all_groups
  nt_pos     <- match(0L, all_groups)

  t0 <- min(cohorts) - 1L
  if (t0 < T_min) {
    stop("The earliest treated cohort adopts in the first observed period, ",
         "so no common untreated base period exists and no placebo ",
         "restriction can be formed.", call. = FALSE)
  }
  t0_idx <- t0 - T_min + 1L

  ## ---- assemble the placebo moment list -----------------------------------
  mom_g <- integer(0)
  mom_t <- integer(0)
  for (g in cohorts) {
    tt <- seq.int(T_min, g - 1L)
    if (identical(type, "local")) tt <- tt[tt >= g - window]
    tt <- tt[tt != t0]
    if (length(tt) == 0L) next
    mom_g <- c(mom_g, rep(g, length(tt)))
    mom_t <- c(mom_t, tt)
  }
  if (length(mom_g) == 0L) {
    stop("No placebo restriction is available under `type = \"", type, "\"`",
         if (identical(type, "local"))
           sprintf(" with `window = %d`. Try a wider window or type = \"full\".",
                   window)
         else ".", call. = FALSE)
  }

  ## ---- drop moments referencing an unobserved group-period ----------------
  g_pos <- match(mom_g, all_groups)
  t_idx <- mom_t - T_min + 1L
  ok <- it$N_gt[cbind(g_pos, t_idx)]      > 0L &
        it$N_gt[cbind(g_pos, t0_idx)]     > 0L &
        it$N_gt[cbind(nt_pos, t_idx)]     > 0L &
        it$N_gt[cbind(nt_pos, t0_idx)]    > 0L
  n_dropped <- sum(!ok)
  if (n_dropped > 0L) {
    mom_g <- mom_g[ok]; mom_t <- mom_t[ok]
    g_pos <- g_pos[ok]; t_idx <- t_idx[ok]
    warning(sprintf(
      "%d placebo moment(s) reference a group-period with no observed unit ",
      n_dropped),
      "and have been dropped from the test.", call. = FALSE)
  }
  if (length(mom_g) == 0L) {
    stop("Every placebo moment references a group-period with no observed ",
         "unit; the test cannot be formed.", call. = FALSE)
  }
  n_mom <- length(mom_g)

  ## ---- moment values and their covariance ---------------------------------
  ym <- it$coh_ymeans
  p_hat <- (ym[cbind(g_pos, t_idx)] - ym[cbind(g_pos, rep(t0_idx, n_mom))]) -
           (ym[cbind(rep(nt_pos, n_mom), t_idx)] -
            ym[cbind(rep(nt_pos, n_mom), rep(t0_idx, n_mom))])

  U_p  <- matrix(0, nrow = n_mom, ncol = length(all_groups) * TT)
  rows <- seq_len(n_mom)
  U_p[cbind(rows, (g_pos  - 1L) * TT + t_idx)]  <- 1
  U_p[cbind(rows, (g_pos  - 1L) * TT + t0_idx)] <- -1
  U_p[cbind(rows, (nt_pos - 1L) * TT + t_idx)]  <-
    U_p[cbind(rows, (nt_pos - 1L) * TT + t_idx)] - 1
  U_p[cbind(rows, (nt_pos - 1L) * TT + t0_idx)] <-
    U_p[cbind(rows, (nt_pos - 1L) * TT + t0_idx)] + 1

  Omega_p <- U_p %*% it$W %*% t(U_p)
  Omega_p <- (Omega_p + t(Omega_p)) / 2

  eg  <- eigen(Omega_p, symmetric = TRUE)
  lam <- eg$values
  lam_max <- max(lam)
  if (!is.finite(lam_max) || lam_max <= 0) {
    stop("The estimated covariance of the placebo moments is not positive ",
         "definite; the test statistic cannot be formed.", call. = FALSE)
  }
  keep <- lam > 1e-8 * lam_max
  df   <- sum(keep)
  if (df < n_mom) {
    warning(sprintf(
      "The covariance of the %d placebo moment(s) is rank %d; the statistic ",
      n_mom, df),
      "uses a truncated inverse and the degrees of freedom are reduced ",
      "accordingly.", call. = FALSE)
  }

  Vk   <- eg$vectors[, keep, drop = FALSE]
  proj <- crossprod(Vk, p_hat)
  stat <- sum(proj^2 / lam[keep])
  pval <- stats::pchisq(stat, df = df, lower.tail = FALSE)

  out <- list(
    statistic   = stat,
    df          = df,
    p_value     = pval,
    type        = type,
    window      = if (identical(type, "local")) window else NA_integer_,
    base_period = t0,
    n_moments   = n_mom,
    n_dropped   = n_dropped,
    moments     = data.frame(g = mom_g, t = mom_t,
                             event_time = mom_t - mom_g,
                             moment = as.numeric(p_hat),
                             row.names = NULL),
    weighting   = object$weighting
  )
  class(out) <- "staggered_gmm_jtest"
  out
}
