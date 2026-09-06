#' staggeredGMM: GMM Estimation of Treatment Effects Under Staggered Adoption
#'
#' Estimates cohort-by-time average treatment effects (CATTs) under staggered
#' treatment adoption by the generalized method of moments, following Arora
#' and Bijani (2026).
#'
#' The single estimation entry point is [gmm_staggered()]. The
#' over-identification test of parallel trends and no anticipation is
#' [gmm_j_test()].
#'
#' @section Data coding conventions:
#' `gmm_staggered()` expects one row per unit-period. The cohort column gives
#' each unit's first treated period, with `0` reserved for units that are
#' never treated. Time must be integer-valued, and the set of periods present
#' in the data (pooled across units) must be consecutive. Individual units may
#' be missing individual periods; see `vignette("staggeredGMM")` for the
#' unbalanced-panel policy. Full details, including every condition that
#' raises an error, are in [gmm_staggered()].
#'
#' @references
#' Arora, P. and Bijani, R. (2026). "Estimating Treatment Effects under
#' Staggered Timing and Non-Spherical Errors."
#' \doi{10.2139/ssrn.6558759}
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom stats complete.cases lm.fit model.frame model.matrix na.pass
#' @importFrom stats pchisq qnorm reformulate residuals setNames toeplitz
#' @importFrom fixest feols
#' @importFrom MASS ginv
#' @importFrom utils head
NULL
