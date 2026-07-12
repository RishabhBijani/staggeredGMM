#' staggeredGMM: GMM Estimation of Treatment Effects under Staggered Adoption
#'
#' Implements the GMM estimator of Arora and Bijani, "Estimating Treatment
#' Effects under Staggered Timing and Non-Spherical Errors", for
#' cohort-and-time-specific average treatment effects (CATTs) under staggered
#' treatment adoption. Three weighting variants are provided:
#' [gmm_staggered_I()] (pooled Toeplitz / "Eff"), [gmm_staggered_II()]
#' (cohort-specific Toeplitz / "HetVar"), and [gmm_staggered_III()] (fully
#' unrestricted within-cohort covariance / "FullCov"). All three support
#' optional baseline-covariate adjustment via `covar`.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom stats complete.cases lm.fit model.matrix residuals setNames toeplitz
#' @importFrom fixest feols
#' @importFrom MASS ginv
#' @import data.table
NULL

utils::globalVariables(c("unit", "time", "Y", "Y_adj", "tau_hat"))
