# S3 methods for the fitted objects.
#
# Each method carries its own documentation block rather than @describeIn,
# because their formal arguments differ from those of gmm_staggered() and
# roxygen would otherwise attach `x`, `parm` and `level` to that function's
# usage section, which R CMD check flags as documented-but-not-in-usage.

.weighting_label <- function(x) {
  switch(x,
         pooled_toeplitz = "pooled Toeplitz (GMM-T)",
         cohort_toeplitz = "cohort-specific Toeplitz (GMM-HT)",
         unrestricted    = "unrestricted within-cohort (GMM-U)",
         x)
}


.fmt_agg <- function(a, label) {
  if (is.na(a$estimate)) {
    sprintf("%-28s      NA  (identified subset: %.4f, se %.4f, %.1f%% of weight)",
            label, a$estimate_identified, a$std_error_identified,
            100 * a$identified_weight_share)
  } else {
    sprintf("%-28s %9.4f  (se %.4f)", label, a$estimate, a$std_error)
  }
}


#' Print a fitted staggered GMM model
#'
#' @param x A `staggered_gmm` object, as returned by [gmm_staggered()].
#' @param ... Ignored.
#' @return `x`, invisibly. Called for the side effect of printing.
#' @examples
#' fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
#'                      idname = "unit_id", gname = "cohort")
#' print(fit)
#' @export
print.staggered_gmm <- function(x, ...) {
  d <- x$design
  cat("Staggered-adoption GMM estimator\n")
  cat("Weighting: ", .weighting_label(x$weighting), "\n\n", sep = "")

  cat(sprintf("  Units                 %8d\n", d$N_units))
  cat(sprintf("  Periods               %8d\n", d$T_periods))
  cat(sprintf("  Treated cohorts       %8d\n", d$n_cohorts))
  cat(sprintf("  CATTs (identified)    %8s\n",
              sprintf("%d (%d)", d$N_beta, d$N_identified)))
  cat(sprintf("  Clean comparisons     %8d\n", d$N_2x2))
  cat(sprintf("  Moment-space rank     %8d\n", d$moment_rank))
  cat(sprintf("  Iterations            %8d\n", x$convergence$n_iter))
  cat(sprintf("  Converged             %8s\n",
              if (x$convergence$converged) "yes" else "no"))
  if (!x$convergence$converged)
    cat(sprintf("  Termination           %8s\n", x$convergence$termination))
  if (!x$validation$balanced)
    cat(sprintf("  Unbalanced            %8s\n",
                sprintf("%d missing cell(s)", x$validation$n_missing_cells)))
  cat("\n")
  cat(.fmt_agg(x$aggregate$CW, "Treated-observation ATT"), "\n")
  cat(.fmt_agg(x$aggregate$EW, "Cohort-equal ATT"), "\n")
  invisible(x)
}


#' Summarise a fitted staggered GMM model
#'
#' Prints everything [print.staggered_gmm()] shows, followed by the full
#' table of cohort-by-time effects and the covariate-adjustment status.
#'
#' @param object A `staggered_gmm` object, as returned by [gmm_staggered()].
#' @param ... Ignored.
#' @return `object`, invisibly. Called for the side effect of printing.
#' @examples
#' fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
#'                      idname = "unit_id", gname = "cohort")
#' summary(fit)
#' @export
summary.staggered_gmm <- function(object, ...) {
  print(object)
  cat("\nCohort-by-time effects:\n")
  tab <- object$catt
  tab$estimate  <- round(tab$estimate, 4)
  tab$std_error <- round(tab$std_error, 4)
  print(tab, row.names = FALSE)
  if (!is.null(object$controls$covar)) {
    cat("\nBaseline covariates: ",
        paste(object$controls$covar, collapse = ", "), "\n", sep = "")
    if (object$controls$cov_n_fallback > 0L)
      cat(sprintf("  %d comparison(s) fell back to the unconditional form.\n",
                  object$controls$cov_n_fallback))
  }
  invisible(object)
}


#' Extract cohort-by-time effect estimates
#'
#' @param object A `staggered_gmm` object, as returned by [gmm_staggered()].
#' @param ... Ignored.
#' @return A named numeric vector of cohort-by-time effects, one per cell,
#'   named `"g<cohort>:t<period>"`. Cells with no clean comparison are `NA`.
#' @examples
#' fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
#'                      idname = "unit_id", gname = "cohort")
#' head(coef(fit))
#' @export
coef.staggered_gmm <- function(object, ...) {
  object$coefficients
}


#' Covariance matrix of the estimated cohort-by-time effects
#'
#' @param object A `staggered_gmm` object, as returned by [gmm_staggered()].
#' @param ... Ignored.
#' @return A square numeric matrix covering the identified cells only, with
#'   row and column names matching the corresponding entries of `coef()`.
#' @examples
#' fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
#'                      idname = "unit_id", gname = "cohort")
#' dim(vcov(fit))
#' @export
vcov.staggered_gmm <- function(object, ...) {
  object$vcov
}


#' Confidence intervals for cohort-by-time effects
#'
#' Pointwise normal-approximation intervals. They are not adjusted for
#' multiplicity, so reading a whole event-study path off them overstates
#' confidence.
#'
#' @param object A `staggered_gmm` object, as returned by [gmm_staggered()].
#' @param parm Optional character vector of coefficient names. Defaults to
#'   all of them.
#' @param level Confidence level, strictly between 0 and 1.
#' @param ... Ignored.
#' @return A two-column matrix of lower and upper limits, with one row per
#'   requested coefficient. Unidentified cells are `NA`.
#' @examples
#' fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
#'                      idname = "unit_id", gname = "cohort")
#' head(confint(fit))
#' @export
confint.staggered_gmm <- function(object, parm, level = 0.95, ...) {
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  est <- object$coefficients
  se  <- stats::setNames(rep(NA_real_, length(est)), names(est))
  se[rownames(object$vcov)] <- sqrt(pmax(0, diag(object$vcov)))

  if (!missing(parm) && !is.null(parm)) {
    unknown <- setdiff(parm, names(est))
    if (length(unknown) > 0L)
      stop("Unknown coefficient name(s): ", .preview(unknown), ".",
           call. = FALSE)
    est <- est[parm]
    se  <- se[parm]
  }

  a <- (1 - level) / 2
  z <- stats::qnorm(1 - a)
  out <- cbind(est - z * se, est + z * se)
  colnames(out) <- sprintf("%.1f %%", 100 * c(a, 1 - a))
  rownames(out) <- names(est)
  out
}


#' Print a specification test result
#'
#' @param x A `staggered_gmm_jtest` object, as returned by [gmm_j_test()].
#' @param ... Ignored.
#' @return `x`, invisibly. Called for the side effect of printing.
#' @examples
#' fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
#'                      idname = "unit_id", gname = "cohort")
#' print(gmm_j_test(fit))
#' @export
print.staggered_gmm_jtest <- function(x, ...) {
  cat("Specification test for parallel trends and no anticipation\n")
  cat("Hansen J on the pre-treatment placebo restrictions\n\n")
  cat(sprintf("  Restriction set   %s\n",
              if (identical(x$type, "local"))
                sprintf("local pre-window (%d period(s) before adoption)",
                        x$window)
              else "all pre-treatment placebo moments"))
  cat(sprintf("  Base period       %d\n", x$base_period))
  cat(sprintf("  Weighting         %s\n", .weighting_label(x$weighting)))
  cat(sprintf("  Moments           %d\n", x$n_moments))
  cat("\n")
  cat(sprintf("  J = %.4f, df = %d, p = %.4f\n",
              x$statistic, x$df, x$p_value))
  invisible(x)
}
