# Input validation and canonical panel preparation.
#
# Design principle: reject, never coerce. Any condition that would silently
# change the user's data is an error. Conditions that only cost information
# are warnings. Every estimator entry point routes through
# validate_staggered_panel() and consumes its return value, so the three
# weighting schemes cannot drift apart in how they read the data.


# TRUE for a numeric vector all of whose finite entries are whole numbers.
# Deliberately FALSE for logical, character and factor: those are rejected
# rather than coerced.
.is_integerish <- function(x) {
  if (!is.numeric(x)) return(FALSE)
  fin <- is.finite(x)
  if (!any(fin)) return(FALSE)
  all(x[fin] == trunc(x[fin]))
}


.describe_class <- function(x) {
  paste(class(x), collapse = "/")
}


# Truncated, comma-separated preview of a vector, for error messages.
.preview <- function(x, n = 5L) {
  x <- unique(x)
  shown <- utils::head(x, n)
  out <- paste(format(shown, trim = TRUE), collapse = ", ")
  if (length(x) > n) out <- paste0(out, ", ... (", length(x), " total)")
  out
}


.check_name_arg <- function(value, arg) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value)) {
    stop(sprintf("`%s` must be a single non-missing column name.", arg),
         call. = FALSE)
  }
  invisible(TRUE)
}


#' Validate and canonicalise a staggered-adoption panel
#'
#' Internal workhorse shared by [gmm_staggered()] and [gmm_j_test()]. Checks
#' every documented input requirement, then returns a canonical
#' representation of the panel that downstream code consumes without further
#' checking.
#'
#' @param data A data frame with one row per unit-period.
#' @param yname,tname,idname,gname Column names, as in [gmm_staggered()].
#' @param covar Optional character vector of baseline covariate columns.
#' @param never_treated `NULL`, `TRUE` or `FALSE`, as in [gmm_staggered()].
#'
#' @return A list with the canonical panel and its metadata. See the source
#'   for the full set of components.
#'
#' @noRd
validate_staggered_panel <- function(data, yname, tname, idname, gname,
                                     covar = NULL, never_treated = NULL) {

  ## ---- 1. argument types --------------------------------------------------
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame (or an object inheriting from one); ",
         sprintf("got an object of class %s.", .describe_class(data)),
         call. = FALSE)
  }
  data <- as.data.frame(data, stringsAsFactors = FALSE)

  .check_name_arg(yname,  "yname")
  .check_name_arg(tname,  "tname")
  .check_name_arg(idname, "idname")
  .check_name_arg(gname,  "gname")

  if (!is.null(never_treated)) {
    if (!is.logical(never_treated) || length(never_treated) != 1L ||
        is.na(never_treated)) {
      stop("`never_treated` must be NULL, TRUE or FALSE.", call. = FALSE)
    }
  }

  if (!is.null(covar)) {
    if (!is.character(covar)) {
      stop("`covar` must be a character vector of column names; ",
           sprintf("got an object of class %s.", .describe_class(covar)),
           call. = FALSE)
    }
    covar <- covar[nzchar(covar) & !is.na(covar)]
    if (length(covar) == 0L) covar <- NULL
  }

  ## ---- 2. columns exist ---------------------------------------------------
  needed <- c(yname, tname, idname, gname, covar)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0L) {
    stop(sprintf("Column(s) not found in `data`: %s.",
                 .preview(missing_cols, 10L)), call. = FALSE)
  }
  if (any(duplicated(c(yname, tname, idname, gname)))) {
    stop("`yname`, `tname`, `idname` and `gname` must name four distinct ",
         "columns.", call. = FALSE)
  }
  if (!is.null(covar) && length(intersect(covar,
                                          c(yname, tname, idname, gname)))) {
    stop("`covar` must not name the outcome, time, id or cohort column.",
         call. = FALSE)
  }

  n_rows <- nrow(data)
  if (n_rows == 0L) stop("`data` has no rows.", call. = FALSE)

  y_raw  <- data[[yname]]
  t_raw  <- data[[tname]]
  id_raw <- data[[idname]]
  g_raw  <- data[[gname]]

  ## ---- 3. column types ----------------------------------------------------
  if (is.factor(y_raw) || is.character(y_raw) || is.logical(y_raw) ||
      !is.numeric(y_raw)) {
    stop(sprintf("`%s` (the outcome) must be numeric; got %s. ", yname,
                 .describe_class(y_raw)),
         "Convert it explicitly rather than relying on automatic coercion, ",
         "which turns a factor into its integer level codes.", call. = FALSE)
  }
  y_raw <- as.numeric(y_raw)

  if (!.is_integerish(t_raw)) {
    stop(sprintf("`%s` (time) must be an integer-valued numeric column; got %s.",
                 tname, .describe_class(t_raw)),
         if (inherits(t_raw, "Date"))
           " Convert dates to a consecutive integer period index first."
         else
           " Non-integer values would be silently truncated.",
         call. = FALSE)
  }
  if (!.is_integerish(g_raw)) {
    stop(sprintf("`%s` (cohort) must be an integer-valued numeric column; got %s.",
                 gname, .describe_class(g_raw)),
         " Code never-treated units as 0.", call. = FALSE)
  }
  if (is.list(id_raw)) {
    stop(sprintf("`%s` (unit id) must be an atomic vector.", idname),
         call. = FALSE)
  }

  ## ---- 4. missing keys ----------------------------------------------------
  for (nm in list(c(idname, "unit id"), c(tname, "time"), c(gname, "cohort"))) {
    col <- data[[nm[1L]]]
    if (anyNA(col)) {
      stop(sprintf("`%s` (%s) contains %d missing value(s); rows: %s.",
                   nm[1L], nm[2L], sum(is.na(col)),
                   .preview(which(is.na(col)))),
           call. = FALSE)
    }
  }
  # Missing outcomes are permitted: they are the unbalanced-panel case.

  t_int <- as.integer(t_raw)
  g_int <- as.integer(g_raw)

  ## ---- 5. unit-period uniqueness ------------------------------------------
  key <- paste(as.character(id_raw), t_int, sep = "\r")
  if (anyDuplicated(key)) {
    dup <- which(duplicated(key))
    stop("Each unit-period combination must appear exactly once; ",
         sprintf("found %d duplicate row(s). Offending rows: %s.",
                 length(dup), .preview(dup)),
         call. = FALSE)
  }

  ## ---- 6. canonical ordering and index maps -------------------------------
  unit_ids <- sort(unique(id_raw))
  N_units  <- length(unit_ids)
  unit_row <- match(id_raw, unit_ids)

  T_min <- min(t_int)
  T_max <- max(t_int)
  TT    <- T_max - T_min + 1L
  if (TT < 2L) {
    stop("At least two distinct time periods are required; `", tname,
         "` takes only one value.", call. = FALSE)
  }
  time_idx <- t_int - T_min + 1L

  ord      <- order(unit_row, time_idx)
  unit_row <- unit_row[ord]
  time_idx <- time_idx[ord]
  y_ord    <- y_raw[ord]
  g_ord    <- g_int[ord]

  ## ---- 7. cohort constant within unit -------------------------------------
  n_distinct_g <- tapply(g_ord, unit_row, function(z) length(unique(z)))
  bad_units    <- which(as.integer(n_distinct_g) != 1L)
  if (length(bad_units) > 0L) {
    stop(sprintf("`%s` (cohort) must be constant within unit; ", gname),
         sprintf("%d unit(s) have more than one cohort value. Affected: %s.",
                 length(bad_units), .preview(unit_ids[bad_units])),
         call. = FALSE)
  }
  cohort_of <- as.integer(tapply(g_ord, unit_row, function(z) z[1L]))

  ## ---- 8. consecutive time support ----------------------------------------
  present <- sort(unique(t_int))
  if (!identical(present, seq.int(T_min, T_max))) {
    gaps <- setdiff(seq.int(T_min, T_max), present)
    stop("The set of periods present in `data` must be consecutive: no ",
         "period may be absent for every unit. ",
         sprintf("Period(s) absent from the whole panel: %s. ", .preview(gaps)),
         "The stationary weighting indexes autocovariances by calendar ",
         "distance, so an empty period would create lag distances that ",
         "correspond to nothing in the data. Individual units may be ",
         "missing individual periods; that is the unbalanced case and is ",
         "supported.", call. = FALSE)
  }

  ## ---- 9. cohort value ranges ---------------------------------------------
  ucoh <- sort(unique(cohort_of))

  neg <- ucoh[ucoh < 0L]
  if (length(neg) > 0L) {
    stop(sprintf("`%s` (cohort) contains negative value(s): %s. ", gname,
                 .preview(neg)),
         "Cohorts must be 0 (never treated) or a treated period.",
         call. = FALSE)
  }

  always <- ucoh[ucoh > 0L & ucoh <= T_min]
  if (length(always) > 0L) {
    n_always <- sum(cohort_of %in% always)
    stop(sprintf("%d unit(s) in cohort(s) %s are treated at or before the ",
                 n_always, .preview(always)),
         sprintf("first observed period (%d). ", T_min),
         "Such units have no pre-treatment period, so none of their ",
         "cohort-by-time effects is identified, and including them would ",
         "silently contaminate the aggregate. Drop them, extend the panel ",
         "backwards, or recode them as never-treated (0) only if that is ",
         "substantively correct.", call. = FALSE)
  }

  treated_cohorts <- ucoh[ucoh > T_min & ucoh <= T_max]
  late_cohorts    <- ucoh[ucoh > T_max]

  if (length(treated_cohorts) == 0L) {
    stop("No treated cohort is identified: `", gname, "` contains no value ",
         sprintf("strictly inside the observed period range (%d, %d].",
                 T_min, T_max), call. = FALSE)
  }

  ## ---- 10. control design -------------------------------------------------
  # Units never treated within the observation window: cohort 0, plus any
  # cohort adopting after the panel ends. Both are untreated at every
  # observed period, so both are valid clean controls throughout.
  is_never <- cohort_of == 0L | cohort_of > T_max
  n_never_zero <- sum(cohort_of == 0L)
  n_never_late <- sum(cohort_of > T_max)
  n_never      <- n_never_zero + n_never_late

  if (isTRUE(never_treated) && n_never == 0L) {
    stop("`never_treated = TRUE` was requested, but no unit is never treated ",
         sprintf("within the observed window: no unit has `%s == 0` or `%s > %d`.",
                 gname, gname, T_max),
         " Use `never_treated = NULL` to detect automatically, or FALSE to ",
         "rely on not-yet-treated controls alone.", call. = FALSE)
  }

  has_never <- if (is.null(never_treated)) n_never > 0L else never_treated
  if (identical(never_treated, FALSE) && n_never > 0L) {
    warning(sprintf(
      "`never_treated = FALSE`: %d never-treated unit(s) are present but will ",
      n_never),
      "not be used as controls. Precision will be lower and some ",
      "cohort-by-time effects may become unidentified.", call. = FALSE)
  }

  # Internal group label: 0 for the never-treated pool, cohort year otherwise.
  grp_of <- cohort_of
  grp_of[is_never] <- 0L
  all_groups <- if (has_never) c(0L, treated_cohorts) else treated_cohorts
  # Units in the never pool are only carried when they serve as controls.
  in_design <- if (has_never) rep(TRUE, N_units) else !is_never

  N_g <- setNames(
    vapply(all_groups, function(g) sum(grp_of == g & in_design), integer(1L)),
    as.character(all_groups)
  )
  empty_groups <- names(N_g)[N_g == 0L]
  if (length(empty_groups) > 0L) {
    stop("Internal inconsistency: group(s) ", .preview(empty_groups),
         " contain no units. Please report this with a reproducible example.",
         call. = FALSE)
  }

  ## ---- 11. outcome matrix and observation pattern -------------------------
  Y_mat <- matrix(NA_real_, nrow = N_units, ncol = TT)
  Y_mat[cbind(unit_row, time_idx)] <- y_ord
  obs_mat <- !is.na(Y_mat)

  n_missing_cells <- sum(!obs_mat)
  balanced <- n_missing_cells == 0L

  # Observed unit counts per group-period; a zero means that cohort-period
  # mean does not exist and every comparison touching it must be dropped.
  n_grp <- length(all_groups)
  N_gt  <- matrix(0L, nrow = n_grp, ncol = TT,
                  dimnames = list(as.character(all_groups), NULL))
  for (gi in seq_len(n_grp)) {
    rows <- which(grp_of == all_groups[gi] & in_design)
    if (length(rows) > 0L)
      N_gt[gi, ] <- as.integer(colSums(obs_mat[rows, , drop = FALSE]))
  }

  ## ---- 12. baseline covariates --------------------------------------------
  X_unit <- NULL
  if (!is.null(covar)) {
    X_unit <- build_baseline_X(data[ord, , drop = FALSE],
                               stats::reformulate(covar),
                               unit_ids, unit_row, time_idx)
  }

  list(
    unit_ids        = unit_ids,
    N_units         = N_units,
    unit_row        = unit_row,
    time_idx        = time_idx,
    y               = y_ord,
    T_min           = T_min,
    T_max           = T_max,
    TT              = TT,
    cohort_of       = cohort_of,
    grp_of          = grp_of,
    in_design       = in_design,
    treated_cohorts = treated_cohorts,
    late_cohorts    = late_cohorts,
    all_groups      = all_groups,
    has_never       = has_never,
    n_never_zero    = n_never_zero,
    n_never_late    = n_never_late,
    N_g             = N_g,
    Y_mat           = Y_mat,
    obs_mat         = obs_mat,
    N_gt            = N_gt,
    balanced        = balanced,
    n_missing_cells = n_missing_cells,
    X_unit          = X_unit,
    covar           = covar,
    names           = list(yname = yname, tname = tname,
                           idname = idname, gname = gname)
  )
}
