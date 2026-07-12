#' Simulated staggered adoption panel
#'
#' A balanced synthetic panel used in examples, the vignette, and tests. It is
#' generated under a staggered adoption design with five treatment cohorts, a
#' never-treated group, cohort-heterogeneous dynamic treatment effects, AR(1)
#' errors, and two baseline covariates. The core data-generating process
#' (fixed effects, AR(1) errors, the cohort-heterogeneous treatment-effect
#' form) is adapted from `data-raw/sim_panel.R`'s DGP2 backbone; the
#' covariates are a package-specific addition. True CATTs are deterministic
#' given (cohort, year). The true aggregate treatment effect is ATT_CW =
#' -16.793750 (weighted by share of treated observations -- what
#' [gmm_staggered_I()] and its siblings report). For reference, the
#' equal-cohort-weighted aggregate is ATT_EW = -15.820203, though this
#' weighting is not itself returned by the package. See `data-raw/sim_panel.R`
#' for the exact data-generating process.
#'
#' @format A data frame with 1980 rows and 6 variables:
#' \describe{
#'   \item{unit_id}{Unit identifier, integer 1--60.}
#'   \item{year}{Time period, integer 1--33.}
#'   \item{cohort}{0 = never-treated; 10, 13, 16, 19, or 22 = first treated
#'     period.}
#'   \item{y}{Observed outcome.}
#'   \item{x1}{A baseline (pre-treatment) continuous covariate, for `covar`
#'     examples.}
#'   \item{x2}{A baseline (pre-treatment) binary covariate, for
#'     multiple-covariate `covar` examples.}
#' }
#' @source Simulated for Arora and Bijani, "Estimating Treatment Effects
#'   under Staggered Timing and Non-Spherical Errors".
"sim_panel"

#' State-level panel: bank deregulation and inequality
#'
#' A panel of US states, 1976--2006, used to study the relationship between
#' intrastate bank-branching deregulation and income inequality.
#'
#' @format A data frame with 1519 rows and 8 variables:
#' \describe{
#'   \item{state}{Two-letter state abbreviation.}
#'   \item{state_name}{State name.}
#'   \item{statefip}{State FIPS code.}
#'   \item{wrkyr}{Calendar year.}
#'   \item{gini}{Income Gini coefficient.}
#'   \item{branch_reform}{Year the state deregulated intrastate bank
#'     branching (1960 for states already deregulated before the sample
#'     begins).}
#'   \item{ln_gini}{Natural log of `gini`.}
#'   \item{D_branch}{Treatment indicator.}
#' }
#' @source Beck, Thorsten, Ross Levine, and Alexey Levkov. "Big Bad Banks?
#'   The Winners and Losers from Bank Deregulation in the United States."
#'   *The Journal of Finance* 65, no. 5 (2010): 1637-67.
#'   <http://www.jstor.org/stable/40864982>.
"beck_banks"

#' District-level panel: NREGS rollout and night lights in India
#'
#' A panel of Indian districts, 2000--2013, used to study the aggregate
#' effects of India's National Rural Employment Guarantee Scheme (NREGS) on
#' local economic activity, proxied by satellite night-light intensity.
#'
#' @format A data frame with 8666 rows and 20 variables. Core columns:
#' \describe{
#'   \item{district_name, state, state_name, st, censuscode, sno}{District
#'     and state identifiers.}
#'   \item{year}{Calendar year.}
#'   \item{nregs}{NREGS rollout indicator for this district-year.}
#'   \item{nr06, nr07, nr08}{Phase-of-rollout indicators (2006/2007/2008).}
#'   \item{avglt, std_lt, dlt00_05, pre_meanlt}{Night-light intensity
#'     measures (mean, standard deviation, pre-period change, and
#'     pre-period mean respectively).}
#'   \item{wage, outputwage}{Wage measures (missing for some districts).}
#'   \item{rggvy, state_frac, dep_samp}{Additional covariates from the
#'     source study (rural-electrification-program indicator, a
#'     state-level fraction, and a sample-inclusion measure respectively;
#'     some missing).}
#' }
#' Column definitions here are inferred from naming convention and are not
#' independently verified against the source paper's own variable
#' definitions table.
#' @source Cook, C. Justin, and Manisha Shah. "Aggregate Effects from Public
#'   Works: Evidence from India." *The Review of Economics and Statistics*
#'   104, no. 4 (2022): 797-806. <https://doi.org/10.1162/rest_a_00993>.
"cook_shah_nregs"
