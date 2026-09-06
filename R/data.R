#' Simulated staggered adoption panel
#'
#' A balanced synthetic panel used in the examples, the vignette and the
#' tests: 60 units observed over 33 periods, with five treatment cohorts
#' (first treated in periods 10, 13, 16, 19 or 22; 10 units each) and a
#' never-treated group of 10 units. True effects are deterministic given
#' (cohort, period), so both aggregate targets are known exactly:
#' `ATT_CW = -16.79375` under treated-observation weighting and
#' `ATT_EW = -15.82020` under cohort-equal weighting. The two differ because
#' earlier cohorts are observed for more post-treatment periods; see
#' [gmm_staggered()] for the definitions.
#'
#' @details The data-generating process, for unit \eqn{i} in cohort
#'   \eqn{g_i} at period \eqn{t}:
#'   \deqn{y_{it} = \alpha_i + \lambda_t + \tau_{it} + x_{1i}\theta_{1t} +
#'   x_{2i}\theta_{2t} + \varepsilon_{it}}
#'   \itemize{
#'     \item \eqn{\alpha_i \sim N(0, 1)}: unit fixed effect.
#'     \item \eqn{\lambda_t \sim N(0, 1)}: period fixed effect.
#'     \item \eqn{\tau_{it} = \beta_{g_i}(1 + r_{g_i})^{t - g_i}} for
#'       \eqn{t \ge g_i} and 0 otherwise, with
#'       \eqn{\beta = (-16, -12, -10, -9, -2)} and
#'       \eqn{r = (0.01, 0.04, 0.08, 0.10, 0.07)} for cohorts
#'       \eqn{(10, 13, 16, 19, 22)}.
#'     \item \eqn{\varepsilon_{it}}: AR(1) within unit, \eqn{\rho = 0.5}.
#'     \item \eqn{x_{1i}}, \eqn{x_{2i}}: baseline covariates with
#'       cohort-correlated means and linear-in-time loadings
#'       (\eqn{\theta_{1t}}, \eqn{\theta_{2t}}), so an unconditional
#'       difference-in-differences is genuinely confounded by them.
#'   }
#'   `set.seed(312844)` throughout. See `data-raw/sim_panel.R` for the exact
#'   implementation.
#'
#'   This is the design of Appendix E of the paper, so the specification test
#'   has 70 full-set and 14 local-window restrictions on it.
#'
#' @format A data frame with 1980 rows and 6 variables:
#' \describe{
#'   \item{unit_id}{Unit identifier, integer 1--60.}
#'   \item{year}{Period, integer 1--33.}
#'   \item{cohort}{0 = never treated; 10, 13, 16, 19 or 22 = first treated
#'     period.}
#'   \item{y}{Observed outcome.}
#'   \item{x1}{Continuous baseline covariate.}
#'   \item{x2}{Binary baseline covariate.}
#' }
#' @source Simulated for Arora, P. and Bijani, R. (2026). "Estimating
#'   Treatment Effects under Staggered Timing and Non-Spherical Errors."
#'   \doi{10.2139/ssrn.6558759}
"sim_panel"

#' State-level panel: bank branch deregulation and income inequality
#'
#' A panel of 49 US states over 1976--2006, used to study intrastate bank
#' branching deregulation and the distribution of income.
#'
#' @section Always-treated states: This dataset is a worked example of a
#'   coding pitfall as much as of a design. Thirteen of the 49 states
#'   deregulated at or before 1976, the first observed year: AK, AZ, CA, DC,
#'   ID, MD, ME, NC, NV, NY, RI, SC and VT, ten of them carrying
#'   `branch_reform == 1960`. These units have no pre-treatment period, so
#'   none of their effects is identified. No state is coded `0`, so the
#'   dataset also contains no never-treated group.
#'
#'   Passed to [gmm_staggered()] as it stands, the always-treated states
#'   raise an error naming them. That is deliberate: silently including them
#'   would return a plausible-looking aggregate contaminated by cells that
#'   are not estimated at all. Restrict to states adopting inside the
#'   window before estimating:
#'
#'   \preformatted{
#'   dat <- beck_banks[beck_banks$branch_reform > 1976, ]
#'   fit <- gmm_staggered(dat, yname = "ln_gini", tname = "wrkyr",
#'                        idname = "state", gname = "branch_reform")
#'   }
#'
#'   This leaves 36 states across cohorts 1977--1999, identified from
#'   not-yet-treated controls alone.
#'
#' @format A data frame with 1519 rows and 8 variables:
#' \describe{
#'   \item{state}{Two-letter state abbreviation.}
#'   \item{state_name}{State name.}
#'   \item{statefip}{State FIPS code.}
#'   \item{wrkyr}{Calendar year, 1976--2006.}
#'   \item{gini}{Income Gini coefficient.}
#'   \item{branch_reform}{Year the state deregulated intrastate bank
#'     branching. `1960` marks states already deregulated before the sample
#'     begins.}
#'   \item{ln_gini}{Natural log of `gini`.}
#'   \item{D_branch}{Treatment indicator.}
#' }
#'
#' @source Beck, T., Levine, R. and Levkov, A. (2010). "Big Bad Banks? The
#'   Winners and Losers from Bank Deregulation in the United States."
#'   \emph{The Journal of Finance} 65(5), 1637--1667.
#'   \doi{10.1111/j.1540-6261.2010.01589.x}
#'
#'   Replication data deposited by the authors at DataverseNL,
#'   \doi{10.34894/B1K9OU}, and made available under the Creative Commons
#'   Attribution 4.0 International licence,
#'   <https://creativecommons.org/licenses/by/4.0/>.
#'
#'   Modified from the deposited data: a subset of columns was selected and
#'   the remainder discarded. No values were altered, recoded or derived. The
#'   full licence note is installed with the package; see
#'   `system.file("LICENSE.note", package = "staggeredGMM")`.
"beck_banks"
