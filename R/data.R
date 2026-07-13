#' Simulated staggered adoption panel
#'
#' A balanced synthetic panel used in examples, the vignette, and tests: 60
#' units observed over 33 periods, with five treatment cohorts (first
#' treated in periods 10, 13, 16, 19, or 22; 10 units each) and a
#' never-treated group (10 units). True CATTs are deterministic given
#' (cohort, year). The true aggregate treatment effect is ATT_CW =
#' -16.793750 (weighted by share of treated observations, which is what
#' [gmm_staggered_I()] and its siblings report).
#'
#' @details The data-generating process, for unit \eqn{i} in cohort
#'   \eqn{g_i} at period \eqn{t}:
#'   \deqn{y_{it} = \alpha_i + \lambda_t + \tau_{it} + x_{1i}\theta_{1t} +
#'   x_{2i}\theta_{2t} + \varepsilon_{it}}
#'   \itemize{
#'     \item \eqn{\alpha_i \sim N(0, 1)}: unit fixed effect.
#'     \item \eqn{\lambda_t \sim N(0, 1)}: period fixed effect.
#'     \item \eqn{\tau_{it} = \beta_{g_i}(1 + r_{g_i})^{t - g_i}} for
#'       \eqn{t \ge g_i} (0 otherwise): the cohort-heterogeneous, dynamic
#'       treatment effect, with \eqn{\beta = (-16, -12, -10, -9, -2)} and
#'       \eqn{r = (0.01, 0.04, 0.08, 0.10, 0.07)} for cohorts
#'       \eqn{(10, 13, 16, 19, 22)} respectively.
#'     \item \eqn{\varepsilon_{it}}: AR(1)-correlated across \eqn{t} within
#'       each unit, with autocorrelation \eqn{\rho = 0.5}.
#'     \item \eqn{x_{1i}}, \eqn{x_{2i}}: the baseline covariates below,
#'       each with a cohort-correlated mean and a linear-in-time loading
#'       (\eqn{\theta_{1t}}, \eqn{\theta_{2t}}), so an unconditional DiD is
#'       genuinely confounded by them.
#'   }
#'   `set.seed(312844)` is used throughout. See `data-raw/sim_panel.R` for
#'   the exact implementation.
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
#' @source Simulated for Arora, P. and Bijani, R. (2026). "Efficient
#'   Estimation of Treatment Effects under Staggered Adoption: GMM
#'   Approach." Ashoka University Economics Discussion Paper 163; and
#'   Arora, P. and Bijani, R. (2026). "Estimating Treatment Effects under
#'   Staggered Timing and Non-Spherical Errors." Available at SSRN:
#'   <https://ssrn.com/abstract=6558759>.
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
#' @format A data frame with 8666 rows and 20 variables:
#' \describe{
#'   \item{sno}{District ID.}
#'   \item{censuscode}{District census ID.}
#'   \item{state_name}{State name.}
#'   \item{district_name}{District name.}
#'   \item{year}{Calendar year.}
#'   \item{nregs}{NREGS rollout indicator for this district-year.}
#'   \item{nr06}{Wave 1 (2006) rollout indicator.}
#'   \item{nr07}{Wave 2 (2007) rollout indicator.}
#'   \item{nr08}{Wave 3 (2008) rollout indicator.}
#'   \item{avglt}{District average nighttime light intensity.}
#'   \item{std_lt}{Standard deviation of nighttime light intensity.}
#'   \item{rggvy}{RGGVY (Rajiv Gandhi Grameen Vidyutikaran Yojana rural
#'     electrification scheme), 10th plan indicator.}
#'   \item{wage}{Agricultural wage (missing for some districts).}
#'   \item{outputwage}{Output per agricultural worker (missing for some
#'     districts).}
#'   \item{state_frac}{SC/ST (Scheduled Caste/Scheduled Tribe) population
#'     fraction.}
#'   \item{dlt00_05}{Mean growth rate of nighttime lights, 2000--2005.}
#'   \item{pre_meanlt}{Mean level of nighttime lights, 2000--2005.}
#'   \item{st}{State census code (ST_CEN_CD).}
#'   \item{state}{State.}
#'   \item{dep_samp}{Marker for deposit sample (missing for some
#'     districts).}
#' }
#' Variable labels are taken directly from the source Stata file's own
#' `describe, fullnames` metadata.
#' @source Cook, C. Justin, and Manisha Shah. "Aggregate Effects from Public
#'   Works: Evidence from India." *The Review of Economics and Statistics*
#'   104, no. 4 (2022): 797-806. <https://doi.org/10.1162/rest_a_00993>.
"cook_shah_nregs"
