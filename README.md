
<!-- README.md is generated from README.Rmd. Please edit that file -->

# staggeredGMM

<!-- badges: start -->

[![R-CMD-check](https://github.com/RishabhBijani/staggeredGMM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/RishabhBijani/staggeredGMM/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://app.codecov.io/gh/RishabhBijani/staggeredGMM/graph/badge.svg)](https://app.codecov.io/gh/RishabhBijani/staggeredGMM)
<!-- badges: end -->

`staggeredGMM` estimates cohort-by-time average treatment effects
(CATTs) under staggered treatment adoption by the generalized method of
moments, implementing the estimator of Arora and Bijani (2026). Only
clean two-by-two difference-in-differences comparisons – against
never-treated or not-yet-treated controls – enter the moment system.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("RishabhBijani/staggeredGMM")
```

## Usage

``` r
library(staggeredGMM)

fit <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                     idname = "unit_id", gname = "cohort")
fit
#> Staggered-adoption GMM estimator
#> Weighting: pooled Toeplitz (GMM-T)
#> 
#>   Units                       60
#>   Periods                     33
#>   Treated cohorts              5
#>   CATTs (identified)     90 (90)
#>   Clean comparisons         1935
#>   Moment-space rank          160
#>   Iterations                   7
#>   Converged                  yes
#> 
#> Treated-observation ATT       -16.6321  (se 0.2522) 
#> Cohort-equal ATT              -15.5701  (se 0.2521)
```

``` r
head(fit$catt)
#>    g  t event_time  estimate std_error identified n_comparisons
#> 1 10 10          0 -16.01541 0.3200093       TRUE            45
#> 2 10 11          1 -15.67617 0.3786223       TRUE            45
#> 3 10 12          2 -15.98558 0.4040038       TRUE            45
#> 4 10 13          3 -16.42134 0.4216487       TRUE            36
#> 5 10 14          4 -16.48369 0.4310857       TRUE            36
#> 6 10 15          5 -16.92851 0.4391031       TRUE            36
```

## Weighting

The covariance model used to form the optimal weight is chosen with
`weighting`:

| `weighting` | Covariance model | Paper |
|----|----|----|
| `"pooled_toeplitz"` (default) | One stationary autocovariance sequence shared by every group | GMM-T |
| `"cohort_toeplitz"` | A separate stationary sequence per cohort | GMM-HT |
| `"unrestricted"` | Full unrestricted within-cohort covariance | GMM-U |

All three target the same CATT vector from the same moment conditions.
They do not generally return the same numbers: the system is
over-identified, so a different weighting matrix gives a different
estimate.

`"unrestricted"` needs the panel to support it. Each cohort’s covariance
is built from that cohort’s residual vectors, so its rank cannot exceed
`min(N_g, T - 1)`, and the weight’s rank is bounded by that summed over
groups. Unless the total comfortably exceeds the number of effects, no
weighted step completes and the estimator warns that it has fallen back
to the identity-weighted seed.

## Testing the identifying assumptions

``` r
gmm_j_test(fit)
#> Specification test for parallel trends and no anticipation
#> Hansen J on the pre-treatment placebo restrictions
#> 
#>   Restriction set   local pre-window (3 period(s) before adoption)
#>   Base period       9
#>   Weighting         pooled Toeplitz (GMM-T)
#>   Moments           14
#> 
#>   J = 9.0396, df = 14, p = 0.8285
#> 
#> A clean result is necessary but not sufficient: post-treatment
#> parallel trends is untestable.
```

## Baseline covariates

`covar` applies the outcome-regression adjustment of Section 4.5 of the
paper, so that parallel trends need hold only conditional on the named
baseline covariates:

``` r
fit_cov <- gmm_staggered(sim_panel, yname = "y", tname = "year",
                         idname = "unit_id", gname = "cohort",
                         covar = c("x1", "x2"))
fit_cov$aggregate$CW$estimate
#> [1] -16.66011
```

## Bundled data

`sim_panel` is a simulated panel used in the examples and tests.
`beck_banks` is a real state-level panel from Beck, Levine and Levkov
(2010), included under CC BY 4.0; see `?beck_banks` and
`inst/LICENSE.note`. It contains thirteen always-treated states and is
documented as a worked example of that pitfall.

The NREGS district panel of Cook and Shah (2022) is not bundled;
`inst/scripts/get_cook_shah_nregs.R` retrieves and prepares it on
request.

## Citation

``` r
citation("staggeredGMM")
```

> Arora, P. and Bijani, R. (2026). “Estimating Treatment Effects under
> Staggered Timing and Non-Spherical Errors.” Available at SSRN:
> <https://doi.org/10.2139/ssrn.6558759>

## License

MIT (c) Rishabh Bijani, Parush Arora. One bundled dataset carries its
own licence; see `inst/LICENSE.note`.
