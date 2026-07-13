
<!-- README.md is generated from README.Rmd. Please edit that file -->

# staggeredGMM

<!-- badges: start -->

[![R-CMD-check](https://github.com/RishabhBijani/staggeredGMM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/RishabhBijani/staggeredGMM/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://app.codecov.io/gh/RishabhBijani/staggeredGMM/graph/badge.svg)](https://app.codecov.io/gh/RishabhBijani/staggeredGMM)
<!-- badges: end -->

`staggeredGMM` estimates cohort x time-specific average treatment
effects (CATTs) under staggered treatment adoption by generalized method
of moments (GMM). Implements the method proposed by Arora and Bijani
(2026). Optional baseline-covariate adjustment via outcome regression is
supported for all three.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("RishabhBijani/staggeredGMM")
```

## Usage

``` r
library(staggeredGMM)

res <- gmm_staggered_I(sim_panel, y_var = "y", t_var = "year",
                        id_var = "unit_id", gname = "cohort",
                        has_nt = TRUE)
res$ATT_CW
#> [1] -16.63214
res$SE_CW
#> [1] 0.2521646
head(res$catt_out)
#>    g  t  beta_hat        SE
#> 1 10 10 -16.01541 0.3200093
#> 2 10 11 -15.67617 0.3786223
#> 3 10 12 -15.98558 0.4040038
#> 4 10 13 -16.42134 0.4216487
#> 5 10 14 -16.48369 0.4310857
#> 6 10 15 -16.92851 0.4391031
```

## The Three Estimators

| Function              | Weighting                                                            |
|-----------------------|----------------------------------------------------------------------|
| `gmm_staggered_I()`   | Pooled Toeplitz, one autocovariance sequence shared by every cohort  |
| `gmm_staggered_II()`  | Cohort-specific Toeplitz, each cohort estimates its own sequence     |
| `gmm_staggered_III()` | Fully unrestricted within-cohort covariance, no stationarity assumed |

## Baseline-Covariate Adjustment

All three estimators accept an optional `covar` argument: a character
vector of baseline (pre-treatment) covariate columns, applying the
outcome-regression adjustment of Section 4.5 of the paper. More than one
covariate can be supplied at once:

``` r
res_cov <- gmm_staggered_I(sim_panel, y_var = "y", t_var = "year",
                           id_var = "unit_id", gname = "cohort",
                           has_nt = TRUE, covar = c("x1", "x2"))
res_cov$ATT_CW
#> [1] -16.66011
```

If a supplied covariate is not actually constant over time for some
units, `staggeredGMM` does not error out: it uses each affected unit’s
first-observed value and raises a warning naming exactly which units
were affected.

## Other Bundled Datasets

Two real-world datasets ship with the package for reference:
`beck_banks`, data from Beck et al. (2010), and `cook_shah_nregs`, data
from Cook and Shah (2022). See `?beck_banks` and `?cook_shah_nregs` for
details.

## Citation

If you use this package, please cite both underlying papers:

> Arora, P. and Bijani, R. (2026). “Efficient Estimation of Treatment
> Effects under Staggered Adoption: GMM Approach.” Ashoka University
> Economics Discussion Paper 163.
> <https://dp.ashoka.edu.in/ash/wpaper/paper163_0.pdf>

> Arora, P. and Bijani, R. (2026). “Estimating Treatment Effects under
> Staggered Timing and Non-Spherical Errors.” Available at SSRN:
> <https://ssrn.com/abstract=6558759> or
> <http://dx.doi.org/10.2139/ssrn.6558759>.

See `citation("staggeredGMM")` or `CITATION.cff` for the full citations.

## License

MIT © Rishabh Bijani, Parush Arora
