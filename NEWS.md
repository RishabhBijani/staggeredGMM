# staggeredGMM 0.1.0

First release.

## Interface

* Single exported estimator `gmm_staggered()`, replacing the earlier
  `gmm_staggered_I()`, `gmm_staggered_II()` and `gmm_staggered_III()`. The
  covariance model is chosen with `weighting`, one of `"pooled_toeplitz"`,
  `"cohort_toeplitz"` or `"unrestricted"` (GMM-T, GMM-HT and GMM-U in the
  paper).
* Argument names follow the `did` package convention: `yname`, `tname`,
  `idname`, `gname`.
* `never_treated = NULL` detects a never-treated group automatically,
  replacing `has_nt = FALSE`.
* Results are an S3 object of class `staggered_gmm` with `print`, `summary`,
  `coef`, `vcov` and `confint` methods.
* New `gmm_j_test()`: the serial-correlation robust over-identification test
  of parallel trends and no anticipation, in full-set and local-pre-window
  forms.

## Estimation

* Both aggregate targets are reported: `CW` (treated-observation weighting,
  as used throughout the paper) and `EW` (cohort-equal weighting).
* Cohort-by-time effects with no clean comparison are reported as `NA` with
  `identified = FALSE`, never as an estimated zero, and are excluded from
  the parameter vector so that partial identification no longer disables the
  efficient weighting for every other cell.
* Unbalanced panels are supported: autocovariances use pairwise-complete
  denominators, and the variance of a group-period mean accounts for the
  overlap between the unit sets contributing at each pair of periods. Both
  reduce exactly to the complete-data formulas when the panel is balanced.
* Convergence and solve success are reported separately, with a termination
  reason and the final parameter change. `max_iter` now defaults to 100: the
  iteration can converge linearly at a slow rate, and `sim_panel` under
  `weighting = "cohort_toeplitz"` needs 86 steps.
* Failure to complete any weighted step always warns, including under
  `weighting = "unrestricted"`, because the returned estimates are then the
  identity-weighted seed rather than a GMM estimate under the requested
  weighting. The warning reports how many independent moment directions the
  estimated weight supported, against the number of effects to identify.
  Merely failing to converge after a weighted step did run continues not to
  warn under `"unrestricted"`, where it is routine.
* Units whose cohort falls after the last observed period are used as clean
  controls rather than discarded.

## Validation

* A shared validator checks every documented requirement before estimation
  and reports the offending rows, units, periods or cohorts. Column types
  are rejected rather than coerced.
* Units treated at or before the first observed period now raise an error.
  They have no pre-treatment period, so none of their effects is identified.

## Data

* `beck_banks` documentation now records its CC BY 4.0 licence and the
  thirteen always-treated states it contains.
