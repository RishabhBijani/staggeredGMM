## Submission

This is an update of staggeredGMM, from 0.1.0 to 0.2.0.

* Two input-validation bug fixes in `gmm_staggered()` (see NEWS.md).
* `gmm_j_test()` output and documentation now state that the test is a test
  of unconditional parallel trends.
* The third-party dataset `beck_banks` and its licence note
  (inst/LICENSE.note) have been removed. The package now bundles only
  simulated data and contains only original, MIT-licensed material.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* local macOS Sequoia 15.6 (aarch64), R 4.5.1, R CMD check --as-cran -- OK
* local Linux (x86_64), R 4.3.3, R CMD check --as-cran -- OK
* win-builder, R-release (4.6.2) -- OK
* win-builder, R-devel (2026-11-01 r90825) -- OK

## Notes

The DOI link https://doi.org/10.2139/ssrn.6558759 (README and vignette)
resolves to SSRN, which returns 403 to automated URL checkers; the link
works in a browser.

## Downstream dependencies

There are currently no downstream dependencies.
