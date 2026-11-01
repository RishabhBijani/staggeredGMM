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

* local Linux (x86_64), R 4.3.3, R CMD check --as-cran -- OK
* win-builder (R-devel, R-release) -- to be run before submission
* R-hub v2: linux, macOS and Windows (R-devel) -- to be run before submission

## Downstream dependencies

There are currently no downstream dependencies.
