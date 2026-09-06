## Submission

First submission of staggeredGMM.

## R CMD check results

0 errors | 0 warnings | 1 note

The note is the standard one for a first submission:

    Maintainer: 'Rishabh Bijani <rishabhbijani@gmail.com>'
    New submission

## Test environments

* local macOS 15.6, aarch64, R 4.5.1
* GitHub Actions: ubuntu-latest (release, devel, oldrel-1), macOS-latest
  (release), windows-latest (release)
* win-builder, R devel
* R-hub

## Notes on bundled data

The package bundles one third-party dataset, `beck_banks`, derived from the
replication data of Beck, Levine and Levkov (2010) deposited at DataverseNL
(doi:10.34894/B1K9OU) under the Creative Commons Attribution 4.0
International licence, which permits redistribution with attribution. The
only change from the deposited data is a selection of columns; no values were
altered, recoded or derived. The attribution, the licence and the change made
are recorded in inst/LICENSE.note and on the dataset's help page. All other
package content is original and MIT-licensed.

## Downstream dependencies

There are currently no downstream dependencies.
