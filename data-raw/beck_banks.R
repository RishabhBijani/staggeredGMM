# Preparation of the bundled `beck_banks` dataset.
#
# Source: replication data for Beck, T., Levine, R. and Levkov, A. (2010),
# "Big Bad Banks? The Winners and Losers from Bank Deregulation in the
# United States", The Journal of Finance 65(5), 1637-1667, deposited by the
# authors at DataverseNL, doi:10.34894/B1K9OU, under the Creative Commons
# Attribution 4.0 International licence.
#
# CHANGES MADE relative to the deposited data, as CC BY 4.0 requires be
# indicated: a subset of columns was selected and the remainder discarded.
# No values were altered, recoded or derived. The retained columns are
# exactly those in panel_A_beck_replication.csv:
#
#   state, state_name, statefip, wrkyr, gini, branch_reform, ln_gini,
#   D_branch
#
# See inst/LICENSE.note and ?beck_banks.
#
# NOTE ON THE DATA. Thirteen of the 49 states carry branch_reform <= 1976,
# the first observed year, ten of them coded 1960. Those units are treated
# throughout the observation window and have no pre-treatment period, so
# gmm_staggered() rejects them by design. Restrict to branch_reform > 1976
# before estimating. See ?beck_banks.

library(readr)   # not a package dependency; data-raw is build-ignored

beck_banks <- read_csv("data-raw/panel_A_beck_replication.csv",
                       show_col_types = FALSE)
beck_banks <- as.data.frame(beck_banks)

stopifnot(
  nrow(beck_banks) == 1519L,
  identical(names(beck_banks),
            c("state", "state_name", "statefip", "wrkyr", "gini",
              "branch_reform", "ln_gini", "D_branch"))
)

usethis::use_data(beck_banks, overwrite = TRUE)
