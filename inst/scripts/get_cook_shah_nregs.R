# Retrieve and prepare the district-level NREGS panel of Cook and Shah
# (2022), used as a second real-world staggered rollout example.
#
# This dataset is NOT bundled with staggeredGMM. No redistribution licence
# has been identified for it, so the package ships this retrieval script
# instead. Running it downloads the file from the author's website; the data
# it retrieves are the authors' and are not covered by this package's
# licence. Please cite:
#
#   Cook, C. J. and Shah, M. (2022). "Aggregate Effects from Public Works:
#   Evidence from India." The Review of Economics and Statistics 104(4),
#   797-806. doi:10.1162/rest_a_00993
#
# Source page: https://cjustincook.weebly.com/research.html
#
# Usage:
#   source(system.file("scripts", "get_cook_shah_nregs.R",
#                      package = "staggeredGMM"))
#   nregs <- get_cook_shah_nregs("~/Downloads/district_light_panel.dta")

get_cook_shah_nregs <- function(path) {
  if (!requireNamespace("haven", quietly = TRUE)) {
    stop("Package 'haven' is required to read the Stata file. ",
         "Install it with install.packages(\"haven\").", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("File not found: ", path, "\n",
         "Download 'district_light_panel.dta' from the replication ",
         "materials at https://cjustincook.weebly.com/research.html and ",
         "pass its path to this function.", call. = FALSE)
  }
  as.data.frame(haven::read_dta(path))
}


# The panel records the NREGS rollout in three waves (2006, 2007, 2008) via
# the indicators nr06, nr07 and nr08. staggeredGMM expects a single cohort
# column giving each district's first treated year, so construct one:
#
#   dat <- get_cook_shah_nregs(path)
#   dat$cohort <- with(dat, ifelse(nr06 == 1, 2006L,
#                           ifelse(nr07 == 1, 2007L,
#                           ifelse(nr08 == 1, 2008L, 0L))))
#
# Note that every district is treated by 2008 in this design, so there is no
# never-treated group and identification rests on not-yet-treated controls
# alone. Check the resulting `identified` column before interpreting the
# aggregate.
