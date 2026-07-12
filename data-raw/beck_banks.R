library(readr)

beck_banks <- read_csv("data-raw/panel_A_beck_replication.csv", show_col_types = FALSE)
beck_banks <- as.data.frame(beck_banks)

usethis::use_data(beck_banks, overwrite = TRUE)
