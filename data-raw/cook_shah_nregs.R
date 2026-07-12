library(haven)

cook_shah_nregs <- as.data.frame(read_dta("data-raw/district_light_panel.dta"))

usethis::use_data(cook_shah_nregs, overwrite = TRUE)
