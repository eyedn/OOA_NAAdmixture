# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# st9_ld.R
# ______________________________________________________________________________


# set up ----
setwd("~/OOA_NAAdmixture/analysis_scripts/")
library(tidyverse)
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
if (
  !exists("ld.summary", envir = .GlobalEnv)
) {
  source("ld_decay.R")
}
chrom.levels <- c("all", as.character(1:22))


# main ----
st9.ld <- ld.summary %>%
  filter(chrom == 1, distance_bin_bp <= 200000) %>%
  select(data.type, pop, distance_bin_bp, mean) %>%
  pivot_wider(names_from = distance_bin_bp, values_from = mean)
readr::write_csv(
  st9.ld,
  file.path(OUTPUT.DIR, "st9_ld.csv")
)
