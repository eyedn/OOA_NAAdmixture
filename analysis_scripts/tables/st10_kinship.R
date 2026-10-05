# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# st10_kinship.R
# ______________________________________________________________________________


# set up ----
setwd("~/OOA_NAAdmixture/analysis_scripts/")
library(tidyverse)
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
if (
  !exists("kinship.summary", envir = .GlobalEnv)
) {
  source("kinship.R")
}
chrom.levels <- c("all", as.character(1:22))


# main ----
st10.kinship <- kinship.summary %>%
  mutate(xmin = round(xmin, 2), xmax = round(xmax, 2)) %>%
  filter(chrom == 1, xmin >= -0.20, xmax <= 0.05) %>%
  select(data.type, pop, xmin, xmax, mean) %>%
  pivot_wider(
    names_from = c(xmin, xmax), values_from = mean, 
    names_glue = "[{xmin}, {xmax})"
      )
readr::write_csv(
  st10.kinship,
  file.path(OUTPUT.DIR, "st10_kinship.csv")
)
