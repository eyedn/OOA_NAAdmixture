# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# st4_ancestry_sd.R
# ______________________________________________________________________________


# set up ----
setwd("~/OOA_NAAdmixture/analysis_scripts/")
library(tidyverse)
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
if (
  !exists("ancestry.bootstrap.summary", envir = .GlobalEnv)
) {
  source("ancestry.R")
}
chrom.levels <- c("all", as.character(1:22))


# main ----
st3.ancestry.sd <- ancestry.bootstrap.summary %>%
  filter(stat == "sd") %>%
  select(data.type, chrom, mean) %>%
  mutate(chrom = factor(chrom, levels = chrom.levels)) %>%
  pivot_wider(names_from = chrom, values_from = mean) %>%
  select(
    any_of(names(.)[!names(.) %in% chrom.levels]), any_of(chrom.levels)
  ) %>%
  rename(`genome-wide` = all)
readr::write_csv(
  st3.ancestry.sd,
  file.path(OUTPUT.DIR, "st3_ancestry_sd.csv")
)