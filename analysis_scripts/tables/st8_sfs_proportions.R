# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# st8_sfs_proportions.R
# ______________________________________________________________________________


# set up ----
setwd("~/OOA_NAAdmixture/analysis_scripts/")
library(tidyverse)
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
if (
  !exists("sfs.summaries", envir = .GlobalEnv)
) {
  source("sfs.R")
}
chrom.levels <- c("all", as.character(1:22))


# main ----
st8.sfs.proportions <- sfs.summaries$proportion %>%
  filter(chrom == 1, minor.allele.count <= 15) %>%
  select(data.type, pop, minor.allele.count, mean) %>%
  pivot_wider(names_from = minor.allele.count, values_from = mean)
readr::write_csv(
  st8.sfs.proportions,
  file.path(OUTPUT.DIR, "st8_sfs_proportions.csv")
)