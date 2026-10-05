# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# st7_sfs_counts.R
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
st7.sfs.counts <- sfs.summaries$count %>%
  filter(chrom == 1, minor.allele.count <= 15) %>%
  select(data.type, pop, minor.allele.count, mean) %>%
  pivot_wider(names_from = minor.allele.count, values_from = mean)
readr::write_csv(
  st7.sfs.counts,
  file.path(OUTPUT.DIR, "st7_sfs_counts.csv")
)