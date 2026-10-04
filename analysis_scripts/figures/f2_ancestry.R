# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# f2_ancestry.R
# ______________________________________________________________________________


# set up ----
setwd("~/OOA_NAAdmixture/analysis_scripts/")
library(tidyverse)
library(patchwork)
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
if (
  !exists("ancestry.bootstrap.tcd.1kg.bar", envir = .GlobalEnv) || 
  !exists("ancestry.bootstrap.tcd.1kg.histograms", envir = .GlobalEnv)
  ) {
  source("ancestry.R")
}


# main ----
f2.ancestry <- (
  ancestry.bootstrap.tcd.1kg.bar + 
    labs(title = NULL, subtitle = NULL)
) / (
  ancestry.bootstrap.tcd.1kg.histograms$chr1 + 
    labs(title = NULL, subtitle = NULL) + 
    theme(legend.position = "none")
) +
  plot_annotation(tag_levels = 'A', tag_suffix = '.)')
dir.create(OUTPUT.DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(f2.ancestry, file.path(OUTPUT.DIR, "f2_ancestry.rds"))
print(f2.ancestry)
