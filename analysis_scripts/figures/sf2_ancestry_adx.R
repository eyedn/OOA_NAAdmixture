# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# sf2_ancestry_adx.R
# ______________________________________________________________________________


# set up ----
setwd("~/OOA_NAAdmixture/analysis_scripts/")
library(tidyverse)
library(patchwork)
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
if (
  !exists("ancestry.bootstrap.all.datatypes.adx.asw.bar", envir = .GlobalEnv) || 
  !exists(
    "ancestry.bootstrap.all.datatypes.adx.asw.histograms", envir = .GlobalEnv
    )
) {
  source("ancestry.R")
}


# main ----
sf2.ancestry.adx <- (
  ancestry.bootstrap.all.datatypes.adx.asw.bar +
    labs(title = NULL, subtitle = NULL)
  ) / (
    ancestry.bootstrap.all.datatypes.adx.asw.histograms$chr1 +
      labs(title = NULL, subtitle = NULL) +
      theme(legend.position = "none")
    ) +
  plot_annotation(tag_levels = 'A', tag_suffix = '.)')
dir.create(OUTPUT.DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(sf2.ancestry.adx, file.path(OUTPUT.DIR, "sf2_ancestry_adx.rds"))
print(sf2.ancestry.adx)
