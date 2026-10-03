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
print(sf2.ancestry.adx)
