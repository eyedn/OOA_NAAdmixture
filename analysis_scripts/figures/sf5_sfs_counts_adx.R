# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# sf5_sfs_counts_adx.R
# ______________________________________________________________________________


# set up ----
setwd("~/OOA_NAAdmixture/analysis_scripts/")
library(tidyverse)
library(patchwork)
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
if (
  !exists("sfs.bootstrap.count.datatype.interval.plots", envir = .GlobalEnv) || 
  !exists("sfs.bootstrap.count.plots", envir = .GlobalEnv)
) {
  source("sfs.R")
}


# main ----
sf5_sfs_counts_adx <- (
  sfs.bootstrap.count.datatype.interval.plots$tcd.1kg +
    labs(title = NULL, subtitle = NULL) +
    theme(axis.title.x = element_blank(), strip.text = element_blank())
  ) / (
    sfs.bootstrap.count.plots$all.datatypes.adx.asw +
      labs(title = NULL, subtitle = NULL) +
      plot_spacer() + 
      plot_layout(widths = c(1, 1))
  ) +
  plot_annotation(tag_levels = 'A', tag_suffix = '.)')
dir.create(OUTPUT.DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(sf5_sfs_counts_adx, file.path(OUTPUT.DIR, "sf5_sfs_counts_adx.rds"))
print(sf5_sfs_counts_adx)
