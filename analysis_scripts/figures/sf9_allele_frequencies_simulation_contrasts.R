# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# sf9_allele_frequencies_simulation_contrasts.R
# ______________________________________________________________________________


# set up ----
setwd("~/OOA_NAAdmixture/analysis_scripts/")
library(tidyverse)
library(patchwork)
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
if (
  !exists("diversity.simulation.contrast.plots", envir = .GlobalEnv) || 
  !exists("sfs.simulation.contrast.plots", envir = .GlobalEnv) ||
  !exists("ld.simulation.contrast.plots", envir = .GlobalEnv)
) {
  source("diversity.R")
  source("sfs.R")
  source("ld_decay.R")
}


# main ----
sf9.allele.frequencies.simulation.contrasts <-
  (
    diversity.simulation.contrast.plots$bonferroni +
      labs(title = NULL, subtitle = NULL) +
      theme(legend.position = "none")
   ) / (
     sfs.simulation.contrast.plots$bonferroni +
       labs(title = NULL, subtitle = NULL) +
       theme(legend.position = "none", strip.text.x = element_blank())
     ) / (
       ld.simulation.contrast.plots$bonferroni +
         labs(title = NULL, subtitle = NULL) +
         facet_grid(
           rows = vars(ld.stat = "r^2"),
           cols = vars(contrast),
           scales = "free_y",
           labeller = labeller(
             .rows = label_parsed,
             contrast = SFS.PLOT.STYLES$contrast.labels
           )
         ) +
         theme(legend.position = "none", strip.text.x = element_blank())
       ) +
  plot_annotation(tag_levels = 'A', tag_suffix = '.)')
dir.create(OUTPUT.DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(
  sf9.allele.frequencies.simulation.contrasts,
  file.path(OUTPUT.DIR, "sf9_allele_frequencies_simulation_contrasts.rds")
)
print(sf9.allele.frequencies.simulation.contrasts)
