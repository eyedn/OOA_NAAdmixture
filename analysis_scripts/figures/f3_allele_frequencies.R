# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# f3_allele_frequencies.R
# ______________________________________________________________________________


# set up ----
setwd("~/OOA_NAAdmixture/analysis_scripts/")
library(tidyverse)
library(patchwork)
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
if (
  !exists("diversity.bootstrap.plots", envir = .GlobalEnv) || 
  !exists(
    "sfs.bootstrap.proportion.datatype.interval.plots", envir = .GlobalEnv
    ) ||
  !exists("ld.bootstrap.tcd.1kg.datatype.interval", envir = .GlobalEnv)
) {
  source("diversity.R")
  source("sfs.R")
  source("ld_decay.R")
}


# main ----
f3.allele.frequncies <- guide_area() / (
  (
    diversity.bootstrap.plots$tc.tcd.1kg +
      labs(title = NULL, subtitle = NULL) +
      facet_grid(
        stat ~ ., scales = "free_y", switch = "y",
        labeller = labeller(
          stat = c(pi = "π", theta = "θ")
        )
      ) +
      guides(
        color = "none",
        fill = guide_legend(order = 1, nrow = 1, byrow = TRUE)
      ) +
      theme(strip.placement = "outside")
    ) + (
      (
        sfs.bootstrap.proportion.datatype.interval.plots$tcd.1kg +
          labs(title = NULL, subtitle = NULL) +
          theme(legend.position = "none", strip.text = element_blank())
        ) / (
          ld.bootstrap.tcd.1kg.datatype.interval +
            labs(title = NULL, subtitle = NULL) +
            theme(legend.position = "none", strip.text = element_blank())
          )
        ) +
      plot_layout(widths = c(5, 11))
    ) + 
  plot_annotation(tag_levels = 'A', tag_suffix = '.)') +
  plot_layout(guides = 'collect', heights = c(0.1, 1))
dir.create(OUTPUT.DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(f3.allele.frequncies, file.path(OUTPUT.DIR, "f3_allele_frequencies.rds"))
print(f3.allele.frequncies)
