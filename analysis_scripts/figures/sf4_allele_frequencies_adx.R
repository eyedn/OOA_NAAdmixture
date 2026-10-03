# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# sf4_allele_frequencies.R
# ______________________________________________________________________________


# set up ----
setwd("~/OOA_NAAdmixture/analysis_scripts/")
library(tidyverse)
library(patchwork)
if (
  !exists("diversity.bootstrap.plots", envir = .GlobalEnv) || 
  !exists("sfs.bootstrap.proportion.plots", envir = .GlobalEnv) ||
  !exists("ld.bootstrap.all.datatypes.adx.asw.all.lines", envir = .GlobalEnv)
) {
  source("diversity.R")
  source("sfs.R")
  source("ld_decay.R")
}


# main ----
sf4.allele.frequencies.adx <- guide_area() / (
  (
    diversity.bootstrap.plots$all.datatypes.adx.asw +
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
      sfs.bootstrap.proportion.plots$all.datatypes.adx.asw +
        labs(title = NULL, subtitle = NULL) +
        theme(legend.position = "none", strip.text = element_blank())
    ) / (
      ld.bootstrap.all.datatypes.adx.asw.all.lines +
        labs(title = NULL, subtitle = NULL) +
        theme(legend.position = "none", strip.text = element_blank())
    )
  ) +
    plot_layout(widths = c(5, 11))
) + 
  plot_annotation(tag_levels = 'A', tag_suffix = '.)') +
  plot_layout(guides = 'collect', heights = c(0.1, 1))
print(sf4.allele.frequencies.adx)
