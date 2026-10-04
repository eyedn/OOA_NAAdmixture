# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# f1_model.R
# ______________________________________________________________________________


# set up ----
setwd("~/OOA_NAAdmixture/analysis_scripts/")
library(tidyverse)
library(patchwork)
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
if (
  !exists("combined.ne.plot", envir = .GlobalEnv) || 
  !exists("admix.plot", envir = .GlobalEnv)
) {
  source("calc_ADX_Ne.R")
}


# main ----
f1.model <- (
  guide_area() / (
    ggplot() + theme_void() + theme(text = element_text(size = 24)) + (
      (
        combined.ne.plot +
          labs(title = NULL, subtitle = NULL) +
          guides(
            color = guide_legend(order = 1, nrow = 1),
            linetype = guide_legend(order = 2, nrow = 1)
          )
      ) / (
        admix.plot +
          labs(title = NULL, subtitle = NULL) +
          guides(color = "none")
      ) +
        plot_layout(heights = c(1, 1))
    )
  ) +
    plot_layout(
      guides = "collect",
      heights = c(0.1, 1)
    ) +
    plot_annotation(
      tag_levels = "A",
      tag_suffix = ".)",
      theme = theme(
        plot.tag = element_text(
          family = "Arial",
          size = 24,
          face = "plain",
          color = "black"
        )
      )
    )
) &
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal"
  )

dir.create(OUTPUT.DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(f1.model, file.path(OUTPUT.DIR, "f1_model.rds"))
print(f1.model)
