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
library(svgparser)
library(grid)
library(patchwork)
if (
  !exists("combined.ne.plot", envir = .GlobalEnv) || 
  !exists("admix.plot", envir = .GlobalEnv)
) {
  source("calc_ADX_Ne.R")
}


# main ----
svg.file <- paste0(
  "/project2/jazlynmo_738/karatas/",
  "OOA_NAAdmixture_data/demography_plot.svg"
)

normalized.svg.file <- paste0(
  "/project2/jazlynmo_738/karatas/",
  "OOA_NAAdmixture_data/demography_plot_normalized.svg"
)

# Rewrite the Matplotlib SVG into a simpler normalized SVG.
rsvg::rsvg_svg(
  svg = svg.file,
  file = normalized.svg.file
)

# Read the normalized file, not the original Matplotlib file.
my.grob <- svgparser::read_svg(normalized.svg.file)

# Test the grob.
grid::grid.newpage()
grid::grid.draw(my.grob)

demesdraw.res <- patchwork::wrap_elements(
  full = my.grob,
  clip = FALSE
  )

# Draw the vector graphic on your screen
right.panels <- (
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

figure.panels <- demesdraw.res | right.panels

f1.model <- (
  guide_area() / figure.panels +
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

print(f1.model)
