# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# diversity.R
# ______________________________________________________________________________


library(tidyverse)
library(nanoparquet)
library(glue)
library(ggridges)


# read in data ----
data.dir <- "~/scratch/OOA_NAAdmixture_small/stats"
div.chr <- map_dfr(
  1:22,
  \(chr) {
    file.path(data.dir, glue("pi_theta_stats.chr{chr}.parquet")) |>
      read_parquet() |>
      mutate(
        value = as.numeric(value)
      )
  }
)

div.genWide <- read_parquet(
  file.path(data.dir, "pi_theta_stats.parquet")
) %>%
  mutate(
    value = as.numeric(value),
    chrom = "all"
  )

div <- bind_rows(div.chr, div.genWide)

# plot pi and theta ----
pop_cols <- c(
  AFR = "#56B4E9",
  ADX = "#4B1FA8",
  EUR = "#fb8072"
)

ggplot(
  div %>%
    # filter(chrom != "all") %>%
    mutate(
      pop = factor(pop, levels = c("AFR", "ADX", "EUR")),
      stat = factor(stat, levels = c("pi", "theta"))
    ),
  aes(
    x = factor(chrom, levels = rev(chrom.levels)),
    y = value,
    color = pop
  )
) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7,
    linewidth = 0.75,
    outlier.shape = NA,
    fill = "white"
  ) +
  facet_wrap(
    ~stat,
    ncol = 1,
    scales = "free_y",
    labeller = as_labeller(
      c(
        pi = "pi",
        theta = "theta[W]"
      ),
      label_parsed
    )
  ) +
  scale_y_continuous(
    labels = scales::label_number()
  ) +
  scale_color_manual(
    values = pop_cols,
    breaks = c("AFR", "ADX", "EUR")
  ) +
  labs(
    x = "Chromosome",
    y = "Statistic value",
    title = expression(
      "Nucleotide Diversity (" * pi * ") and Watterson's " *
        theta[W] * " across Replicates"
    )
  ) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1, "lines")
  )
