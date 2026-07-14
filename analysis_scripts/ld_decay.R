# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# ld_decay.R
# ______________________________________________________________________________


library(tidyverse)
library(nanoparquet)
library(glue)
library(ggridges)


# read in data ----
data.dir <- "~/scratch/OOA_NAAdmixture_small/stats"
ld.chr <- map_dfr(
  1:22,
  \(chr) {
    file.path(data.dir, glue("ld_decay.chr{chr}.parquet")) |>
      read_parquet()
  }
)

ld.genWide <- read_parquet(
  file.path(data.dir, "ld_decay.parquet")
) %>%
  mutate(
    chrom = "all"
  )

anc <- bind_rows(anc.chr, anc.genWide)

# plot ld decay genome wide ----
pop_cols <- c(
  AFR = "#56B4E9",
  ADX = "#4B1FA8",
  EUR = "#fb8072"
)

ld.summary <- ld.genWide %>%
  mutate(
    distance_bin_bp = as.numeric(distance_bin_bp),
    mean_r2 = as.numeric(mean_r2),
    pop = factor(pop, levels = c("AFR", "ADX", "EUR"))
  ) %>%
  group_by(pop, distance_bin_bp) %>%
  summarize(
    mean = mean(mean_r2, na.rm = TRUE),
    sd   = sd(mean_r2, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(
  ld.summary,
  aes(
    x = distance_bin_bp,
    y = mean,
    color = pop,
    fill = pop
  )
) +
  geom_ribbon(
    aes(
      ymin = pmax(0, mean - 2 * sd),
      ymax = mean + 2 * sd
    ),
    alpha = 0.20,
    color = NA
  ) +
  geom_line(
    linewidth = 1.25
  ) +
  scale_color_manual(
    values = pop_cols,
    breaks = c("AFR", "ADX", "EUR")
  ) +
  scale_fill_manual(
    values = pop_cols,
    breaks = c("AFR", "ADX", "EUR")
  ) +
  scale_x_continuous(
    limits = c(0, 200000),
    breaks = seq(0, 200000, by = 25000),
    labels = scales::label_comma()
  ) +
  labs(
    x = "Distance between SNPs (bp)",
    y = expression("Mean " * r^2),
    title = "Genome-wide LD Decay across Replicates",
    subtitle = "Lines represent the mean LD decay across replicates; ribbons denote ±2 SD"
  ) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )
