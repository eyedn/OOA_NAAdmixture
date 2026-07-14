# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# ancestry.R
# ______________________________________________________________________________


library(tidyverse)
library(nanoparquet)
library(glue)
library(ggridges)


# read in data ----
data.dir <- "~/scratch/OOA_NAAdmixture_small/stats"
anc.chr <- map_dfr(
  1:22,
  \(chr) {
    file.path(data.dir, glue("ancestry.chr{chr}.parquet")) |>
      read_parquet() |>
      mutate(
        afr_tspop = as.numeric(afr_tspop),
        eur_tspop = as.numeric(eur_tspop),
        afr_q = as.numeric(afr_q),
        eur_q = as.numeric(eur_q)
      )
  }
)

anc.genWide <- read_parquet(
  file.path(data.dir, "ancestry.parquet")
) %>%
  mutate(
    afr_tspop = as.numeric(afr_tspop),
    eur_tspop = as.numeric(eur_tspop),
    afr_q = as.numeric(afr_q),
    eur_q = as.numeric(eur_q),
    chrom = "all"
  )

anc <- bind_rows(anc.chr, anc.genWide)

# summarize data ----
anc.adx <- anc %>%
  filter(pop == "ADX", chrom != "all")

chrom.levels <- c(as.character(22:1), "all")

anc.summary <- anc.adx %>%
  group_by(chrom, rep) %>%
  summarize(
    mean = mean(afr_tspop, na.rm = TRUE),
    sd   = sd(afr_tspop, na.rm = TRUE),
    .groups = "drop"
  )

# ridge line plots of mean ancestry ----
anc.summary.medians <- anc.summary %>%
  group_by(chrom) %>%
  summarize(
    median = median(mean),
    .groups = "drop"
  ) %>%
  mutate(
    y = match(chrom, chrom.levels)
  )

ggplot(
  anc.summary,
  aes(
    x = mean,
    y = factor(chrom, levels = chrom.levels),
    fill = factor(chrom, levels = chrom.levels)
  )
) +
  geom_density_ridges(
    linewidth = 1,
    scale = 1.5,
    rel_min_height = 0.005,
    alpha = 0.7,
    color = "black"
  ) +
  geom_point(
    data = anc.summary.medians,
    aes(
      x = median,
      y = y + 0.35,
      fill = factor(chrom, levels = chrom.levels)
    ),
    inherit.aes = FALSE,
    shape = 23,
    stroke = 1,
    color = "black",
    size = 3
  ) +
  labs(
    x = "Mean African ancestry",
    y = "Chromosome",
    title = "Mean African Ancestry of ADX across replicates"
  ) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# ridgel line plots of sd ancestry ----
anc.summary.sd.medians <- anc.summary %>%
  group_by(chrom) %>%
  summarize(
    median = median(sd),
    .groups = "drop"
  ) %>%
  mutate(
    y = match(chrom, chrom.levels)
  )

ggplot(
  anc.summary,
  aes(
    x = sd,
    y = factor(chrom, levels = chrom.levels),
    fill = factor(chrom, levels = chrom.levels)
  )
) +
  geom_density_ridges(
    linewidth = 1,
    scale = 1.5,
    rel_min_height = 0.005,
    alpha = 0.7,
    color = "black"
  ) +
  geom_point(
    data = anc.summary.sd.medians,
    aes(
      x = median,
      y = y + 0.35,
      fill = factor(chrom, levels = chrom.levels)
    ),
    inherit.aes = FALSE,
    shape = 23,
    stroke = 1,
    color = "black",
    size = 3
  ) +
  labs(
    x = "Standard deviation of African ancestry",
    y = "Chromosome",
    title = "Standard Deviation of African Ancestry of ADX across Replicates"
  ) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# mean density curves ----
# ==============================================================================
# Mean ancestry histograms across replicates
# ==============================================================================
breaks <- seq(0, 1, by = 0.05)

hist.df <- anc %>%
  filter(pop == "ADX") %>%
  group_by(chrom, rep) %>%
  group_modify(~{
    h <- hist(
      .x$afr_tspop,
      breaks = breaks,
      plot = FALSE
    )
    tibble(
      xmin  = head(h$breaks, -1),
      xmax  = tail(h$breaks, -1),
      xmid  = h$mids,
      count = h$counts,
      frac = h$counts / sum(count)
    )
  }) %>%
  ungroup()

hist.summary <- hist.df %>%
  group_by(chrom, xmin, xmax, xmid) %>%
  summarize(
    mean.frac = mean(frac),
    sd.frac   = sd(frac),
    .groups = "drop"
  )

ggplot(
  hist.summary,
  aes(
    x = xmid,
    y = mean.frac,
    fill = factor(chrom, levels = rev(chrom.levels))
  )
) +
  geom_col(
    width = diff(breaks)[1] * 0.95,
    color = "black",
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(
      ymin = pmax(0, mean.frac - 2 * sd.frac),
      ymax = mean.frac + 2 * sd.frac
    ),
    width = 0,
    linewidth = 0.5,
    color = "black"
  ) +
  facet_wrap(
    ~factor(chrom, levels = rev(chrom.levels)),
    ncol = 6
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  labs(
    x = "African ancestry",
    y = "Mean fraction of individuals per bin",
    title = "African Ancestry Distributions of ADX across chromosomes",
    subtitle = "Bars represent the mean histogram across replicates; error bars denote ±2 SD"
  ) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1, "lines")
  )
