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
library(scales)

# read in data ----
data.dir <- "~/scratch/OOA_NAAdmixture_small/stats"
kin.chr <- map_dfr(
  1:22,
  \(chr) {
    file.path(data.dir, glue("kinship_unrelated.chr{chr}.parquet")) |>
      read_parquet() 
  }
)

kin.genWide <- read_parquet(
  file.path(data.dir, "kinship.parquet")
) %>%
  mutate(
    chrom = "all"
  )

kin <- bind_rows(kin.chr, kin.genWide) %>%
  mutate(
    kinship = as.numeric(kinship),
    pop = factor(pop, levels = c("AFR", "ADX", "EUR"))
  )

# summarize kinship structure ----
kin.summary <- kin %>%
  group_by(rep, chrom, pop) %>%
  summarize(
    mean = mean(kinship),
    sd = sd(kinship),
    .groups = "drop"
  )

pop_cols <- c(
  AFR = "#56B4E9",
  ADX = "#4B1FA8",
  EUR = "#fb8072"
)

chrom.levels <- c(as.character(22:1), "all")

# ridgeline of mean pairwise kinship ----
kin.summary.mean.medians <- kin.summary %>%
  group_by(chrom, pop) %>%
  summarize(
    median = median(mean),
    .groups = "drop"
  ) %>%
  mutate(
    y = match(chrom, chrom.levels)
  )

ggplot(
  kin.summary,
  aes(
    x = mean,
    y = factor(chrom, levels = chrom.levels),
    fill = pop
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
    data = kin.summary.mean.medians,
    aes(
      x = median,
      y = y + 0.35,
      fill = pop
    ),
    inherit.aes = FALSE,
    shape = 23,
    stroke = 1,
    color = "black",
    size = 3
  ) +
  scale_fill_manual(
    values = pop_cols,
    breaks = c("AFR", "ADX", "EUR")
  ) +
  labs(
    x = "Mean Pairwise Kinship",
    y = "Chromosome",
    title = "Mean Pairwise Kinship across replicates"
  ) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# ridgeline of sd pairwise kinship ----
kin.summary.sd.medians <- kin.summary %>%
  group_by(chrom, pop) %>%
  summarize(
    median = median(sd),
    .groups = "drop"
  ) %>%
  mutate(
    y = match(chrom, chrom.levels)
  )

ggplot(
  kin.summary,
  aes(
    x = sd,
    y = factor(chrom, levels = chrom.levels),
    fill = pop
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
    data = kin.summary.sd.medians,
    aes(
      x = median,
      y = y + 0.35,
      fill = pop
    ),
    inherit.aes = FALSE,
    shape = 23,
    stroke = 1,
    color = "black",
    size = 3
  ) +
  scale_fill_manual(
    values = pop_cols,
    breaks = c("AFR", "ADX", "EUR")
  ) +
  labs(
    x = "Standard deviation of Pairwise Kinship",
    y = "Chromosome",
    title = "Standard deviation of Pairwise Kinship across replicates"
  ) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# mean histogram of pairwise kinship ----
break_width <- 0.01

breaks <- seq(
  floor(min(kin$kinship, na.rm = TRUE) / break_width) * break_width,
  ceiling(max(kin$kinship, na.rm = TRUE) / break_width) * break_width,
  by = break_width
)

kin.hist.df <- kin %>%
  group_by(rep, chrom, pop) %>%
  group_modify(~{
    h <- hist(
      .x$kinship,
      breaks = breaks,
      plot = FALSE
    )
    tibble(
      xmin  = head(h$breaks, -1),
      xmax  = tail(h$breaks, -1),
      xmid  = h$mids,
      count = h$counts,
      frac = h$counts / sum(h$counts)
    )
  }) %>%
  ungroup()

kin.hist.summary <- kin.hist.df %>%
  group_by(pop, chrom, xmin, xmax, xmid) %>%
  summarize(
    mean.frac = mean(frac),
    sd.frac   = sd(frac),
    .groups = "drop"
  )

ggplot(
  kin.hist.summary,
  aes(
    x = xmid,
    y = mean.frac,
    fill = pop
  )
) +
  geom_col(
    width = diff(breaks)[1] * 0.95,
    color = "black",
    linewidth = 0.5
  ) +
  # geom_errorbar(
  #   aes(
  #     ymin = pmax(0, mean.frac - 2 * sd.frac),
  #     ymax = mean.frac + 2 * sd.frac
  #   ),
  #   width = 0,
  #   linewidth = 0.5,
  #   color = "black"
  # ) +
  facet_wrap(
    ~factor(chrom, levels = rev(chrom.levels)),
    ncol = 6
  ) +
  scale_fill_manual(
    values = pop_cols,
    breaks = c("AFR", "ADX", "EUR")
  ) +
  labs(
    x = "Pairwise Kinship",
    y = "Mean fraction of individuals per bin",
    title = "Pairwise kinship distributions of ADX across chromosomes",
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

# per-individual average kinship -----
kin.long <- bind_rows(
  kin %>%
    transmute(
      rep,
      chrom,
      pop,
      sample = sample1,
      partner = sample2,
      kinship
    ),
  
  kin %>%
    transmute(
      rep,
      chrom,
      pop,
      sample = sample2,
      partner = sample1,
      kinship
    )
)

individual.kin <- kin.long %>%
  group_by(rep, chrom, pop, sample) %>%
  summarize(
    mean.kinship = mean(kinship),
    .groups = "drop"
  )

ggplot(
  individual.kin,
  aes(
    pop,
    mean.kinship,
    fill = pop,
    color = pop
  )
) +
  geom_violin(
    linewidth = 0.75,
    alpha = 0.35,
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.18,
    linewidth = 0.75,
    outlier.shape = NA,
    fill = "white"
  ) +
  geom_jitter(
    width = 0.08,
    size = 0.7,
    alpha = 0.2
  ) +
  scale_fill_manual(
    values = pop_cols
  ) +
  scale_color_manual(
    values = pop_cols
  )
