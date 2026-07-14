# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# sfs.R
# ______________________________________________________________________________


library(tidyverse)
library(nanoparquet)
library(glue)
library(ggridges)
library(scales)

# read in data ----
data.dir <- "~/scratch/OOA_NAAdmixture_small/stats"
sfs1D.chr <- map_dfr(
  1:22,
  \(chr) {
    file.path(data.dir, glue("sfs.chr{chr}.parquet")) |>
      read_parquet() 
  }
)

sfs1D.genWide <- read_parquet(
  file.path(data.dir, "sfs.parquet")
) %>%
  mutate(
    chrom = "all"
  )

sfs1D <- bind_rows(sfs1D.chr, sfs1D.genWide) %>%
  filter(derived_allele_count != "0")

# summarize sfs ----
sfs1D.summary <- sfs1D %>%
  mutate(
    derived_allele_count = as.numeric(derived_allele_count),
    count = as.numeric(count),
    pop = factor(pop, levels = c("AFR", "ADX", "EUR"))
  ) %>%
  group_by(rep, chrom, pop) %>%
  mutate(
    fraction = count / sum(count)
  ) %>%
  ungroup() %>%
  group_by(chrom, pop, derived_allele_count) %>%
  summarize(
    mean.frac = mean(fraction),
    sd.frac = sd(fraction),
    mean.count = mean(count),
    sd.count = sd(count),
    .groups = "drop"
  ) 

# plot pseudo log raw counts ----
ggplot(
  sfs1D.summary %>% filter(chrom == "1", derived_allele_count <= 10),
  aes(
    x = derived_allele_count,
    y = mean.count,
    color = pop,
    fill = pop
  )
) +
  geom_ribbon(
    aes(
      ymin = pmax(0, mean.count - 2 * sd.count),
      ymax = mean.count + 2 * sd.count
    ),
    alpha = 0.20,
    color = NA
  ) +
  geom_line(
    linewidth = 1
  ) +
  facet_wrap(
    ~factor(chrom, levels = rev(chrom.levels)),
    ncol = 6
  ) +
  scale_color_manual(
    values = pop_cols
  ) +
  scale_fill_manual(
    values = pop_cols
  ) +
  scale_y_continuous(
    trans = pseudo_log_trans(base = 10)
  ) +
  labs(
    x = "Derived allele count",
    y = "Mean number of sites",
    title = "One-dimensional Site Frequency Spectrum",
    subtitle = "Lines represent the mean count across replicates; shaded ribbons denote ±2 SD"
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

# plot pseudo log frequency ----
ggplot(
  sfs1D.summary %>% filter(chrom == "1", derived_allele_count <= 10),
  aes(
    x = derived_allele_count,
    y = mean.frac,
    color = pop,
    fill = pop
  )
) +
  geom_ribbon(
    aes(
      ymin = pmax(0, mean.frac - 2 * sd.frac),
      ymax = mean.frac + 2 * sd.frac
    ),
    alpha = 0.20,
    color = NA
  ) +
  geom_line(
    linewidth = 1
  ) +
  facet_wrap(
    ~factor(chrom, levels = rev(chrom.levels)),
    ncol = 6
  ) +
  scale_color_manual(
    values = pop_cols
  ) +
  scale_fill_manual(
    values = pop_cols
  ) +
  scale_y_continuous(
    trans = pseudo_log_trans(base = 10)
  ) +
  labs(
    x = "Derived allele count",
    y = "Mean fraction of segregating sites",
    title = "Normalized One-dimensional Site Frequency Spectrum",
    subtitle = "Lines represent the mean fraction across replicates; shaded ribbons denote ±2 SD"
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

# plot number of segregating sites per chromosome ----
segsites.summary <- sfs1D %>%
  filter(chrom != "all") %>%
  mutate(
    count = as.numeric(count),
    pop = factor(pop, levels = c("AFR", "ADX", "EUR")),
    chrom = factor(as.integer(chrom), levels = 1:22)
  ) %>%
  group_by(rep, chrom, pop) %>%
  summarize(
    segregating_sites = sum(count),
    .groups = "drop"
  ) %>%
  group_by(chrom, pop) %>%
  summarize(
    mean.segsites = mean(segregating_sites),
    sd.segsites = sd(segregating_sites),
    .groups = "drop"
  )

ggplot(
  segsites.summary,
  aes(
    x = chrom,
    y = mean.segsites,
    fill = pop
  )
) +
  geom_col(
    position = position_dodge(width = 0.9),
    width = 0.8
  ) +
  geom_errorbar(
    aes(
      ymin = pmax(0, mean.segsites - 2 * sd.segsites),
      ymax = mean.segsites + 2 * sd.segsites
    ),
    position = position_dodge(width = 0.9),
    width = 0.25,
    linewidth = 0.6
  ) +
  scale_fill_manual(
    values = pop_cols
  ) +
  labs(
    x = "Chromosome",
    y = "Mean number of segregating sites",
    title = "Segregating sites by chromosome and population",
    subtitle = "Bars represent the mean across replicates; error bars denote ±2 SD"
  ) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )
