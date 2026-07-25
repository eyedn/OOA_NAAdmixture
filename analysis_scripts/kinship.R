# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# kinship.R
# ______________________________________________________________________________

# pattern: Imperative Shell


library(tidyverse)
library(scales)

source("analysis_scripts/analysis_utils.R")


# constants ----
KIN_CUTOFF <- 0.0442


# read and harmonize unrelated-kinship tables ----
kin.schema <- c(
  rep = "integer",
  chrom = "character",
  pop = "character",
  sample1 = "character",
  sample2 = "character",
  kinship = "double"
)

kin <- read.harmonized.table(
  sim.table = "kinship_unrelated",
  empirical.table = "kinship_unrelated",
  schema = kin.schema,
  sim.required = c(
    "rep", "pop", "sample1", "sample2", "kinship"
  ),
  empirical.required = c(
    "rep", "pop", "sample1", "sample2", "kinship"
  ),
  empirical.rename = c(id1 = "sample1", id2 = "sample2")
)


# pairwise and per-individual summaries ----
kin.summary <- kin |>
  group_by(
    data_source, analysis_level, chrom, rep, pop, pop_role
  ) |>
  summarize(
    mean = mean(kinship, na.rm = TRUE),
    median = median(kinship, na.rm = TRUE),
    sd = sd(kinship, na.rm = TRUE),
    q05 = quantile(kinship, 0.05, na.rm = TRUE),
    q95 = quantile(kinship, 0.95, na.rm = TRUE),
    negative_fraction = mean(kinship < 0, na.rm = TRUE),
    above_cutoff_fraction = mean(kinship > KIN_CUTOFF, na.rm = TRUE),
    .groups = "drop"
  )

kin.long <- bind_rows(
  kin |>
    transmute(
      rep, chrom, pop, pop_role, data_source, analysis_level,
      sample = sample1, partner = sample2, kinship
    ),
  kin |>
    transmute(
      rep, chrom, pop, pop_role, data_source, analysis_level,
      sample = sample2, partner = sample1, kinship
    )
)

individual.kin <- kin.long |>
  group_by(
    rep, chrom, pop, pop_role, data_source, analysis_level, sample
  ) |>
  summarize(mean_kinship = mean(kinship, na.rm = TRUE), .groups = "drop")

kin.contrasts <- population.contrasts(
  kin.summary,
  value.column = "mean",
  group.columns = c(
    "rep", "chrom", "data_source", "analysis_level"
  ),
  metric.type = "kinship"
)


# kinship figures ----
kin.distribution.plot <- ggplot(
  kin |> filter(analysis_level == "genome"),
  aes(
    x = kinship,
    color = data_source,
    group = interaction(data_source, rep)
  )
) +
  stat_ecdf(linewidth = 1) +
  facet_wrap(~pop_role, scales = "free_x") +
  labs(
    x = "Pairwise kinship coefficient",
    y = "Cumulative fraction",
    title = "Unrelated-sample kinship distributions"
  ) +
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank())

kin.individual.plot <- ggplot(
  individual.kin |> filter(analysis_level == "genome"),
  aes(x = pop_role, y = mean_kinship, color = data_source)
) +
  geom_violin(
    data = \(data) filter(data, data_source == "simulation"),
    aes(fill = pop_role),
    alpha = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    data = \(data) filter(data, data_source == "empirical"),
    position = position_jitter(width = 0.08),
    alpha = 0.45
  ) +
  scale_fill_manual(values = POP.COLORS) +
  labs(
    x = "Population role",
    y = "Per-individual mean kinship",
    title = "Per-individual kinship distributions"
  ) +
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank())

kin.tail.plot <- ggplot(
  kin |> filter(analysis_level == "genome", kinship >= 0),
  aes(
    x = kinship,
    color = data_source,
    group = interaction(data_source, rep)
  )
) +
  stat_ecdf(aes(y = after_stat(1 - y)), linewidth = 1) +
  geom_vline(xintercept = KIN_CUTOFF, linetype = "dashed") +
  facet_wrap(~pop_role) +
  scale_y_log10(labels = label_percent()) +
  labs(
    x = "Pairwise kinship coefficient",
    y = "Upper-tail fraction",
    title = "Unrelated-kinship upper tail"
  ) +
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank())


# interactive results ----
# Pairwise rows are dependent, so summaries are descriptive and no pair-level
# inferential p-values are calculated.
print(kin.summary)
print(kin.contrasts)
print(kin.distribution.plot)
print(kin.individual.plot)
print(kin.tail.plot)
