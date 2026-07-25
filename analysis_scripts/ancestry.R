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
library(scales)

source("analysis_scripts/analysis_utils.R")


# read and harmonize ancestry tables ----
anc.schema <- c(
  rep = "integer",
  chrom = "character",
  pop = "character",
  sample_id = "character",
  vcf_sample_id = "character",
  afr_tspop = "double",
  eur_tspop = "double",
  afr_q = "double",
  eur_q = "double",
  afr_unsupervised_q = "double",
  eur_unsupervised_q = "double",
  span = "double"
)

anc <- read.harmonized.table(
  sim.table = "ancestry",
  empirical.table = "ancestry",
  schema = anc.schema,
  sim.required = c(
    "rep", "pop", "sample_id", "vcf_sample_id",
    "afr_tspop", "eur_tspop", "afr_q", "eur_q", "span"
  ),
  empirical.required = c(
    "rep", "pop", "sample_id", "vcf_sample_id", "afr_q", "eur_q",
    "afr_unsupervised_q", "eur_unsupervised_q"
  )
) |>
  mutate(
    comparison_afr = if_else(
      data_source == "simulation",
      afr_tspop,
      afr_q
    )
  )


# ancestry distribution summaries ----
anc.admixed <- anc |>
  filter(pop_role == "ADMIXED")

anc.sample.summary <- anc.admixed |>
  group_by(data_source, analysis_level, chrom, rep) |>
  summarize(
    mean = mean(comparison_afr, na.rm = TRUE),
    sd = sd(comparison_afr, na.rm = TRUE),
    median = median(comparison_afr, na.rm = TRUE),
    iqr = IQR(comparison_afr, na.rm = TRUE),
    lower_tail = quantile(
      comparison_afr,
      0.025,
      na.rm = TRUE,
      names = FALSE
    ),
    upper_tail = quantile(
      comparison_afr,
      0.975,
      na.rm = TRUE,
      names = FALSE
    ),
    .groups = "drop"
  )

anc.predictive.summary <- predictive.summary(
  anc.sample.summary,
  value.column = "mean",
  group.columns = c("analysis_level", "chrom")
)

anc.chromosome.deviation <- anc.sample.summary |>
  select(data_source, rep, chrom, analysis_level, mean) |>
  group_by(data_source, rep) |>
  mutate(genome_mean = mean[chrom == "all"][1]) |>
  ungroup() |>
  filter(analysis_level == "chromosome") |>
  mutate(chromosome_to_genome_deviation = mean - genome_mean)


# ancestry figures ----
anc.distribution.plot <- ggplot(
  anc.admixed |> filter(analysis_level == "genome"),
  aes(
    x = comparison_afr,
    color = data_source,
    linetype = data_source,
    group = interaction(data_source, rep)
  )
) +
  stat_ecdf(linewidth = 1.1) +
  geom_rug(
    data = \(data) filter(data, data_source == "empirical"),
    aes(color = data_source),
    sides = "b",
    alpha = 0.45
  ) +
  scale_color_manual(
    values = c(simulation = "#4B1FA8", empirical = "black")
  ) +
  labs(
    x = "African ancestry",
    y = "Cumulative fraction",
    title = "Admixed ancestry: simulation and empirical observation"
  ) +
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank())

anc.chromosome.envelope <- anc.sample.summary |>
  filter(data_source == "simulation", analysis_level == "chromosome") |>
  pivot_longer(c(mean, sd), names_to = "summary", values_to = "value") |>
  group_by(chrom, summary) |>
  summarize(
    median = median(value, na.rm = TRUE),
    lower = quantile(value, 0.025, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

anc.chromosome.empirical <- anc.sample.summary |>
  filter(data_source == "empirical", analysis_level == "chromosome") |>
  pivot_longer(c(mean, sd), names_to = "summary", values_to = "value")

anc.chromosome.plot <- ggplot(
  anc.chromosome.envelope,
  aes(x = as.integer(chrom), y = median)
) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25) +
  geom_line(linewidth = 1) +
  geom_point(
    data = anc.chromosome.empirical,
    aes(y = value),
    color = "#D55E00",
    size = 2
  ) +
  facet_wrap(~summary, scales = "free_y") +
  scale_x_continuous(
    breaks = as.integer(unique(anc.chromosome.envelope$chrom))
  ) +
  labs(
    x = "Chromosome",
    y = "African ancestry summary",
    title = "Chromosome ancestry predictive envelope"
  ) +
  theme_bw(base_size = 16)

anc.calibration.data <- bind_rows(
  anc.admixed |>
    filter(data_source == "simulation") |>
    transmute(
      data_source,
      chrom,
      reference = afr_tspop,
      estimate = afr_q
    ),
  anc.admixed |>
    filter(data_source == "empirical") |>
    transmute(
      data_source,
      chrom,
      reference = afr_unsupervised_q,
      estimate = afr_q
    )
)

anc.calibration.plot <- ggplot(
  anc.calibration.data,
  aes(x = reference, y = estimate, color = data_source)
) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~data_source) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    x = "Tree-sequence truth or unsupervised estimate",
    y = "Supervised estimate",
    title = "Ancestry calibration"
  ) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none")


# interactive results ----
print(anc.sample.summary)
print(anc.predictive.summary)
print(anc.chromosome.deviation)
print(anc.distribution.plot)
print(anc.chromosome.plot)
print(anc.calibration.plot)
