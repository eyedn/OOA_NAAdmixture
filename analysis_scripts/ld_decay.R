# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# ld_decay.R
# ______________________________________________________________________________

# pattern: Imperative Shell


library(tidyverse)
library(scales)

source("analysis_scripts/analysis_utils.R")


# read and harmonize LD-decay tables ----
ld.schema <- c(
  rep = "integer",
  chrom = "character",
  pop = "character",
  window_start = "double",
  window_end = "double",
  distance_bin_bp = "double",
  mean_r2 = "double",
  sum_r2 = "double",
  n_pairs = "double"
)

ld <- read.harmonized.table(
  sim.table = "ld_decay",
  empirical.table = "ld_decay",
  schema = ld.schema,
  sim.required = c(
    "rep", "pop", "distance_bin_bp", "sum_r2", "n_pairs"
  ),
  empirical.required = c(
    "rep", "pop", "distance_bin_bp", "sum_r2", "n_pairs"
  )
)

ld.weighted <- aggregate.ld(
  ld,
  c(
    "rep", "chrom", "pop", "pop_role", "data_source",
    "analysis_level", "distance_bin_bp"
  )
) |>
  filter(n_pairs > 0, is.finite(mean_r2))

ld.series.counts <- ld.weighted |>
  distinct(
    analysis_level, chrom, data_source, rep, pop_role
  ) |>
  count(analysis_level, chrom, name = "expected_series")

ld.shared.bins <- ld.weighted |>
  distinct(
    analysis_level, chrom, distance_bin_bp,
    data_source, rep, pop_role
  ) |>
  count(
    analysis_level,
    chrom,
    distance_bin_bp,
    name = "observed_series"
  ) |>
  left_join(
    ld.series.counts,
    by = c("analysis_level", "chrom")
  ) |>
  filter(observed_series == expected_series)

ld.comparable <- ld.weighted |>
  semi_join(
    ld.shared.bins,
    by = c("analysis_level", "chrom", "distance_bin_bp")
  )
if (nrow(ld.comparable) == 0) {
  stop("LD series have no shared distance bins", call. = FALSE)
}


# predictive decay curves and diagnostics ----
ld.predictive <- predictive.summary(
  ld.comparable,
  value.column = "mean_r2",
  group.columns = c(
    "analysis_level", "chrom", "pop_role", "distance_bin_bp"
  )
)

ld.residuals <- ld.predictive |>
  mutate(residual = empirical_value - simulation_median)

ld.contrasts <- population.contrasts(
  ld.comparable,
  value.column = "mean_r2",
  group.columns = c(
    "rep", "chrom", "data_source", "analysis_level",
    "distance_bin_bp"
  ),
  metric.type = "positive"
)

ld.ranges <- ld.comparable |>
  group_by(
    analysis_level, chrom, data_source, rep, pop_role
  ) |>
  summarize(
    series_min = min(distance_bin_bp),
    series_max = max(distance_bin_bp),
    .groups = "drop"
  ) |>
  group_by(analysis_level, chrom) |>
  summarize(
    common_min = max(series_min),
    common_max = min(series_max),
    .groups = "drop"
  )

if (any(ld.ranges$common_min > ld.ranges$common_max)) {
  stop("LD series do not share a common distance range", call. = FALSE)
}

ld.common <- ld.comparable |>
  left_join(ld.ranges, by = c("analysis_level", "chrom")) |>
  filter(
    distance_bin_bp >= common_min,
    distance_bin_bp <= common_max
  )

trapezoid.area <- function(distance, value) {
  ordering <- order(distance)
  distance <- distance[ordering]
  value <- value[ordering]
  sum(diff(distance) * (head(value, -1) + tail(value, -1)) / 2)
}

ld.metrics <- ld.common |>
  group_by(
    rep, chrom, pop, pop_role, data_source, analysis_level
  ) |>
  arrange(distance_bin_bp, .by_group = TRUE) |>
  summarize(
    first_bin_ld = first(mean_r2),
    half_decay_distance = {
      threshold <- first(mean_r2) / 2
      crossing <- distance_bin_bp[mean_r2 <= threshold]
      if (length(crossing) == 0) NA_real_ else min(crossing)
    },
    decay_auc = trapezoid.area(distance_bin_bp, mean_r2),
    short_long_ratio = first(mean_r2) / last(mean_r2),
    .groups = "drop"
  )


# LD figures ----
ld.decay.plot <- ggplot(
  ld.predictive |> filter(analysis_level == "genome"),
  aes(
    x = distance_bin_bp,
    y = simulation_median,
    color = pop_role,
    fill = pop_role
  )
) +
  geom_ribbon(
    aes(ymin = simulation_lower, ymax = simulation_upper),
    alpha = 0.15,
    color = NA
  ) +
  geom_line(linewidth = 1) +
  geom_line(aes(y = empirical_value), linetype = "dashed") +
  scale_color_manual(values = POP.COLORS) +
  scale_fill_manual(values = POP.COLORS) +
  scale_x_continuous(labels = label_comma()) +
  labs(
    x = "Distance between SNPs (bp)",
    y = expression("Pair-weighted mean " * r^2),
    title = "LD decay predictive comparison"
  ) +
  theme_bw(base_size = 16)

ld.residual.plot <- ggplot(
  ld.residuals |> filter(analysis_level == "genome"),
  aes(x = distance_bin_bp, y = residual, color = pop_role)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 1) +
  scale_color_manual(values = POP.COLORS) +
  scale_x_continuous(labels = label_comma()) +
  labs(
    x = "Distance between SNPs (bp)",
    y = "Empirical minus simulation median",
    title = "LD decay residuals"
  ) +
  theme_bw(base_size = 16)

ld.contrast.plot <- ggplot(
  ld.contrasts |> filter(analysis_level == "genome"),
  aes(x = distance_bin_bp, y = value, color = contrast)
) +
  geom_line(
    data = \(data) filter(data, data_source == "simulation"),
    aes(group = interaction(rep, contrast)),
    alpha = 0.15
  ) +
  geom_line(
    data = \(data) filter(data, data_source == "empirical"),
    linewidth = 1.1
  ) +
  scale_x_continuous(labels = label_comma()) +
  labs(
    x = "Distance between SNPs (bp)",
    y = "Population log ratio",
    title = "LD population-role contrasts"
  ) +
  theme_bw(base_size = 16)


# interactive results ----
print(ld.metrics)
print(ld.contrasts)
print(ld.decay.plot)
print(ld.residual.plot)
print(ld.contrast.plot)
