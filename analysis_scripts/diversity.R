# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# diversity.R
# ______________________________________________________________________________

# pattern: Imperative Shell


library(tidyverse)
library(scales)

source("analysis_scripts/analysis_utils.R")


# read and harmonize diversity tables ----
div.schema <- c(
  rep = "integer",
  chrom = "character",
  pop = "character",
  stat = "character",
  value = "double",
  ne_value = "double",
  mutation_rate = "double",
  span = "double",
  segregating_sites = "double",
  wattersons_const = "double"
)

read.diversity.scope <- function(empirical.table, empirical.scope) {
  read.harmonized.table(
    sim.table = "pi_theta_stats",
    empirical.table = empirical.table,
    schema = div.schema,
    sim.required = c(
      "rep", "pop", "stat", "value", "mutation_rate"
    ),
    empirical.required = c(
      "rep", "pop", "stat", "value", "mutation_rate"
    )
  ) |>
    mutate(
      sequence_scope = if_else(
        data_source == "simulation",
        "simulated_full",
        empirical.scope
      )
    )
}

div.intergenic <- read.diversity.scope(
  "pi_theta_stats_intergenic",
  "onekg_intergenic_callable"
)
div.full.callable <- read.diversity.scope(
  "pi_theta_stats_full_callable_chrom",
  "onekg_full_callable"
)

div <- bind_rows(
  filter(div.intergenic, data_source == "simulation"),
  filter(div.intergenic, data_source == "empirical"),
  filter(div.full.callable, data_source == "empirical")
)
div.views <- make.diversity.views(div)


# descriptive diversity metrics ----
div.metrics <- div |>
  select(
    rep, chrom, pop, pop_role, data_source, analysis_level,
    sequence_scope, stat, value
  ) |>
  pivot_wider(names_from = stat, values_from = value) |>
  mutate(
    pi_theta_ratio = pi / theta,
    pi_minus_theta = pi - theta
  ) |>
  pivot_longer(
    c(pi, theta, pi_theta_ratio, pi_minus_theta),
    names_to = "metric",
    values_to = "value"
  )

div.predictive.input <- bind_rows(lapply(
  names(div.views),
  function(scope) {
    div.views[[scope]] |>
      mutate(empirical_scope = scope)
  }
))

div.metrics.comparison <- bind_rows(lapply(
  names(div.views),
  function(scope) {
    div.metrics |>
      filter(
        sequence_scope %in% c("simulated_full", scope)
      ) |>
      mutate(empirical_scope = scope)
  }
))

div.stat.predictive <- predictive.summary(
  div.predictive.input,
  value.column = "value",
  group.columns = c(
    "analysis_level", "chrom", "pop_role", "stat", "empirical_scope"
  )
)

div.predictive.summary <- predictive.summary(
  div.metrics.comparison,
  value.column = "value",
  group.columns = c(
    "analysis_level", "chrom", "pop_role", "metric",
    "empirical_scope"
  )
)

div.contrasts <- population.contrasts(
  div.metrics.comparison |>
    filter(metric %in% c("pi", "theta", "pi_theta_ratio")),
  value.column = "value",
  group.columns = c(
    "rep", "chrom", "data_source", "analysis_level",
    "empirical_scope", "metric"
  ),
  metric.type = "positive"
)


# diversity figures ----
div.predictive.plot <- ggplot(
  div.predictive.input |> filter(analysis_level == "genome"),
  aes(x = pop_role, y = value, color = data_source)
) +
  geom_boxplot(
    data = \(data) filter(data, data_source == "simulation"),
    outlier.shape = NA
  ) +
  geom_point(
    data = \(data) filter(data, data_source == "empirical"),
    size = 2.5
  ) +
  facet_grid(empirical_scope ~ stat, scales = "free_y") +
  scale_y_continuous(labels = label_number()) +
  labs(
    x = "Population role",
    y = "Statistic value",
    title = "Diversity predictive comparisons"
  ) +
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank())

div.chromosome.envelope <- div.stat.predictive |>
  filter(analysis_level == "chromosome")

div.chromosome.plot <- ggplot(
  div.chromosome.envelope,
  aes(
    x = as.integer(chrom),
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
  geom_line(linewidth = 0.9) +
  geom_point(aes(y = empirical_value), shape = 21, color = "black") +
  facet_grid(empirical_scope ~ stat, scales = "free_y") +
  scale_color_manual(values = POP.COLORS) +
  scale_fill_manual(values = POP.COLORS) +
  labs(
    x = "Chromosome",
    y = "Statistic value",
    title = "Chromosome diversity predictive profiles"
  ) +
  theme_bw(base_size = 16)

div.contrast.plot <- ggplot(
  div.contrasts |> filter(analysis_level == "genome"),
  aes(x = contrast, y = value, color = data_source)
) +
  geom_boxplot(
    data = \(data) filter(data, data_source == "simulation"),
    outlier.shape = NA
  ) +
  geom_point(
    data = \(data) filter(data, data_source == "empirical"),
    size = 2.5
  ) +
  facet_grid(empirical_scope ~ metric, scales = "free_y") +
  labs(
    x = "Population contrast",
    y = "Log ratio",
    title = "Relative diversity differences"
  ) +
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank())


# interactive results ----
print(div.metrics)
print(div.predictive.summary)
print(div.contrasts)
print(div.predictive.plot)
print(div.chromosome.plot)
print(div.contrast.plot)
