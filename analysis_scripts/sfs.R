# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# sfs.R
# ______________________________________________________________________________

# pattern: Imperative Shell


library(tidyverse)
library(scales)

source("analysis_scripts/analysis_utils.R")


# constants ----
RARE.MAF.CUTOFF <- 0.05


# read and harmonize SFS tables ----
sfs.schema <- c(
  rep = "integer",
  chrom = "character",
  pop = "character",
  derived_allele_count = "integer",
  minor_allele_count = "integer",
  count = "double"
)

sfs <- read.harmonized.table(
  sim.table = "sfs",
  empirical.table = "sfs",
  schema = sfs.schema,
  sim.required = c(
    "rep", "pop", "derived_allele_count", "count"
  ),
  empirical.required = c(
    "rep", "pop", "minor_allele_count", "count"
  )
)

sfs.chroms <- discover.chromosomes(ONEKG.STATS.DIR, "sfs")
variant.qc.chroms <- discover.chromosomes(ONEKG.STATS.DIR, "variant_qc")
if (!identical(sfs.chroms, variant.qc.chroms)) {
  stop(
    "Empirical SFS and variant-QC chromosome outputs do not match",
    call. = FALSE
  )
}

# empirical target haplotype counts ----
variant.qc.path <- require.genome.file(ONEKG.STATS.DIR, "variant_qc")
variant.qc.paths <- c(
  file.path(
    ONEKG.STATS.DIR,
    paste0("variant_qc.chr", variant.qc.chroms, ".parquet")
  ),
  variant.qc.path
)
variant.qc.tables <- read.parquet.paths(variant.qc.paths)
variant.qc.chrom.labels <- c(variant.qc.chroms, "all")
target.columns <- c(
  AFR = "retained_YRI_samples",
  EUR = "retained_CEU_samples",
  ADMIXED = "retained_ASW_samples"
)
target.haplotype.map <- map2_dfr(
  variant.qc.tables,
  variant.qc.chrom.labels,
  function(variant.qc, chrom) {
    missing.target.columns <- setdiff(target.columns, names(variant.qc))
    if (length(missing.target.columns) > 0) {
      stop(
        "Variant-QC table is missing retained sample count(s): ",
        paste(missing.target.columns, collapse = ", "),
        call. = FALSE
      )
    }
    tibble(
      chrom = chrom,
      pop_role = names(target.columns),
      num_haplotypes = 2L * vapply(
        target.columns,
        function(column) {
          values <- unique(as.integer(variant.qc[[column]]))
          values <- values[!is.na(values)]
          if (length(values) != 1 || values <= 0) {
            stop(
              "Variant-QC retained sample count must be one ",
              "positive value: ",
              column,
              call. = FALSE
            )
          }
          values
        },
        integer(1)
      )
    )
  }
)


# callable spans ----
span.schema <- c(
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

span.tables <- read.harmonized.table(
  sim.table = "pi_theta_stats",
  empirical.table = "pi_theta_stats_full_callable_chrom",
  schema = span.schema,
  sim.required = c("rep", "pop", "stat"),
  empirical.required = c("rep", "pop", "stat")
) |>
  filter(analysis_level == "chromosome", stat == "pi") |>
  distinct(data_source, rep, pop_role, chrom, span)

if (any(is.na(span.tables$span))) {
  stop("Chromosome pi/theta span metadata is missing", call. = FALSE)
}

genome.spans <- span.tables |>
  group_by(data_source, rep, pop_role) |>
  summarize(span = sum(span), .groups = "drop") |>
  mutate(chrom = "all")

span.map <- bind_rows(span.tables, genome.spans)


# project simulation spectra and normalize both sources ----
project.one.spectrum <- function(data, chrom, pop.role) {
  source.haplotypes <- infer.source.haplotypes(data)
  target <- target.haplotype.map |>
    filter(
      .data$chrom == .env$chrom,
      .data$pop_role == .env$pop.role
    ) |>
    pull(num_haplotypes)
  if (length(target) != 1) {
    stop("Missing unique target haplotype count", call. = FALSE)
  }
  projected <- project.unfolded.sfs(
    data,
    source.haplotypes,
    target
  )
  fold.projected.sfs(projected, target) |>
    mutate(num_haplotypes = target)
}

sim.sfs <- sfs |>
  filter(data_source == "simulation") |>
  group_by(
    rep, chrom, pop, pop_role, data_source, analysis_level
  ) |>
  group_modify(
    ~project.one.spectrum(
      .x,
      .y$chrom[[1]],
      .y$pop_role[[1]]
    )
  ) |>
  ungroup()

empirical.sfs <- sfs |>
  filter(data_source == "empirical") |>
  transmute(
    rep, chrom, pop, pop_role, data_source, analysis_level,
    minor_allele_count,
    count
  ) |>
  left_join(
    target.haplotype.map,
    by = c("chrom", "pop_role")
  )
if (any(is.na(empirical.sfs$num_haplotypes))) {
  stop("Empirical SFS target haplotype metadata is missing", call. = FALSE)
}

sfs.comparable <- bind_rows(sim.sfs, empirical.sfs) |>
  left_join(
    span.map,
    by = c("data_source", "rep", "pop_role", "chrom")
  )
if (any(is.na(sfs.comparable$span))) {
  stop("Required SFS callable-span metadata is missing", call. = FALSE)
}

sfs.normalized <- sfs.comparable |>
  group_by(
    rep, chrom, pop, pop_role, data_source,
    analysis_level, num_haplotypes, span
  ) |>
  group_modify(
    ~normalize.folded.sfs(
      .x,
      .y$num_haplotypes[[1]],
      .y$span[[1]]
    )
  ) |>
  ungroup()


# folded-spectrum summaries and distances ----
sfs.metrics <- sfs.normalized |>
  group_by(
    rep, chrom, pop, pop_role, data_source, analysis_level
  ) |>
  summarize(
    singleton_fraction = sum(
      fraction[minor_allele_count == 1],
      na.rm = TRUE
    ),
    doubleton_fraction = sum(
      fraction[minor_allele_count == 2],
      na.rm = TRUE
    ),
    mac_le_5_fraction = sum(
      fraction[minor_allele_count <= 5],
      na.rm = TRUE
    ),
    rare_maf_fraction = sum(
      fraction[minor_allele_frequency <= RARE.MAF.CUTOFF],
      na.rm = TRUE
    ),
    .groups = "drop"
  )

sfs.metric.predictive <- sfs.metrics |>
  pivot_longer(
    ends_with("_fraction"),
    names_to = "metric",
    values_to = "value"
  ) |>
  predictive.summary(
    value.column = "value",
    group.columns = c(
      "analysis_level", "chrom", "pop_role", "metric"
    )
  )

empirical.spectra <- sfs.normalized |>
  filter(data_source == "empirical") |>
  select(chrom, pop_role, minor_allele_count, fraction)

sfs.distance <- sfs.normalized |>
  filter(data_source == "simulation") |>
  group_by(rep, chrom, pop_role, analysis_level) |>
  group_modify(
    ~{
      empirical <- empirical.spectra |>
        filter(
          chrom == .y$chrom[[1]],
          pop_role == .y$pop_role[[1]]
        )
      tibble(
        hellinger_distance = spectrum.hellinger(.x, empirical)
      )
    }
  ) |>
  ungroup()

sfs.rare.contrasts <- population.contrasts(
  sfs.metrics,
  value.column = "rare_maf_fraction",
  group.columns = c(
    "rep", "chrom", "data_source", "analysis_level"
  ),
  metric.type = "positive"
)


# folded-SFS figures ----
sfs.envelope <- sfs.normalized |>
  filter(data_source == "simulation") |>
  group_by(
    analysis_level, chrom, pop_role, minor_allele_frequency
  ) |>
  summarize(
    median = median(fraction, na.rm = TRUE),
    lower = quantile(fraction, 0.025, na.rm = TRUE),
    upper = quantile(fraction, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

sfs.empirical <- sfs.normalized |>
  filter(data_source == "empirical")

sfs.shape.plot <- ggplot(
  sfs.envelope |> filter(analysis_level == "genome"),
  aes(
    x = minor_allele_frequency,
    y = median,
    color = pop_role,
    fill = pop_role
  )
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.15,
    color = NA
  ) +
  geom_line(linewidth = 1) +
  geom_point(
    data = sfs.empirical |> filter(analysis_level == "genome"),
    aes(y = fraction),
    size = 1.6
  ) +
  scale_color_manual(values = POP.COLORS) +
  scale_fill_manual(values = POP.COLORS) +
  scale_y_continuous(trans = pseudo_log_trans()) +
  labs(
    x = "Minor allele frequency",
    y = "Fraction of segregating sites",
    title = "Projected and folded 1D SFS"
  ) +
  theme_bw(base_size = 16)

sfs.density.plot <- ggplot(
  sfs.normalized |> filter(analysis_level == "genome"),
  aes(
    x = minor_allele_frequency,
    y = sites_per_mb,
    color = data_source
  )
) +
  geom_line(aes(group = interaction(rep, data_source)), alpha = 0.35) +
  facet_wrap(~pop_role, scales = "free_y") +
  scale_y_continuous(trans = pseudo_log_trans()) +
  labs(
    x = "Minor allele frequency",
    y = "Sites per callable Mb",
    title = "Folded-SFS site density"
  ) +
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank())

sfs.rare.contrast.plot <- ggplot(
  sfs.rare.contrasts |> filter(analysis_level == "genome"),
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
  labs(
    x = "Population contrast",
    y = "Log ratio of rare-variant mass",
    title = "Rare-variant population contrasts"
  ) +
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank())


# interactive results ----
# Raw counts enter comparisons only after simulation projection and folding.
print(sfs.metrics)
print(sfs.metric.predictive)
print(sfs.distance)
print(sfs.rare.contrasts)
print(sfs.shape.plot)
print(sfs.density.plot)
print(sfs.rare.contrast.plot)
