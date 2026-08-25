# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# ld_decay.R
# ______________________________________________________________________________


# set up ----
library(tidyverse)
library(nanoparquet)


SIM.SMALL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_small/stats"
SIM.LARGE.DATA.DIR <- "~/scratch/OOA_NAAdmixture_large/stats"
EMPIRICAL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_1kG/stats"
CHROMOSOMES <- as.character(1:22)
SELECTED.CHROMOSOMES <- c("1", "5", "10", "14", "18", "22")
LD.X.LOWER <- 0
LD.X.UPPER <- 250000
LD.X.BREAKS <- seq(LD.X.LOWER, LD.X.UPPER, by = 50000)
PLOT.BASE.SIZE <- 24
PLOT.STYLES <- list(
  population.colors = c(
    AFR = "#56B4E9", ADX = "#4B1FA8", EUR = "#fb8072",
    YRI = "#eec4dc", ASW = "#e44b8d", CEU = "#bb437e"
  ),
  series.linetypes = c(
    Simulation_small = "dashed",
    Simulation_large = "dotdash",
    Empirical = "solid"
  ),
  series.labels = c(
    Simulation_small = "Simulation small",
    Simulation_large = "Simulation large",
    Empirical = "Empirical"
  )
)


# internal functions ----


# map source population labels to shared population roles
add.population.roles <- function(data) {
  data <- data %>%
    mutate(role = case_when(
      pop %in% c("AFR", "YRI") ~ "AFR",
      pop %in% c("ADX", "ASW") ~ "ADX",
      pop %in% c("EUR", "CEU") ~ "EUR",
      TRUE ~ NA_character_
    ))
  if (any(is.na(data$role))) {
    stop("LD data contain an unsupported population label")
  }

  return(data)
}


# standardize one LD table to the shared analysis schema
normalize.ld.table <- function(data, data.type.input) {
  data <- data %>%
    mutate(
      rep = as.numeric(rep),
      chrom = as.character(chrom),
      pop = as.character(pop),
      distance_bin_bp = as.numeric(distance_bin_bp),
      mean_r2 = as.numeric(mean_r2),
      sum_r2 = as.numeric(sum_r2),
      n_pairs = as.numeric(n_pairs),
      data.type = data.type.input
    ) %>%
    add.population.roles()

  return(data)
}


# read chromosome-labelled LD files for one source
read.ld.chromosomes <- function(
    data.directory, chromosomes, data.type.input
) {
  paths <- file.path(
    path.expand(data.directory),
    paste0("ld_decay.chr", chromosomes, ".parquet")
  )
  data <- map2_dfr(paths, chromosomes, function(path, chrom) {
    table <- read_parquet(path)
    table$chrom <- chrom
    return(normalize.ld.table(table, data.type.input))
  })

  return(data)
}


# read the empirical genome-wide LD table
read.empirical.ld.genome <- function(data.directory) {
  data <- read_parquet(file.path(
    path.expand(data.directory), "ld_decay.parquet"
  ))
  data$chrom <- "all"
  data <- normalize.ld.table(data, "Empirical")

  return(data)
}


# reconstruct chromosome or genome curves from pooled sums and pair counts
pool.ld.curves <- function(data, include.chromosome) {
  if (!"role" %in% names(data)) data <- add.population.roles(data)
  grouping.columns <- c(
    "data.type", "rep", "pop", "role", "distance_bin_bp"
  )
  if (include.chromosome) {
    grouping.columns <- c(grouping.columns, "chrom")
  }
  pooled <- data %>%
    group_by(across(all_of(grouping.columns))) %>%
    summarise(
      sum.r2 = sum(sum_r2),
      n.pairs = sum(n_pairs),
      chromosome.count = if (include.chromosome) {
        NA_integer_
      } else {
        n_distinct(chrom)
      },
      .groups = "drop"
    ) %>%
    mutate(
      mean.r2 = if_else(n.pairs > 0, sum.r2 / n.pairs, NA_real_)
    )
  if (!include.chromosome) pooled$chrom <- "all"

  return(pooled)
}


# summarize replicate simulation curves and fixed empirical curves
summarize.ld.curves <- function(data) {
  simulation <- data %>%
    filter(data.type != "Empirical") %>%
    group_by(data.type, pop, role, chrom, distance_bin_bp) %>%
    summarise(
      mean = mean(mean.r2, na.rm = TRUE),
      sd = sd(mean.r2, na.rm = TRUE),
      replicate.count = n_distinct(rep),
      chromosome.count = if (all(is.na(chromosome.count))) {
        NA_integer_
      } else {
        max(chromosome.count, na.rm = TRUE)
      },
      .groups = "drop"
    )
  empirical <- data %>%
    filter(data.type == "Empirical") %>%
    transmute(
      data.type, pop, role, chrom, distance_bin_bp,
      mean = mean.r2, sd = NA_real_, replicate.count = 1L,
      chromosome.count
    ) %>%
    distinct()
  summary <- bind_rows(simulation, empirical) %>%
    mutate(
      role = factor(role, levels = c("AFR", "ADX", "EUR")),
      chrom = factor(
        chrom, levels = c(SELECTED.CHROMOSOMES, "all")
      ),
      data.type = factor(
        data.type,
        levels = c(
          "Simulation_small", "Simulation_large", "Empirical"
        )
      )
    ) %>%
    filter(!is.na(mean))

  return(summary)
}


# add shared scales, labels, guides, and theme to one LD plot
style.ld.plot <- function(plot, title, styles) {
  plot <- plot +
    scale_color_manual(values = styles$population.colors) +
    scale_fill_manual(values = styles$population.colors) +
    scale_linetype_manual(
      values = styles$series.linetypes,
      labels = styles$series.labels
    ) +
    scale_x_continuous(
      limits = c(LD.X.LOWER, LD.X.UPPER), breaks = LD.X.BREAKS
    ) +
    labs(
      x = "Distance Between SNPs (bp)",
      y = expression("Mean " * r^2),
      title = title, color = NULL, fill = NULL, linetype = NULL
    ) +
    guides(
      color = guide_legend(
        order = 1, override.aes = list(linetype = "solid")
      ),
      fill = "none",
      linetype = guide_legend(
        order = 2, override.aes = list(color = "black")
      )
    ) +
    theme_bw(base_size = PLOT.BASE.SIZE) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_blank(),
      legend.key.width = unit(42, "pt"),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )

  return(plot)
}


# add simulation uncertainty ribbons and all source lines
add.ld.geometries <- function(plot, data) {
  simulation <- data %>% filter(data.type != "Empirical")
  plot <- plot +
    geom_ribbon(
      data = simulation,
      aes(
        ymin = pmax(0, mean - 2 * sd),
        ymax = mean + 2 * sd,
        group = interaction(data.type, pop)
      ),
      alpha = 0.2, color = NA
    ) +
    geom_line(
      data = data,
      aes(
        linetype = data.type,
        group = interaction(data.type, pop)
      ),
      linewidth = 1.25
    )

  return(plot)
}


# construct the three genome views and selected-chromosome view
make.ld.plots <- function(data, styles) {
  genome <- data %>% filter(chrom == "all")
  selected <- data %>%
    filter(chrom %in% SELECTED.CHROMOSOMES) %>%
    mutate(
      chrom = factor(
        as.character(chrom), levels = SELECTED.CHROMOSOMES
      )
    )
  genome.base <- ggplot(
    genome,
    aes(
      x = distance_bin_bp, y = mean,
      color = pop, fill = pop
    )
  )
  selected.base <- ggplot(
    selected,
    aes(
      x = distance_bin_bp, y = mean,
      color = pop, fill = pop
    )
  )

  genome.by.data.set <- add.ld.geometries(genome.base, genome) +
    facet_wrap(~data.type)
  genome.by.data.set <- style.ld.plot(
    genome.by.data.set, "Genome-wide LD Decay by Data Set", styles
  )
  genome.by.role <- add.ld.geometries(genome.base, genome) +
    facet_wrap(~role)
  genome.by.role <- style.ld.plot(
    genome.by.role, "Genome-wide LD Decay by Population Role", styles
  )
  genome.populations <- genome.base %>%
    add.ld.geometries(genome) %>%
    style.ld.plot(
      "Genome-wide LD Decay Across Populations", styles
    )
  selected.chromosomes <- add.ld.geometries(
    selected.base, selected
  ) +
    facet_wrap(~chrom, ncol = 3, drop = FALSE)
  selected.chromosomes <- style.ld.plot(
    selected.chromosomes,
    "LD Decay Across Selected Chromosomes", styles
  )
  plots <- list(
    genome.by.data.set = genome.by.data.set,
    genome.by.role = genome.by.role,
    genome.populations = genome.populations,
    selected.chromosomes = selected.chromosomes
  )

  return(plots)
}


# analysis ----


# read all autosomes for both simulation sizes
sim.small.ld.chromosomes <- read.ld.chromosomes(
  SIM.SMALL.DATA.DIR, CHROMOSOMES, "Simulation_small"
)
sim.large.ld.chromosomes <- read.ld.chromosomes(
  SIM.LARGE.DATA.DIR, CHROMOSOMES, "Simulation_large"
)
simulation.ld.chromosomes <- bind_rows(
  sim.small.ld.chromosomes, sim.large.ld.chromosomes
)

# reconstruct chromosome and genome simulation curves from pooled values
simulation.ld.selected <- simulation.ld.chromosomes %>%
  filter(chrom %in% SELECTED.CHROMOSOMES) %>%
  pool.ld.curves(include.chromosome = TRUE)
simulation.ld.genome <- simulation.ld.chromosomes %>%
  pool.ld.curves(include.chromosome = FALSE)

# read and pool the empirical selected-chromosome and genome curves
empirical.ld.selected <- read.ld.chromosomes(
  EMPIRICAL.DATA.DIR, SELECTED.CHROMOSOMES, "Empirical"
) %>%
  pool.ld.curves(include.chromosome = TRUE)
empirical.ld.genome <- read.empirical.ld.genome(
  EMPIRICAL.DATA.DIR
) %>%
  pool.ld.curves(include.chromosome = TRUE)

# summarize curves and construct all four retained LD views
ld.summary <- bind_rows(
  simulation.ld.selected, simulation.ld.genome,
  empirical.ld.selected, empirical.ld.genome
) %>%
  summarize.ld.curves()
ld.plots <- make.ld.plots(ld.summary, PLOT.STYLES)
walk(ld.plots, print)
