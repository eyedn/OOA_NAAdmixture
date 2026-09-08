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
library(glue)
library(nanoparquet)


SIM.SMALL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_small/stats"
SIMDOWN.SMALL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_smallOnekgDownsample/stats"
SIM.LARGE.DATA.DIR <- "~/scratch/OOA_NAAdmixture_large/stats"
SIMDOWN.LARGE.DATA.DIR <- "~/scratch/OOA_NAAdmixture_largeOnekgDownsample/stats"
EMPIRICAL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_1kG/stats"
CHROMOSOMES <- as.character(1:22)
PLOT.CHROMOSOME <- "1"
SOURCE.LEVELS <- c(
  "Simulation_small", "Simulation_large", "Simulation_small_simDown",
  "Simulation_large_simDown", "Empirical"
)
LD.X.LOWER <- 0
LD.X.UPPER <- 250000
LD.X.BREAKS <- seq(LD.X.LOWER, LD.X.UPPER, by = 50000)
PLOT.BASE.SIZE <- 24
PLOT.STYLES <- list(
  population.colors = c(
    AFR = "#56B4E9", ADX = "#4B1FA8", EUR = "#fb8072",
    YRI = "#eec4dc", ASW = "#e44b8d", CEU = "#bb437e"
  ),
  source.colors = c(
    Simulation_small = "#9A83CE", Simulation_large = "#32146F",
    Simulation_small_simDown = "#6F55B5",
    Simulation_large_simDown = "#4B1FA8"
  ),
  series.labels = c(
    Simulation_small = "Sm. Sim.",
    Simulation_large = "Lg. Sim.",
    Simulation_small_simDown = "Sm. D. Sim.",
    Simulation_large_simDown = "Lg. D. Sim.",
    Empirical = "Emp."
  )
)


# internal functions ----


# describe the selected LD estimator, chromosome, and uncertainty
ld.plot.subtitle <- function(chromosome) {
  subtitle <- paste0(
    "Rogers–Huff r²; chrom.", chromosome
  )
  return(subtitle)
}


# retain source-specific populations before pooling sufficient statistics
apply.ld.source.contract <- function(data) {
  retained <- data %>%
    filter(
      (data.type %in% c(
        "Simulation_small", "Simulation_small_simDown"
      ) & pop %in% c("AFR", "ADX", "EUR")) |
        (data.type %in% c(
          "Simulation_large", "Simulation_large_simDown"
        ) & pop == "ADX") |
        (data.type == "Empirical" & pop %in% c("AFR", "ADX", "EUR",
                                                "YRI", "ASW", "CEU"))
    ) %>%
    mutate(data.type = factor(data.type, levels = SOURCE.LEVELS))

  return(retained)
}


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
    apply.ld.source.contract() %>%
    add.population.roles()

  return(data)
}


# read chromosome-labelled LD files for one source
read.ld.chromosomes <- function(
    data.directory, chromosomes, data.type.input
) {
  paths <- file.path(
    path.expand(data.directory),
    glue("ld_decay.chr{chromosomes}.parquet")
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
summarize.ld.curves <- function(data, chromosome) {
  chromosome <- as.character(chromosome)
  if (length(chromosome) != 1L) {
    stop("LD summarization requires exactly one chromosome")
  }
  scoped <- data %>% filter(as.character(chrom) == chromosome)
  if (!nrow(scoped)) {
    stop("LD data do not contain the requested chromosome")
  }
  if (!identical(unique(as.character(scoped$chrom)), chromosome)) {
    stop("LD summary contains data outside the requested chromosome")
  }
  simulation <- scoped %>%
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
  empirical <- scoped %>%
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
      chrom = as.character(chrom),
      data.type = factor(
        data.type,
        levels = SOURCE.LEVELS
      )
    ) %>%
    filter(!is.na(mean))

  return(summary)
}


# add shared scales, labels, guides, and theme to one LD plot
style.ld.plot <- function(plot, title, subtitle, styles, view) {
  plot <- plot +
    scale_color_manual(values = if (view == "simulated") {
      styles$source.colors
    } else {
      styles$population.colors
    }) +
    scale_fill_manual(values = if (view == "simulated") {
      styles$source.colors
    } else {
      styles$population.colors
    }) +
    scale_x_continuous(
      limits = c(LD.X.LOWER, LD.X.UPPER), breaks = LD.X.BREAKS
    ) +
    labs(
      x = "Distance Between SNPs (bp)",
      y = expression("Mean " * r^2),
      title = title, subtitle = subtitle,
      color = NULL, fill = NULL
    ) +
    guides(
      color = guide_legend(order = 1),
      fill = "none"
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
        group = interaction(data.type, pop)
      ),
      linewidth = 1.25
    )

  return(plot)
}


# construct the source-group view for exactly one chromosome
filter.plot.view <- function(data, view) {
  sources <- switch(
    view,
    small_empirical = c(
      "Simulation_small", "Simulation_small_simDown", "Empirical"
    ),
    simulated = c(
      "Simulation_small", "Simulation_large", "Simulation_small_simDown",
      "Simulation_large_simDown"
    ),
    stop("Unsupported LD plot view: ", view)
  )
  filtered <- data %>%
    filter(as.character(data.type) %in% sources) %>%
    filter(view != "simulated" | pop == "ADX") %>%
    mutate(data.type = factor(as.character(data.type), levels = sources)) %>%
    arrange(data.type)
  return(filtered)
}


# construct one scoped LD view for exactly one chromosome
make.ld.plot <- function(data, chromosome, styles, view) {
  if (length(chromosome) != 1L) {
    stop("LD plotting requires exactly one chromosome")
  }
  chromosome <- as.character(chromosome)
  if (!chromosome %in% CHROMOSOMES) {
    stop("LD plotting requires exactly one autosomal chromosome")
  }
  plot.data <- filter.plot.view(data, view) %>%
    filter(as.character(chrom) == chromosome) %>%
    mutate(
      plot.key = if (view == "simulated") {
        factor(as.character(data.type), levels = c(
          "Simulation_small", "Simulation_large",
          "Simulation_small_simDown", "Simulation_large_simDown"
        ))
      } else {
        factor(as.character(pop), levels = c(
          "AFR", "ADX", "EUR", "YRI", "ASW", "CEU"
        ))
      }
    )
  if (!nrow(plot.data)) {
    stop("LD data do not contain the selected chromosome")
  }
  plot <- ggplot(
    plot.data,
    aes(
      x = distance_bin_bp, y = mean,
      color = plot.key, fill = plot.key
    )
  )
  plot <- add.ld.geometries(plot, plot.data) +
    facet_wrap(~data.type, nrow = 1, drop = FALSE,
               labeller = labeller(data.type = styles$series.labels))
  plot <- style.ld.plot(
    plot, if (view == "simulated") {
      "LD Decay: Simulated ADX"
    } else {
      "LD Decay: Small Simulation and Empirical"
    }, ld.plot.subtitle(chromosome), styles, view
  )

  return(plot)
}


# analysis ----


# read the selected chromosome for all four simulation sources
sim.small.ld.chromosomes <- read.ld.chromosomes(
  SIM.SMALL.DATA.DIR, PLOT.CHROMOSOME, "Simulation_small"
)
simDown.small.ld.chromosomes <- read.ld.chromosomes(
  SIMDOWN.SMALL.DATA.DIR, PLOT.CHROMOSOME,
  "Simulation_small_simDown"
)
sim.large.ld.chromosomes <- read.ld.chromosomes(
  SIM.LARGE.DATA.DIR, PLOT.CHROMOSOME, "Simulation_large"
)
simDown.large.ld.chromosomes <- read.ld.chromosomes(
  SIMDOWN.LARGE.DATA.DIR, PLOT.CHROMOSOME,
  "Simulation_large_simDown"
)
simulation.ld.chromosomes <- bind_rows(
  sim.small.ld.chromosomes, simDown.small.ld.chromosomes,
  sim.large.ld.chromosomes, simDown.large.ld.chromosomes
)

# pool every simulation chromosome curve from producer sufficient statistics
simulation.ld.selected <- simulation.ld.chromosomes %>%
  pool.ld.curves(include.chromosome = TRUE)

# read and pool the empirical selected chromosome
empirical.ld.selected <- read.ld.chromosomes(
  EMPIRICAL.DATA.DIR, PLOT.CHROMOSOME, "Empirical"
) %>%
  pool.ld.curves(include.chromosome = TRUE)

# summarize curves and construct the grouped source view
ld.summary <- bind_rows(
  simulation.ld.selected, empirical.ld.selected
) %>%
  summarize.ld.curves(PLOT.CHROMOSOME)
ld.small.empirical.plot <- make.ld.plot(
  ld.summary, PLOT.CHROMOSOME, PLOT.STYLES, "small_empirical"
)
ld.simulated.plot <- make.ld.plot(
  ld.summary, PLOT.CHROMOSOME, PLOT.STYLES, "simulated"
)
print(ld.small.empirical.plot)
print(ld.simulated.plot)
