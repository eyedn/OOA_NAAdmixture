# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# diversity.R
# ______________________________________________________________________________


# set up ----
library(tidyverse)
library(nanoparquet)


SIM.SMALL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_small/stats"
SIM.LARGE.DATA.DIR <- "~/scratch/OOA_NAAdmixture_large/stats"
EMPIRICAL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_1kG/stats"
CHROMOSOMES <- as.character(1:22)
SELECTED.CHROMOSOMES <- c("1", "5", "10", "14", "18", "22")
PLOT.BASE.SIZE <- 24
PLOT.STYLES <- list(
  population.colors = c(
    AFR = "#56B4E9", ADX = "#4B1FA8", EUR = "#fb8072"
  ),
  series.shapes = c(
    Simulation_small = 21,
    Simulation_large = 22,
    Empirical = 24
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
    stop("Diversity data contain an unsupported population label")
  }

  return(data)
}


# standardize one diversity table to the shared analysis schema
normalize.diversity.table <- function(
    data, data.type.input, mask.input = NA_character_
) {
  data <- data %>%
    mutate(
      rep = as.numeric(rep),
      chrom = as.character(chrom),
      pop = as.character(pop),
      stat = as.character(stat),
      value = as.numeric(value),
      data.type = data.type.input,
      mask = mask.input
    ) %>%
    add.population.roles()

  return(data)
}


# read chromosome-labelled diversity files for one source
read.diversity.chromosomes <- function(
    data.directory, file.family, chromosomes, data.type.input,
    mask.input = NA_character_
) {
  paths <- file.path(
    path.expand(data.directory),
    str_replace(file.family, fixed("{chrom}"), chromosomes)
  )
  data <- map2_dfr(paths, chromosomes, function(path, chrom) {
    table <- read_parquet(path)
    table$chrom <- chrom
    return(normalize.diversity.table(
      table, data.type.input, mask.input
    ))
  })

  return(data)
}


# read one genome-wide empirical diversity file
read.diversity.genome <- function(
    data.directory, file.name, mask.input
) {
  data <- read_parquet(file.path(
    path.expand(data.directory), file.name
  ))
  data$chrom <- "all"
  data <- normalize.diversity.table(
    data, "Empirical", mask.input
  )

  return(data)
}


# calculate replicate means and standard deviations for simulations
summarize.simulation.diversity <- function(data) {
  if (!"role" %in% names(data)) data <- add.population.roles(data)
  summary <- data %>%
    filter(stat %in% c("pi", "theta")) %>%
    group_by(data.type, role, stat, chrom) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      replicate.count = n_distinct(rep),
      .groups = "drop"
    )

  return(summary)
}


# duplicate simulation summaries across empirical mask comparisons
duplicate.simulation.masks <- function(data) {
  duplicated <- crossing(
    data,
    mask = c("Intergenic", "Full callable")
  )

  return(duplicated)
}


# assemble selected-chromosome points and empirical genome references
build.diversity.plot.data <- function(
    simulation.summary, empirical.chromosome, empirical.genome,
    chromosomes
) {
  if (!"role" %in% names(empirical.chromosome)) {
    empirical.chromosome <- add.population.roles(empirical.chromosome)
  }
  if (!"role" %in% names(empirical.genome)) {
    empirical.genome <- add.population.roles(empirical.genome)
  }
  if (!"data.type" %in% names(empirical.chromosome)) {
    empirical.chromosome$data.type <- "Empirical"
  }
  simulation.points <- simulation.summary %>%
    filter(chrom %in% chromosomes) %>%
    duplicate.simulation.masks() %>%
    transmute(
      data.type, role, stat, chrom, mask,
      estimate = mean, sd
    )
  empirical.points <- empirical.chromosome %>%
    filter(chrom %in% chromosomes, stat %in% c("pi", "theta")) %>%
    transmute(
      data.type, role, stat, chrom, mask,
      estimate = value, sd = NA_real_
    )
  points <- bind_rows(simulation.points, empirical.points) %>%
    mutate(
      chrom = factor(chrom, levels = chromosomes),
      role = factor(role, levels = c("AFR", "ADX", "EUR")),
      data.type = factor(
        data.type,
        levels = c(
          "Simulation_small", "Simulation_large", "Empirical"
        )
      ),
      mask = factor(mask, levels = c("Intergenic", "Full callable")),
      stat = factor(stat, levels = c("pi", "theta"))
    )
  genome.lines <- empirical.genome %>%
    filter(chrom == "all", stat %in% c("pi", "theta")) %>%
    transmute(role, stat, mask, estimate = value) %>%
    mutate(
      role = factor(role, levels = c("AFR", "ADX", "EUR")),
      mask = factor(mask, levels = c("Intergenic", "Full callable")),
      stat = factor(stat, levels = c("pi", "theta"))
    )

  return(list(points = points, genome.lines = genome.lines))
}


# construct the selected-chromosome diversity plot
make.diversity.plot <- function(points, genome.lines, styles) {
  dodge <- position_dodge(width = 0.75)
  plot <- ggplot(
    points,
    aes(
      x = chrom, y = estimate, color = role, fill = role,
      shape = data.type, group = interaction(role, data.type)
    )
  ) +
    geom_hline(
      data = genome.lines,
      aes(yintercept = estimate, color = role),
      linetype = "dotted", linewidth = 0.9
    ) +
    geom_errorbar(
      data = points,
      aes(ymin = estimate - 2 * sd, ymax = estimate + 2 * sd),
      position = dodge, width = 0.15, linewidth = 0.8,
      na.rm = TRUE
    ) +
    geom_point(
      position = dodge, color = "black", stroke = 1, size = 3.5,
      aes(fill = role)
    ) +
    facet_grid(
      stat ~ mask, scales = "free_y",
      labeller = labeller(
        stat = c(pi = "π", theta = "θ[w]")
      )
    ) +
    scale_color_manual(values = styles$population.colors) +
    scale_fill_manual(values = styles$population.colors) +
    scale_shape_manual(
      values = styles$series.shapes,
      labels = styles$series.labels
    ) +
    labs(
      x = "Chromosome", y = NULL,
      title = "Genetic Diversity Across Selected Chromosomes",
      color = NULL, fill = NULL, shape = NULL
    ) +
    guides(
      color = guide_legend(
        order = 1, override.aes = list(shape = 21, size = 3.5)
      ),
      fill = "none",
      shape = guide_legend(
        order = 2,
        override.aes = list(fill = "white", color = "black")
      )
    ) +
    theme_bw(base_size = PLOT.BASE.SIZE) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )

  return(plot)
}


# analysis ----


# read selected-chromosome simulation diversity estimates
sim.small.diversity <- read.diversity.chromosomes(
  SIM.SMALL.DATA.DIR, "pi_theta_stats.chr{chrom}.parquet",
  SELECTED.CHROMOSOMES, "Simulation_small"
)
sim.large.diversity <- read.diversity.chromosomes(
  SIM.LARGE.DATA.DIR, "pi_theta_stats.chr{chrom}.parquet",
  SELECTED.CHROMOSOMES, "Simulation_large"
)

# read selected-chromosome empirical diversity estimates
emp.intergenic.chromosome <- read.diversity.chromosomes(
  EMPIRICAL.DATA.DIR,
  "pi_theta_stats_intergenic.chr{chrom}.parquet",
  SELECTED.CHROMOSOMES, "Empirical", "Intergenic"
)
emp.full.callable.chromosome <- read.diversity.chromosomes(
  EMPIRICAL.DATA.DIR,
  "pi_theta_stats_full_callable_chrom.chr{chrom}.parquet",
  SELECTED.CHROMOSOMES, "Empirical", "Full callable"
)

# read genome-wide empirical diversity references
emp.intergenic.genome <- read.diversity.genome(
  EMPIRICAL.DATA.DIR, "pi_theta_stats_intergenic.parquet",
  "Intergenic"
)
emp.full.callable.genome <- read.diversity.genome(
  EMPIRICAL.DATA.DIR,
  "pi_theta_stats_full_callable_chrom.parquet", "Full callable"
)

# summarize sources and construct the primary diversity plot
simulation.diversity.summary <- bind_rows(
  sim.small.diversity, sim.large.diversity
) %>%
  summarize.simulation.diversity()
diversity.plot.data <- build.diversity.plot.data(
  simulation.diversity.summary,
  bind_rows(
    emp.intergenic.chromosome, emp.full.callable.chromosome
  ),
  bind_rows(emp.intergenic.genome, emp.full.callable.genome),
  SELECTED.CHROMOSOMES
)
diversity.plot <- make.diversity.plot(
  diversity.plot.data$points,
  diversity.plot.data$genome.lines,
  PLOT.STYLES
)
print(diversity.plot)
