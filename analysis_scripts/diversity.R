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
SIMDOWN.SMALL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_smallOnekgDownsample/stats"
SIM.LARGE.DATA.DIR <- "~/scratch/OOA_NAAdmixture_large/stats"
SIMDOWN.LARGE.DATA.DIR <- "~/scratch/OOA_NAAdmixture_largeOnekgDownsample/stats"
EMPIRICAL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_1kG/stats"
CHROMOSOMES <- as.character(1:22)
SELECTED.CHROMOSOMES <- c("1")
SOURCE.LEVELS <- c(
  "Simulation_small", "Simulation_large", "Simulation_small_simDown",
  "Simulation_large_simDown", "Empirical"
)
PLOT.BASE.SIZE <- 24
PLOT.STYLES <- list(
  population.colors = c(
    AFR = "#56B4E9", ADX = "#4B1FA8", EUR = "#fb8072",
    YRI = "#eec4dc", ASW = "#e44b8d", CEU = "#bb437e"
  ),
  fill.colors = c(
    "small AFR" = "#9BD5F2", "small ADX" = "#9A83CE",
    "small EUR" = "#FBB4AE", "large ADX" = "#32146F",
    "small simDown AFR" = "#56B4E9",
    "small simDown ADX" = "#6F55B5",
    "small simDown EUR" = "#FB8072",
    "large simDown ADX" = "#4B1FA8",
    "empirical YRI" = "#EEC4DC", "empirical ASW" = "#E44B8D",
    "empirical CEU" = "#BB437E"
  ),
  fill.labels = c(
    "small AFR" = "Sm. Sim. AFR",
    "small ADX" = "Sm. Sim. ADX",
    "small EUR" = "Sm. Sim. EUR",
    "large ADX" = "Lg. Sim. ADX",
    "small simDown AFR" = "Sm. D. Sim. AFR",
    "small simDown ADX" = "Sm. D. Sim. ADX",
    "small simDown EUR" = "Sm. D. Sim. EUR",
    "large simDown ADX" = "Lg. D. Sim. ADX",
    "empirical YRI" = "Emp. YRI",
    "empirical ASW" = "Emp. ASW",
    "empirical CEU" = "Emp. CEU"
  ),
  empirical.colors = c(
    YRI = "#EEC4DC", ASW = "#E44B8D", CEU = "#BB437E"
  )
)


# internal functions ----


# describe the diversity estimators, scope, and uncertainty concisely
diversity.plot.subtitle <- function(chromosomes) {
  subtitle <- paste0(
    "π and Watterson’s θ · chromosomes ",
    paste(chromosomes, collapse = ", "), " · ±2 SD simulations"
  )
  return(subtitle)
}


# retain source-specific populations before any replicate summaries
apply.diversity.source.contract <- function(data) {
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
    apply.diversity.source.contract() %>%
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
    group_by(data.type, pop, role, stat, chrom, mask) %>%
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
  original <- data %>%
    filter(!str_detect(as.character(data.type), "simDown$")) %>%
    select(-mask) %>%
    crossing(mask = c("Intergenic", "Full callable"))
  simDown <- data %>%
    filter(str_detect(as.character(data.type), "simDown$"))
  duplicated <- bind_rows(original, simDown)

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
      data.type, pop, role, stat, chrom, mask,
      estimate = mean, sd, replicate.count
    )
  empirical.points <- empirical.chromosome %>%
    filter(chrom %in% chromosomes, stat %in% c("pi", "theta")) %>%
    transmute(
      data.type, pop, role, stat, chrom, mask,
      estimate = value, sd = NA_real_, replicate.count = 1L
    )
  points <- bind_rows(simulation.points, empirical.points) %>%
    mutate(
      chrom = factor(chrom, levels = chromosomes),
      pop = factor(pop, levels = c("AFR", "ADX", "EUR", "YRI", "ASW", "CEU")),
      data.type = factor(
        data.type, levels = SOURCE.LEVELS
      ),
      mask = factor(mask, levels = c("Intergenic", "Full callable")),
      stat = factor(stat, levels = c("pi", "theta")),
      fill.key = factor(case_when(
        data.type == "Simulation_small" ~ paste("small", pop),
        data.type == "Simulation_large" ~ "large ADX",
        data.type == "Simulation_small_simDown" ~
          paste("small simDown", pop),
        data.type == "Simulation_large_simDown" ~ "large simDown ADX",
        data.type == "Empirical" ~ paste("empirical", pop)
      ), levels = names(PLOT.STYLES$fill.colors))
    )
  genome.lines <- empirical.genome %>%
    filter(chrom == "all", stat %in% c("pi", "theta")) %>%
    transmute(pop, role, stat, mask, estimate = value) %>%
    mutate(
      pop = factor(pop, levels = c("AFR", "ADX", "EUR", "YRI", "ASW", "CEU")),
      mask = factor(mask, levels = c("Intergenic", "Full callable")),
      stat = factor(stat, levels = c("pi", "theta"))
    )

  return(list(points = points, genome.lines = genome.lines))
}


# retain one diversity view and its eligible empirical references
filter.diversity.plot.view <- function(points, genome.lines, view) {
  sources <- switch(
    view,
    small_empirical = c(
      "Simulation_small", "Simulation_small_simDown", "Empirical"
    ),
    simulated = c(
      "Simulation_small", "Simulation_large", "Simulation_small_simDown",
      "Simulation_large_simDown"
    ),
    stop("Unsupported diversity plot view: ", view)
  )
  view.points <- points %>%
    filter(as.character(data.type) %in% sources) %>%
    filter(view != "simulated" | pop == "ADX") %>%
    mutate(
      data.type = factor(as.character(data.type), levels = sources),
      fill.key = factor(
        as.character(fill.key),
        levels = names(PLOT.STYLES$fill.colors)[
          names(PLOT.STYLES$fill.colors) %in% fill.key
        ]
      )
    ) %>%
    arrange(data.type, fill.key) %>%
    droplevels()
  view.lines <- if (view == "small_empirical") genome.lines else {
    genome.lines[0, , drop = FALSE]
  }
  return(list(points = view.points, genome.lines = view.lines))
}


# construct the selected-chromosome diversity plot
make.diversity.plot <- function(points, genome.lines, styles) {
  chromosome.scope <- points$chrom %>%
    as.character() %>%
    unique()
  subtitle <- diversity.plot.subtitle(chromosome.scope)
  dodge <- position_dodge(width = 0.75)
  fill.keys <- names(styles$fill.colors)[
    names(styles$fill.colors) %in% as.character(points$fill.key)
  ]
  plot <- ggplot(
    points,
    aes(
      x = chrom, y = estimate, fill = fill.key,
      group = interaction(pop, data.type)
    )
  ) +
    geom_hline(
      data = genome.lines,
      aes(yintercept = estimate, color = pop),
      linetype = "dotted", linewidth = 0.9
    ) +
    geom_col(
      position = dodge, width = 0.7,
      color = "black", linewidth = 0.15
    ) +
    geom_errorbar(
      data = points,
      aes(ymin = estimate - 2 * sd, ymax = estimate + 2 * sd),
      position = dodge, width = 0.15, linewidth = 0.8,
      na.rm = TRUE
    ) +
    facet_grid(
      stat ~ mask, scales = "free_y",
      labeller = labeller(
        stat = c(pi = "π", theta = "θ[w]")
      )
    ) +
    scale_color_manual(values = styles$empirical.colors) +
    scale_fill_manual(
      values = styles$fill.colors[fill.keys],
      breaks = fill.keys,
      labels = styles$fill.labels[fill.keys],
      limits = fill.keys,
      drop = FALSE
    ) +
    scale_y_continuous(
      labels = scales::label_number(accuracy = 0.00001)
    ) +
    labs(
      x = "Chromosome", y = NULL,
      title = "Genetic Diversity Across Selected Chromosomes",
      subtitle = subtitle,
      color = NULL, fill = NULL, shape = NULL
    ) +
    guides(
      color = "none",
      fill = guide_legend(order = 1)
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
simDown.small.intergenic.diversity <- read.diversity.chromosomes(
  SIMDOWN.SMALL.DATA.DIR,
  "pi_theta_stats_intergenic.chr{chrom}.parquet",
  SELECTED.CHROMOSOMES, "Simulation_small_simDown", "Intergenic"
)
simDown.small.full.callable.diversity <- read.diversity.chromosomes(
  SIMDOWN.SMALL.DATA.DIR,
  "pi_theta_stats_full_callable_chrom.chr{chrom}.parquet",
  SELECTED.CHROMOSOMES, "Simulation_small_simDown", "Full callable"
)
sim.large.diversity <- read.diversity.chromosomes(
  SIM.LARGE.DATA.DIR, "pi_theta_stats.chr{chrom}.parquet",
  SELECTED.CHROMOSOMES, "Simulation_large"
)
simDown.large.intergenic.diversity <- read.diversity.chromosomes(
  SIMDOWN.LARGE.DATA.DIR,
  "pi_theta_stats_intergenic.chr{chrom}.parquet",
  SELECTED.CHROMOSOMES, "Simulation_large_simDown", "Intergenic"
)
simDown.large.full.callable.diversity <- read.diversity.chromosomes(
  SIMDOWN.LARGE.DATA.DIR,
  "pi_theta_stats_full_callable_chrom.chr{chrom}.parquet",
  SELECTED.CHROMOSOMES, "Simulation_large_simDown", "Full callable"
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

# summarize sources and construct the two diversity views
simulation.diversity.summary <- bind_rows(
  sim.small.diversity,
  simDown.small.intergenic.diversity,
  simDown.small.full.callable.diversity,
  sim.large.diversity,
  simDown.large.intergenic.diversity,
  simDown.large.full.callable.diversity
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
diversity.small.empirical.data <- filter.diversity.plot.view(
  diversity.plot.data$points, diversity.plot.data$genome.lines,
  "small_empirical"
)
diversity.simulated.data <- filter.diversity.plot.view(
  diversity.plot.data$points, diversity.plot.data$genome.lines, "simulated"
)
diversity.small.empirical.plot <- make.diversity.plot(
  diversity.small.empirical.data$points,
  diversity.small.empirical.data$genome.lines,
  PLOT.STYLES
)
diversity.simulated.plot <- make.diversity.plot(
  diversity.simulated.data$points,
  diversity.simulated.data$genome.lines,
  PLOT.STYLES
)
print(diversity.small.empirical.plot)
print(diversity.simulated.plot)
