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


SIM.TC.DATA.DIR <- "~/scratch/OOA_NAAdmixture_2T12Consistent/stats"
SIMDOWN.TC.DATA.DIR <- "~/scratch/OOA_NAAdmixture_2T12ConsistentOnekgDownsample/stats"
SIM.LG.DATA.DIR <- "~/scratch/OOA_NAAdmixture_largeGrowth/stats"
SIMDOWN.LG.DATA.DIR <- "~/scratch/OOA_NAAdmixture_largeGrowthOnekgDownsample/stats"
EMPIRICAL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_1kG/stats"
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
CHROMOSOMES <- as.character(1:22)
SELECTED.CHROMOSOMES <- c("1", "10", "20")
SOURCE.LEVELS <- c(
  "Simulation_2T12Consistent", "Simulation_2T12Consistent_simDown",
  "Simulation_largeGrowth", "Simulation_largeGrowth_simDown",
  "Empirical"
  )
SOURCE.DISPLAY.LEVELS <- c("T.C.", "T.C.D.", "L.G.", "L.G.D.", "Emp.")
SOURCE.LABELS <- setNames(SOURCE.DISPLAY.LEVELS, SOURCE.LEVELS)
POPULATION.LEVELS <- c("AFR", "ADX", "EUR", "YRI", "ASW", "CEU")
PLOT.CONFIGS <- list(
  TC.1kG = SOURCE.LEVELS[c(1, 5)],
  TC.TCD = SOURCE.LEVELS[c(1, 2)],
  TCD.1kG = SOURCE.LEVELS[c(2, 5)],
  onlyADX = SOURCE.LEVELS[1:4],
  tcd.1kg = SOURCE.LEVELS[c(2, 5)],
  tc.tcd.1kg = SOURCE.LEVELS[c(1, 2, 5)],
  all.datatypes.adx.asw = SOURCE.LEVELS
  )
RANDOM.SEED <- 123L
BOOTSTRAP.REPLICATES <- 1000L
PLOT.BASE.SIZE <- 24
CATEGORICAL.BAR.DODGE <- 0.9
CATEGORICAL.BAR.WIDTH <- 0.8
CATEGORICAL.BAR.LINEWIDTH <- 1
DENSE.BAR.WIDTH.MULTIPLIER <- 0.8
DENSE.BAR.LINEWIDTH <- 0.75
PLOT.STYLES <- list(
  population.colors = c(
    AFR = "#56B4E9", ADX = "#4B1FA8", EUR = "#fb8072",
    YRI = "#eec4dc", ASW = "#e44b8d", CEU = "#bb437e"
    ),
  fill.colors = c(
    "TC AFR" = "#9BD5F2", "TC ADX" = "#9A83CE",
    "TC EUR" = "#FBB4AE",
    "TC D. AFR" = "#56B4E9",
    "TC D. ADX" = "#6F55B5",
    "TC D. EUR" = "#FB8072",
    "LG ADX" = "#32146F",
    "LG D. ADX" = "#4B1FA8",
    "empirical YRI" = "#EEC4DC", "empirical ASW" = "#E44B8D",
    "empirical CEU" = "#BB437E"
    ),
  fill.labels = c(
    "TC AFR" = "T.C. AFR",
    "TC ADX" = "T.C. ADX",
    "TC EUR" = "T.C. EUR",
    "TC D. AFR" = "T.C.D. AFR",
    "TC D. ADX" = "T.C.D. ADX",
    "TC D. EUR" = "T.C.D. EUR",
    "LG ADX" = "L.G. ADX",
    "LG D. ADX" = "L.G.D. ADX",
    "empirical YRI" = "YRI",
    "empirical ASW" = "ASW",
    "empirical CEU" = "CEU"
    ),
  empirical.colors = c(
    YRI = "#EEC4DC", ASW = "#E44B8D", CEU = "#BB437E"
    ),
  series.labels = SOURCE.LABELS
  )


# internal functions ----


# return the canonical levels represented by a filtered plot view
order.active.levels <- function(values, canonical.levels) {
  active.levels <- canonical.levels[
    canonical.levels %in% as.character(values)
    ]
  return(active.levels)
  }


# summarize complete simulation replicates with direct percentile intervals
summarize.bootstrap.interval <- function(
    data, grouping.columns, value.column
  ) {
  if (any(!is.finite(data[[value.column]]))) {
    stop("Simulation values must be finite")
    }
  summary <- data %>%
    group_by(across(all_of(setdiff(grouping.columns, "rep")))) %>%
    summarise(
      replicate.count = n_distinct(rep),
      mean = mean(.data[[value.column]]),
      lower = quantile(.data[[value.column]], 0.025, names = FALSE),
      upper = quantile(.data[[value.column]], 0.975, names = FALSE),
      .groups = "drop"
      )
  if (any(summary$replicate.count != 50L)) {
    stop("Simulation summaries require exactly 50 complete replicates")
    }
  return(summary)
  }


# retain source-specific populations before any replicate summaries
apply.diversity.source.contract <- function(data) {
  retained <- data %>%
    filter(
      (data.type %in% c(
        "Simulation_2T12Consistent", "Simulation_2T12Consistent_simDown"
        ) & pop %in% c("AFR", "ADX", "EUR")) |
        (data.type %in% c(
          "Simulation_largeGrowth", "Simulation_largeGrowth_simDown"
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
  if (!"span" %in% names(data)) {
    stop("Diversity data are missing the required span column")
    }
  data <- data %>%
    mutate(
      rep = as.numeric(rep),
      chrom = as.character(chrom),
      pop = as.character(pop),
      stat = as.character(stat),
      value = as.numeric(value),
      span = as.numeric(span),
      data.type = data.type.input,
      mask = mask.input
      ) %>%
    apply.diversity.source.contract() %>%
    add.population.roles()
  if (any(!is.finite(data$span) | data$span <= 0)) {
    stop("Diversity spans must be positive and finite")
    }

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
  missing <- !file.exists(paths)
  if (any(missing)) {
    warning(
      paste0(
        data.type.input, " ", file.family,
        " is unavailable for chromosomes: ",
        paste(chromosomes[missing], collapse = ", ")
        ),
      call. = FALSE
      )
    }
  paths <- paths[!missing]
  chromosomes <- chromosomes[!missing]
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
  replicate.summary <- data %>%
    filter(stat %in% c("pi", "theta")) %>%
    group_by(data.type, rep, pop, role, stat, chrom, mask, span) %>%
    summarise(
      value = mean(value),
      .groups = "drop"
      )
  summary <- summarize.bootstrap.interval(
    replicate.summary,
    c(
      "data.type", "rep", "pop", "role", "stat", "chrom", "mask",
      "span"
      ),
    "value"
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
      data.type, pop, role, stat, chrom, mask, span,
      estimate = mean, lower, upper, replicate.count
      )
  empirical.points <- empirical.chromosome %>%
    filter(chrom %in% chromosomes, stat %in% c("pi", "theta")) %>%
    transmute(
      data.type, pop, role, stat, chrom, mask, span,
      estimate = value, lower = NA_real_, upper = NA_real_, replicate.count = 1L
      )
  points <- bind_rows(simulation.points, empirical.points) %>%
    mutate(
      chrom = factor(chrom, levels = chromosomes),
      pop = factor(pop, levels = POPULATION.LEVELS),
      data.type = factor(
        data.type, levels = SOURCE.LEVELS
        ),
      mask = factor(mask, levels = c("Intergenic", "Full callable")),
      stat = factor(stat, levels = c("pi", "theta")),
      fill.key = factor(case_when(
        data.type == "Simulation_2T12Consistent" ~ paste("TC", pop),
        data.type == "Simulation_largeGrowth" ~ "LG ADX",
        data.type == "Simulation_2T12Consistent_simDown" ~
          paste("TC D.", pop),
        data.type == "Simulation_largeGrowth_simDown" ~ "LG D. ADX",
        data.type == "Empirical" ~ paste("empirical", pop)
        ), levels = names(PLOT.STYLES$fill.colors))
      )
  genome.lines <- empirical.genome %>%
    filter(chrom == "all", stat %in% c("pi", "theta")) %>%
    transmute(pop, role, stat, mask, span, estimate = value) %>%
    mutate(
      pop = factor(pop, levels = POPULATION.LEVELS),
      mask = factor(mask, levels = c("Intergenic", "Full callable")),
      stat = factor(stat, levels = c("pi", "theta"))
      )

  return(list(points = points, genome.lines = genome.lines))
  }


# retain normalized values or scale all diversity estimates by their span
rescale.diversity.plot.values <- function(data, value.view) {
  if (!value.view %in% c("normalized", "unscaled")) {
    stop("Unsupported diversity value view: ", value.view)
    }
  if (!"span" %in% names(data)) {
    stop("Diversity plot data are missing span")
    }
  if (any(!is.finite(data$span) | data$span <= 0)) {
    stop("Diversity plot spans must be positive and finite")
    }
  if (value.view == "normalized") {
    return(data)
    }
  scaled <- data %>% mutate(estimate = estimate * span)
  if (all(c("lower", "upper") %in% names(scaled))) {
    scaled <- scaled %>%
      mutate(lower = lower * span, upper = upper * span)
    }
  return(scaled)
  }


# retain one configured diversity view and its eligible empirical references
filter.diversity.plot.view <- function(
    points, genome.lines, data.types, tag
  ) {
  if (!tag %in% names(PLOT.CONFIGS)) {
    stop("Unsupported diversity plot tag: ", tag)
    }
  if (!identical(data.types, PLOT.CONFIGS[[tag]])) {
    stop("Diversity data types do not match the configured tag")
    }
  view.points <- points %>%
    filter(as.character(data.type) %in% data.types) %>%
    filter(tag != "onlyADX" | pop == "ADX") %>%
    filter(
      tag != "all.datatypes.adx.asw" |
        (data.type != "Empirical" & pop == "ADX") |
        (data.type == "Empirical" & pop == "ASW")
      )
  active.sources <- order.active.levels(
    view.points$data.type, SOURCE.LEVELS
    )
  active.populations <- order.active.levels(
    view.points$pop, POPULATION.LEVELS
    )
  view.points <- view.points %>%
    mutate(
      data.type = factor(as.character(data.type), levels = active.sources),
      pop = factor(as.character(pop), levels = active.populations),
      fill.key = factor(
        as.character(fill.key),
        levels = order.active.levels(
          fill.key, names(PLOT.STYLES$fill.colors)
          )
        )
      ) %>%
    arrange(data.type, fill.key) %>%
    droplevels()
  view.lines <- if ("Empirical" %in% active.sources) genome.lines else {
    genome.lines[0, , drop = FALSE]
    }
  view.lines <- view.lines %>%
    filter(as.character(pop) %in% active.populations) %>%
    mutate(pop = factor(as.character(pop), levels = active.populations))
  return(list(points = view.points, genome.lines = view.lines))
  }


# construct one configured selected-chromosome diversity plot
make.diversity.plot <- function(
    points, genome.lines, styles, data.types, tag, value.view = "normalized"
  ) {
  view <- filter.diversity.plot.view(
    points, genome.lines, data.types, tag
    )
  points <- rescale.diversity.plot.values(view$points, value.view)
  genome.lines <- rescale.diversity.plot.values(
    view$genome.lines, value.view
    )
  points <- filter(points, mask == "Intergenic")
  genome.lines <- filter(genome.lines, mask == "Intergenic")
  dodge <- position_dodge(width = CATEGORICAL.BAR.DODGE)
  fill.keys <- levels(points$fill.key)
  plot <- ggplot(
    points,
    aes(
      x = chrom, y = estimate, fill = fill.key,
      group = interaction(pop, data.type)
      )
    ) +
    ggfx::with_outer_glow(
      geom_hline(
        data = genome.lines,
        aes(yintercept = estimate, color = pop),
        linetype = "longdash", linewidth = CATEGORICAL.BAR.LINEWIDTH
        ),
      colour = "black", sigma = 0, expand = 3
      ) +
    geom_col(
      position = dodge, width = CATEGORICAL.BAR.WIDTH,
      color = "black", linewidth = CATEGORICAL.BAR.LINEWIDTH
      ) +
    geom_errorbar(
      data = points,
      aes(ymin = lower, ymax = upper),
      position = dodge, width = 0.15,
      linewidth = CATEGORICAL.BAR.LINEWIDTH,
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
      drop = TRUE
      ) +
    scale_y_continuous(
      labels = scales::label_number(accuracy = 0.00001)
      ) +
    labs(
      x = "Chromosome", y = NULL,
      title = if (value.view == "normalized") {
        "Genetic Diversity Across Selected Chromosomes"
        } else {
        "Span-Unscaled Genetic Diversity Across Selected Chromosomes"
        },
      color = NULL, fill = NULL, shape = NULL
      ) +
    guides(
      color = "none",
      fill = guide_legend(order = 1, nrow = 2, byrow = TRUE)
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


# read all available chromosome-level simulation diversity estimates
sim.tc.diversity <- read.diversity.chromosomes(
  SIM.TC.DATA.DIR, "pi_theta_stats.chr{chrom}.parquet",
  CHROMOSOMES, "Simulation_2T12Consistent"
  )
simDown.tc.intergenic.diversity <- read.diversity.chromosomes(
  SIMDOWN.TC.DATA.DIR,
  "pi_theta_stats_intergenic.chr{chrom}.parquet",
  CHROMOSOMES, "Simulation_2T12Consistent_simDown", "Intergenic"
  )
simDown.tc.full.callable.diversity <- read.diversity.chromosomes(
  SIMDOWN.TC.DATA.DIR,
  "pi_theta_stats_full_callable_chrom.chr{chrom}.parquet",
  CHROMOSOMES, "Simulation_2T12Consistent_simDown", "Full callable"
  )
sim.lg.diversity <- read.diversity.chromosomes(
  SIM.LG.DATA.DIR, "pi_theta_stats.chr{chrom}.parquet",
  CHROMOSOMES, "Simulation_largeGrowth"
  )
simDown.lg.intergenic.diversity <- read.diversity.chromosomes(
  SIMDOWN.LG.DATA.DIR,
  "pi_theta_stats_intergenic.chr{chrom}.parquet",
  CHROMOSOMES, "Simulation_largeGrowth_simDown", "Intergenic"
  )
simDown.lg.full.callable.diversity <- read.diversity.chromosomes(
  SIMDOWN.LG.DATA.DIR,
  "pi_theta_stats_full_callable_chrom.chr{chrom}.parquet",
  CHROMOSOMES, "Simulation_largeGrowth_simDown", "Full callable"
  )

# read all available chromosome-level empirical diversity estimates
emp.intergenic.chromosome <- read.diversity.chromosomes(
  EMPIRICAL.DATA.DIR,
  "pi_theta_stats_intergenic.chr{chrom}.parquet",
  CHROMOSOMES, "Empirical", "Intergenic"
  )
emp.full.callable.chromosome <- read.diversity.chromosomes(
  EMPIRICAL.DATA.DIR,
  "pi_theta_stats_full_callable_chrom.chr{chrom}.parquet",
  CHROMOSOMES, "Empirical", "Full callable"
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

# summarize sources and construct all configured diversity views
simulation.diversity.summary <- bind_rows(
  sim.tc.diversity,
  simDown.tc.intergenic.diversity,
  simDown.tc.full.callable.diversity,
  sim.lg.diversity,
  simDown.lg.intergenic.diversity,
  simDown.lg.full.callable.diversity
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
diversity.bootstrap.plots <- imap(list(
  tc.tcd.1kg = PLOT.CONFIGS$tc.tcd.1kg,
  all.datatypes.adx.asw = PLOT.CONFIGS$all.datatypes.adx.asw
  ), function(data.types, tag) {
  return(make.diversity.plot(
    diversity.plot.data$points,
    diversity.plot.data$genome.lines,
    PLOT.STYLES, data.types,
    tag
    ))
  })
diversity.bootstrap.unscaled.plots <- imap(list(
  tc.tcd.1kg = PLOT.CONFIGS$tc.tcd.1kg,
  all.datatypes.adx.asw = PLOT.CONFIGS$all.datatypes.adx.asw
  ), function(data.types, tag) {
  return(make.diversity.plot(
    diversity.plot.data$points,
    diversity.plot.data$genome.lines,
    PLOT.STYLES, data.types, tag,
    value.view = "unscaled"
    ))
  })

# persist every plot before printing figures at the end of the script
dir.create(OUTPUT.DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(diversity.bootstrap.plots$tc.tcd.1kg, file.path(
  OUTPUT.DIR, "diversity.bootstrap.tc.tcd.1kg.rds"
  ))
saveRDS(diversity.bootstrap.plots$all.datatypes.adx.asw, file.path(
  OUTPUT.DIR, "diversity.bootstrap.all.datatypes.adx.asw.rds"
  ))
saveRDS(diversity.bootstrap.unscaled.plots$tc.tcd.1kg, file.path(
  OUTPUT.DIR, "diversity.bootstrap.tc.tcd.1kg.unscaled.rds"
  ))
saveRDS(diversity.bootstrap.unscaled.plots$all.datatypes.adx.asw, file.path(
  OUTPUT.DIR, "diversity.bootstrap.all.datatypes.adx.asw.unscaled.rds"
  ))

print(diversity.bootstrap.plots$tc.tcd.1kg)
print(diversity.bootstrap.plots$all.datatypes.adx.asw)
print(diversity.bootstrap.unscaled.plots$tc.tcd.1kg)
print(diversity.bootstrap.unscaled.plots$all.datatypes.adx.asw)
