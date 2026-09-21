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


SIM.TC.DATA.DIR <- "~/scratch/OOA_NAAdmixture_2T12Consistent/stats"
SIMDOWN.TC.DATA.DIR <- "~/scratch/OOA_NAAdmixture_2T12ConsistentOnekgDownsample/stats"
SIM.LG.DATA.DIR <- "~/scratch/OOA_NAAdmixture_largeGrowth/stats"
SIMDOWN.LG.DATA.DIR <- "~/scratch/OOA_NAAdmixture_largeGrowthOnekgDownsample/stats"
EMPIRICAL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_1kG/stats"
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
CHROMOSOMES <- as.character(1:22)
SELECTED.CHROMOSOMES <- c("1")
SOURCE.LEVELS <- c(
  "Simulation_2T12Consistent", "Simulation_2T12Consistent_simDown",
  "Simulation_largeGrowth", "Simulation_largeGrowth_simDown",
  "Empirical"
  )
PLOT.CONFIGS <- list(
  TC.1kG = SOURCE.LEVELS[c(1, 5)],
  TC.TCD = SOURCE.LEVELS[c(1, 2)],
  TCD.1kG = SOURCE.LEVELS[c(2, 5)],
  onlyADX = SOURCE.LEVELS[1:4]
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
    Simulation_2T12Consistent = "#9A83CE",
    Simulation_2T12Consistent_simDown = "#6F55B5",
    Simulation_largeGrowth = "#32146F",
    Simulation_largeGrowth_simDown = "#4B1FA8"
    ),
  series.labels = c(
    Simulation_2T12Consistent = "T.C.",
    Simulation_2T12Consistent_simDown = "T.C.D.",
    Simulation_largeGrowth = "L.G.",
    Simulation_largeGrowth_simDown = "L.G.D.",
    Empirical = "Emp."
    )
  )


# internal functions ----


# summarize complete pooled replicate curves with percentile intervals
summarize.simulation.interval <- function(
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


# retain source-specific populations before pooling sufficient statistics
apply.ld.source.contract <- function(data) {
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
    glue::glue("ld_decay.chr{chromosomes}.parquet")
    )
  missing <- !file.exists(paths)
  if (any(missing)) {
    warning(
      paste0(
        data.type.input,
        " LD files are unavailable for chromosomes: ",
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
summarize.ld.curves <- function(data, chromosomes) {
  chromosomes <- as.character(chromosomes)
  if (!length(chromosomes) || any(!chromosomes %in% CHROMOSOMES)) {
    stop("LD summarization requires autosomal chromosomes")
    }
  scoped <- data %>% filter(as.character(chrom) %in% chromosomes)
  if (!nrow(scoped)) {
    stop("LD data do not contain any requested chromosomes")
    }
  if (any(!as.character(scoped$chrom) %in% chromosomes)) {
    stop("LD summary contains data outside the requested chromosomes")
    }
  simulation <- scoped %>%
    filter(data.type != "Empirical") %>%
    summarize.simulation.interval(
      c("data.type", "rep", "pop", "role", "chrom", "distance_bin_bp"),
      "mean.r2"
      )
  empirical <- scoped %>%
    filter(data.type == "Empirical") %>%
    transmute(
      data.type, pop, role, chrom, distance_bin_bp,
      mean = mean.r2, lower = NA_real_, upper = NA_real_, replicate.count = 1L,
      chromosome.count
      ) %>%
    distinct()
  empirical.genome <- data %>%
    filter(data.type == "Empirical", as.character(chrom) == "all") %>%
    transmute(
      data.type, pop, role, chrom, distance_bin_bp,
      mean = mean.r2, lower = NA_real_, upper = NA_real_, replicate.count = 1L,
      chromosome.count
      ) %>%
    distinct()
  summary <- bind_rows(simulation, empirical, empirical.genome) %>%
    mutate(
      role = factor(role, levels = c("AFR", "ADX", "EUR")),
      chrom = factor(as.character(chrom), levels = c(chromosomes, "all")),
      data.type = factor(
        data.type,
        levels = SOURCE.LEVELS
        )
      ) %>%
    filter(!is.na(mean))

  return(summary)
  }


# add shared scales, labels, guides, and theme to one LD plot
style.ld.plot <- function(plot, title, subtitle, styles, tag) {
  source.view <- tag == "onlyADX"
  plot <- plot +
    scale_color_manual(
      values = if (source.view) {
        styles$source.colors
        } else {
        styles$population.colors
        },
      breaks = if (source.view) {
        PLOT.CONFIGS[[tag]]
        } else {
        names(styles$population.colors)
        },
      labels = if (source.view) {
        styles$series.labels[PLOT.CONFIGS[[tag]]]
        } else {
        waiver()
        }
      ) +
    scale_fill_manual(values = if (source.view) {
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
      color = guide_legend(order = 1, nrow = 1, byrow = TRUE),
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
        ymin = lower,
        ymax = upper,
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


# retain one configured source view
filter.plot.view <- function(data, data.types, tag) {
  if (!tag %in% names(PLOT.CONFIGS)) {
    stop("Unsupported LD plot tag: ", tag)
    }
  if (!identical(data.types, PLOT.CONFIGS[[tag]])) {
    stop("LD data types do not match the configured tag")
    }
  filtered <- data %>%
    filter(as.character(data.type) %in% data.types) %>%
    filter(tag != "onlyADX" | pop == "ADX") %>%
    mutate(data.type = factor(as.character(data.type), levels = data.types)) %>%
    arrange(data.type) %>%
    droplevels()
  return(filtered)
  }


# construct one scoped LD view for selected chromosomes
make.ld.plot <- function(
    data, chromosomes, styles, data.types, tag, show.all = FALSE
  ) {
  chromosomes <- as.character(chromosomes)
  if (!length(chromosomes) || any(!chromosomes %in% CHROMOSOMES)) {
    stop("LD plotting requires autosomal chromosomes")
    }
  source.view <- tag == "onlyADX"
  plot.data <- filter.plot.view(data, data.types, tag) %>%
    filter(
      as.character(chrom) %in% chromosomes |
        (show.all & data.type == "Empirical" & chrom == "all")
      ) %>%
    mutate(
      chrom = factor(
        as.character(chrom),
        levels = c(chromosomes, if (show.all) "all")
        ),
      plot.key = if (source.view) {
        factor(as.character(data.type), levels = data.types)
        } else {
        factor(as.character(pop), levels = c(
          "AFR", "ADX", "EUR", "YRI", "ASW", "CEU"
          ))
        }
      )
  if (!nrow(plot.data)) {
    stop("LD data do not contain any selected chromosomes")
    }
  plot <- ggplot(
    plot.data,
    aes(
      x = distance_bin_bp, y = mean,
      color = plot.key, fill = plot.key
      )
    )
  plot <- add.ld.geometries(plot, plot.data) +
    facet_grid(
      chrom ~ data.type, drop = TRUE,
      scales = "free_y",
      labeller = labeller(data.type = styles$series.labels)
      )
  plot <- style.ld.plot(
    plot, paste("LD Decay:", tag),
    paste0(
      "Rogers–Huff r²; chrom.", paste(chromosomes, collapse = ", ")
      ),
    styles, tag
    )

  return(plot)
  }


# build one chromosome-1 bootstrap LD view over the requested distance range
make.bootstrap.ld.plot <- function(data, data.types, view) {
  source.view <- grepl("all.datatypes.adx.asw", view)
  plotted <- data %>%
    filter(
      as.character(chrom) == "1",
      as.character(data.type) %in% data.types,
      between(distance_bin_bp, 5000, 250000)
      ) %>%
    filter(
      !source.view |
        (data.type != "Empirical" & pop == "ADX") |
        (data.type == "Empirical" & pop == "ASW")
      ) %>%
    mutate(
      plot.key = case_when(
        source.view & data.type == "Empirical" ~ "ASW",
        source.view ~ as.character(data.type),
        TRUE ~ as.character(pop)
        )
      )
  plot <- ggplot(plotted, aes(distance_bin_bp, mean, color = plot.key,
    fill = plot.key, group = interaction(data.type, pop))) +
    scale_x_continuous(limits = c(5000, 250000)) +
    labs(x = "Distance between SNPs (bp)", y = expression("Mean " * r^2),
      color = NULL, fill = NULL) +
    theme_bw(base_size = PLOT.BASE.SIZE) +
    theme(legend.position = "top", panel.grid.minor = element_blank())
  if (!grepl("all.lines$", view)) {
    plot <- plot + geom_ribbon(
      data = filter(plotted, data.type != "Empirical"),
      aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA
      )
    }
  source.colors <- c(
    PLOT.STYLES$source.colors,
    ASW = PLOT.STYLES$population.colors[["ASW"]]
    )
  source.labels <- c(PLOT.STYLES$series.labels, ASW = "ASW")
  plot <- plot + geom_line(linewidth = 1) +
    scale_color_manual(
      values = if (source.view) {
        source.colors
        } else {
        PLOT.STYLES$population.colors
        },
      breaks = if (source.view) names(source.colors) else NULL,
      labels = if (source.view) source.labels[names(source.colors)] else NULL
      ) +
    scale_fill_manual(
      values = if (source.view) {
        source.colors
        } else {
        PLOT.STYLES$population.colors
        },
      breaks = if (source.view) names(source.colors) else NULL,
      labels = if (source.view) source.labels[names(source.colors)] else NULL
      )
  if (view == "role.interval") plot <- plot + facet_wrap(~role)
  if (grepl("datatype.interval$", view)) {
    plot <- plot + facet_wrap(
      ~data.type,
      labeller = labeller(data.type = PLOT.STYLES$series.labels)
      )
    }
  return(plot)
  }


# analysis ----


# read all available chromosomes for all four simulation sources
sim.tc.ld.chromosomes <- read.ld.chromosomes(
  SIM.TC.DATA.DIR, CHROMOSOMES, "Simulation_2T12Consistent"
  )
simDown.tc.ld.chromosomes <- read.ld.chromosomes(
  SIMDOWN.TC.DATA.DIR, CHROMOSOMES,
  "Simulation_2T12Consistent_simDown"
  )
sim.lg.ld.chromosomes <- read.ld.chromosomes(
  SIM.LG.DATA.DIR, CHROMOSOMES, "Simulation_largeGrowth"
  )
simDown.lg.ld.chromosomes <- read.ld.chromosomes(
  SIMDOWN.LG.DATA.DIR, CHROMOSOMES,
  "Simulation_largeGrowth_simDown"
  )
simulation.ld.chromosomes <- bind_rows(
  sim.tc.ld.chromosomes, simDown.tc.ld.chromosomes,
  sim.lg.ld.chromosomes, simDown.lg.ld.chromosomes
  )

# pool every simulation chromosome curve from producer sufficient statistics
simulation.ld.selected <- simulation.ld.chromosomes %>%
  pool.ld.curves(include.chromosome = TRUE)

# read and pool all available empirical chromosomes
empirical.ld.selected <- read.ld.chromosomes(
  EMPIRICAL.DATA.DIR, CHROMOSOMES, "Empirical"
  ) %>%
  pool.ld.curves(include.chromosome = TRUE)
empirical.ld.genome <- read.empirical.ld.genome(EMPIRICAL.DATA.DIR) %>%
  pool.ld.curves(include.chromosome = TRUE)

# summarize curves and construct all configured source views
ld.summary <- bind_rows(
  simulation.ld.selected, empirical.ld.selected, empirical.ld.genome
  ) %>%
  summarize.ld.curves(SELECTED.CHROMOSOMES)
ld.plots <- imap(PLOT.CONFIGS, function(data.types, tag) {
  return(make.ld.plot(
    ld.summary, SELECTED.CHROMOSOMES, PLOT.STYLES, data.types, tag,
    show.all = FALSE
    ))
  })

# persist every plot before printing figures at the end of the script
# Legacy plot writes are retained as inactive reference code.
if (FALSE) {
  dir.create(OUTPUT.DIR, recursive = TRUE, showWarnings = FALSE)
  iwalk(ld.plots, function(plot, tag) {
    saveRDS(plot, file.path(
      OUTPUT.DIR,
      str_replace("ld.decay.{tag}.rds", fixed("{tag}"), tag)
      ))
    })
  print(ld.plots$TC.1kG)
  print(ld.plots$TC.TCD)
  print(ld.plots$TCD.1kG)
  print(ld.plots$onlyADX)
  }

# save TCD/1kG and all-datatype ADX/ASW chromosome-1 bootstrap LD views
bootstrap.tcd.1kg.ld.all.lines <- make.bootstrap.ld.plot(
  ld.summary, c("Simulation_2T12Consistent_simDown", "Empirical"),
  "all.lines"
  )
bootstrap.tcd.1kg.ld.role.interval <- make.bootstrap.ld.plot(
  ld.summary, c("Simulation_2T12Consistent_simDown", "Empirical"),
  "role.interval"
  )
bootstrap.tcd.1kg.ld.datatype.interval <- make.bootstrap.ld.plot(
  ld.summary, c("Simulation_2T12Consistent_simDown", "Empirical"),
  "datatype.interval"
  )
bootstrap.all.datatypes.adx.asw.ld.all.lines <- make.bootstrap.ld.plot(
  ld.summary, SOURCE.LEVELS, "all.datatypes.adx.asw.all.lines"
  )
bootstrap.all.datatypes.adx.asw.ld.datatype.interval <- make.bootstrap.ld.plot(
  ld.summary, SOURCE.LEVELS, "all.datatypes.adx.asw.datatype.interval"
  )
dir.create(OUTPUT.DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(bootstrap.tcd.1kg.ld.all.lines, file.path(
  OUTPUT.DIR, "ld.bootstrap.tcd.1kg.all.lines.rds"
  ))
saveRDS(bootstrap.tcd.1kg.ld.role.interval, file.path(
  OUTPUT.DIR, "ld.bootstrap.tcd.1kg.role.interval.rds"
  ))
saveRDS(bootstrap.tcd.1kg.ld.datatype.interval, file.path(
  OUTPUT.DIR, "ld.bootstrap.tcd.1kg.datatype.interval.rds"
  ))
saveRDS(bootstrap.all.datatypes.adx.asw.ld.all.lines, file.path(
  OUTPUT.DIR, "ld.bootstrap.all.datatypes.adx.asw.all.lines.rds"
  ))
saveRDS(bootstrap.all.datatypes.adx.asw.ld.datatype.interval, file.path(
  OUTPUT.DIR, "ld.bootstrap.all.datatypes.adx.asw.datatype.interval.rds"
  ))

print(bootstrap.tcd.1kg.ld.all.lines)
print(bootstrap.tcd.1kg.ld.role.interval)
print(bootstrap.tcd.1kg.ld.datatype.interval)
print(bootstrap.all.datatypes.adx.asw.ld.all.lines)
print(bootstrap.all.datatypes.adx.asw.ld.datatype.interval)
