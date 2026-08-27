# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# sfs.R
# ______________________________________________________________________________


# pattern: Mixed (unavoidable)
# Reason: self-contained helpers accompany sequential I/O orchestration.

# set up ----
library(tidyverse)
library(glue)
library(nanoparquet)
library(scales)


SIM.SMALL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_small/stats"
SIM.LARGE.DATA.DIR <- "~/scratch/OOA_NAAdmixture_large/stats"
EMPIRICAL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_1kG/stats"
CHROMOSOMES <- as.character(1:22)
SELECTED.CHROMOSOMES <- c("1", "5", "10", "14", "18", "22")
DISPLAY.BIN.MAX <- 10
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


# require a producer table to contain the expected schema
check.required.columns <- function(data, required.columns, table.name) {
  missing.columns <- setdiff(required.columns, names(data))
  if (length(missing.columns) > 0L) {
    stop(glue(
      "{table.name} is missing required columns: ",
      "{glue_collapse(missing.columns, sep = ', ')}"
    ))
  }

  return(invisible(NULL))
}


# map source population labels to shared roles
add.population.roles <- function(data) {
  data <- data %>%
    mutate(role = case_when(
      pop %in% c("AFR", "YRI") ~ "AFR",
      pop %in% c("ADX", "ASW") ~ "ADX",
      pop %in% c("EUR", "CEU") ~ "EUR",
      TRUE ~ NA_character_
    ))
  if (any(is.na(data$role))) {
    stop("SFS data contain an unsupported population label")
  }

  return(data)
}


# validate numeric SFS counts and the analysis-owned allele-count field
validate.sfs.values <- function(data, table.name) {
  invalid <- !is.finite(data$allele.count) |
    data$allele.count < 0 |
    data$allele.count != floor(data$allele.count) |
    !is.finite(data$count) |
    data$count < 0
  if (any(invalid)) {
    stop(glue("{table.name} contains invalid allele or site counts"))
  }

  return(data)
}


# normalize simulation fields without renaming producer-owned columns
normalize.simulation.sfs <- function(data, data.type.input) {
  check.required.columns(
    data,
    c("rep", "chrom", "pop", "derived_allele_count", "count"),
    "Simulation SFS"
  )
  data <- data %>%
    mutate(
      rep = as.numeric(rep),
      chrom = as.character(chrom),
      pop = as.character(pop),
      derived_allele_count = as.numeric(derived_allele_count),
      count = as.numeric(count),
      allele.count = derived_allele_count,
      data.type = data.type.input
    ) %>%
    add.population.roles() %>%
    validate.sfs.values("Simulation SFS")

  return(data)
}


# normalize folded empirical fields without renaming producer-owned columns
normalize.empirical.sfs <- function(data) {
  check.required.columns(
    data,
    c("rep", "chrom", "pop", "minor_allele_count", "count"),
    "Empirical SFS"
  )
  data <- data %>%
    mutate(
      rep = as.numeric(rep),
      chrom = as.character(chrom),
      pop = as.character(pop),
      minor_allele_count = as.numeric(minor_allele_count),
      count = as.numeric(count),
      allele.count = minor_allele_count,
      data.type = "Empirical"
    ) %>%
    add.population.roles() %>%
    validate.sfs.values("Empirical SFS")

  return(data)
}


# infer and validate complete unfolded bins from zero through N
validate.simulation.bins <- function(data) {
  grouping.columns <- c(
    "data.type", "rep", "pop", "role", "chrom"
  )
  groups <- data %>%
    group_by(across(all_of(grouping.columns))) %>%
    group_split(.keep = TRUE)
  validated <- map_dfr(groups, function(group) {
    bins <- sort(group$allele.count)
    sample.size <- max(bins)
    expected <- seq.int(0, sample.size)
    if (
      sample.size <= 0 |
        anyDuplicated(bins) > 0L |
        !identical(as.numeric(bins), as.numeric(expected))
    ) {
      stop(
        "Simulation SFS groups require exactly one bin from 0 through N"
      )
    }
    group$sample.size <- sample.size
    return(group)
  })

  return(validated)
}


# derive role-specific empirical haplotype targets from retained samples
derive.empirical.target.sizes <- function(variant.qc) {
  qc.columns <- c(
    AFR = "retained_YRI_samples",
    ADX = "retained_ASW_samples",
    EUR = "retained_CEU_samples"
  )
  check.required.columns(
    variant.qc, unname(qc.columns), "Empirical variant QC"
  )
  target.sizes <- imap_dfr(qc.columns, function(column, role.input) {
    values <- as.numeric(variant.qc[[column]])
    values <- unique(values[is.finite(values)])
    if (
      length(values) != 1L ||
        values <= 0 ||
        values != floor(values)
    ) {
      stop(glue(
        "Empirical variant QC requires one positive sample count for ",
        "{role.input}"
      ))
    }
    return(tibble(
      role = role.input,
      target.size = 2 * values
    ))
  })

  return(target.sizes)
}


# reconstruct simulation genome spectra by summing chromosome counts
aggregate.simulation.genome <- function(data) {
  genome <- data %>%
    group_by(data.type, rep, pop, role, allele.count) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    mutate(chrom = "all") %>%
    select(
      data.type, rep, chrom, pop, role, allele.count, count
    ) %>%
    validate.simulation.bins()

  return(genome)
}


# validate already-folded empirical bins and attach haplotype targets
prepare.empirical.comparable <- function(data, target.sizes) {
  prepared <- data %>%
    left_join(target.sizes, by = "role")
  if (
    any(is.na(prepared$target.size)) |
      any(prepared$target.size <= 0) |
      any(prepared$allele.count > floor(prepared$target.size / 2))
  ) {
    stop("Empirical SFS contains bins outside its folded sample size")
  }
  has.duplicates <- prepared %>%
    group_by(data.type, rep, pop, role, chrom) %>%
    summarise(
      duplicate.bin = anyDuplicated(allele.count) > 0L,
      .groups = "drop"
    ) %>%
    pull(duplicate.bin) %>%
    any()
  if (has.duplicates) {
    stop("Empirical SFS groups require unique folded bins")
  }
  prepared <- prepared %>%
    mutate(sample.size = target.size) %>%
    select(-target.size)

  return(prepared)
}


# calculate exact expected projection for one unfolded spectrum
project.one.spectrum <- function(group) {
  source.size <- unique(group$sample.size)
  target.size <- unique(group$target.size)
  if (
    length(source.size) != 1L ||
      length(target.size) != 1L ||
      target.size <= 0 ||
      target.size > source.size
  ) {
    stop("Projection target must be positive and no larger than source N")
  }
  projected.bins <- seq.int(0, target.size)
  projected.counts <- map_dbl(projected.bins, function(target.bin) {
    probabilities <- dhyper(
      target.bin,
      group$allele.count,
      source.size - group$allele.count,
      target.size
    )
    return(sum(group$count * probabilities))
  })
  if (!isTRUE(all.equal(
    sum(projected.counts), sum(group$count), tolerance = 1e-10
  ))) {
    stop("Hypergeometric projection did not conserve spectrum mass")
  }
  projected <- tibble(
    allele.count = projected.bins,
    count = projected.counts,
    sample.size = target.size
  )

  return(projected)
}


# project complete simulation spectra to empirical role targets
project.simulation.sfs <- function(data, target.sizes) {
  projected <- data %>%
    validate.simulation.bins() %>%
    left_join(target.sizes, by = "role")
  if (any(is.na(projected$target.size))) {
    stop("Every simulation population role requires a projection target")
  }
  projected <- projected %>%
    group_by(data.type, rep, pop, role, chrom) %>%
    group_modify(~project.one.spectrum(.x)) %>%
    ungroup()

  return(projected)
}


# fold unfolded spectra with the midpoint counted exactly once
fold.simulation.sfs <- function(data) {
  folded <- data %>%
    validate.simulation.bins() %>%
    mutate(
      allele.count = pmin(allele.count, sample.size - allele.count)
    ) %>%
    group_by(
      data.type, rep, pop, role, chrom, sample.size, allele.count
    ) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    arrange(data.type, rep, pop, chrom, allele.count)

  return(folded)
}


# retain segregating bins and remove zero or fixed-bin mass
filter.segregating.bins <- function(data) {
  filtered <- data %>%
    filter(allele.count > 0, allele.count < sample.size)

  return(filtered)
}


# normalize segregating-site mass before display-bin restriction
calculate.sfs.proportions <- function(data) {
  proportions <- data %>%
    filter.segregating.bins() %>%
    group_by(data.type, rep, pop, role, chrom) %>%
    mutate(
      segregating.mass = sum(count),
      proportion = count / segregating.mass
    ) %>%
    ungroup()
  if (
    any(!is.finite(proportions$segregating.mass)) |
      any(proportions$segregating.mass <= 0) |
      any(!is.finite(proportions$proportion))
  ) {
    stop("SFS proportions require positive segregating-site mass")
  }

  return(proportions)
}


# restrict displays only after full-spectrum calculations
filter.display.bins <- function(data, maximum.bin = DISPLAY.BIN.MAX) {
  filtered <- data %>%
    filter(allele.count >= 1, allele.count <= maximum.bin)

  return(filtered)
}


# summarize simulation replicates and retain fixed empirical curves
summarize.sfs.curves <- function(data, value.type) {
  if (!value.type %in% c("count", "proportion")) {
    stop("SFS summary value type must be count or proportion")
  }
  if (!"value" %in% names(data)) {
    value.column <- if (value.type == "count") "count" else "proportion"
    data$value <- data[[value.column]]
  }
  if (!"sample.size" %in% names(data)) {
    data$sample.size <- NA_real_
  }
  grouping.columns <- c(
    "data.type", "pop", "role", "chrom", "sample.size",
    "allele.count"
  )
  simulation <- data %>%
    filter(data.type != "Empirical") %>%
    group_by(across(all_of(grouping.columns))) %>%
    summarise(
      mean = mean(value),
      sd = sd(value),
      replicate.count = n_distinct(rep),
      .groups = "drop"
    )
  empirical <- data %>%
    filter(data.type == "Empirical") %>%
    group_by(across(all_of(grouping.columns))) %>%
    summarise(
      mean = mean(value),
      sd = NA_real_,
      replicate.count = n_distinct(rep),
      .groups = "drop"
    )
  summary <- bind_rows(simulation, empirical) %>%
    mutate(
      ribbon.lower = if_else(
        data.type == "Empirical", NA_real_, pmax(0, mean - 2 * sd)
      ),
      ribbon.upper = case_when(
        data.type == "Empirical" ~ NA_real_,
        value.type == "proportion" ~ pmin(1, mean + 2 * sd),
        TRUE ~ mean + 2 * sd
      ),
      role = factor(role, levels = c("AFR", "ADX", "EUR")),
      data.type = factor(
        data.type,
        levels = c(
          "Simulation_small", "Simulation_large", "Empirical"
        )
      )
    )

  return(summary)
}


# read all independent Parquet inputs in one future/furrr batch
read.sfs.inputs <- function(
    sim.small.directory, sim.large.directory, empirical.directory,
    chromosomes
) {
  chromosome.tasks <- crossing(
    data.type = c(
      "Simulation_small", "Simulation_large", "Empirical"
    ),
    chrom = chromosomes
  ) %>%
    mutate(
      directory = case_when(
        data.type == "Simulation_small" ~ sim.small.directory,
        data.type == "Simulation_large" ~ sim.large.directory,
        TRUE ~ empirical.directory
      ),
      kind = "chromosome",
      path = file.path(
        path.expand(directory), glue("sfs.chr{chrom}.parquet")
      )
    )
  extra.tasks <- tibble(
    data.type = c("Empirical", "Empirical"),
    chrom = c("all", "all"),
    directory = empirical.directory,
    kind = c("genome", "variant.qc"),
    path = file.path(
      path.expand(empirical.directory),
      c("sfs.parquet", "variant_qc.parquet")
    )
  )
  tasks <- bind_rows(chromosome.tasks, extra.tasks)
  previous.plan <- future::plan()
  on.exit(future::plan(previous.plan), add = TRUE)
  future::plan(future::multisession)
  tables <- furrr::future_map(
    tasks$path,
    nanoparquet::read_parquet,
    .options = furrr::furrr_options(seed = TRUE)
  )
  chromosome.rows <- which(tasks$kind == "chromosome")
  chromosome.tables <- map2_dfr(
    tables[chromosome.rows], chromosome.rows,
    function(table, index) {
      table$chrom <- tasks$chrom[index]
      table$data.type.input <- tasks$data.type[index]
      return(table)
    }
  )
  simulation.small <- chromosome.tables %>%
    filter(data.type.input == "Simulation_small") %>%
    select(-data.type.input) %>%
    normalize.simulation.sfs("Simulation_small")
  simulation.large <- chromosome.tables %>%
    filter(data.type.input == "Simulation_large") %>%
    select(-data.type.input) %>%
    normalize.simulation.sfs("Simulation_large")
  empirical.chromosomes <- chromosome.tables %>%
    filter(data.type.input == "Empirical") %>%
    select(-data.type.input) %>%
    normalize.empirical.sfs()
  genome.index <- which(tasks$kind == "genome")
  empirical.genome <- tables[[genome.index]]
  empirical.genome$chrom <- "all"
  empirical.genome <- normalize.empirical.sfs(empirical.genome)
  qc.index <- which(tasks$kind == "variant.qc")
  inputs <- list(
    simulation.small = simulation.small,
    simulation.large = simulation.large,
    empirical.chromosomes = empirical.chromosomes,
    empirical.genome = empirical.genome,
    variant.qc = tables[[qc.index]]
  )

  return(inputs)
}


# apply shared colors, labels, scales, guides, and theme
style.sfs.plot <- function(
    plot, title, subtitle, x.label, y.label, styles,
    include.linetype = FALSE
) {
  plot <- plot +
    scale_color_manual(values = styles$population.colors) +
    scale_fill_manual(values = styles$population.colors) +
    scale_x_continuous(
      limits = c(1, DISPLAY.BIN.MAX),
      breaks = seq_len(DISPLAY.BIN.MAX)
    ) +
    scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
    labs(
      x = x.label, y = y.label, title = title,
      subtitle = subtitle, color = NULL, fill = NULL,
      linetype = NULL
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
  if (include.linetype) {
    plot <- plot +
      scale_linetype_manual(
        values = styles$series.linetypes,
        labels = styles$series.labels
      ) +
      guides(
        color = guide_legend(
          order = 1,
          override.aes = list(linetype = "solid")
        ),
        fill = "none",
        linetype = guide_legend(
          order = 2,
          override.aes = list(color = "black")
        )
      )
  } else {
    plot <- plot + guides(fill = "none")
  }

  return(plot)
}


# add simulation ribbons and source lines
add.sfs.geometries <- function(plot, data, include.linetype) {
  simulation <- data %>% filter(data.type != "Empirical")
  if (include.linetype) {
    plot <- plot +
      geom_ribbon(
        data = simulation,
        aes(
          ymin = ribbon.lower,
          ymax = ribbon.upper,
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
  } else {
    plot <- plot +
      geom_ribbon(
        data = simulation,
        aes(
          ymin = ribbon.lower,
          ymax = ribbon.upper,
          group = pop
        ),
        alpha = 0.2, color = NA
      ) +
      geom_line(
        data = data,
        aes(group = pop),
        linewidth = 1.25
      )
  }

  return(plot)
}


# construct one native-axis standalone SFS plot
make.one.standalone.sfs.plot <- function(
    data, data.type.input, value.type, styles
) {
  selected <- data %>%
    filter(
      data.type == data.type.input,
      chrom %in% c("all", SELECTED.CHROMOSOMES)
    ) %>%
    filter.display.bins() %>%
    mutate(
      chrom = factor(
        as.character(chrom),
        levels = c("all", SELECTED.CHROMOSOMES)
      )
    )
  is.empirical <- data.type.input == "Empirical"
  x.label <- if (is.empirical) {
    "Minor allele count"
  } else {
    "Derived allele count"
  }
  y.label <- if (value.type == "count") {
    "Mean number of sites"
  } else {
    "Mean proportion of segregating sites"
  }
  value.label <- if (value.type == "count") {
    "Raw Counts"
  } else {
    "Proportions"
  }
  data.label <- styles$series.labels[[data.type.input]]
  replicate.count <- max(selected$replicate.count)
  subtitle <- if (is.empirical) {
    "Fixed folded empirical spectrum · Genome and selected chromosomes"
  } else {
    glue(
      "Unfolded simulation spectrum · Replicates: ",
      "{replicate.count} · Genome and selected chromosomes"
    )
  }
  plot <- ggplot(
    selected,
    aes(
      x = allele.count, y = mean,
      color = pop, fill = pop
    )
  ) %>%
    add.sfs.geometries(selected, include.linetype = FALSE) +
    facet_wrap(~chrom, ncol = 4, drop = FALSE)
  plot <- style.sfs.plot(
    plot,
    glue("{data.label} SFS {value.label}"),
    subtitle, x.label, y.label, styles
  )

  return(plot)
}


# construct the six named standalone plots
make.standalone.sfs.plots <- function(
    raw.summary, proportion.summary, styles
) {
  plots <- list(
    Simulation_small.raw.counts = make.one.standalone.sfs.plot(
      raw.summary, "Simulation_small", "count", styles
    ),
    Simulation_small.proportions = make.one.standalone.sfs.plot(
      proportion.summary, "Simulation_small", "proportion", styles
    ),
    Simulation_large.raw.counts = make.one.standalone.sfs.plot(
      raw.summary, "Simulation_large", "count", styles
    ),
    Simulation_large.proportions = make.one.standalone.sfs.plot(
      proportion.summary, "Simulation_large", "proportion", styles
    ),
    Empirical.raw.counts = make.one.standalone.sfs.plot(
      raw.summary, "Empirical", "count", styles
    ),
    Empirical.proportions = make.one.standalone.sfs.plot(
      proportion.summary, "Empirical", "proportion", styles
    )
  )

  return(plots)
}


# construct the four combined comparison views
make.combined.sfs.plots <- function(data, value.type, styles) {
  data <- data %>% filter.display.bins()
  genome <- data %>% filter(chrom == "all")
  selected <- data %>%
    filter(chrom %in% SELECTED.CHROMOSOMES) %>%
    mutate(
      chrom = factor(
        as.character(chrom), levels = SELECTED.CHROMOSOMES
      )
    )
  replicate.count <- max(
    data$replicate.count[data$data.type != "Empirical"]
  )
  comparison.label <- if (value.type == "count") {
    "Projected and folded counts"
  } else {
    "Native-size folded proportions"
  }
  y.label <- if (value.type == "count") {
    "Mean number of sites"
  } else {
    "Mean proportion of segregating sites"
  }
  genome.subtitle <- glue(
    "{comparison.label} · Simulation replicates: {replicate.count} · ",
    "Scope: genome-wide"
  )
  selected.subtitle <- glue(
    "{comparison.label} · Simulation replicates: {replicate.count} · ",
    "Scope: selected chromosomes"
  )
  genome.base <- ggplot(
    genome,
    aes(
      x = allele.count, y = mean,
      color = pop, fill = pop
    )
  )
  selected.base <- ggplot(
    selected,
    aes(
      x = allele.count, y = mean,
      color = pop, fill = pop
    )
  )
  genome.by.data.set <- add.sfs.geometries(
    genome.base, genome, include.linetype = TRUE
  ) +
    facet_wrap(~data.type)
  genome.by.data.set <- style.sfs.plot(
    genome.by.data.set,
    "Genome-wide SFS by Data Set",
    genome.subtitle,
    "Minor allele count", y.label, styles,
    include.linetype = TRUE
  )
  genome.by.role <- add.sfs.geometries(
    genome.base, genome, include.linetype = TRUE
  ) +
    facet_wrap(~role)
  genome.by.role <- style.sfs.plot(
    genome.by.role,
    "Genome-wide SFS by Population Role",
    genome.subtitle,
    "Minor allele count", y.label, styles,
    include.linetype = TRUE
  )
  genome.populations <- genome.base %>%
    add.sfs.geometries(genome, include.linetype = TRUE) %>%
    style.sfs.plot(
      "Genome-wide SFS Across Populations",
      genome.subtitle,
      "Minor allele count", y.label, styles,
      include.linetype = TRUE
    )
  selected.chromosomes <- add.sfs.geometries(
    selected.base, selected, include.linetype = TRUE
  ) +
    facet_wrap(~chrom, ncol = 3, drop = FALSE)
  selected.chromosomes <- style.sfs.plot(
    selected.chromosomes,
    "SFS Across Selected Chromosomes",
    selected.subtitle,
    "Minor allele count", y.label, styles,
    include.linetype = TRUE
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


# read chromosome SFS, empirical genome SFS, and empirical QC in parallel
sfs.inputs <- read.sfs.inputs(
  SIM.SMALL.DATA.DIR,
  SIM.LARGE.DATA.DIR,
  EMPIRICAL.DATA.DIR,
  CHROMOSOMES
)

# validate native spectra and reconstruct simulation genome spectra
simulation.small.chromosomes <- sfs.inputs$simulation.small %>%
  validate.simulation.bins()
simulation.large.chromosomes <- sfs.inputs$simulation.large %>%
  validate.simulation.bins()
simulation.chromosomes <- bind_rows(
  simulation.small.chromosomes,
  simulation.large.chromosomes
)
simulation.genome <- aggregate.simulation.genome(
  simulation.chromosomes
)
simulation.native <- bind_rows(
  simulation.chromosomes,
  simulation.genome
)

# attach empirical targets without folding empirical spectra again
empirical.target.sizes <- derive.empirical.target.sizes(
  sfs.inputs$variant.qc
)
empirical.native <- bind_rows(
  sfs.inputs$empirical.chromosomes,
  sfs.inputs$empirical.genome
) %>%
  prepare.empirical.comparable(empirical.target.sizes)

# construct native raw counts and native-axis proportions
native.raw <- bind_rows(
  simulation.native %>% filter.segregating.bins(),
  empirical.native %>% filter.segregating.bins()
)
native.proportions <- bind_rows(
  simulation.native %>% calculate.sfs.proportions(),
  empirical.native %>% calculate.sfs.proportions()
)
native.raw.summary <- summarize.sfs.curves(native.raw, "count")
native.proportion.summary <- summarize.sfs.curves(
  native.proportions, "proportion"
)

# project simulations to empirical sizes and fold comparable counts
projected.simulation <- project.simulation.sfs(
  simulation.native,
  empirical.target.sizes
) %>%
  fold.simulation.sfs() %>%
  filter.segregating.bins()
projected.counts <- bind_rows(
  projected.simulation,
  empirical.native %>% filter.segregating.bins()
)
projected.count.summary <- summarize.sfs.curves(
  projected.counts, "count"
)

# fold simulations at native sizes and normalize full segregating spectra
comparable.proportions <- bind_rows(
  simulation.native %>%
    fold.simulation.sfs() %>%
    calculate.sfs.proportions(),
  empirical.native %>% calculate.sfs.proportions()
)
comparable.proportion.summary <- summarize.sfs.curves(
  comparable.proportions, "proportion"
)

# construct and print all retained plot collections sequentially
standalone.sfs.plots <- make.standalone.sfs.plots(
  native.raw.summary,
  native.proportion.summary,
  PLOT.STYLES
)
combined.projected.count.plots <- make.combined.sfs.plots(
  projected.count.summary,
  "count",
  PLOT.STYLES
)
combined.proportion.plots <- make.combined.sfs.plots(
  comparable.proportion.summary,
  "proportion",
  PLOT.STYLES
)
walk(standalone.sfs.plots, print)
walk(combined.projected.count.plots, print)
walk(combined.proportion.plots, print)
