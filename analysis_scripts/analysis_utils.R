# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# analysis_utils.R
# ______________________________________________________________________________


# constants ----
SIM.STATS.DIR <- path.expand("~/scratch/OOA_NAAdmixture_small/stats")
ONEKG.STATS.DIR <- path.expand("~/scratch/OOA_NAAdmixture_1kG/stats")

POP.COLORS <- c(
  AFR = "#56B4E9",
  EUR = "#fb8072",
  ADMIXED = "#4B1FA8"
)

POPULATION.CROSSWALK <- data.frame(
  data_source = c(
    rep("simulation", 3),
    rep("empirical", 3)
  ),
  pop = c("AFR", "EUR", "ADX", "YRI", "CEU", "ASW"),
  pop_role = c("AFR", "EUR", "ADMIXED", "AFR", "EUR", "ADMIXED"),
  stringsAsFactors = FALSE
)


# schema validation and harmonization ----
typed.na <- function(type) {
  switch(
    type,
    character = NA_character_,
    integer = NA_integer_,
    double = NA_real_,
    logical = NA,
    stop("Unsupported canonical type: ", type, call. = FALSE)
  )
}

cast.character <- function(value, column) {
  as.character(value)
}

cast.double <- function(value, column) {
  if (is.factor(value)) {
    value <- as.character(value)
  }
  converted <- suppressWarnings(as.numeric(value))
  invalid <- !is.na(value) & is.na(converted)
  if (any(invalid)) {
    stop(
      "Unsafe double conversion for column ",
      column,
      call. = FALSE
    )
  }
  converted
}

cast.integer <- function(value, column) {
  numeric.value <- cast.double(value, column)
  unsafe <- !is.na(numeric.value) & (
    !is.finite(numeric.value) |
      abs(numeric.value - round(numeric.value)) >
        sqrt(.Machine$double.eps) |
      numeric.value > .Machine$integer.max |
      numeric.value < -.Machine$integer.max
  )
  if (any(unsafe)) {
    stop(
      "Unsafe integer conversion for column ",
      column,
      call. = FALSE
    )
  }
  as.integer(numeric.value)
}

cast.column <- function(value, type, column) {
  switch(
    type,
    character = cast.character(value, column),
    integer = cast.integer(value, column),
    double = cast.double(value, column),
    logical = as.logical(value),
    stop("Unsupported canonical type: ", type, call. = FALSE)
  )
}

rename.columns <- function(data, rename.map = NULL) {
  if (is.null(rename.map) || length(rename.map) == 0) {
    return(data)
  }
  old.names <- names(rename.map)
  present <- old.names %in% names(data)
  names(data)[match(old.names[present], names(data))] <-
    unname(rename.map[present])
  data
}

add.population.role <- function(data, data.source) {
  crosswalk <- POPULATION.CROSSWALK[
    POPULATION.CROSSWALK$data_source == data.source,
    c("pop", "pop_role")
  ]
  unknown <- setdiff(unique(data$pop), crosswalk$pop)
  if (length(unknown) > 0) {
    stop(
      "Unknown ",
      data.source,
      " population label(s): ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  data$pop_role <- crosswalk$pop_role[match(data$pop, crosswalk$pop)]
  data
}

harmonize.table <- function(
    data,
    schema,
    required.columns,
    data.source,
    analysis.level,
    rename.map = NULL) {
  data <- rename.columns(as.data.frame(data), rename.map)
  missing.required <- setdiff(required.columns, names(data))
  if (length(missing.required) > 0) {
    stop(
      "Missing required column(s): ",
      paste(missing.required, collapse = ", "),
      call. = FALSE
    )
  }

  for (column in names(schema)) {
    if (!column %in% names(data)) {
      data[[column]] <- rep(typed.na(schema[[column]]), nrow(data))
    }
    data[[column]] <- cast.column(
      data[[column]],
      schema[[column]],
      column
    )
  }
  data <- data[names(schema)]
  if (analysis.level == "genome") {
    data$chrom <- rep("all", nrow(data))
  }
  data$data_source <- rep(data.source, nrow(data))
  data$analysis_level <- rep(analysis.level, nrow(data))
  add.population.role(data, data.source)
}


# common file discovery and input ----
discover.chromosomes <- function(stats.dir, table.name) {
  if (!dir.exists(stats.dir)) {
    stop("Missing statistics directory: ", stats.dir, call. = FALSE)
  }
  pattern <- paste0("^", table.name, "\\.chr([0-9]+)\\.parquet$")
  files <- list.files(stats.dir, pattern = pattern, full.names = FALSE)
  if (length(files) == 0) {
    stop(
      "No chromosome Parquet files found for ",
      table.name,
      " in ",
      stats.dir,
      call. = FALSE
    )
  }
  chroms <- sub(pattern, "\\1", files)
  as.character(sort(unique(as.integer(chroms))))
}

require.genome.file <- function(stats.dir, table.name) {
  path <- file.path(stats.dir, paste0(table.name, ".parquet"))
  if (!file.exists(path)) {
    stop("Missing genome Parquet file: ", path, call. = FALSE)
  }
  path
}

matching.chromosomes <- function(
    sim.table,
    empirical.table,
    sim.stats.dir = SIM.STATS.DIR,
    onekg.stats.dir = ONEKG.STATS.DIR) {
  sim.chroms <- discover.chromosomes(sim.stats.dir, sim.table)
  empirical.chroms <- discover.chromosomes(
    onekg.stats.dir,
    empirical.table
  )
  if (!identical(sim.chroms, empirical.chroms)) {
    stop(
      "Unmatched chromosome outputs for ",
      sim.table,
      " and ",
      empirical.table,
      ": simulation={",
      paste(sim.chroms, collapse = ","),
      "}; empirical={",
      paste(empirical.chroms, collapse = ","),
      "}",
      call. = FALSE
    )
  }
  require.genome.file(sim.stats.dir, sim.table)
  require.genome.file(onekg.stats.dir, empirical.table)
  sim.chroms
}

read.parquet.paths <- function(paths) {
  required.packages <- c("future", "furrr", "nanoparquet", "vctrs")
  unavailable <- required.packages[
    !vapply(required.packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(unavailable) > 0) {
    stop(
      "Missing required R package(s): ",
      paste(unavailable, collapse = ", "),
      call. = FALSE
    )
  }

  previous.plan <- future::plan()
  on.exit(future::plan(previous.plan), add = TRUE)
  workers <- min(length(paths), future::availableCores())
  future::plan(future::multisession, workers = max(1L, workers))
  furrr::future_map(
    paths,
    nanoparquet::read_parquet,
    .options = furrr::furrr_options(seed = TRUE)
  )
}

read.source.table <- function(
    stats.dir,
    table.name,
    chroms,
    schema,
    required.columns,
    data.source,
    rename.map = NULL) {
  chr.paths <- file.path(
    stats.dir,
    paste0(table.name, ".chr", chroms, ".parquet")
  )
  genome.path <- require.genome.file(stats.dir, table.name)
  missing.paths <- chr.paths[!file.exists(chr.paths)]
  if (length(missing.paths) > 0) {
    stop(
      "Missing chromosome Parquet file(s): ",
      paste(missing.paths, collapse = ", "),
      call. = FALSE
    )
  }

  chr.tables <- read.parquet.paths(chr.paths)
  chr.tables <- Map(
    function(data, chrom) {
      data$chrom <- rep(chrom, nrow(data))
      harmonize.table(
        data,
        schema,
        required.columns,
        data.source,
        "chromosome",
        rename.map
      )
    },
    chr.tables,
    chroms
  )
  genome.table <- read.parquet.paths(genome.path)[[1]]
  genome.table <- harmonize.table(
    genome.table,
    schema,
    required.columns,
    data.source,
    "genome",
    rename.map
  )
  vctrs::vec_rbind(!!!c(chr.tables, list(genome.table)))
}

read.harmonized.table <- function(
    sim.table,
    empirical.table,
    schema,
    sim.required,
    empirical.required,
    sim.rename = NULL,
    empirical.rename = NULL,
    sim.stats.dir = SIM.STATS.DIR,
    onekg.stats.dir = ONEKG.STATS.DIR) {
  chroms <- matching.chromosomes(
    sim.table,
    empirical.table,
    sim.stats.dir,
    onekg.stats.dir
  )
  sim <- read.source.table(
    sim.stats.dir,
    sim.table,
    chroms,
    schema,
    sim.required,
    "simulation",
    sim.rename
  )
  empirical <- read.source.table(
    onekg.stats.dir,
    empirical.table,
    chroms,
    schema,
    empirical.required,
    "empirical",
    empirical.rename
  )
  vctrs::vec_rbind(sim, empirical)
}


# predictive summaries ----
predictive.summary <- function(
    data,
    value.column,
    group.columns = character()) {
  required <- c("data_source", "rep", value.column, group.columns)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(
      "Predictive summary missing column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  sim <- data[data$data_source == "simulation", , drop = FALSE]
  empirical <- data[data$data_source == "empirical", , drop = FALSE]
  if (nrow(sim) == 0 || nrow(empirical) == 0) {
    stop(
      "Predictive summaries require simulation and empirical rows",
      call. = FALSE
    )
  }
  if (length(group.columns) > 0) {
    sim.keys <- dplyr::distinct(
      sim,
      dplyr::across(dplyr::all_of(group.columns))
    )
    empirical.keys <- dplyr::distinct(
      empirical,
      dplyr::across(dplyr::all_of(group.columns))
    )
    if (
      nrow(dplyr::anti_join(
        sim.keys,
        empirical.keys,
        by = group.columns
      )) > 0 ||
        nrow(dplyr::anti_join(
          empirical.keys,
          sim.keys,
          by = group.columns
        )) > 0
    ) {
      stop(
        "Simulation and empirical predictive groups do not match",
        call. = FALSE
      )
    }
  }

  value.symbol <- rlang::sym(value.column)
  sim.by.rep <- sim |>
    dplyr::group_by(dplyr::across(dplyr::all_of(
      c(group.columns, "rep")
    ))) |>
    dplyr::summarise(
      .value = mean(!!value.symbol, na.rm = TRUE),
      .groups = "drop"
    )
  sim.summary <- sim.by.rep |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group.columns))) |>
    dplyr::summarise(
      simulation_median = stats::median(.value, na.rm = TRUE),
      simulation_lower = stats::quantile(
        .value,
        0.025,
        na.rm = TRUE,
        names = FALSE
      ),
      simulation_upper = stats::quantile(
        .value,
        0.975,
        na.rm = TRUE,
        names = FALSE
      ),
      simulation_sd = stats::sd(.value, na.rm = TRUE),
      .simulation_values = list(.value),
      .groups = "drop"
    )
  empirical.summary <- empirical |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group.columns))) |>
    dplyr::summarise(
      empirical_value = mean(!!value.symbol, na.rm = TRUE),
      .groups = "drop"
    )
  joined <- dplyr::left_join(
    sim.summary,
    empirical.summary,
    by = group.columns
  )
  joined |>
    dplyr::rowwise() |>
    dplyr::mutate(
      empirical_percentile = mean(
        .simulation_values <= empirical_value,
        na.rm = TRUE
      ),
      absolute_discrepancy = abs(
        empirical_value - simulation_median
      ),
      relative_discrepancy = dplyr::if_else(
        simulation_median > 0,
        absolute_discrepancy / simulation_median,
        NA_real_
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-.simulation_values)
}

make.diversity.views <- function(data) {
  empirical.scopes <- c(
    "onekg_intergenic_callable",
    "onekg_full_callable"
  )
  stats::setNames(
    lapply(
      empirical.scopes,
      function(scope) {
        data[
          data$sequence_scope %in% c("simulated_full", scope),
          ,
          drop = FALSE
        ]
      }
    ),
    empirical.scopes
  )
}


# population contrasts ----
population.contrasts <- function(
    data,
    value.column,
    group.columns,
    metric.type = c("positive", "kinship")) {
  metric.type <- match.arg(metric.type)
  required <- c(group.columns, "pop_role", value.column)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(
      "Population contrast missing column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  value.symbol <- rlang::sym(value.column)
  wide <- data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(
      c(group.columns, "pop_role")
    ))) |>
    dplyr::summarise(
      .value = mean(!!value.symbol, na.rm = TRUE),
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(names_from = pop_role, values_from = .value)
  required.roles <- c("AFR", "EUR", "ADMIXED")
  if (!all(required.roles %in% names(wide))) {
    stop(
      "Population contrasts require AFR, EUR, and ADMIXED roles",
      call. = FALSE
    )
  }
  if (any(is.na(wide[required.roles]))) {
    stop(
      "Every population contrast group requires all population roles",
      call. = FALSE
    )
  }

  if (metric.type == "positive") {
    if (any(wide[required.roles] <= 0, na.rm = TRUE)) {
      stop(
        "Positive-metric population contrasts require positive values",
        call. = FALSE
      )
    }
    expressions <- list(
      "AFR/EUR" = log(wide$AFR / wide$EUR),
      "ADMIXED/EUR" = log(wide$ADMIXED / wide$EUR),
      "ADMIXED/AFR" = log(wide$ADMIXED / wide$AFR)
    )
  } else {
    expressions <- list(
      "AFR-EUR" = wide$AFR - wide$EUR,
      "ADMIXED-EUR" = wide$ADMIXED - wide$EUR,
      "ADMIXED-AFR" = wide$ADMIXED - wide$AFR
    )
  }

  dplyr::bind_rows(lapply(
    names(expressions),
    function(name) {
      dplyr::bind_cols(
        wide[group.columns],
        data.frame(
          contrast = name,
          value = expressions[[name]],
          stringsAsFactors = FALSE
        )
      )
    }
  ))
}


# kinship and LD helpers ----
normalize.kinship.identifiers <- function(data) {
  names(data)[names(data) == "id1"] <- "sample1"
  names(data)[names(data) == "id2"] <- "sample2"
  if (!all(c("sample1", "sample2") %in% names(data))) {
    stop(
      "Kinship table requires sample1/sample2 or id1/id2",
      call. = FALSE
    )
  }
  data
}

pair.weighted.ld <- function(data) {
  sum.r2 <- sum(data$sum_r2, na.rm = TRUE)
  n.pairs <- sum(data$n_pairs, na.rm = TRUE)
  if (n.pairs <= 0) {
    return(NA_real_)
  }
  sum.r2 / n.pairs
}

aggregate.ld <- function(data, group.columns) {
  data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group.columns))) |>
    dplyr::summarise(
      sum_r2 = sum(sum_r2, na.rm = TRUE),
      n_pairs = sum(n_pairs, na.rm = TRUE),
      mean_r2 = dplyr::if_else(
        n_pairs > 0,
        sum_r2 / n_pairs,
        NA_real_
      ),
      .groups = "drop"
    )
}


# folded SFS helpers ----
infer.source.haplotypes <- function(data) {
  bins <- sort(unique(data$derived_allele_count))
  expected <- 0:max(bins)
  if (!identical(as.integer(bins), as.integer(expected))) {
    stop(
      "Simulation SFS must contain a complete unfolded spectrum",
      call. = FALSE
    )
  }
  as.integer(max(bins))
}

project.unfolded.sfs <- function(
    data,
    source.haplotypes,
    target.haplotypes) {
  if (target.haplotypes > source.haplotypes) {
    stop(
      "Empirical target sample size exceeds simulation sample size",
      call. = FALSE
    )
  }
  required.bins <- 0:source.haplotypes
  if (!identical(
    sort(unique(as.integer(data$derived_allele_count))),
    required.bins
  )) {
    stop(
      "Simulation SFS must contain every unfolded allele-count bin",
      call. = FALSE
    )
  }

  projected <- numeric(target.haplotypes + 1L)
  for (row in seq_len(nrow(data))) {
    source.count <- as.integer(data$derived_allele_count[[row]])
    probabilities <- stats::dhyper(
      0:target.haplotypes,
      source.count,
      source.haplotypes - source.count,
      target.haplotypes
    )
    projected <- projected + data$count[[row]] * probabilities
  }
  data.frame(
    derived_allele_count = 0:target.haplotypes,
    count = projected
  )
}

fold.projected.sfs <- function(data, num.haplotypes) {
  midpoint <- floor(num.haplotypes / 2)
  counts <- numeric(midpoint + 1L)
  for (minor.count in 0:midpoint) {
    complement <- num.haplotypes - minor.count
    counts[[minor.count + 1L]] <- data$count[
      data$derived_allele_count == minor.count
    ]
    if (complement != minor.count) {
      counts[[minor.count + 1L]] <- counts[[minor.count + 1L]] +
        data$count[data$derived_allele_count == complement]
    }
  }
  data.frame(
    minor_allele_count = 0:midpoint,
    count = counts
  )
}

normalize.folded.sfs <- function(data, num.haplotypes, span) {
  if (
    length(num.haplotypes) != 1 ||
      is.na(num.haplotypes) ||
      !is.finite(num.haplotypes) ||
      num.haplotypes < 2 ||
      num.haplotypes != as.integer(num.haplotypes)
  ) {
    stop(
      "SFS normalization requires one integer haplotype count",
      call. = FALSE
    )
  }
  if (!is.finite(span) || span <= 0) {
    stop("SFS normalization requires a positive span", call. = FALSE)
  }
  minor.counts <- data$minor_allele_count
  invalid.counts <- (
    is.na(minor.counts) |
      minor.counts < 0 |
      minor.counts > floor(num.haplotypes / 2)
  )
  if (
    any(invalid.counts) ||
      anyDuplicated(minor.counts) ||
      any(!is.finite(data$count)) ||
      any(data$count < 0)
  ) {
    stop(
      "Folded SFS contains invalid bins or counts",
      call. = FALSE
    )
  }
  segregating <- data[data$minor_allele_count > 0, , drop = FALSE]
  total <- sum(segregating$count, na.rm = TRUE)
  if (total <= 0) {
    stop(
      "SFS normalization requires positive segregating-site mass",
      call. = FALSE
    )
  }
  segregating$num_haplotypes <- as.integer(num.haplotypes)
  segregating$minor_allele_frequency <-
    segregating$minor_allele_count / num.haplotypes
  segregating$fraction <- segregating$count / total
  segregating$sites_per_mb <- segregating$count / (span / 1e6)
  segregating
}

spectrum.hellinger <- function(simulation, empirical) {
  joined <- dplyr::full_join(
    simulation[c("minor_allele_count", "fraction")],
    empirical[c("minor_allele_count", "fraction")],
    by = "minor_allele_count",
    suffix = c("_simulation", "_empirical")
  ) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::starts_with("fraction_"),
        ~tidyr::replace_na(.x, 0)
      )
    )
  sqrt(
    sum(
      (
        sqrt(joined$fraction_simulation) -
          sqrt(joined$fraction_empirical)
      )^2
    ) / 2
  )
}
