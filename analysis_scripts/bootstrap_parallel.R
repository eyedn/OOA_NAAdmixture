# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# bootstrap_parallel.R
# ______________________________________________________________________________


# require the packages used to schedule reproducible bootstrap work.
if (!requireNamespace("future", quietly = TRUE) ||
    !requireNamespace("furrr", quietly = TRUE)) {
  stop("Parallel bootstraps require the future and furrr packages")
  }


# return the configured multisession worker count.
bootstrap.parallel.workers <- function() {
  configured <- Sys.getenv("OOA_BOOTSTRAP_WORKERS", unset = "")
  if (nzchar(configured)) {
    workers <- suppressWarnings(as.integer(configured))
    if (is.na(workers) || workers < 1L) {
      stop("OOA_BOOTSTRAP_WORKERS must be a positive integer")
      }
    return(workers)
    }
  available <- parallelly::availableCores()
  if (is.na(available) || available < 1L) {
    available <- parallel::detectCores()
    }
  if (is.na(available) || available < 1L) {
    available <- 1L
    }
  return(max(1L, as.integer(available) - 1L))
  }


# run ordered work-table rows with L'Ecuyer streams and restore the caller plan.
bootstrap.parallel.pmap <- function(work, worker, seed, workers = NULL) {
  if (!is.data.frame(work)) {
    stop("Parallel bootstrap work must be a data frame")
    }
  if (!is.function(worker)) {
    stop("Parallel bootstrap worker must be a function")
    }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
    stop("Parallel bootstrap seed must be one finite number")
    }
  if (is.null(workers)) {
    workers <- bootstrap.parallel.workers()
    }
  if (!is.numeric(workers) || length(workers) != 1L ||
      !is.finite(workers) || workers < 1L || workers %% 1L != 0) {
    stop("Parallel bootstrap workers must be a positive integer")
    }
  previous.plan <- future::plan()
  on.exit(future::plan(previous.plan), add = TRUE)
  future::plan(future::multisession, workers = as.integer(workers))
  return(furrr::future_pmap(
    work, worker,
    .options = furrr::furrr_options(seed = as.integer(seed))
    ))
  }
