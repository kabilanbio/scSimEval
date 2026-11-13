#' Normalize Time Unit String
#' @keywords internal
.normalize_time_unit <- function(time_unit = "seconds") {
  if (is.null(time_unit) || length(time_unit) == 0 || is.na(time_unit[1])) {
    return("seconds")
  }
  u <- tolower(trimws(as.character(time_unit[1])))
  if (u %in% c("s", "sec", "secs", "second", "seconds")) {
    return("seconds")
  }
  if (u %in% c("m", "min", "mins", "minute", "minutes")) {
    return("minutes")
  }
  stop(sprintf("Unsupported `time_unit`: '%s'. Use 'seconds' (or 's') or 'minutes' (or 'min').", time_unit[1]), call. = FALSE)
}

#' Normalize Memory Unit String
#' @keywords internal
.normalize_memory_unit <- function(memory_unit = "MiB") {
  if (is.null(memory_unit) || length(memory_unit) == 0 || is.na(memory_unit[1])) {
    return("MiB")
  }
  u <- tolower(trimws(as.character(memory_unit[1])))
  if (u %in% c("mib", "mebibyte", "mebibytes")) {
    return("MiB")
  }
  if (u %in% c("mb", "megabyte", "megabytes")) {
    return("MB")
  }
  stop(sprintf("Unsupported `memory_unit`: '%s'. Use 'MiB' or 'MB'.", memory_unit[1]), call. = FALSE)
}

#' Convert Elapsed Runtime to Canonical Seconds
#' @keywords internal
.convert_time_to_seconds <- function(elapsed_time, time_unit = "seconds") {
  if (is.null(elapsed_time) || length(elapsed_time) == 0) return(NULL)
  val <- as.numeric(elapsed_time[1])
  if (is.na(val)) return(NA_real_)
  u <- .normalize_time_unit(time_unit)
  if (identical(u, "minutes")) val * 60 else val
}

#' Convert Canonical Elapsed Seconds to Requested Display Unit
#' @keywords internal
.convert_time_from_seconds <- function(elapsed_seconds, time_unit = "seconds") {
  if (is.null(elapsed_seconds) || length(elapsed_seconds) == 0) return(elapsed_seconds)
  u <- .normalize_time_unit(time_unit)
  val <- as.numeric(elapsed_seconds)
  if (identical(u, "minutes")) val / 60 else val
}

#' Convert Memory Usage to Canonical MiB (Binary Mebibytes, 1024^2 bytes)
#' @keywords internal
.convert_memory_to_mib <- function(memory_val, memory_unit = "MiB") {
  if (is.null(memory_val) || length(memory_val) == 0) return(NULL)
  val <- as.numeric(memory_val[1])
  if (is.na(val)) return(NA_real_)
  u <- .normalize_memory_unit(memory_unit)
  if (identical(u, "MB")) val * (1000^2 / 1024^2) else val
}

#' Convert Canonical MiB to Requested Display Memory Unit
#' @keywords internal
.convert_memory_from_mib <- function(memory_mib, memory_unit = "MiB") {
  if (is.null(memory_mib) || length(memory_mib) == 0) return(memory_mib)
  u <- .normalize_memory_unit(memory_unit)
  val <- as.numeric(memory_mib)
  if (identical(u, "MB")) val * (1024^2 / 1000^2) else val
}

#' Validate and Normalize Simulation CPU Core Count
#' @keywords internal
.normalize_sim_cores <- function(sim_cores = NULL) {
  if (is.null(sim_cores) || length(sim_cores) == 0 || is.na(sim_cores[1])) {
    return(NULL)
  }
  val <- suppressWarnings(as.integer(sim_cores[1]))
  if (is.na(val) || val < 1L) {
    stop("`sim_cores` must be a positive integer >= 1.", call. = FALSE)
  }
  val
}

#' Compute Core-Seconds Runtime Metric
#'
#' Quantifies cumulative CPU work across all allocated processing cores
#' (\eqn{\text{Elapsed Time (s)} \times N_{\text{cores}}}; Lannelongue et al., 2021;
#' Wratten et al., 2021; Grealey et al., 2022).
#'
#' @param elapsed_time Numeric elapsed wall-clock runtime.
#' @param time_unit Time unit of \code{elapsed_time}: \code{"seconds"} (default) or \code{"minutes"}.
#' @param sim_cores Positive integer specifying the number of CPU cores allocated (default \code{1L}).
#' @return Numeric Core-Seconds (\eqn{\text{core}\cdot\text{s}}).
#' @references
#' Lannelongue, L., Grealey, J., & Inouye, M. (2021). Green Algorithms: Quantifying the carbon footprint
#' of computation. \emph{Advanced Science}, 8(12), 2100707.
#'
#' Wratten, L., Wilm, A., & Göke, J. (2021). Reproducible, scalable, and shareable analysis pipelines
#' with bioinformatics workflow managers. \emph{Nature Methods}, 18, 1161–1168.
#'
#' Grealey, J., et al. (2022). The Carbon Footprint of Bioinformatics. \emph{Molecular Biology and Evolution},
#' 39(3), msac034.
#' @export
compute_core_seconds <- function(elapsed_time, time_unit = "seconds", sim_cores = 1L) {
  sec <- .convert_time_to_seconds(elapsed_time, time_unit)
  if (is.null(sec) || is.na(sec)) return(NA_real_)
  cores <- .normalize_sim_cores(sim_cores)
  if (is.null(cores)) cores <- 1L
  as.numeric(sec * cores)
}

#' Compute Peak RAM per Allocated CPU Core
#'
#' Normalizes peak memory consumption by the number of allocated processor cores
#' (\eqn{\text{Peak RAM} / N_{\text{cores}}}; Holtgrewe et al., 2015; Boisvert et al., 2012;
#' Chen et al., 2022).
#'
#' @param memory_val Numeric peak memory consumption.
#' @param memory_unit Memory unit of \code{memory_val}: \code{"MiB"} (default) or \code{"MB"}.
#' @param sim_cores Positive integer specifying the number of CPU cores allocated (default \code{1L}).
#' @param target_unit Target memory unit: \code{"MiB"} (MiB/core, default), \code{"GiB"} (GiB/core),
#'   \code{"MB"} (MB/core), or \code{"GB"} (GB/core).
#' @return Numeric Peak RAM per allocated core in requested target unit.
#' @references
#' Holtgrewe, M., Kuchenbecker, L., & Reinert, K. (2015). Methods for the detection and assembly of
#' novel sequence in high-throughput sequencing data. \emph{Bioinformatics}, 31(12), 1904–1912.
#'
#' Boisvert, S., Raymond, F., Godzaridis, E., Laviolette, F., & Corbeil, J. (2012). Ray Meta: scalable
#' de novo metagenome assembly and profiling. \emph{Genome Biology}, 13(12), R122.
#'
#' Chen, W., et al. (2022). STEPS 4.0: Fast and memory-efficient molecular simulations of neurons
#' at the nanoscale. \emph{Frontiers in Neuroinformatics}, 16, 883742.
#' @export
compute_memory_per_core <- function(memory_val, memory_unit = "MiB", sim_cores = 1L,
                                    target_unit = c("MiB", "GiB", "MB", "GB")) {
  mib <- .convert_memory_to_mib(memory_val, memory_unit)
  if (is.null(mib) || is.na(mib)) return(NA_real_)
  cores <- .normalize_sim_cores(sim_cores)
  if (is.null(cores)) cores <- 1L
  mib_per_core <- as.numeric(mib / cores)

  t_u <- tolower(trimws(as.character(target_unit[1])))
  if (t_u %in% c("gib", "gib/core", "gb_bin")) {
    mib_per_core / 1024
  } else if (t_u %in% c("mb", "mb/core", "megabyte")) {
    mib_per_core * (1024^2 / 1000^2)
  } else if (t_u %in% c("gb", "gb/core", "gigabyte")) {
    (mib_per_core * (1024^2 / 1000^2)) / 1000
  } else {
    mib_per_core
  }
}

#' Benchmark Execution Runtime and Memory Usage
#'
#' Times an expression and measures change in memory allocation, computing both raw
#' runtime/memory metrics and core-normalized Category VIII scalability indicators
#' (\code{core_seconds} and \code{memory_per_core_mib}).
#'
#' @param expr Expression or function to benchmark.
#' @param sim_cores Optional positive integer specifying CPU cores used during execution (default \code{1L}).
#' @return A list with elapsed seconds/minutes, peak memory allocated (MiB and MB),
#'   \code{core_seconds}, \code{core_hours}, \code{memory_per_core_mib}, \code{memory_per_core_gib},
#'   and \code{sim_cores}.
#' @references
#' Lannelongue, L., et al. (2021). Green Algorithms. \emph{Advanced Science}, 8(12), 2100707.
#' Holtgrewe, M., et al. (2015). \emph{Bioinformatics}, 31(12), 1904–1912.
#' @examples
#' benchmark_resource_usage({ Sys.sleep(0.01); 1 + 1 })
#' @export
benchmark_resource_usage <- function(expr, sim_cores = NULL) {
  sim_cores_norm <- .normalize_sim_cores(sim_cores)
  cores_active <- if (!is.null(sim_cores_norm)) sim_cores_norm else 1L
  gc(verbose = FALSE, reset = TRUE)
  start_time <- proc.time()
  
  res <- tryCatch(
    eval.parent(substitute(expr)),
    error = function(e) {
      message("Expression evaluation failed: ", e$message)
      NULL
    }
  )
  
  elapsed_time <- (proc.time() - start_time)[["elapsed"]]
  mem_info <- gc()
  
  # Approximate MiB used in Vcells + Ncells (1024^2 bytes)
  mem_mib <- sum(mem_info[, 2]) * (8 / (1024^2))
  mem_mb_si <- .convert_memory_from_mib(mem_mib, "MB")

  core_sec <- as.numeric(elapsed_time * cores_active)
  mem_per_core_mib <- as.numeric(mem_mib / cores_active)
  
  out <- list(
    result = res,
    core_seconds = core_sec,
    core_hours = core_sec / 3600,
    memory_per_core_mib = mem_per_core_mib,
    memory_per_core_gib = mem_per_core_mib / 1024,
    memory_per_core_mb = mem_mb_si / cores_active,
    memory_per_core_gb = (mem_mb_si / cores_active) / 1000,
    elapsed_seconds = elapsed_time,
    elapsed_minutes = elapsed_time / 60,
    memory_mb = mem_mib,
    memory_mib = mem_mib,
    memory_mb_si = mem_mb_si
  )
  if (!is.null(sim_cores_norm)) {
    out$sim_cores <- sim_cores_norm
  }
  out
}
