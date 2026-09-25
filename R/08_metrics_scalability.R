#' Benchmark Execution Runtime and Memory Usage
#'
#' Times an expression and measures change in memory allocation.
#'
#' @param expr Expression or function to benchmark.
#' @return A list with elapsed seconds and peak memory allocated (MB).
#' @export
benchmark_resource_usage <- function(expr) {
  gc(verbose = FALSE, reset = TRUE)
  start_time <- proc.time()
  
  res <- tryCatch(
    eval.parent(substitute(expr)),
    error = function(e) {
      warning("Expression failed with error: ", e$message)
      NULL
    }
  )
  
  elapsed_time <- (proc.time() - start_time)[["elapsed"]]
  mem_info <- gc()
  
  # Approximate MB used in Vcells + Ncells
  mem_mb <- sum(mem_info[, 2]) * (8 / (1024^2))
  
  list(
    result = res,
    elapsed_seconds = elapsed_time,
    memory_mb = mem_mb
  )
}
