#' @title Launch Interactive scSimEval Benchmarking Studio
#' @description Launches the Shiny web application embedded in the \pkg{scSimEval} package.
#'   Provides an intuitive graphical interface to ingest single-cell reference and simulated datasets,
#'   configure ground-truth-free evaluation pipelines, interactively inspect the flagship 62-measure
#'   bubble matrix, adjust category weighting for custom method rankings, and export
#'   publication-ready figures and tables.
#'
#' @param port Optional port number for the local web server. Default is \code{NULL} (random open port).
#' @param host Character string specifying the IP address to listen on. Defaults to \code{"127.0.0.1"}.
#' @param launch.browser Logical, whether to automatically launch the default web browser.
#'   Defaults to \code{TRUE} in interactive sessions.
#'
#' @return Invisibly returns the Shiny app process object.
#' @export
#'
#' @examples
#' \dontrun{
#' library(scSimEval)
#' # Launch the interactive studio
#' launch_scSimEval_app()
#'
#' }
launch_scSimEval_app <- function(port = NULL, host = "127.0.0.1", launch.browser = interactive()) {
  # Verify suggested packages required for the interactive UI
  required_pkgs <- c("shiny", "bslib", "DT", "ggplot2", "Matrix")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  
  if (length(missing_pkgs) > 0) {
    stop(
      sprintf(
        "The following suggested packages are required to run the scSimEval Shiny app: %s.\nPlease install them via:\n  install.packages(c(%s))",
        paste(paste0("'", missing_pkgs, "'"), collapse = ", "),
        paste(paste0("\"", missing_pkgs, "\""), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  app_dir <- system.file("shiny", "scSimEvalApp", package = "scSimEval")
  if (app_dir == "" || !dir.exists(app_dir)) {
    # Fallback check for local development checkout
    fallback_dir <- file.path("inst", "shiny", "scSimEvalApp")
    if (dir.exists(fallback_dir)) {
      app_dir <- fallback_dir
    } else {
      stop("Could not find the 'scSimEvalApp' directory in the installed package. Please ensure scSimEval is properly installed.", call. = FALSE)
    }
  }
  
  message("Starting scSimEval Benchmarking Studio...")
  message(sprintf("Listening on http://%s:%s", host, ifelse(is.null(port), "random", port)))
  
  shiny::runApp(
    appDir = app_dir,
    port = port,
    host = host,
    launch.browser = launch.browser,
    display.mode = "normal"
  )
}
