find_test_app_dir <- function() {
  local_cand <- c(
    file.path("..", "..", "inst", "shiny", "scSimEvalApp"),
    file.path("inst", "shiny", "scSimEvalApp")
  )
  for (cand in local_cand) {
    if (dir.exists(cand)) return(normalizePath(cand, winslash = "/"))
  }
  system.file("shiny", "scSimEvalApp", package = "scSimEval")
}

test_that("Shiny application files exist and structure is valid", {
  app_dir <- find_test_app_dir()
  expect_true(dir.exists(app_dir))
  expect_true(file.exists(file.path(app_dir, "app.R")))
  
  demo_file <- file.path(app_dir, "data", "demo_benchmark_data.rds")
  expect_true(file.exists(demo_file))
  
  demo_obj <- readRDS(demo_file)
  expect_true(is.list(demo_obj))
  expect_true("benchmark_summary_table" %in% names(demo_obj))
  expect_gte(length(unique(demo_obj$benchmark_summary_table$Metric)), 40)
  expect_gte(length(unique(demo_obj$benchmark_summary_table$Category)), 6)
})

test_that("Shiny app object instantiates cleanly", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  
  app_dir <- find_test_app_dir()
  app <- shiny::shinyAppDir(app_dir)
  expect_s3_class(app, "shiny.appobj")
})

test_that("Shiny app and package launcher configure maxRequestSize to 10 GB", {
  app_dir <- find_test_app_dir()
  app_lines <- readLines(file.path(app_dir, "app.R"))
  expect_true(any(grepl("shiny\\.maxRequestSize\\s*=\\s*10\\s*\\*\\s*1024\\^3", app_lines)))
  
  # Check default argument in launch_scSimEval_app
  expect_equal(formals(launch_scSimEval_app)$maxRequestSize, quote(10 * 1024^3))
})
