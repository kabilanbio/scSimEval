test_that("Shiny application files exist and structure is valid", {
  app_dir <- system.file("shiny", "scSimEvalApp", package = "scSimEval")
  if (app_dir == "") {
    app_dir <- file.path("..", "..", "inst", "shiny", "scSimEvalApp")
  }
  
  expect_true(dir.exists(app_dir))
  expect_true(file.exists(file.path(app_dir, "app.R")))
  
  demo_file <- file.path(app_dir, "data", "demo_benchmark_data.rds")
  expect_true(file.exists(demo_file))
  
  demo_obj <- readRDS(demo_file)
  expect_true(is.list(demo_obj))
  expect_true("benchmark_summary_table" %in% names(demo_obj))
  expect_equal(length(unique(demo_obj$benchmark_summary_table$Metric)), 62)
  expect_equal(length(unique(demo_obj$benchmark_summary_table$Category)), 8)
})

test_that("Shiny app object instantiates cleanly", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  
  app_dir <- system.file("shiny", "scSimEvalApp", package = "scSimEval")
  if (app_dir == "") {
    app_dir <- file.path("..", "..", "inst", "shiny", "scSimEvalApp")
  }
  
  app <- shiny::shinyAppDir(app_dir)
  expect_s3_class(app, "shiny.appobj")
})
