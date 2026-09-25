test_that("evaluate_multiple_datasets handles flat prefixed multiomics and unimodal datasets", {
  data("example_scrna", package = "scSimEval", envir = environment())
  data("example_scatac", package = "scSimEval", envir = environment())

  # Construct the user-specified collection:
  # Method 1 (paired multiomics), Method 2 (paired multiomics), Method 3 (unimodal ATAC), Method 4 (unimodal RNA)
  datasets <- list(
    "Method 1-scRNA-seq" = list(ref = example_scrna$ref, sim = example_scrna$sim),
    "Method 1-scATAC-seq" = list(ref = example_scatac$ref, sim = example_scatac$sim),
    "Method 2-scRNA-seq" = list(ref = example_scrna$ref, sim = example_scrna$sim),
    "Method 2-scATAC-seq" = list(ref = example_scatac$ref, sim = example_scatac$sim),
    "Method 3-scATAC-seq" = list(ref = example_scatac$ref, sim = example_scatac$sim),
    "Method 4-scRNA-seq" = list(ref = example_scrna$ref, sim = example_scrna$sim)
  )

  res <- evaluate_multiple_datasets(
    datasets,
    pair_by_prefix = TRUE,
    compute_bivariate = FALSE,
    verbose = FALSE
  )

  expect_s3_class(res, "scSimEval_consolidated")
  expect_true(is.data.frame(res$consolidated_summary_table))
  expect_true(is.matrix(res$consolidated_score_matrix))
  expect_true(is.data.frame(res$dataset_overview))

  # Check columns in consolidated summary table
  expected_cols <- c("Dataset", "Data_Name", "Modality", "Category", "Property", "Metric", "Value")
  expect_true(all(expected_cols %in% colnames(res$consolidated_summary_table)))

  # Check all datasets are represented
  unique_dsets <- unique(res$consolidated_summary_table$Dataset)
  expect_true(all(c("Method 1", "Method 2", "Method 3", "Method 4") %in% unique_dsets))

  # Check modalities
  unique_mods <- unique(res$consolidated_summary_table$Modality)
  expect_true("scRNA-seq" %in% unique_mods)
  expect_true("scATAC-seq" %in% unique_mods)
  expect_true("Joint (Cross-Modal)" %in% unique_mods)

  # Check print method
  output <- capture.output(print(res))
  expect_true(any(grepl("scSimEval Consolidated Multi-Dataset Benchmark Results", output)))

  # Test consolidated visualization
  p <- plot_consolidated_summary(res, category = "(I) Distributional Properties", metric = "KS")
  expect_s3_class(p, "ggplot")
})

test_that("evaluate_multiple_datasets handles grouped dataset structure", {
  data("example_scrna", package = "scSimEval", envir = environment())
  data("example_scatac", package = "scSimEval", envir = environment())

  grouped_datasets <- list(
    "Experiment_A" = list(
      ref = list("scRNA-seq" = example_scrna$ref, "scATAC-seq" = example_scatac$ref),
      sim = list("scRNA-seq" = example_scrna$sim, "scATAC-seq" = example_scatac$sim)
    ),
    "Experiment_B" = list(
      modality = "scRNA-seq",
      ref = example_scrna$ref,
      sim = example_scrna$sim
    )
  )

  res2 <- evaluate_multiple_datasets(
    grouped_datasets,
    pair_by_prefix = FALSE,
    compute_bivariate = FALSE,
    verbose = FALSE
  )

  expect_s3_class(res2, "scSimEval_consolidated")
  expect_true(all(c("Experiment_A", "Experiment_B") %in% res2$consolidated_summary_table$Dataset))

  p2 <- plot_consolidated_summary(res2, category = "(I) Distributional Properties", metric = "Wasserstein")
  expect_s3_class(p2, "ggplot")
})

