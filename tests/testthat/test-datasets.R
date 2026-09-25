test_that("example_scrna dataset has valid structure and dimensions", {
  data("example_scrna", package = "scSimEval", envir = environment())
  expect_true(exists("example_scrna"))
  expect_type(example_scrna, "list")
  expect_named(example_scrna, c("ref", "sim", "cell_types", "batch_info"), ignore.order = TRUE)
  
  # Matrix dimensions
  expect_equal(dim(example_scrna$ref), c(60, 80))
  expect_equal(dim(example_scrna$sim), c(60, 80))
  
  # Cell annotations
  expect_length(example_scrna$cell_types, 80)
  expect_length(example_scrna$batch_info, 80)
  expect_true(is.factor(example_scrna$cell_types))
  expect_true(is.factor(example_scrna$batch_info))
  expect_equal(levels(example_scrna$cell_types), c("TypeA", "TypeB"))
  expect_equal(levels(example_scrna$batch_info), c("Batch1", "Batch2"))
})

test_that("example_scatac dataset has valid structure and dimensions", {
  data("example_scatac", package = "scSimEval", envir = environment())
  expect_true(exists("example_scatac"))
  expect_type(example_scatac, "list")
  expect_named(example_scatac, c("ref", "sim", "cell_types", "batch_info"), ignore.order = TRUE)
  
  # Matrix dimensions
  expect_equal(dim(example_scatac$ref), c(60, 80))
  expect_equal(dim(example_scatac$sim), c(60, 80))
  
  # Cell annotations
  expect_length(example_scatac$cell_types, 80)
  expect_length(example_scatac$batch_info, 80)
  expect_true(is.factor(example_scatac$cell_types))
  expect_true(is.factor(example_scatac$batch_info))
})

test_that("example_multiomics dataset has paired modalities and resource stats", {
  data("example_multiomics", package = "scSimEval", envir = environment())
  expect_true(exists("example_multiomics"))
  expect_type(example_multiomics, "list")
  expect_true(all(c("ref_multi", "sim_multi", "cell_types", "batch_info", "resource_stats") %in% names(example_multiomics)))
  
  expect_named(example_multiomics$ref_multi, c("rna", "atac"), ignore.order = TRUE)
  expect_named(example_multiomics$sim_multi, c("rna", "atac"), ignore.order = TRUE)
  
  # Cross-modality cell alignment
  expect_equal(ncol(example_multiomics$ref_multi$rna), ncol(example_multiomics$ref_multi$atac))
  expect_equal(colnames(example_multiomics$ref_multi$rna), colnames(example_multiomics$ref_multi$atac))
  
  # Resource stats
  expect_true(is.numeric(example_multiomics$resource_stats$cpu_time))
  expect_true(is.numeric(example_multiomics$resource_stats$memory_mb))
})
