test_that("evaluate_simulation_accuracy executes correctly on example_scrna", {
  data("example_scrna", package = "scSimEval", envir = environment())
  
  res <- evaluate_simulation_accuracy(
    ref_data = example_scrna$ref,
    sim_data = example_scrna$sim,
    compute_bivariate = FALSE,
    verbose = FALSE
  )
  
  expect_type(res, "list")
  expect_true(!is.null(res$metrics_summary_table))
  expect_s3_class(res$metrics_summary_table, "data.frame")
  expect_true(nrow(res$metrics_summary_table) > 0)
})

test_that("evaluate_clustering_metrics executes correctly with cell types", {
  data("example_scrna", package = "scSimEval", envir = environment())
  
  clust_res <- evaluate_clustering_metrics(
    data = example_scrna$sim,
    cell_types = example_scrna$cell_types,
    ref_data = example_scrna$ref
  )
  
  expect_type(clust_res, "list")
  expect_true("silhouette" %in% names(clust_res) || "silhouette_sim" %in% names(clust_res))
  expect_true("davies_bouldin" %in% names(clust_res) || "davies_bouldin_sim" %in% names(clust_res))
  expect_true("ARI" %in% names(clust_res) || "ari" %in% names(clust_res))
  expect_true("NMI" %in% names(clust_res) || "nmi" %in% names(clust_res))
})

test_that("evaluate_batch_metrics executes correctly with batch info", {
  data("example_scrna", package = "scSimEval", envir = environment())
  
  batch_res <- evaluate_batch_metrics(
    data = example_scrna$sim,
    batch_info = example_scrna$batch_info,
    cell_types = example_scrna$cell_types
  )
  
  expect_type(batch_res, "list")
  expect_true("shannon_entropy" %in% names(batch_res))
  expect_true("pcr_r2" %in% names(batch_res))
})

test_that("trajectory pseudotime is inferred correctly from scRNA-seq", {
  data("example_scrna", package = "scSimEval", envir = environment())
  
  pt <- infer_scrna_pseudotime(example_scrna$ref)
  expect_type(pt, "double")
  expect_length(pt, ncol(example_scrna$ref))
})

test_that("evaluate_multiomics_accuracy executes end-to-end on example_multiomics", {
  data("example_multiomics", package = "scSimEval", envir = environment())
  
  master_res <- evaluate_multiomics_accuracy(
    ref_multi = example_multiomics$ref_multi,
    sim_multi = example_multiomics$sim_multi,
    cell_types = example_multiomics$cell_types,
    batch_info = example_multiomics$batch_info,
    memory_mb = example_multiomics$resource_stats$memory_mb,
    elapsed_time = example_multiomics$resource_stats$elapsed_time,
    verbose = FALSE
  )
  
  expect_type(master_res, "list")
  expect_true("benchmark_summary_table" %in% names(master_res))
  expect_s3_class(master_res$benchmark_summary_table, "data.frame")
  expect_true(nrow(master_res$benchmark_summary_table) > 100)
  
  # Check representation across categories
  cats <- unique(master_res$benchmark_summary_table$Category)
  expect_true("Distributional Properties" %in% cats)
  expect_true("Cellular Structure & Mixing" %in% cats)
  expect_true("Cross-Modal Relationships" %in% cats)
  expect_true("Computational Scalability" %in% cats)
  
  # Check batch mixing property
  props <- unique(master_res$benchmark_summary_table$Property)
  expect_true("batch_mixing" %in% props)
})
