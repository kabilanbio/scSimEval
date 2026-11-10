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
  expect_true("neighborhood_purity" %in% names(clust_res))
  expect_true("ARI" %in% names(clust_res) || "ari" %in% names(clust_res))
  expect_true("NMI" %in% names(clust_res) || "nmi" %in% names(clust_res))
  expect_true("cdi" %in% names(clust_res) && is.finite(clust_res$cdi))
  expect_true("CDI_AIC" %in% names(clust_res) && is.finite(clust_res$CDI_AIC))
  expect_true("CDI_BIC" %in% names(clust_res) && is.finite(clust_res$CDI_BIC))
})

test_that("calc_cdi uses Bioconductor CDI accurately on sparse and stratified matrices", {
  set.seed(123)
  # Cluster 1 expresses genes 1:20 high / 21:40 low; Cluster 2 expresses 1:20 low / 21:40 high
  mat_a <- rbind(
    matrix(stats::rnbinom(20 * 60, size = 5, mu = 10), nrow = 20, ncol = 60),
    matrix(stats::rnbinom(20 * 60, size = 5, mu = 1),  nrow = 20, ncol = 60)
  )
  mat_b <- rbind(
    matrix(stats::rnbinom(20 * 60, size = 5, mu = 1),  nrow = 20, ncol = 60),
    matrix(stats::rnbinom(20 * 60, size = 5, mu = 10), nrow = 20, ncol = 60)
  )
  mat_full <- methods::as(cbind(mat_a, mat_b), "dgCMatrix")
  true_cl <- factor(rep(c("Cluster1", "Cluster2"), each = 60))
  rand_cl <- sample(true_cl)

  cdi_true <- calc_cdi(mat_full, true_cl, max_cells = 120L)
  cdi_rand <- calc_cdi(mat_full, rand_cl, max_cells = 120L)

  expect_true(is.finite(cdi_true$cdi_aic) && is.finite(cdi_true$cdi_bic))
  expect_lt(cdi_true$cdi_aic, cdi_rand$cdi_aic)
  expect_lt(cdi_true$cdi_bic, cdi_rand$cdi_bic)

  # Also verify stratified subsampling executes cleanly
  cdi_sub <- calc_cdi(mat_full, true_cl, max_cells = 80L)
  expect_true(is.finite(cdi_sub$cdi_aic) && cdi_sub$n_cells <= 80L)
})

test_that("evaluate_batch_metrics executes correctly with batch info", {
  data("example_scrna", package = "scSimEval", envir = environment())
  
  batch_res <- evaluate_batch_metrics(
    data = example_scrna$sim,
    batch_info = example_scrna$batch_info,
    cell_types = example_scrna$cell_types
  )
  
  expect_type(batch_res, "list")
  expect_setequal(names(batch_res), c("batch_silhouette", "pcr_r2", "cms", "isi", "cross_batch_F1"))
  expect_true(is.finite(batch_res$pcr_r2))
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
  expect_true("Cellular Manifold & Clustering Concordance" %in% cats)
  # expect_true("Cross-Modal Relationships" %in% cats) # Cleared: pending new metrics
  expect_true("(VIII) Computational Scalability & Efficiency" %in% cats)
  
  # Check batch mixing property
  props <- unique(master_res$benchmark_summary_table$Property)
  expect_true("batch_mixing" %in% props)
})


test_that("evaluate_simulation_accuracy enforces Option B (strictly requires Seurat or SingleCellExperiment)", {
  raw_mat1 <- matrix(rpois(100, 2), 10, 10)
  raw_mat2 <- matrix(rpois(100, 2), 10, 10)
  
  # Plain matrices must throw informative error
  expect_error(
    evaluate_simulation_accuracy(raw_mat1, raw_mat2, verbose = FALSE),
    "scSimEval strictly requires a 'Seurat' or 'SingleCellExperiment' object"
  )
})
