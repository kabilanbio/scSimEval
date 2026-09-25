test_that("evaluate_deg_fidelity computes multi-framework DEG metrics correctly", {
  set.seed(42)
  n_genes <- 150
  n_cells <- 80
  
  ref_mat <- matrix(rpois(n_genes * n_cells, lambda = 4), nrow = n_genes, ncol = n_cells)
  rownames(ref_mat) <- paste0("Gene_", seq_len(n_genes))
  colnames(ref_mat) <- paste0("Cell_", seq_len(n_cells))
  cell_types <- factor(rep(c("Type_A", "Type_B"), each = 40))
  ref_mat[1:25, 41:80] <- ref_mat[1:25, 41:80] + rpois(25 * 40, lambda = 6)
  
  sim_mat <- matrix(rpois(n_genes * n_cells, lambda = 4), nrow = n_genes, ncol = n_cells)
  rownames(sim_mat) <- rownames(ref_mat)
  colnames(sim_mat) <- colnames(ref_mat)
  sim_mat[1:20, 41:80] <- sim_mat[1:20, 41:80] + rpois(20 * 40, lambda = 5.5)
  
  res <- evaluate_deg_fidelity(
    ref_data = ref_mat,
    sim_data = sim_mat,
    ref_celltypes = cell_types,
    sim_celltypes = cell_types,
    fdr_cutoff = 0.05,
    logfc_cutoff = 0.5,
    top_n_de = 50,
    classifier = "knn",
    run_pvalue_uniformity = TRUE
  )
  
  expect_true(is.list(res))
  expect_true("deg_summary_table" %in% names(res))
  expect_equal(nrow(res$deg_summary_table), 15)
  expect_true(all(c("Framework", "Metric", "Value", "Direction") %in% colnames(res$deg_summary_table)))
  
  # Check key metric names
  metrics <- res$deg_summary_table$Metric
  expect_true("DEG_Ratio" %in% metrics)
  expect_true("SimBench_SMAPE" %in% metrics)
  expect_true("Log2FC_Pearson_Corr" %in% metrics)
  expect_true("Silhouette_Sim" %in% metrics)
  expect_true("Classifier_Accuracy" %in% metrics)
})
