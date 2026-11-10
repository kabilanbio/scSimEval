test_that("run_edger_qlf_detrate identifies DEGs with cellular detection rate covariate", {
  set.seed(42)
  n_genes <- 100
  n_cells <- 60
  counts <- matrix(rpois(n_genes * n_cells, lambda = 3), nrow = n_genes, ncol = n_cells)
  rownames(counts) <- paste0("Gene_", seq_len(n_genes))
  colnames(counts) <- paste0("Cell_", seq_len(n_cells))
  groups <- factor(rep(c("Group_A", "Group_B"), each = 30))
  
  # Inject clear DE signal in first 15 genes
  counts[1:15, 31:60] <- counts[1:15, 31:60] + rpois(15 * 30, lambda = 8)
  
  res <- run_edger_qlf_detrate(
    counts = counts,
    group = groups,
    fdr_cutoff = 0.05,
    logfc_cutoff = 0.5
  )
  
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), n_genes)
  expect_true(all(c("gene", "logFC", "logCPM", "PValue", "FDR", "is_deg") %in% colnames(res)))
  
  # Top injected genes should be detected as DEGs
  sig_genes <- res$gene[res$is_deg]
  expect_true(length(sig_genes) >= 10)
  expect_true(any(paste0("Gene_", 1:15) %in% sig_genes))
})

test_that("calc_simpipe_deg_ratio computes Approach 1 true DEG recovery P and DEG ratio", {
  set.seed(42)
  n_genes <- 120
  n_cells <- 80
  
  ref_mat <- matrix(rpois(n_genes * n_cells, lambda = 3), nrow = n_genes, ncol = n_cells)
  rownames(ref_mat) <- paste0("Gene_", seq_len(n_genes))
  colnames(ref_mat) <- paste0("Cell_", seq_len(n_cells))
  cell_types_2g <- factor(rep(c("Type1", "Type2"), each = 40))
  ref_mat[1:20, 41:80] <- ref_mat[1:20, 41:80] + rpois(20 * 40, lambda = 7)
  
  sim_mat <- matrix(rpois(n_genes * n_cells, lambda = 3), nrow = n_genes, ncol = n_cells)
  rownames(sim_mat) <- rownames(ref_mat)
  colnames(sim_mat) <- colnames(ref_mat)
  sim_mat[1:18, 41:80] <- sim_mat[1:18, 41:80] + rpois(18 * 40, lambda = 6.5)
  
  # Test n = 2 groups
  res_2g <- calc_simpipe_deg_ratio(
    ref_data = ref_mat,
    sim_data = sim_mat,
    ref_celltypes = cell_types_2g,
    sim_celltypes = cell_types_2g,
    fdr_cutoff = 0.05,
    logfc_cutoff = 0.5
  )
  
  expect_type(res_2g, "list")
  expect_true(is.numeric(res_2g$P) && res_2g$P >= 0 && res_2g$P <= 1)
  expect_true(is.numeric(res_2g$deg_ratio) && res_2g$deg_ratio >= 0)
  expect_equal(nrow(res_2g$pairwise_results), 1)
  expect_equal(res_2g$pairwise_results$Weight[1], 1.0)
  
  # Test n > 2 groups (3 cell types -> 3 pairwise contrasts weighted)
  cell_types_3g <- factor(rep(c("TypeA", "TypeB", "TypeC"), length.out = n_cells))
  res_3g <- calc_simpipe_deg_ratio(
    ref_data = ref_mat,
    sim_data = sim_mat,
    ref_celltypes = cell_types_3g,
    sim_celltypes = cell_types_3g,
    fdr_cutoff = 0.05,
    logfc_cutoff = 0.5
  )
  
  expect_type(res_3g, "list")
  expect_equal(nrow(res_3g$pairwise_results), 3)
  expect_equal(sum(res_3g$pairwise_results$Weight), 1.0, tolerance = 1e-6)
  expect_true(is.numeric(res_3g$P) && res_3g$P >= 0 && res_3g$P <= 1)
})

test_that("calc_simpipe_distribution_score computes Approach 2 p-value uniformity score", {
  set.seed(42)
  n_genes <- 120
  n_cells <- 60
  
  sim_mat <- matrix(rpois(n_genes * n_cells, lambda = 4), nrow = n_genes, ncol = n_cells)
  rownames(sim_mat) <- paste0("Gene_", seq_len(n_genes))
  colnames(sim_mat) <- paste0("Cell_", seq_len(n_cells))
  cell_types <- factor(rep(c("TypeA", "TypeB"), each = 30))
  sim_mat[1:15, 31:60] <- sim_mat[1:15, 31:60] + rpois(15 * 30, lambda = 6)
  
  res_dist <- calc_simpipe_distribution_score(
    sim_data = sim_mat,
    sim_celltypes = cell_types,
    n_bins = 10,
    alpha = 0.05
  )
  
  expect_type(res_dist, "list")
  expect_true(!is.na(res_dist$distribution_score))
  expect_true(res_dist$distribution_score %in% c(0, 1))
  expect_equal(nrow(res_dist$pairwise_results), 1)
  expect_true(all(c("Chisq_Stat", "Chisq_Pval", "Accepted_Null", "S_i") %in% colnames(res_dist$pairwise_results)))
})

test_that("evaluate_deg_fidelity strictly requires biological replicates and evaluates pseudobulk DE fidelity", {
  set.seed(42)
  n_genes <- 100
  n_cells <- 60
  
  ref_mat <- matrix(rpois(n_genes * n_cells, lambda = 4), nrow = n_genes, ncol = n_cells)
  rownames(ref_mat) <- paste0("Gene_", seq_len(n_genes))
  colnames(ref_mat) <- paste0("Cell_", seq_len(n_cells))
  cell_types <- factor(rep(c("Type_A", "Type_B"), each = 30))
  batch_labels <- factor(rep(c("Rep1", "Rep2"), length.out = n_cells))
  ref_mat[1:20, 31:60] <- ref_mat[1:20, 31:60] + rpois(20 * 30, lambda = 6)
  
  sim_mat <- matrix(rpois(n_genes * n_cells, lambda = 4), nrow = n_genes, ncol = n_cells)
  rownames(sim_mat) <- rownames(ref_mat)
  colnames(sim_mat) <- colnames(ref_mat)
  sim_mat[1:18, 31:60] <- sim_mat[1:18, 31:60] + rpois(18 * 30, lambda = 5.5)
  
  # Without biological replicates: gracefully skipped
  res_skip <- evaluate_deg_fidelity(
    ref_data = ref_mat,
    sim_data = sim_mat,
    ref_celltypes = cell_types,
    sim_celltypes = cell_types,
    batch_ref = NULL,
    batch_sim = NULL
  )
  expect_equal(nrow(res_skip$deg_summary_table), 0)
  
  # With biological replicates: computes pseudobulk DE fidelity
  res <- evaluate_deg_fidelity(
    ref_data = ref_mat,
    sim_data = sim_mat,
    ref_celltypes = cell_types,
    sim_celltypes = cell_types,
    batch_ref = batch_labels,
    batch_sim = batch_labels,
    fdr_cutoff = 0.05,
    logfc_cutoff = 0,
    classifier = "knn"
  )
  
  expect_type(res, "list")
  expect_true("deg_summary_table" %in% names(res))
  expect_s3_class(res$deg_summary_table, "data.frame")
  expect_equal(nrow(res$deg_summary_table), 6)
  
  metrics <- res$deg_summary_table$Metric
  expected_6 <- c("pval_dist_score", "accuracy", "precision", "recall", "f1_score", "auc")
  expect_true(all(expected_6 %in% metrics))
  
  # Check directions and subcategory
  expect_true(all(c("SubCategory", "Framework", "Metric", "Value", "Direction") %in% colnames(res$deg_summary_table)))
  expect_true(all(res$deg_summary_table$SubCategory == "scRNA-seq Differential Expression"))
  expect_true(all(!is.na(res$deg_summary_table$Value)))
  expect_true(all(res$deg_summary_table$Value >= 0 & res$deg_summary_table$Value <= 1))
})
