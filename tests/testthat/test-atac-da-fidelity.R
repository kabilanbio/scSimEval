test_that("binarize_matrix binarizes dense and sparse matrices", {
  # Dense matrix
  m_dense <- matrix(c(0, 2, 1, 0, 5, 0), nrow = 2, ncol = 3)
  b_dense <- binarize_matrix(m_dense)
  expect_true(is.matrix(b_dense))
  expect_equal(as.vector(b_dense), c(0, 1, 1, 0, 1, 0))
  
  # Sparse matrix
  m_sparse <- Matrix::Matrix(m_dense, sparse = TRUE)
  b_sparse <- binarize_matrix(m_sparse)
  expect_s4_class(b_sparse, "dgCMatrix")
  expect_equal(as.vector(as.matrix(b_sparse)), c(0, 1, 1, 0, 1, 0))
})

test_that("to_pseudobulk constructs pseudobulk counts per cell type and replicate", {
  set.seed(42)
  n_peaks <- 50
  n_cells <- 40
  counts <- matrix(rpois(n_peaks * n_cells, lambda = 2), nrow = n_peaks, ncol = n_cells)
  rownames(counts) <- paste0("Peak_", seq_len(n_peaks))
  colnames(counts) <- paste0("Cell_", seq_len(n_cells))
  
  cell_types <- factor(rep(c("CT_A", "CT_B"), each = 20))
  reps <- factor(rep(c("Rep1", "Rep2"), length.out = n_cells))
  
  # Missing replicates should stop
  expect_error(to_pseudobulk(counts, cell_types = cell_types, replicates = NULL))
  
  pb <- to_pseudobulk(
    input = counts,
    cell_types = cell_types,
    replicates = reps,
    binarize = TRUE,
    min_cells = 3,
    min_reps = 2
  )
  
  expect_type(pb, "list")
  expect_true(all(c("CT_A", "CT_B") %in% names(pb)))
  expect_equal(ncol(pb$CT_A), 2)
  expect_equal(ncol(pb$CT_B), 2)
  expect_equal(nrow(pb$CT_A), n_peaks)
})

test_that("run_edger_qlf_da detects differential accessibility between groups", {
  set.seed(42)
  n_peaks <- 80
  mat <- matrix(rpois(n_peaks * 4, lambda = 8), nrow = n_peaks, ncol = 4)
  rownames(mat) <- paste0("Peak_", seq_len(n_peaks))
  colnames(mat) <- c("Rep1:A", "Rep2:A", "Rep1:B", "Rep2:B")
  
  # Inject DA in first 10 peaks
  mat[1:10, 3:4] <- mat[1:10, 3:4] + rpois(10 * 2, lambda = 25)
  groups <- c("A", "A", "B", "B")
  
  res <- run_edger_qlf_da(mat, groups = groups, fdr_cutoff = 0.05, logfc_cutoff = 0)
  
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), n_peaks)
  expect_true(all(c("feature", "logFC", "logCPM", "PValue", "FDR", "is_da") %in% colnames(res)))
  expect_true(any(res$is_da))
})

test_that("evaluate_atac_da_fidelity strictly enforces biological replicates and evaluates DAR metrics", {
  data(example_scatac, package = "scSimEval")
  
  # Without biological replicates: should gracefully skip Category V
  res_skip <- evaluate_atac_da_fidelity(
    ref_data = example_scatac$ref,
    sim_data = example_scatac$sim,
    ref_celltypes = example_scatac$cell_types,
    sim_celltypes = example_scatac$cell_types,
    batch_ref = NULL,
    batch_sim = NULL
  )
  expect_equal(nrow(res_skip$da_summary_table), 0)
  
  # With biological replicates: computes the 4 Category V DAR metrics
  res_full <- evaluate_atac_da_fidelity(
    ref_data = example_scatac$ref,
    sim_data = example_scatac$sim,
    ref_celltypes = example_scatac$cell_types,
    sim_celltypes = example_scatac$cell_types,
    batch_ref = example_scatac$batch_info,
    batch_sim = example_scatac$batch_info
  )
  
  expect_type(res_full, "list")
  expect_true("da_summary_table" %in% names(res_full))
  expect_equal(nrow(res_full$da_summary_table), 6)
  
  metrics <- res_full$da_summary_table$Metric
  expected_6 <- c("pval_dist_score", "accuracy", "precision", "recall", "f1_score", "auc")
  expect_true(all(expected_6 %in% metrics))
  
  expect_true(all(c("SubCategory", "Framework", "Metric", "Value", "Direction") %in% colnames(res_full$da_summary_table)))
  expect_true(all(res_full$da_summary_table$SubCategory == "scATAC-seq Differential Accessibility"))
  expect_true(all(!is.na(res_full$da_summary_table$Value)))
  expect_true(all(res_full$da_summary_table$Value >= 0 & res_full$da_summary_table$Value <= 1))
})

