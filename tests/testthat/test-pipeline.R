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


test_that("evaluate_simulation_accuracy strictly requires Seurat or SingleCellExperiment", {
  raw_mat1 <- matrix(rpois(100, 2), 10, 10)
  raw_mat2 <- matrix(rpois(100, 2), 10, 10)
  
  # Plain matrices must throw informative error
  expect_error(
    evaluate_simulation_accuracy(raw_mat1, raw_mat2, verbose = FALSE),
    "scSimEval strictly requires a 'Seurat' or 'SingleCellExperiment' object"
  )
})

test_that("scalability inputs support minutes, MB vs MiB conversion, and sim_cores", {
  data("example_scrna", package = "scSimEval", envir = environment())

  res <- evaluate_simulation_accuracy(
    ref_data = example_scrna$ref,
    sim_data = example_scrna$sim,
    compute_bivariate = FALSE,
    memory_mb = 500,
    elapsed_time = 2.5,
    time_unit = "minutes",
    memory_unit = "MB",
    sim_cores = 8L,
    verbose = FALSE
  )

  tbl <- res$metrics_summary_table
  t_row <- tbl[tbl$Metric == "core_seconds", ]
  m_row <- tbl[tbl$Metric == "memory_per_core_mib", ]

  expect_equal(nrow(t_row), 1)
  expect_equal(t_row$Value, 150 * 8, tolerance = 1e-6)
  expect_equal(nrow(m_row), 1)
  expect_equal(m_row$Value, (500 * (1000^2 / 1024^2)) / 8, tolerance = 1e-4)
  expect_true("Sim_Cores" %in% colnames(tbl))
  expect_equal(unique(tbl$Sim_Cores), 8L)
  expect_equal(res$resource_usage$time_unit, "minutes")
  expect_equal(res$resource_usage$memory_unit, "MB")
  expect_equal(res$resource_usage$sim_cores, 8L)
  expect_equal(res$resource_usage$core_seconds, 1200, tolerance = 1e-6)

  data("example_multiomics", package = "scSimEval", envir = environment())
  res_m <- evaluate_multiomics_accuracy(
    ref_multi = example_multiomics$ref_multi,
    sim_multi = example_multiomics$sim_multi,
    cell_types = example_multiomics$cell_types,
    batch_info = example_multiomics$batch_info,
    memory_mb = 1024,
    elapsed_time = 1.5,
    time_unit = "min",
    memory_unit = "MB",
    sim_cores = 4L,
    verbose = FALSE
  )
  mt <- res_m$benchmark_summary_table
  expect_equal(mt$Value[mt$Metric == "core_seconds"], 90 * 4, tolerance = 1e-6)
  expect_equal(mt$Value[mt$Metric == "memory_per_core_mib"], (1024 * (1000^2 / 1024^2)) / 4, tolerance = 1e-4)
  expect_equal(unique(mt$Sim_Cores), 4L)
})

test_that("multiomics supports both combined and separate RNA/ATAC scalability in paired and unpaired modes", {
  data("example_multiomics", package = "scSimEval", envir = environment())

  # 1. Paired mode with Separate RNA & ATAC scalability inputs
  res_sep <- evaluate_multiomics_accuracy(
    ref_multi = example_multiomics$ref_multi,
    sim_multi = example_multiomics$sim_multi,
    cell_types = example_multiomics$cell_types,
    batch_info = example_multiomics$batch_info,
    pairing = "paired",
    elapsed_time_rna = 2.0,
    elapsed_time_atac = 30.0,
    memory_mb_rna = 1000,
    memory_mb_atac = 1200,
    time_unit_rna = "minutes",
    time_unit_atac = "seconds",
    memory_unit_rna = "MB",
    memory_unit_atac = "MiB",
    sim_cores_rna = 4L,
    sim_cores_atac = 8L,
    verbose = FALSE
  )

  mt_sep <- res_sep$benchmark_summary_table
  t_rows <- mt_sep[mt_sep$Metric == "core_seconds", ]
  m_rows <- mt_sep[mt_sep$Metric == "memory_per_core_mib", ]
  expect_equal(nrow(t_rows), 1L)
  expect_equal(nrow(m_rows), 1L)
  expect_equal(t_rows$Modality, "Computation")
  expect_equal(m_rows$Modality, "Computation")
  expect_equal(t_rows$Value, 120 * 4 + 30 * 8, tolerance = 1e-6)
  expect_equal(m_rows$Value, (1000 * (1000^2 / 1024^2)) / 4, tolerance = 1e-4)
  expect_equal(unique(mt_sep$Sim_Cores), 8L)
  expect_equal(res_sep$resource_usage$mode, "separate")
  expect_equal(res_sep$resource_usage$rna$elapsed_time, 120, tolerance = 1e-6)
  expect_equal(res_sep$resource_usage$atac$elapsed_time, 30, tolerance = 1e-6)
  expect_equal(res_sep$resource_usage$rna$memory_mb, 1000 * (1000^2 / 1024^2), tolerance = 1e-4)
  expect_equal(res_sep$resource_usage$atac$memory_mb, 1200, tolerance = 1e-6)

  # 2. Unpaired mode with Combined (RNA + ATAC together) scalability inputs
  res_unp <- evaluate_multiomics_accuracy(
    ref_multi = example_multiomics$ref_multi,
    sim_multi = example_multiomics$sim_multi,
    cell_types = example_multiomics$cell_types,
    batch_info = example_multiomics$batch_info,
    pairing = "unpaired",
    elapsed_time = 2.5,
    memory_mb = 2048,
    time_unit = "minutes",
    memory_unit = "MiB",
    sim_cores = 6L,
    verbose = FALSE
  )

  mt_unp <- res_unp$benchmark_summary_table
  expect_equal(mt_unp$Value[mt_unp$Metric == "core_seconds"], 150 * 6, tolerance = 1e-6)
  expect_equal(mt_unp$Value[mt_unp$Metric == "memory_per_core_mib"], 2048 / 6, tolerance = 1e-6)
  expect_equal(unique(mt_unp$Sim_Cores), 6L)
  expect_equal(res_unp$resource_usage$mode, "combined")

  # 3. Layer-specific scalability plotting works from resource_usage[[modality]]
  p_rna <- plot_scalability_benchmark(
    list(PairedSep = res_sep, UnpairedComb = res_unp),
    modality = "rna"
  )
  expect_true(inherits(p_rna, "gg"))
})

test_that("Category VII Saelens et al. (2019) HVG preprocessing and dynverse metrics execute cleanly", {
  data("example_scrna", package = "scSimEval", envir = environment())

  # 1. HVG preprocessing pipeline (3-MAD QC, low-expression filter, scran HVG selection, convergence)
  prep <- preprocess_trajectory_hvgs(
    ref_data = example_scrna$ref,
    sim_data = example_scrna$sim,
    cell_types_ref = example_scrna$cell_types,
    cell_types_sim = example_scrna$cell_types
  )
  expect_type(prep, "list")
  expect_true(length(prep$hvgs) >= 20L)
  expect_equal(nrow(prep$ref_expr), length(prep$hvgs))
  expect_equal(nrow(prep$sim_expr), length(prep$hvgs))
  expect_true(all(Matrix::colSums(prep$ref_expr) > 0))
  expect_true(all(Matrix::colSums(prep$sim_expr) > 0))

  # 2. Synthetic branching 3-lineage trajectory self-comparison (identical ref vs ref)
  set.seed(7)
  n_g <- 120L
  n_c <- 180L
  t_seq <- seq(0, 1, length.out = 60L)
  branch_lbl <- rep(c("Stem", "Lineage1", "Lineage2"), each = 60L)
  base_mat <- matrix(stats::rpois(n_g * n_c, lambda = 2), nrow = n_g, ncol = n_c)
  rownames(base_mat) <- paste0("Gene_", seq_len(n_g))
  colnames(base_mat) <- paste0("Cell_", seq_len(n_c))
  # Program strong branch-specific HVGs
  for (g in 1:15) base_mat[g, 61:120] <- stats::rpois(60L, lambda = 18 * t_seq + 2)
  for (g in 16:30) base_mat[g, 121:180] <- stats::rpois(60L, lambda = 18 * t_seq + 2)
  for (g in 31:45) base_mat[g, 1:60] <- stats::rpois(60L, lambda = 12 * (1 - t_seq) + 2)

  self_res <- evaluate_trajectory_metrics(
    ref_data = base_mat,
    sim_data = base_mat,
    cell_types_ref = branch_lbl,
    cell_types_sim = branch_lbl,
    num_trees = 500L
  )
  expect_setequal(
    c("him", "edge_flip", "geodesic_correlation", "F1_branches", "wcor_features"),
    names(self_res)[seq_len(5)]
  )
  expect_equal(self_res$him, 1, tolerance = 1e-6)
  expect_equal(self_res$edge_flip, 1, tolerance = 1e-6)
  expect_gt(self_res$geodesic_correlation, 0.95)
  expect_gt(self_res$F1_branches, 0.95)
  expect_gt(self_res$wcor_features, 0.85)

  # 3. Cross-dataset unsupervised trajectory evaluation (no cell types supplied)
  unsup_res <- evaluate_trajectory_metrics(
    ref_data = example_scrna$ref,
    sim_data = example_scrna$sim,
    num_trees = 500L
  )
  for (m in c("him", "edge_flip", "geodesic_correlation", "F1_branches", "wcor_features")) {
    expect_true(is.finite(unsup_res[[m]]))
    expect_gte(unsup_res[[m]], 0)
    expect_lte(unsup_res[[m]], 1)
  }

  # 4. Gold-standard filter preserves all 5 Category VII dynverse metrics
  eval_uni <- evaluate_simulation_accuracy(
    ref_data = example_scrna$ref,
    sim_data = example_scrna$sim,
    compute_bivariate = FALSE,
    trajectory_num_trees = 500L,
    verbose = FALSE
  )
  gs_tbl <- filter_gold_standard_metrics(eval_uni$metrics_summary_table)
  cat7_tbl <- gs_tbl[gs_tbl$Category == "(VII) Trajectory & Lineage Topology", ]
  expect_setequal(
    cat7_tbl$Metric,
    c("him", "edge_flip", "geodesic_correlation", "F1_branches", "wcor_features")
  )

  # 5. Category VII requires scRNA-seq only and is skipped for scATAC-seq
  data("example_scatac", package = "scSimEval", envir = environment())
  eval_atac <- evaluate_simulation_accuracy(
    ref_data = example_scatac$ref,
    sim_data = example_scatac$sim,
    modality = "scATAC-seq",
    compute_bivariate = FALSE,
    trajectory_num_trees = 500L,
    verbose = FALSE
  )
  atac_cat7 <- eval_atac$metrics_summary_table[
    eval_atac$metrics_summary_table$Category == "(VII) Trajectory & Lineage Topology",
  ]
  expect_equal(nrow(atac_cat7), 0L)

  # 6. Multiomics evaluates Category VII once on scRNA-seq only (never ATAC, never duplicated)
  data("example_multiomics", package = "scSimEval", envir = environment())
  eval_mo <- evaluate_multiomics_accuracy(
    ref_multi = example_multiomics$ref_multi,
    sim_multi = example_multiomics$sim_multi,
    cell_types = example_multiomics$cell_types,
    compute_bivariate = FALSE,
    trajectory_num_trees = 500L,
    verbose = FALSE
  )
  mo_cat7 <- eval_mo$benchmark_summary_table[
    eval_mo$benchmark_summary_table$Category == "(VII) Trajectory & Lineage Topology",
  ]
  expect_equal(nrow(mo_cat7), 5L)
  expect_equal(unique(mo_cat7$Modality), "rna")
  expect_setequal(
    mo_cat7$Metric,
    c("him", "edge_flip", "geodesic_correlation", "F1_branches", "wcor_features")
  )
})
