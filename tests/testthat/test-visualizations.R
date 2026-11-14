test_that("plot_distribution_qc returns valid ggplot with single and multiple simulators", {
  data("example_scrna", package = "scSimEval", envir = environment())
  
  p1 <- plot_distribution_qc(example_scrna$ref, example_scrna$sim)
  expect_s3_class(p1, "ggplot")

  sim_list <- list(
    "SimA" = example_scrna$sim,
    "SimB" = example_scrna$sim
  )
  p2 <- plot_distribution_qc(example_scrna$ref, sim_list)
  expect_s3_class(p2, "ggplot")

  p_qq <- plot_distribution_qc(example_scrna$ref, example_scrna$sim, layout = "qq")
  expect_s3_class(p_qq, "ggplot")

  p_qq_multi <- plot_distribution_qc(example_scrna$ref, sim_list, layout = "qq")
  expect_s3_class(p_qq_multi, "ggplot")

  p_qq_cell <- plot_distribution_qc(example_scrna$ref, example_scrna$sim, layout = "qq", properties = "cell")
  expect_s3_class(p_qq_cell, "ggplot")

  p_qq_feat <- plot_distribution_qc(example_scrna$ref, example_scrna$sim, layout = "qq", properties = "feature")
  expect_s3_class(p_qq_feat, "ggplot")

  p_qq_biv <- plot_distribution_qc(example_scrna$ref, example_scrna$sim, layout = "qq", properties = "bivariate")
  expect_s3_class(p_qq_biv, "ggplot")

  p_qq_single_prop <- plot_distribution_qc(example_scrna$ref, example_scrna$sim, layout = "qq", properties = "mean_vs_dispersion")
  expect_s3_class(p_qq_single_prop, "ggplot")

  p_dens_multi_prop <- plot_distribution_qc(example_scrna$ref, example_scrna$sim, layout = "density", properties = "mean_vs_dispersion")
  expect_s3_class(p_dens_multi_prop, "ggplot")

  p_dens <- plot_distribution_qc(example_scrna$ref, example_scrna$sim, layout = "density")
  expect_s3_class(p_dens, "ggplot")

  p_dens_cell <- plot_distribution_qc(example_scrna$ref, sim_list, layout = "density", properties = "cell")
  expect_s3_class(p_dens_cell, "ggplot")
})

test_that("plot_metric_boxplots returns a valid ggplot object", {
  demo_df <- data.frame(
    Method = rep(c("Splatter", "scDesign3", "SymSim"), each = 4),
    Category = rep(c("(I) Distributional Properties", "(III) Cellular Structure & Concordance"), 6),
    Metric = rep(c("KS", "Wasserstein", "ARI", "NMI"), 3),
    Score = runif(12, 0.4, 0.95)
  )
  p <- plot_metric_boxplots(demo_df)
  expect_s3_class(p, "ggplot")
})

test_that("plot_category_metric_bars returns a valid ggplot object", {
  demo_df <- data.frame(
    Method = rep(c("Splatter", "scDesign3", "SymSim"), each = 2),
    Category = "(I) Distributional Properties",
    Metric = rep(c("KS", "Wasserstein"), 3),
    Score = runif(6, 0.4, 0.95)
  )
  p_norm <- plot_category_metric_bars(demo_df, category = "(I) Distributional Properties", score_type = "normalized")
  expect_s3_class(p_norm, "ggplot")

  p_raw <- plot_category_metric_bars(demo_df, category = "(I) Distributional Properties", score_type = "raw")
  expect_s3_class(p_raw, "ggplot")
})

test_that("compute_method_leaderboard computes valid ranking table", {
  demo_df <- data.frame(
    Method = rep(c("Splatter", "scDesign3", "SymSim"), each = 2),
    Category = "(I) Distributional Properties",
    Metric = rep(c("KS", "Wasserstein"), 3),
    Score = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
  )
  lb <- compute_method_leaderboard(demo_df)
  expect_s3_class(lb, "data.frame")
  expect_named(lb, c("Overall_Rank", "Method", "Average_Fidelity", "Fidelity_Score"))
  expect_equal(nrow(lb), 3)
  expect_equal(lb$Overall_Rank, c(1, 2, 3))
  expect_true(all(lb$Fidelity_Score >= 0 & lb$Fidelity_Score <= 1))
})

test_that("plot_scalability_benchmark returns valid ggplot and patchwork objects", {
  demo_sc <- data.frame(
    Method = rep(c("Splatter", "scDesign3", "SymSim"), each = 2),
    Category = "(VIII) Computational Scalability & Efficiency",
    Metric = rep(c("elapsed_time_seconds", "peak_memory_mb"), 3),
    Score = c(30, 700, 65, 1300, 45, 800),
    Sim_Cores = rep(c(1L, 4L, 8L), each = 2)
  )
  p_run <- plot_scalability_benchmark(demo_sc, type = "runtime")
  expect_s3_class(p_run, "ggplot")

  p_run_min <- plot_scalability_benchmark(demo_sc, type = "runtime", time_unit = "minutes")
  expect_s3_class(p_run_min, "ggplot")
  expect_true(grepl("Minutes", p_run_min$labels$x, fixed = TRUE))

  p_mem <- plot_scalability_benchmark(demo_sc, type = "memory")
  expect_s3_class(p_mem, "ggplot")

  p_mem_mb <- plot_scalability_benchmark(demo_sc, type = "memory", memory_unit = "MB")
  expect_s3_class(p_mem_mb, "ggplot")
  expect_true(grepl("Peak RAM.*MB", p_mem_mb$labels$x))

  # Aliased inputs in minutes and decimal MB ingest cleanly
  alias_sc <- data.frame(
    Method = rep(c("SimA", "SimB"), each = 2),
    Category = "(VIII) Computational Scalability & Efficiency",
    Metric = rep(c("Runtime_Minutes", "PeakRAM_MB"), 2),
    Score = c(2, 1000, 3, 2000),
    Sim_Cores = rep(c(2L, 4L), each = 2)
  )
  filtered_alias <- filter_gold_standard_metrics(alias_sc)
  expect_equal(filtered_alias$Score[filtered_alias$Metric == "elapsed_time_seconds"], c(120, 180))
  expect_equal(filtered_alias$Score[filtered_alias$Metric == "peak_memory_mb"],
               c(1000, 2000) * (1000^2 / 1024^2), tolerance = 1e-4)

  p_comp <- plot_scalability_benchmark(demo_sc, type = "composite", layout = "4panel",
                                       time_unit = "minutes", memory_unit = "MB")
  expect_true(inherits(p_comp, "ggplot") || inherits(p_comp, "patchwork"))
})

test_that("plot_metric_mds and plot_metric_pca return valid ggplot objects", {
  demo_df <- data.frame(
    Method = rep(c("Splatter", "scDesign3", "SymSim", "dyngen"), each = 3),
    Category = "(I) Distributional Properties",
    Metric = rep(c("KS", "Wasserstein", "MAD"), 4),
    Score = runif(12, 0.5, 0.95)
  )
  p_mds <- plot_metric_mds(demo_df)
  expect_s3_class(p_mds, "ggplot")

  p_pca <- plot_metric_pca(demo_df)
  expect_s3_class(p_pca, "ggplot")
})

test_that("plot_benchmark_summary returns a valid ggplot object", {
  data("example_multiomics", package = "scSimEval", envir = environment())
  
  res <- evaluate_multiomics_accuracy(
    ref_multi = example_multiomics$ref_multi,
    sim_multi = example_multiomics$sim_multi,
    verbose = FALSE
  )
  
  p <- plot_benchmark_summary(res$benchmark_summary_table, metric_type = "KS")
  expect_s3_class(p, "ggplot")
})

test_that("plot_metric_heatmap returns a valid ggplot object", {
  data("example_multiomics", package = "scSimEval", envir = environment())
  
  res <- evaluate_multiomics_accuracy(
    ref_multi = example_multiomics$ref_multi,
    sim_multi = example_multiomics$sim_multi,
    verbose = FALSE
  )
  
  sim_list <- list(
    Sim_Tool_A = res$benchmark_summary_table,
    Sim_Tool_B = res$benchmark_summary_table
  )
  
  p <- plot_metric_heatmap(sim_list, metric_name = "KS", top_n_properties = 8)
  expect_s3_class(p, "ggplot")
})

test_that("plot_benchmark_bubble_matrix returns a valid ggplot object with data.frame and list inputs", {
  demo_df <- data.frame(
    Method = rep(c("Splat", "scDesign3", "SymSim", "dyngen"), each = 4),
    Category = rep(c("Accuracy", "Accuracy", "Cellular Structure", "Scalability"), 4),
    Metric = rep(c("KS Distance", "Wasserstein", "Silhouette ASW", "CPU Time"), 4),
    Score = c(0.8, 0.7, 0.9, 0.6, 0.5, 0.8, 0.6, 0.5, 0.8, 0.7, 0.6, 0.5, 0.7, 0.6, 0.9, 0.5)
  )
  
  classes <- list(
    "Class 1" = c("Splat", "scDesign3"),
    "Class 2" = c("SymSim", "dyngen")
  )
  
  p1 <- plot_benchmark_bubble_matrix(demo_df, method_classes = classes)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(ggplot2::ggplotGrob(p1), "gtable")
  expect_equal(p1$theme$legend.box, "horizontal")

  p1_nobars <- plot_benchmark_bubble_matrix(demo_df, method_classes = classes, show_category_bars = FALSE)
  expect_s3_class(p1_nobars, "ggplot")
  expect_s3_class(ggplot2::ggplotGrob(p1_nobars), "gtable")
  expect_true(length(p1$layers) > length(p1_nobars$layers))
  
  p2 <- plot_bubble_matrix(demo_df)
  expect_s3_class(p2, "ggplot")

  bm_mod <- data.frame(
    Method = rep(c("SimA", "SimB"), each = 5),
    Category = rep(c(
      "(I) Distributional Properties",
      "(II) Sparsity & Zero-Inflation",
      "(III) Cellular Structure & Concordance",
      "(VII) Trajectory & Lineage Topology",
      "(VI) Multiomics Integration & Alignment"
    ), 2),
    Metric = rep(c("KS", "curve_rmse", "ari", "geodesic_correlation", "FOSCTTM"), 2),
    Modality = rep(c("RNA", "ATAC", "RNA", "RNA", "Joint (Cross-Modal)"), 2),
    Score = c(0.1, 0.2, 0.8, 0.5, 0.15, 0.3, 0.1, 0.6, 0.4, 0.25)
  )
  p_mod <- plot_benchmark_bubble_matrix(bm_mod)
  expect_s3_class(p_mod, "ggplot")
  expect_false(any(grepl("\\[(scRNA|scATAC|Joint)\\]", levels(p_mod$data$Display_Metric))))
  expect_false("Modality Scope" %in% levels(p_mod$data$Method))
  expect_true(all(c("scRNA-seq only", "scATAC-seq only", "Multiomics (both scRNA-seq + scATAC-seq)") %in% as.character(p_mod$data$Modality_Scope)))
  expect_s3_class(ggplot2::ggplotGrob(p_mod), "gtable")

  bm_filt_in <- data.frame(
    Method = "SimA",
    Category = c(
      "(III) Cellular Structure & Concordance",
      "(V) Biological Signal & Downstream Fidelity",
      "(VII) Cross-Modal Coupling & Modularity",
      "(VII) Cross-Modal Coupling & Modularity",
      "(VII) Cross-Modal Coupling & Modularity"
    ),
    Property = c("clustering", "simbench_signal", "cell_alignment", "chromatin_coaccessibility", "coexpression_modularity"),
    Metric = c("davies_bouldin_sim", "DE (Mean Shift)_abs_error", "FOSCTTM", "rv_coefficient", "module_correlation_r"),
    Modality = c("RNA", "RNA", "Joint (Cross-Modal)", "ATAC", "RNA"),
    Score = c(0.6, 0.7, 0.2, 0.8, 0.9),
    stringsAsFactors = FALSE
  )
  bm_filt_out <- filter_gold_standard_metrics(bm_filt_in)
  expect_equal(nrow(bm_filt_out), 5L)
  expect_true("(III) Cellular Manifold & Clustering Concordance" %in% bm_filt_out$Category)
  expect_true("(V) Biological Signal & Differential Fidelity" %in% bm_filt_out$Category)
  expect_true("(VI) Multiomics Integration & Alignment" %in% bm_filt_out$Category)

  # Additional file 1 complete 55-measure alignment check
  af1_df <- data.frame(
    Method = rep(c("SimA", "SimB"), each = 10),
    Category = rep(c(
      "(I) Distributional Properties", "(I) Distributional Properties", "(I) Distributional Properties",
      "(II) Sparsity & Zero-Inflation Dynamics",
      "(III) Cellular Manifold & Clustering Concordance",
      "(V) Biological Signal & Differential Fidelity", "(V) Biological Signal & Differential Fidelity",
      "(VI) Multiomics Integration & Alignment", "(VI) Multiomics Integration & Alignment",
      "(VI) Multiomics Integration & Alignment"
    ), 2),
    Property = rep(c(
      "pca_dist", "lib_vs_cellzero", "scaled_variance",
      "accessibility_sparsity_curve",
      "clustering",
      "scRNA_differential_expression", "scATAC_differential_accessibility",
      "multiomics_integration", "multiomics_integration", "multiomics_integration"
    ), 2),
    Metric = rep(c(
      "KDE_ISE", "EMD_2D", "Runs_Statistic",
      "curve_rmse",
      "davies_bouldin",
      "pval_dist_score", "auc",
      "batch_asw", "ari", "nmi"
    ), 2),
    Modality = rep(c("RNA", "ATAC", "RNA", "ATAC", "RNA", "RNA", "ATAC", "Joint (Cross-Modal)", "Joint (Cross-Modal)", "Joint (Cross-Modal)"), 2),
    Score = c(
      0.05, 0.2, -0.5, 0.04, 0.9, 0.92, 0.88, 0.91, 0.85, 0.82,
      0.40, 1.5, -25.0, 0.25, 2.1, 0.45, 0.55, 0.60, 0.40, 0.35
    ),
    stringsAsFactors = FALSE
  )
  af1_filt <- filter_gold_standard_metrics(af1_df)
  expect_equal(nrow(af1_filt), 20L)
  p_af1 <- plot_benchmark_bubble_matrix(af1_filt)
  expect_s3_class(p_af1, "ggplot")
  expect_true(all(c("KDE ISE", "2D EMD", "Runs Statistic", "Accessibility Sparsity RMSE",
                    "DE/DA Null P-Value Uniformity", "DE/DA ROC-AUC", "Batch ASW",
                    "Integration ARI", "Integration NMI") %in% as.character(p_af1$data$Display_Metric)))
  runs_a <- p_af1$data$Score_Norm[p_af1$data$Metric == "Runs_Statistic" & p_af1$data$Method == "SimA"]
  runs_b <- p_af1$data$Score_Norm[p_af1$data$Metric == "Runs_Statistic" & p_af1$data$Method == "SimB"]
  expect_gt(runs_a, runs_b)
  auc_a <- p_af1$data$Score_Norm[p_af1$data$Metric == "auc" & p_af1$data$Method == "SimA"]
  auc_b <- p_af1$data$Score_Norm[p_af1$data$Metric == "auc" & p_af1$data$Method == "SimB"]
  expect_gt(auc_a, auc_b)

  p_rna <- plot_benchmark_bubble_matrix(af1_filt, assay_modality = "RNA")
  expect_true("DE Null P-Value Uniformity" %in% as.character(p_rna$data$Display_Metric))
  expect_true("DE ROC-AUC" %in% as.character(p_rna$data$Display_Metric))

  p_atac <- plot_benchmark_bubble_matrix(af1_filt, assay_modality = "ATAC")
  expect_true("DA Null P-Value Uniformity" %in% as.character(p_atac$data$Display_Metric))
  expect_true("DA ROC-AUC" %in% as.character(p_atac$data$Display_Metric))
})

test_that("compute_dataset_embeddings, plot_dataset_embeddings and compute_embedding_quality_metrics execute properly", {
  set.seed(42)
  ref <- matrix(stats::rpois(2000, 3), nrow = 40, ncol = 50)
  sim1 <- matrix(stats::rpois(2000, 3), nrow = 40, ncol = 50)
  sim2 <- matrix(stats::rpois(2000, 4), nrow = 40, ncol = 50)
  ct <- sample(c("TypeA", "TypeB"), 50, replace = TRUE)
  batch <- sample(c("Batch1", "Batch2"), 50, replace = TRUE)
  
  # Compute embeddings
  emb <- compute_dataset_embeddings(
    reference = ref,
    simulated = list(Splatter = sim1, scDesign3 = sim2),
    reduction = "pca",
    n_pcs = 5,
    cell_types = ct,
    batch = batch
  )
  
  expect_s3_class(emb, "data.frame")
  expect_true(all(c("Dataset", "Role", "Dim1", "Dim2", "Cell_Type", "Batch") %in% colnames(emb)))
  expect_equal(nrow(emb), 150)
  
  # Plot faceted grid
  p_facet <- plot_dataset_embeddings(emb, reduction = "pca", layout = "facet", color_by = "cell_type")
  expect_s3_class(p_facet, "ggplot")
  
  # Plot side-by-side
  p_side <- plot_dataset_embeddings(emb, reduction = "pca", layout = "side_by_side", selected_methods = "Splatter", color_by = "cluster")
  expect_s3_class(p_side, "ggplot")
  
  # Compute quality metrics table
  qm <- compute_embedding_quality_metrics(emb)
  expect_s3_class(qm, "data.frame")
  expect_equal(nrow(qm), 3)
  expect_true(all(c("Dataset", "Role", "Mean Silhouette", "ARI (Cluster Fidelity)", "Library Size Diff (%)") %in% colnames(qm)))
})

test_that(".normalize_bubble_scores keeps 2D KS distances on a proper min-max scale", {
  # Realistic [0, 1] Fasano-Franceschini 2D KS test statistics
  df_real <- data.frame(
    Method = c("Method_A", "Method_B"),
    Metric = rep("Fasano_Franceschini_2D_KS", 2),
    Score_Raw = c(0.18, 0.06),
    stringsAsFactors = FALSE
  )
  out_real <- .normalize_bubble_scores(df_real)
  expect_true(all(out_real$Score_Norm >= 0 & out_real$Score_Norm <= 1))
  expect_equal(out_real$Score_Norm, c(0, 1), tolerance = 1e-8)

  # Arbitrary / large raw test statistics
  df_large <- data.frame(
    Method = c("Method_A", "Method_B"),
    Metric = rep("Fasano_Franceschini_2D_KS", 2),
    Score_Raw = c(850, 100),
    stringsAsFactors = FALSE
  )
  out_large <- .normalize_bubble_scores(df_large)
  expect_true(all(out_large$Score_Norm >= 0 & out_large$Score_Norm <= 1))
  expect_equal(out_large$Score_Norm, c(0, 1), tolerance = 1e-8)
  expect_false(all(out_large$Score_Norm == 0))
})
