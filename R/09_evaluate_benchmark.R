
#' Evaluate Comprehensive Simulation Accuracy (Unimodal Omics Layer)
#'
#' Computes distribution distances, cellular and feature-level properties,
#' zero-inflation curves, and manifold distances between reference and simulated datasets.
#'
#' \if{html}{\figure{workflow_diagram.png}{options: width="100\%" alt="scSimEval Unified Benchmarking Workflow Diagram"}}
#'
#' @param ref_data Reference count matrix (features x cells).
#' @param sim_data Simulated count matrix (features x cells).
#' @param compute_bivariate Logical, whether to compute 2D bivariate tests. Default is TRUE.
#' @param threads Number of threads for parallel computation. Default is 1.
#' @param cpu_time Optional numeric value of CPU execution time in seconds.
#' @param memory_mb Optional numeric value of peak memory allocation in MB.
#' @param system_time Optional numeric value of system CPU time in seconds.
#' @param elapsed_time Optional numeric value of wall-clock elapsed time in seconds.
#' @param verbose Logical, whether to print execution messages. Default is TRUE.
#'
#' @return A list containing tidy summary tables and detailed metric lists.
#' @export
evaluate_simulation_accuracy <- function(
  ref_data,
  sim_data,
  compute_bivariate = TRUE,
  threads = 1,
  cpu_time = NULL,
  memory_mb = NULL,
  system_time = NULL,
  elapsed_time = NULL,
  verbose = TRUE
) {
  # Step 1: Extract properties
  if (verbose) message("[1/5] Extracting cell-level properties...")
  ref_cell <- extract_cell_properties(ref_data, verbose = FALSE)
  sim_cell <- extract_cell_properties(sim_data, verbose = FALSE)
  
  if (verbose) message("[2/5] Extracting feature-level properties...")
  ref_feat <- extract_feature_properties(ref_data, verbose = FALSE)
  sim_feat <- extract_feature_properties(sim_data, verbose = FALSE)
  
  # Step 2: Compute Univariate Accuracy Metrics (Cell Level)
  if (verbose) message("[3/5] Computing univariate accuracy metrics...")
  cell_props_to_eval <- c("library_size", "log_library_size", "zero_fraction_cell", "cell_cor")
  cell_results <- list()
  for (p in cell_props_to_eval) {
    if (!is.null(ref_cell[[p]]) && !is.null(sim_cell[[p]]) && length(ref_cell[[p]]) > 0) {
      cell_results[[p]] <- calc_all_univariate_metrics(ref_cell[[p]], sim_cell[[p]], p)
    }
  }
  
  # Step 3: Compute Univariate Accuracy Metrics (Feature Level)
  feat_props_to_eval <- c("mean_expression", "variance", "sd", "cv", "zero_fraction_feature", "dispersion", "gene_cor")
  feat_results <- list()
  for (p in feat_props_to_eval) {
    if (!is.null(ref_feat[[p]]) && !is.null(sim_feat[[p]]) && length(ref_feat[[p]]) > 0) {
      feat_results[[p]] <- calc_all_univariate_metrics(ref_feat[[p]], sim_feat[[p]], p)
    }
  }
  
  # Step 4: Compute Bivariate Accuracy Metrics
  bivariate_results <- list()
  if (compute_bivariate) {
    if (verbose) message("[4/5] Computing bivariate metrics (Fasano-Franceschini, Peacock, KDE zstat, 2D EMD)...")
    
    # 1) Library size vs Cell zero fraction
    ref_lib_zero <- cbind(log10(ref_cell$library_size + 1), ref_cell$zero_fraction_cell)
    sim_lib_zero <- cbind(log10(sim_cell$library_size + 1), sim_cell$zero_fraction_cell)
    bivariate_results$lib_vs_cellzero <- calc_all_bivariate_metrics(
      ref_lib_zero, sim_lib_zero, "lib_vs_cellzero", threads = threads
    )
    
    # 2) Mean vs SD
    ref_mean_sd <- cbind(ref_feat$mean_expression, ref_feat$sd)
    sim_mean_sd <- cbind(sim_feat$mean_expression, sim_feat$sd)
    bivariate_results$mean_vs_sd <- calc_all_bivariate_metrics(
      ref_mean_sd, sim_mean_sd, "mean_vs_sd", threads = threads
    )
    
    # 3) Mean vs Feature Zero Fraction (Mean-Dropout)
    ref_mean_zero <- cbind(ref_feat$mean_expression, ref_feat$zero_fraction_feature)
    sim_mean_zero <- cbind(sim_feat$mean_expression, sim_feat$zero_fraction_feature)
    bivariate_results$mean_vs_featzero <- calc_all_bivariate_metrics(
      ref_mean_zero, sim_mean_zero, "mean_vs_featzero", threads = threads
    )
    
    # 4) Mean vs Dispersion (Mean-Variance)
    ref_mean_disp <- cbind(ref_feat$mean_expression, ref_feat$dispersion)
    sim_mean_disp <- cbind(sim_feat$mean_expression, sim_feat$dispersion)
    bivariate_results$mean_vs_dispersion <- calc_all_bivariate_metrics(
      ref_mean_disp, sim_mean_disp, "mean_vs_dispersion", threads = threads
    )
  }
  
  # Step 5: Dropout, Zero-Inflation, & Generative Manifold
  if (verbose) message("[5/5] Computing zero-probability and manifold distances...")
  zero_prob_eval <- evaluate_zero_probability_curve(ref_data, sim_data)
  ref_excess_w <- calc_excess_zero_weights(ref_data)
  sim_excess_w <- calc_excess_zero_weights(sim_data)
  
  mmd_eval <- calc_mmd(ref_data, sim_data)
  fsd_eval <- calc_frechet_singlecell_distance(ref_data, sim_data)
  gen_pr <- calc_generative_precision_recall(ref_data, sim_data)
  
  # Build tidy summary table
  summary_rows <- list()
  all_univ_metric_names <- c(
    "MAD", "KS", "MAE", "RMSE", "OV", "Bhattacharyya", "Wasserstein",
    "ECDF_DiffArea", "Runs_Statistic", "Runs_PValue", "NN_Mismatch",
    "Between_Dataset_Silh_Global", "Between_Dataset_Silh_Local"
  )
  
  for (prop in names(cell_results)) {
    for (m in all_univ_metric_names) {
      key <- paste0(prop, "_", m)
      if (!is.null(cell_results[[prop]][[key]])) {
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          Category = "Distributional Properties",
          Property = prop,
          Metric = m,
          Value = as.numeric(cell_results[[prop]][[key]]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  for (prop in names(feat_results)) {
    for (m in all_univ_metric_names) {
      key <- paste0(prop, "_", m)
      if (!is.null(feat_results[[prop]][[key]])) {
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          Category = "Distributional Properties",
          Property = prop,
          Metric = m,
          Value = as.numeric(feat_results[[prop]][[key]]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (compute_bivariate) {
    all_biv_metric_names <- c(
      "Fasano_Franceschini_2D_KS", "Peacock_2D_KS", "KDE_Bivariate_zstat",
      "EMD_2D", "NN_Mismatch_2D", "Between_Dataset_Silh_2D_Global",
      "Between_Dataset_Silh_2D_Local"
    )
    for (b_name in names(bivariate_results)) {
      for (m in all_biv_metric_names) {
        key <- paste0(b_name, "_", m)
        if (!is.null(bivariate_results[[b_name]][[key]])) {
          summary_rows[[length(summary_rows) + 1]] <- data.frame(
            Category = "Correlation & Dependencies",
            Property = b_name,
            Metric = m,
            Value = as.numeric(bivariate_results[[b_name]][[key]]),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  for (m in names(zero_prob_eval)) {
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      Category = "Distributional Properties",
      Property = "dropout_curve",
      Metric = m,
      Value = as.numeric(zero_prob_eval[[m]]),
      stringsAsFactors = FALSE
    )
  }
  
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    Category = "Distributional Properties",
    Property = "excess_zero_weights",
    Metric = "excess_zero_weight_discrepancy",
    Value = abs(ref_excess_w$mean_excess_zero_weight - sim_excess_w$mean_excess_zero_weight),
    stringsAsFactors = FALSE
  )
  
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    Category = "Cellular Structure & Mixing",
    Property = "manifold",
    Metric = "MMD",
    Value = as.numeric(mmd_eval$mmd),
    stringsAsFactors = FALSE
  )
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    Category = "Cellular Structure & Mixing",
    Property = "manifold",
    Metric = "Frechet_SingleCell_Distance",
    Value = as.numeric(fsd_eval$fsd),
    stringsAsFactors = FALSE
  )
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    Category = "Cellular Structure & Mixing",
    Property = "manifold",
    Metric = "Generative_Precision",
    Value = as.numeric(gen_pr$generative_precision),
    stringsAsFactors = FALSE
  )
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    Category = "Cellular Structure & Mixing",
    Property = "manifold",
    Metric = "Generative_Recall",
    Value = as.numeric(gen_pr$generative_recall),
    stringsAsFactors = FALSE
  )
  
  # Step 6: Computational Scalability (Category VIII)
  if (!is.null(elapsed_time) || !is.null(memory_mb)) {
    if (verbose) message("[6/6] Logging computational scalability metrics...")
    if (!is.null(elapsed_time)) {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        Category = "Computational Scalability",
        Property = "resource_usage",
        Metric = "elapsed_time_seconds",
        Value = as.numeric(elapsed_time),
        stringsAsFactors = FALSE
      )
    }
    if (!is.null(memory_mb)) {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        Category = "Computational Scalability",
        Property = "resource_usage",
        Metric = "peak_memory_mb",
        Value = as.numeric(memory_mb),
        stringsAsFactors = FALSE
      )
    }
  }

  summary_table <- do.call(rbind, summary_rows)
  if (verbose) message("Unimodal accuracy evaluation complete.")
  
  list(
    metrics_summary_table = summary_table,
    cell_metrics = cell_results,
    feature_metrics = feat_results,
    bivariate_metrics = bivariate_results,
    zero_probability_curve = zero_prob_eval,
    excess_zero_weights = list(ref = ref_excess_w, sim = sim_excess_w),
    generative_manifold = list(
      mmd = mmd_eval,
      frechet_distance = fsd_eval,
      precision_recall = gen_pr
    ),
    ref_properties = list(cell = ref_cell, feature = ref_feat),
    sim_properties = list(cell = sim_cell, feature = sim_feat),
    resource_usage = list(cpu_time = cpu_time, memory_mb = memory_mb, system_time = system_time, elapsed_time = elapsed_time)
  )
}

#' Master Multiomics Simulation Benchmarking Pipeline
#'
#' Orchestrates comprehensive benchmarking across all 7 evaluation categories
#' using strictly the user's provided inputs:
#' \itemize{
#'   \item \code{ref_multi}: List of reference omics matrices (e.g. \code{list(rna = ..., atac = ...)}).
#'   \item \code{sim_multi}: List of simulated omics matrices (e.g. \code{list(rna = ..., atac = ...)}).
#'   \item \code{cell_types}: Vector/factor of cell-type identities.
#'   \item \code{batch_info}: Vector/factor of batch or donor annotations.
#'   \item Computational resource usage: \code{cpu_time}, \code{memory_mb}, \code{system_time}.
#'   \item Trajectory features inferred directly from scRNA-seq counts.
#' }
#'
#' \if{html}{\figure{workflow_diagram.png}{options: width="100\%" alt="scSimEval Unified Benchmarking Workflow Diagram"}}
#'
#' @param ref_multi Named list of reference matrices for Modality 1 and 2.
#' @param sim_multi Named list of simulated matrices for Modality 1 and 2.
#' @param cell_types Optional factor or vector of cell type labels for cells.
#' @param batch_info Optional factor or vector of batch labels for cells.
#' @param cpu_time Optional numeric value of CPU execution time in seconds.
#' @param memory_mb Optional numeric value of peak memory allocation in MB.
#' @param system_time Optional numeric value of system CPU time in seconds.
#' @param elapsed_time Optional numeric value of wall-clock elapsed time in seconds.
#' @param feature_pairs Optional 2-column data.frame of linked feature pairs.
#' @param compute_bivariate Logical, whether to compute 2D bivariate tests. Default FALSE for speed.
#' @param threads Number of CPU threads. Default is 1.
#' @param verbose Logical, whether to print execution progress. Default is TRUE.
#'
#' @return A list containing detailed results and a unified tidy master summary table.
#' @export
evaluate_multiomics_accuracy <- function(
  ref_multi,
  sim_multi,
  cell_types = NULL,
  batch_info = NULL,
  cpu_time = NULL,
  memory_mb = NULL,
  system_time = NULL,
  elapsed_time = NULL,
  feature_pairs = NULL,
  compute_bivariate = FALSE,
  threads = 1,
  verbose = TRUE
) {
  mod_names <- names(ref_multi)
  if (is.null(mod_names) || length(mod_names) < 2) {
    mod_names <- c("Modality_1", "Modality_2")
  }
  m1_name <- mod_names[1]
  m2_name <- mod_names[2]
  
  master_rows <- list()
  
  # ============================================================================
  # 1. Unimodal Evaluation for Modality 1 and Modality 2 (Categories 1 & 2)
  # ============================================================================
  if (verbose) message(sprintf("=== [Step 1] Evaluating %s Layer Accuracy ===", m1_name))
  m1_eval <- evaluate_simulation_accuracy(
    ref_data = ref_multi[[1]],
    sim_data = sim_multi[[1]],
    compute_bivariate = compute_bivariate,
    threads = threads,
    verbose = verbose
  )
  t1 <- m1_eval$metrics_summary_table
  t1$Modality <- m1_name
  master_rows[[length(master_rows) + 1]] <- t1
  
  if (verbose) message(sprintf("=== [Step 2] Evaluating %s Layer Accuracy ===", m2_name))
  m2_eval <- evaluate_simulation_accuracy(
    ref_data = ref_multi[[2]],
    sim_data = sim_multi[[2]],
    compute_bivariate = compute_bivariate,
    threads = threads,
    verbose = verbose
  )
  t2 <- m2_eval$metrics_summary_table
  t2$Modality <- m2_name
  master_rows[[length(master_rows) + 1]] <- t2
  
  # ============================================================================
  # 2. Cellular Structure & Clustering Accuracy (Category 3)
  # ============================================================================
  clustering_res <- list()
  if (!is.null(cell_types)) {
    if (verbose) message("=== [Step 3] Evaluating Cellular Structure & Clustering Preservation ===")
    clu_eval <- evaluate_clustering_metrics(
      data = sim_multi[[1]],
      cell_types = cell_types,
      ref_data = ref_multi[[1]]
    )
    clustering_res <- clu_eval
    
    clu_metrics_to_record <- c(
      "silhouette_sim", "silhouette_ref", "dunn_sim", "davies_bouldin_sim",
      "calinski_harabasz_sim", "ari", "nmi", "ami", "fmi", "homogeneity",
      "completeness", "v_measure", "clustering_accuracy", "neighborhood_purity", "cdi"
    )
    for (cm in clu_metrics_to_record) {
      if (!is.null(clu_eval[[cm]])) {
        master_rows[[length(master_rows) + 1]] <- data.frame(
          Category = "Cellular Structure & Mixing",
          Property = "clustering",
          Metric = cm,
          Value = as.numeric(clu_eval[[cm]]),
          Modality = m1_name,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # ============================================================================
  # 3. Batch Effect & Technical Confounder Mixing (Category 4)
  # ============================================================================
  batch_res <- list()
  if (!is.null(batch_info)) {
    if (verbose) message("=== [Step 4] Evaluating Batch Integration & Technical Confounder Mixing ===")
    b_eval <- evaluate_batch_metrics(
      data = sim_multi[[1]],
      batch_info = batch_info,
      pre_data = ref_multi[[1]],
      cell_types = cell_types,
      verbose = FALSE
    )
    batch_res <- b_eval
    
    batch_metrics_to_record <- c(
      "batch_silhouette", "shannon_entropy", "pcr_r2", "cms", "isi",
      "seurat_mixing_metric", "ldf_diff", "local_structure",
      "cross_batch_accuracy", "cross_batch_F1"
    )
    for (bm in batch_metrics_to_record) {
      if (!is.null(b_eval[[bm]])) {
        master_rows[[length(master_rows) + 1]] <- data.frame(
          Category = "Cellular Structure & Mixing",
          Property = "batch_mixing",
          Metric = bm,
          Value = as.numeric(b_eval[[bm]]),
          Modality = m1_name,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # ============================================================================
  # 4. Biological Signal & Downstream Performance (Category 5)
  # ============================================================================
  signal_res <- list()
  if (!is.null(cell_types)) {
    if (verbose) message("=== [Step 5] Evaluating Biological Signal Fidelity (SimBench DE/DV/DD/DP/BD) ===")
    sig_eval <- evaluate_simbench_signals(
      ref_mat = ref_multi[[1]],
      sim_mat = sim_multi[[1]],
      ref_celltypes = cell_types,
      sim_celltypes = cell_types
    )
    signal_res$simbench_signals <- sig_eval
    
    for (i in seq_len(nrow(sig_eval))) {
      master_rows[[length(master_rows) + 1]] <- data.frame(
        Category = "Biological Signal & Downstream",
        Property = "simbench_signal",
        Metric = paste0(sig_eval$Signal_Type[i], "_abs_error"),
        Value = as.numeric(sig_eval$Absolute_Error[i]),
        Modality = m1_name,
        stringsAsFactors = FALSE
      )
    }
    
    # Predictive cell identity classifier
    pred_model <- evaluate_predictive_de_model(
      data = sim_multi[[1]],
      group = cell_types
    )
    signal_res$predictive_model <- pred_model
    master_rows[[length(master_rows) + 1]] <- data.frame(
      Category = "Biological Signal & Downstream",
      Property = "cell_identity_prediction",
      Metric = "predictive_accuracy",
      Value = as.numeric(pred_model$accuracy),
      Modality = m1_name,
      stringsAsFactors = FALSE
    )
    master_rows[[length(master_rows) + 1]] <- data.frame(
      Category = "Biological Signal & Downstream",
      Property = "cell_identity_prediction",
      Metric = "predictive_macro_F1",
      Value = as.numeric(pred_model$F1),
      Modality = m1_name,
      stringsAsFactors = FALSE
    )
    
    # Deconvolution mixture accuracy
    deconv_res <- calc_deconvolution_accuracy(
      true_proportions = cell_types,
      estimated_proportions = cell_types,
      batch_ref = batch_info,
      batch_sim = batch_info
    )
    signal_res$deconvolution <- deconv_res
    master_rows[[length(master_rows) + 1]] <- data.frame(
      Category = "Biological Signal & Downstream",
      Property = "cell_mixture",
      Metric = "deconvolution_RMSE",
      Value = as.numeric(deconv_res$deconvolution_rmse),
      Modality = m1_name,
      stringsAsFactors = FALSE
    )
    master_rows[[length(master_rows) + 1]] <- data.frame(
      Category = "Biological Signal & Downstream",
      Property = "cell_mixture",
      Metric = "deconvolution_JSD",
      Value = as.numeric(deconv_res$deconvolution_jsd),
      Modality = m1_name,
      stringsAsFactors = FALSE
    )
    
    # Delta variance if batch_info available
    if (!is.null(batch_info)) {
      dv_ref <- calc_delta_variance(ref_multi[[1]], replicates = batch_info, conditions = cell_types)
      dv_sim <- calc_delta_variance(sim_multi[[1]], replicates = batch_info, conditions = cell_types)
      dv_mae <- mean(abs(dv_ref$delta_variance - dv_sim$delta_variance), na.rm = TRUE)
      master_rows[[length(master_rows) + 1]] <- data.frame(
        Category = "Biological Signal & Downstream",
        Property = "pseudoreplication",
        Metric = "delta_variance_MAE",
        Value = as.numeric(dv_mae),
        Modality = m1_name,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # ============================================================================
  # 5. Trajectory & Differentiation Lineage Dynamics (Category 6)
  # ============================================================================
  if (verbose) message("=== [Step 6] Evaluating Trajectory & Lineage Dynamics Inferred from scRNA-seq ===")
  traj_eval <- evaluate_trajectory_metrics(
    ref_data = ref_multi[[1]],
    sim_data = sim_multi[[1]],
    cell_types_ref = cell_types,
    cell_types_sim = cell_types
  )
  master_rows[[length(master_rows) + 1]] <- data.frame(
    Category = "Trajectory Dynamics",
    Property = "pseudotime",
    Metric = "pseudotime_correlation",
    Value = as.numeric(traj_eval$pseudotime_correlation),
    Modality = m1_name,
    stringsAsFactors = FALSE
  )
  if (!is.null(traj_eval$tree_height_rmse) && !is.na(traj_eval$tree_height_rmse)) {
    master_rows[[length(master_rows) + 1]] <- data.frame(
      Category = "Trajectory Dynamics",
      Property = "lineage_tree",
      Metric = "tree_height_rmse",
      Value = as.numeric(traj_eval$tree_height_rmse),
      Modality = m1_name,
      stringsAsFactors = FALSE
    )
  }
  
  # ============================================================================
  # 6. Cross-Modality Coupling & Modularity (Category 7)
  # ============================================================================
  if (verbose) message("=== [Step 7] Evaluating Cross-Modality Regulatory Linkage & Alignment ===")
  cross_cor_res <- calc_cross_modality_correlation(
    ref_mod1 = ref_multi[[1]],
    ref_mod2 = ref_multi[[2]],
    sim_mod1 = sim_multi[[1]],
    sim_mod2 = sim_multi[[2]],
    feature_pairs = feature_pairs
  )
  
  cross_cor_metrics <- c("MAD", "KS", "MAE", "RMSE", "OV", "Bhattacharyya", "Wasserstein")
  for (m in cross_cor_metrics) {
    key <- paste0("cross_modality_cor_", m)
    if (!is.null(cross_cor_res[[key]])) {
      master_rows[[length(master_rows) + 1]] <- data.frame(
        Category = "Cross-Modal Relationships",
        Property = "cross_feature_correlation",
        Metric = m,
        Value = as.numeric(cross_cor_res[[key]]),
        Modality = "Joint (Cross-Modal)",
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Cross-modal label transfer
  if (!is.null(cell_types)) {
    cpred <- evaluate_cross_modal_prediction(sim_multi[[2]], cell_types)
    master_rows[[length(master_rows) + 1]] <- data.frame(
      Category = "Cross-Modal Relationships",
      Property = "label_transfer",
      Metric = "cross_modal_accuracy",
      Value = as.numeric(cpred$cross_modal_accuracy),
      Modality = "Joint (Cross-Modal)",
      stringsAsFactors = FALSE
    )
    master_rows[[length(master_rows) + 1]] <- data.frame(
      Category = "Cross-Modal Relationships",
      Property = "label_transfer",
      Metric = "cross_modal_F1",
      Value = as.numeric(cpred$cross_modal_F1),
      Modality = "Joint (Cross-Modal)",
      stringsAsFactors = FALSE
    )
  }
  
  # FOSCTTM cell alignment in joint latent space
  sub_m1 <- as.matrix(sim_multi[[1]])
  sub_m2 <- as.matrix(sim_multi[[2]])
  min_c <- min(ncol(sub_m1), ncol(sub_m2))
  p1 <- stats::prcomp(t(sub_m1[, seq_len(min_c)]), rank. = 5)$x
  p2 <- stats::prcomp(t(sub_m2[, seq_len(min_c)]), rank. = 5)$x
  fos_res <- calc_foscttm(p1, p2)
  master_rows[[length(master_rows) + 1]] <- data.frame(
    Category = "Cross-Modal Relationships",
    Property = "cell_alignment",
    Metric = "FOSCTTM",
    Value = as.numeric(fos_res$foscttm),
    Modality = "Joint (Cross-Modal)",
    stringsAsFactors = FALSE
  )
  master_rows[[length(master_rows) + 1]] <- data.frame(
    Category = "Cross-Modal Relationships",
    Property = "cell_alignment",
    Metric = "match_at_1",
    Value = as.numeric(fos_res$match_at_1),
    Modality = "Joint (Cross-Modal)",
    stringsAsFactors = FALSE
  )
  
  # In silico generation fidelity
  gen_fid <- calc_cross_modal_generation(ref_multi[[2]], sim_multi[[2]])
  master_rows[[length(master_rows) + 1]] <- data.frame(
    Category = "Cross-Modal Relationships",
    Property = "in_silico_generation",
    Metric = "cell_pearson_cor",
    Value = as.numeric(gen_fid$mean_cell_pcc),
    Modality = m2_name,
    stringsAsFactors = FALSE
  )
  master_rows[[length(master_rows) + 1]] <- data.frame(
    Category = "Cross-Modal Relationships",
    Property = "in_silico_generation",
    Metric = "feature_pearson_cor",
    Value = as.numeric(gen_fid$mean_feat_pcc),
    Modality = m2_name,
    stringsAsFactors = FALSE
  )
  
  # Regulatory coupling
  coupling_res <- calc_atac_rna_coupling(sim_multi[[2]], sim_multi[[1]])
  master_rows[[length(master_rows) + 1]] <- data.frame(
    Category = "Cross-Modal Relationships",
    Property = "regulatory_coupling",
    Metric = "mean_coupling_cor",
    Value = as.numeric(coupling_res$mean_coupling_cor),
    Modality = "Joint (Cross-Modal)",
    stringsAsFactors = FALSE
  )
  
  # Peak co-accessibility fidelity
  peak_coacc <- calc_peak_coaccessibility_fidelity(ref_multi[[2]], sim_multi[[2]])
  master_rows[[length(master_rows) + 1]] <- data.frame(
    Category = "Cross-Modal Relationships",
    Property = "chromatin_coaccessibility",
    Metric = "rv_coefficient",
    Value = as.numeric(peak_coacc$rv_coefficient),
    Modality = m2_name,
    stringsAsFactors = FALSE
  )
  
  # Co-expression module fidelity
  coexpr_fid <- calc_coexpression_module_fidelity(ref_multi[[1]], sim_multi[[1]])
  master_rows[[length(master_rows) + 1]] <- data.frame(
    Category = "Cross-Modal Relationships",
    Property = "coexpression_modularity",
    Metric = "module_correlation_r",
    Value = as.numeric(coexpr_fid$module_correlation_r),
    Modality = m1_name,
    stringsAsFactors = FALSE
  )
  
  # Accessibility profile concordance
  acc_prof <- calc_accessibility_profile_concordance(ref_multi[[2]], sim_multi[[2]], cell_types = cell_types)
  master_rows[[length(master_rows) + 1]] <- data.frame(
    Category = "Cross-Modal Relationships",
    Property = "accessibility_profile",
    Metric = "global_PCC",
    Value = as.numeric(acc_prof$global_pcc),
    Modality = m2_name,
    stringsAsFactors = FALSE
  )
  
  # ============================================================================
  # 7. Computational Scalability & Usability (Category 8)
  # ============================================================================
  if (!is.null(elapsed_time) || !is.null(memory_mb)) {
    if (verbose) message("=== [Step 8] Logging Computational Scalability & Usability ===")
    if (!is.null(elapsed_time)) {
      master_rows[[length(master_rows) + 1]] <- data.frame(
        Category = "Computational Scalability",
        Property = "resource_usage",
        Metric = "elapsed_time_seconds",
        Value = as.numeric(elapsed_time),
        Modality = "Computation",
        stringsAsFactors = FALSE
      )
    }
    if (!is.null(memory_mb)) {
      master_rows[[length(master_rows) + 1]] <- data.frame(
        Category = "Computational Scalability",
        Property = "resource_usage",
        Metric = "peak_memory_mb",
        Value = as.numeric(memory_mb),
        Modality = "Computation",
        stringsAsFactors = FALSE
      )
    }
  }
  
  master_table <- do.call(rbind, master_rows)
  if (verbose) message("Master multiomics benchmarking complete.")
  
  list(
    benchmark_summary_table = master_table,
    mod1_unimodal = m1_eval,
    mod2_unimodal = m2_eval,
    clustering = clustering_res,
    batch = batch_res,
    signal = signal_res,
    trajectory = traj_eval,
    cross_modality = list(
      correlation = cross_cor_res,
      foscttm = fos_res,
      generation = gen_fid,
      coupling = coupling_res,
      coaccessibility = peak_coacc,
      coexpression = coexpr_fid,
      profile = acc_prof
    )
  )
}


#' Evaluate Multiple Single-Cell Datasets Simultaneously
#'
#' Orchestrates simultaneous benchmarking across a collection of single-cell datasets,
#' including paired multiomics (e.g., scRNA-seq and scATAC-seq from the same biological experiment)
#' and unimodal datasets (e.g., standalone scRNA-seq or scATAC-seq). Produces a unified,
#' consolidated summary table, comparative ranking matrix, and dataset overview.
#'
#' @param datasets A named list of datasets to evaluate. Can be provided in two formats:
#'   \itemize{
#'     \item \strong{Format A (Grouped Multiomics & Unimodal):} A list of dataset specifications, e.g.:
#'       \code{list("Method 1" = list(ref = list(rna = ..., atac = ...), sim = list(rna = ..., atac = ...)), "Method 3" = list(ref = ..., sim = ..., modality = "scATAC-seq"))}
#'     \item \strong{Format B (Flat Named Datasets):} A list with named modality datasets, e.g.:
#'       \code{list("Method 1-scRNA-seq" = list(ref = ..., sim = ...), "Method 1-scATAC-seq" = list(ref = ..., sim = ...), "Method 3-scATAC-seq" = list(ref = ..., sim = ...), "Method 4-scRNA-seq" = list(ref = ..., sim = ...))}
#'   }
#' @param pair_by_prefix Logical; if \code{TRUE} (default), datasets sharing a common root prefix
#'   (such as "Method 1-scRNA-seq" and "Method 1-scATAC-seq") are automatically paired as multiomics
#'   for joint cross-modal evaluation. If \code{FALSE}, each dataset entry is evaluated independently.
#' @param compute_bivariate Logical; whether to compute 2D bivariate tests (default FALSE for speed).
#' @param threads Number of CPU threads for parallel metric computation. Default is 1.
#' @param verbose Logical; whether to print progress messages. Default is TRUE.
#'
#' @return An S3 object of class \code{scSimEval_consolidated} containing:
#'   \itemize{
#'     \item \code{consolidated_summary_table}: Master tidy data.frame with columns \code{Dataset}, \code{Data_Name}, \code{Modality}, \code{Category}, \code{Property}, \code{Metric}, and \code{Value}.
#'     \item \code{consolidated_score_matrix}: Wide comparison matrix comparing all datasets across metrics.
#'     \item \code{dataset_overview}: Aggregate summary table with total metrics and mean scores per dataset.
#'     \item \code{dataset_results}: Named list containing detailed evaluation results for each individual dataset.
#'   }
#' @export
evaluate_multiple_datasets <- function(
  datasets,
  pair_by_prefix = TRUE,
  compute_bivariate = FALSE,
  threads = 1,
  verbose = TRUE
) {
  if (!is.list(datasets) || length(datasets) == 0) {
    stop("`datasets` must be a non-empty named list of datasets.")
  }
  if (is.null(names(datasets)) || any(names(datasets) == "")) {
    names(datasets) <- paste0("Dataset_", seq_along(datasets))
  }

  parsed_experiments <- list()

  if (pair_by_prefix) {
    mod_pattern <- "[-_ ]*(scRNA[-_]?seq|scATAC[-_]?seq|scRNA|scATAC|RNA|ATAC|CITE[-_]?seq|CITE|Methyl)"
    prefixes <- sub(paste0(mod_pattern, ".*$"), "", names(datasets), ignore.case = TRUE)
    prefixes <- sub("[-_ ]+$", "", prefixes)
    prefixes[prefixes == ""] <- names(datasets)[prefixes == ""]

    unique_prefixes <- unique(prefixes)

    for (pfx in unique_prefixes) {
      matching_idx <- which(prefixes == pfx)
      if (length(matching_idx) == 1) {
        idx <- matching_idx[1]
        item <- datasets[[idx]]
        item_name <- names(datasets)[idx]

        ref_mat <- if (!is.null(item$ref)) item$ref else if (!is.null(item$ref_multi)) item$ref_multi else item$ref_data
        sim_mat <- if (!is.null(item$sim)) item$sim else if (!is.null(item$sim_multi)) item$sim_multi else item$sim_data

        if (is.list(ref_mat) && length(ref_mat) >= 2 && !is.data.frame(ref_mat)) {
          parsed_experiments[[pfx]] <- list(
            type = "multiomics",
            name = pfx,
            data_names = if (!is.null(names(ref_mat))) paste0(pfx, "-", names(ref_mat)) else paste0(pfx, c("-Mod1", "-Mod2")),
            ref = ref_mat,
            sim = sim_mat,
            cell_types = item$cell_types,
            batch_info = item$batch_info,
            cpu_time = item$cpu_time,
            memory_mb = item$memory_mb,
            system_time = item$system_time,
            elapsed_time = item$elapsed_time,
            feature_pairs = item$feature_pairs
          )
        } else {
          inferred_mod <- if (!is.null(item$modality)) item$modality else {
            if (grepl("ATAC", item_name, ignore.case = TRUE)) "scATAC-seq"
            else if (grepl("RNA", item_name, ignore.case = TRUE)) "scRNA-seq"
            else "Single-Cell"
          }
          parsed_experiments[[item_name]] <- list(
            type = "unimodal",
            name = pfx,
            data_name = item_name,
            modality = inferred_mod,
            ref = ref_mat,
            sim = sim_mat,
            cell_types = item$cell_types,
            batch_info = item$batch_info
          )
        }
      } else {
        sub_items <- datasets[matching_idx]
        sub_names <- names(datasets)[matching_idx]

        is_atac <- grepl("ATAC", sub_names, ignore.case = TRUE)
        is_rna <- grepl("RNA", sub_names, ignore.case = TRUE)

        if (any(is_rna) && any(is_atac)) {
          rna_idx <- which(is_rna)[1]
          atac_idx <- which(is_atac)[1]

          rna_ref <- if (!is.null(sub_items[[rna_idx]]$ref)) sub_items[[rna_idx]]$ref else sub_items[[rna_idx]]$ref_data
          rna_sim <- if (!is.null(sub_items[[rna_idx]]$sim)) sub_items[[rna_idx]]$sim else sub_items[[rna_idx]]$sim_data
          atac_ref <- if (!is.null(sub_items[[atac_idx]]$ref)) sub_items[[atac_idx]]$ref else sub_items[[atac_idx]]$ref_data
          atac_sim <- if (!is.null(sub_items[[atac_idx]]$sim)) sub_items[[atac_idx]]$sim else sub_items[[atac_idx]]$sim_data

          ref_list <- list("scRNA-seq" = rna_ref, "scATAC-seq" = atac_ref)
          sim_list <- list("scRNA-seq" = rna_sim, "scATAC-seq" = atac_sim)

          ct <- if (!is.null(sub_items[[rna_idx]]$cell_types)) sub_items[[rna_idx]]$cell_types else sub_items[[atac_idx]]$cell_types
          bi <- if (!is.null(sub_items[[rna_idx]]$batch_info)) sub_items[[rna_idx]]$batch_info else sub_items[[atac_idx]]$batch_info

          parsed_experiments[[pfx]] <- list(
            type = "multiomics",
            name = pfx,
            data_names = c(sub_names[rna_idx], sub_names[atac_idx]),
            ref = ref_list,
            sim = sim_list,
            cell_types = ct,
            batch_info = bi,
            cpu_time = sub_items[[rna_idx]]$cpu_time,
            memory_mb = sub_items[[rna_idx]]$memory_mb,
            system_time = sub_items[[rna_idx]]$system_time,
            elapsed_time = sub_items[[rna_idx]]$elapsed_time,
            feature_pairs = sub_items[[rna_idx]]$feature_pairs
          )
        } else {
          for (k in seq_along(matching_idx)) {
            idx_k <- matching_idx[k]
            it_name <- names(datasets)[idx_k]
            it <- datasets[[idx_k]]
            ref_m <- if (!is.null(it$ref)) it$ref else it$ref_data
            sim_m <- if (!is.null(it$sim)) it$sim else it$sim_data
            inferred_m <- if (!is.null(it$modality)) it$modality else {
              if (grepl("ATAC", it_name, ignore.case = TRUE)) "scATAC-seq"
              else if (grepl("RNA", it_name, ignore.case = TRUE)) "scRNA-seq"
              else "Single-Cell"
            }
            parsed_experiments[[it_name]] <- list(
              type = "unimodal",
              name = pfx,
              data_name = it_name,
              modality = inferred_m,
              ref = ref_m,
              sim = sim_m,
              cell_types = it$cell_types,
              batch_info = it$batch_info
            )
          }
        }
      }
    }
  } else {
    for (d_name in names(datasets)) {
      item <- datasets[[d_name]]
      ref_mat <- if (!is.null(item$ref)) item$ref else if (!is.null(item$ref_multi)) item$ref_multi else item$ref_data
      sim_mat <- if (!is.null(item$sim)) item$sim else if (!is.null(item$sim_multi)) item$sim_multi else item$sim_data

      if (is.list(ref_mat) && length(ref_mat) >= 2 && !is.data.frame(ref_mat)) {
        parsed_experiments[[d_name]] <- list(
          type = "multiomics",
          name = d_name,
          data_names = if (!is.null(names(ref_mat))) paste0(d_name, "-", names(ref_mat)) else paste0(d_name, c("-Mod1", "-Mod2")),
          ref = ref_mat,
          sim = sim_mat,
          cell_types = item$cell_types,
          batch_info = item$batch_info,
          cpu_time = item$cpu_time,
          memory_mb = item$memory_mb,
          system_time = item$system_time,
          elapsed_time = item$elapsed_time,
          feature_pairs = item$feature_pairs
        )
      } else {
        inferred_mod <- if (!is.null(item$modality)) item$modality else {
          if (grepl("ATAC", d_name, ignore.case = TRUE)) "scATAC-seq"
          else if (grepl("RNA", d_name, ignore.case = TRUE)) "scRNA-seq"
          else "Single-Cell"
        }
        parsed_experiments[[d_name]] <- list(
          type = "unimodal",
          name = d_name,
          data_name = d_name,
          modality = inferred_mod,
          ref = ref_mat,
          sim = sim_mat,
          cell_types = item$cell_types,
          batch_info = item$batch_info
        )
      }
    }
  }

  all_master_rows <- list()
  dataset_results <- list()
  n_exp <- length(parsed_experiments)

  if (verbose) {
    message(sprintf(">>> Initializing scSimEval multi-dataset benchmarking for %d experiment(s)...", n_exp))
  }

  for (i in seq_along(parsed_experiments)) {
    exp <- parsed_experiments[[i]]
    exp_name <- exp$name
    if (verbose) {
      message(sprintf("\n[%d/%d] Benchmarking dataset: '%s' (%s)...", i, n_exp, exp_name, exp$type))
    }

    if (exp$type == "multiomics") {
      res <- evaluate_multiomics_accuracy(
        ref_multi = exp$ref,
        sim_multi = exp$sim,
        cell_types = exp$cell_types,
        batch_info = exp$batch_info,
        cpu_time = exp$cpu_time,
        memory_mb = exp$memory_mb,
        system_time = exp$system_time,
        elapsed_time = exp$elapsed_time,
        feature_pairs = exp$feature_pairs,
        compute_bivariate = compute_bivariate,
        threads = threads,
        verbose = verbose
      )

      df <- res$benchmark_summary_table
      df$Dataset <- exp_name

      m1_tag <- names(exp$ref)[1]
      m2_tag <- names(exp$ref)[2]
      if (is.null(m1_tag)) m1_tag <- "scRNA-seq"
      if (is.null(m2_tag)) m2_tag <- "scATAC-seq"

      dn1 <- if (length(exp$data_names) >= 1) exp$data_names[1] else paste0(exp_name, "-", m1_tag)
      dn2 <- if (length(exp$data_names) >= 2) exp$data_names[2] else paste0(exp_name, "-", m2_tag)

      data_name_col <- rep("", nrow(df))
      for (r in seq_len(nrow(df))) {
        mod_r <- df$Modality[r]
        if (mod_r == m1_tag || grepl("RNA", mod_r, ignore.case = TRUE)) {
          data_name_col[r] <- dn1
        } else if (mod_r == m2_tag || grepl("ATAC", mod_r, ignore.case = TRUE)) {
          data_name_col[r] <- dn2
        } else if (grepl("Cross|Joint", mod_r, ignore.case = TRUE)) {
          data_name_col[r] <- paste0(exp_name, "-Joint")
        } else {
          data_name_col[r] <- paste0(exp_name, "-", mod_r)
        }
      }
      df$Data_Name <- data_name_col
      df <- df[, c("Dataset", "Data_Name", "Modality", "Category", "Property", "Metric", "Value")]
      all_master_rows[[length(all_master_rows) + 1]] <- df
      dataset_results[[exp_name]] <- res

    } else {
      res <- evaluate_simulation_accuracy(
        ref_data = exp$ref,
        sim_data = exp$sim,
        compute_bivariate = compute_bivariate,
        threads = threads,
        verbose = verbose
      )

      df <- res$metrics_summary_table
      df$Dataset <- exp_name
      df$Data_Name <- exp$data_name
      df$Modality <- exp$modality
      df <- df[, c("Dataset", "Data_Name", "Modality", "Category", "Property", "Metric", "Value")]
      all_master_rows[[length(all_master_rows) + 1]] <- df
      dataset_results[[exp$data_name]] <- res
    }
  }

  consolidated_table <- do.call(rbind, all_master_rows)
  rownames(consolidated_table) <- NULL

  # Build consolidated score matrix (Metric x Data_Name)
  metric_keys <- paste0(consolidated_table$Category, " | ", consolidated_table$Property, " [", consolidated_table$Metric, "]")
  unique_keys <- unique(metric_keys)
  unique_data_names <- unique(consolidated_table$Data_Name)

  score_mat <- matrix(NA_real_, nrow = length(unique_keys), ncol = length(unique_data_names),
                      dimnames = list(unique_keys, unique_data_names))
  for (r in seq_len(nrow(consolidated_table))) {
    k <- metric_keys[r]
    dn <- consolidated_table$Data_Name[r]
    score_mat[k, dn] <- consolidated_table$Value[r]
  }

  # Build dataset overview summary
  overview_list <- list()
  for (dn in unique_data_names) {
    sub_df <- consolidated_table[consolidated_table$Data_Name == dn, ]
    ks_vals <- sub_df$Value[sub_df$Metric == "KS"]
    wass_vals <- sub_df$Value[sub_df$Metric == "Wasserstein"]
    rmse_vals <- sub_df$Value[sub_df$Metric == "RMSE"]

    overview_list[[length(overview_list) + 1]] <- data.frame(
      Dataset = sub_df$Dataset[1],
      Data_Name = dn,
      Modality = sub_df$Modality[1],
      Total_Metrics = nrow(sub_df),
      Mean_KS_Distance = if (length(ks_vals) > 0) round(mean(ks_vals, na.rm = TRUE), 4) else NA_real_,
      Mean_Wasserstein = if (length(wass_vals) > 0) round(mean(wass_vals, na.rm = TRUE), 4) else NA_real_,
      Mean_RMSE = if (length(rmse_vals) > 0) round(mean(rmse_vals, na.rm = TRUE), 4) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  dataset_overview <- do.call(rbind, overview_list)

  out <- list(
    consolidated_summary_table = consolidated_table,
    consolidated_score_matrix = score_mat,
    dataset_overview = dataset_overview,
    dataset_results = dataset_results
  )
  class(out) <- c("scSimEval_consolidated", "list")

  if (verbose) {
    message(sprintf("\n>>> Consolidated benchmark complete across %d dataset entries (%d total metric instances).",
                    length(unique_data_names), nrow(consolidated_table)))
  }

  out
}

#' @export
print.scSimEval_consolidated <- function(x, ...) {
  cat("=====================================================================\n")
  cat("  scSimEval Consolidated Multi-Dataset Benchmark Results\n")
  cat("=====================================================================\n")
  cat(sprintf("Total Evaluated Metrics : %d across %d dataset entries\n",
              nrow(x$consolidated_summary_table),
              length(unique(x$consolidated_summary_table$Data_Name))))
  cat(sprintf("Unique Datasets        : %s\n", paste(unique(x$consolidated_summary_table$Dataset), collapse = ", ")))
  cat(sprintf("Unique Modalities      : %s\n", paste(unique(x$consolidated_summary_table$Modality), collapse = ", ")))
  cat("\nDataset Overview Summary:\n")
  print(x$dataset_overview, row.names = FALSE)
  cat("=====================================================================\n")
  invisible(x)
}

