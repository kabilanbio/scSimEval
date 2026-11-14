
#' Append Category VII dynverse Trajectory Metrics to Summary Rows
#' @keywords internal
#' @noRd
.append_trajectory_rows <- function(traj_res, modality = NULL) {
  if (is.null(traj_res)) return(list())
  spec <- list(
    list(field = "him", prop = "trajectory_topology"),
    list(field = "edge_flip", prop = "trajectory_topology"),
    list(field = "geodesic_correlation", prop = "geodesic_ordering"),
    list(field = "F1_branches", prop = "branch_concordance"),
    list(field = "wcor_features", prop = "dynamic_features")
  )
  rows <- list()
  for (s in spec) {
    val <- traj_res[[s$field]]
    if (!is.null(val) && length(val) == 1L && is.finite(val)) {
      df <- data.frame(
        Category = "(VII) Trajectory & Lineage Topology",
        Property = s$prop,
        Metric = s$field,
        Value = as.numeric(val),
        stringsAsFactors = FALSE
      )
      if (!is.null(modality)) df$Modality <- modality
      rows[[length(rows) + 1L]] <- df
    }
  }
  rows
}

#' Evaluate Comprehensive Simulation Accuracy (Unimodal Omics Layer)
#'
#' Computes distribution distances, cellular and feature-level properties,
#' zero-inflation curves, and manifold distances between reference and simulated datasets.
#'
#' \if{html}{\figure{workflow_diagram.png}{options: width="800" alt="scSimEval Unified Benchmarking Workflow Diagram"}}
#'
#' @details
#' \strong{Metric Measurement Capacity by Data & Metadata Mode (56 Evaluation Measures Suite):}
#' \itemize{
#'   \item \strong{Count Matrix Alone (No Metadata):} Evaluates Category I (Distributional
#'     Properties: 10 univariate metrics---Wasserstein, KS, MAD, MAE, RMSE, Bhattacharyya,
#'     KDE-ISE, Runs statistic, NN mismatch, Inter-dataset silhouette---and 4 bivariate metrics---
#'     Fasano-Franceschini 2D KS, 2D EMD, 2D NN mismatch, 2D Inter-dataset silhouette),
#'     Category II (Sparsity & Zero-Inflation Dynamics: logistic dropout curve midpoint/slope/R2
#'     or accessibility-sparsity curve RMSE/R2, and excess zero-inflation weight discrepancy),
#'     Category III manifold metrics (MMD, Generative Precision/Recall), and Category VII
#'     (Trajectory & Lineage Topology via the HVG-filtered dynverse workflow, evaluated for
#'     \code{scRNA-seq} only and skipped for \code{scATAC-seq}: HIM topology similarity,
#'     edge flip similarity, geodesic distance correlation, branch assignment F1,
#'     and Random Forest feature importance weighted correlation wcor_features).
#'   \item \strong{+ Cell Types (Optional):} Unlocks Category III subpopulation metrics
#'     (Silhouette ASW, CDI, Neighborhood Purity, ARI, NMI, and auxiliary Davies-Bouldin Index).
#'   \item \strong{+ Batch / Biological Replicate Annotations (Optional):} Unlocks Category IV
#'     (Technical Variation & Batch Mixing Dynamics: Batch Silhouette, PCR R2, CMS score,
#'     LISI/ISI, Cross-Batch F1) and Category V (Biological Signal & Differential Fidelity:
#'     Libra pseudobulk edgeR-QLF null P-value uniformity, accuracy, precision, recall,
#'     F1-score, and ROC-AUC for scRNA-seq DE or scATAC-seq DA).
#'   \item \strong{Full Metadata (+ Runtime & Memory):} Unlocks Category VIII scalability tracking
#'     (\code{core_seconds}, \code{memory_per_core_mib}) for the complete 56-measure suite.
#' }
#'
#' @param ref_data Reference count matrix (features x cells).
#' @param sim_data Simulated count matrix (features x cells).
#' @param compute_bivariate Logical, whether to compute 2D bivariate tests. Default is TRUE.
#' @param cell_types Cell type metadata details
#' @param cell_types_ref Cell type details of reference dataset
#' @param cell_types_sim Cell type details of simulated dataset
#' @param batch_info Batch details
#' @param batch_ref Batch details of reference dataset
#' @param batch_sim Batch details of simulated dataset
#' @param modality Name of modality
#' @param threads Number of threads for parallel computation. Default is 1.
#' @param memory_mb Optional numeric value of peak memory allocation (in \code{memory_unit}).
#' @param elapsed_time Optional numeric value of wall-clock elapsed time (in \code{time_unit}).
#' @param verbose Logical, whether to print execution messages. Default is TRUE.
#' @param time_unit Character string specifying the unit of \code{elapsed_time}:
#'   \code{"seconds"} (default, or \code{"s"}) or \code{"minutes"} (or \code{"min"}).
#' @param memory_unit Character string specifying the unit of \code{memory_mb}:
#'   \code{"MiB"} (default, binary mebibytes) or \code{"MB"} (decimal megabytes, converted to MiB).
#' @param sim_cores Optional positive integer specifying the number of CPU cores used
#'   by the simulator to generate \code{sim_data}. Default is NULL.
#' @param trajectory_num_trees Number of Random Forest trees for Category VII \code{wcor_features}
#'   (default \code{10000L}, per Saelens et al. 2019).
#' @param include_trajectory Logical; whether to evaluate Category VII trajectory & lineage topology
#'   metrics when \code{modality} is scRNA-seq (automatically skipped for \code{scATAC-seq}). Default is \code{TRUE}.
#'
#' @return A list containing tidy summary tables and detailed metric lists.
#' @examples
#' data(example_scrna, package = "scSimEval")
#' res <- evaluate_simulation_accuracy(example_scrna$ref, example_scrna$sim)
#' @export
evaluate_simulation_accuracy <- function(
  ref_data,
  sim_data,
  cell_types = NULL,
  cell_types_ref = NULL,
  cell_types_sim = NULL,
  batch_info = NULL,
  batch_ref = NULL,
  batch_sim = NULL,
  modality = NULL,
  compute_bivariate = TRUE,
  threads = 1,
  memory_mb = NULL,
  elapsed_time = NULL,
  verbose = TRUE,
  time_unit = c("seconds", "minutes"),
  memory_unit = c("MiB", "MB"),
  sim_cores = NULL,
  trajectory_num_trees = 10000L,
  include_trajectory = TRUE
) {
  # Strict standardization check (Seurat / SingleCellExperiment)
  val_ref <- validate_sc_object(ref_data, "ref_data")
  val_sim <- validate_sc_object(sim_data, "sim_data")
  
  ref_data <- val_ref$counts
  sim_data <- val_sim$counts
  
  if (is.null(cell_types_ref) && !is.null(cell_types)) cell_types_ref <- cell_types
  if (is.null(cell_types_sim) && !is.null(cell_types)) cell_types_sim <- cell_types
  if (is.null(cell_types_ref) && !is.null(val_ref$cell_types)) cell_types_ref <- val_ref$cell_types
  if (is.null(cell_types_sim) && !is.null(val_sim$cell_types)) cell_types_sim <- val_sim$cell_types
  
  if (is.null(batch_ref) && !is.null(batch_info)) batch_ref <- batch_info
  if (is.null(batch_sim) && !is.null(batch_info)) batch_sim <- batch_info
  if (is.null(batch_ref) && !is.null(val_ref$batch_info)) batch_ref <- val_ref$batch_info
  if (is.null(batch_sim) && !is.null(val_sim$batch_info)) batch_sim <- val_sim$batch_info
  
  # SimBench standard: filter out rare cell types with < 10 cells if cell types are provided
  if (!is.null(cell_types_ref) && length(cell_types_ref) == ncol(ref_data)) {
    tab_ref <- table(cell_types_ref)
    rare_types_ref <- names(tab_ref)[tab_ref < 10]
    if (length(rare_types_ref) > 0 && length(tab_ref) > length(rare_types_ref)) {
      if (verbose) message(sprintf("SimBench standard: skipping cell type(s) with < 10 cells in reference: %s", paste(rare_types_ref, collapse = ", ")))
      keep_ref <- which(!cell_types_ref %in% rare_types_ref)
      if (length(keep_ref) >= 10) {
        ref_data <- ref_data[, keep_ref, drop = FALSE]
        cell_types_ref <- cell_types_ref[keep_ref]
        if (!is.null(batch_ref) && length(batch_ref) > length(keep_ref)) {
          batch_ref <- batch_ref[keep_ref]
        }
      }
    }
  }
  if (!is.null(cell_types_sim) && length(cell_types_sim) == ncol(sim_data)) {
    tab_sim <- table(cell_types_sim)
    rare_types_sim <- names(tab_sim)[tab_sim < 10]
    if (length(rare_types_sim) > 0 && length(tab_sim) > length(rare_types_sim)) {
      if (verbose) message(sprintf("SimBench standard: skipping cell type(s) with < 10 cells in simulation: %s", paste(rare_types_sim, collapse = ", ")))
      keep_sim <- which(!cell_types_sim %in% rare_types_sim)
      if (length(keep_sim) >= 10) {
        sim_data <- sim_data[, keep_sim, drop = FALSE]
        cell_types_sim <- cell_types_sim[keep_sim]
        if (!is.null(batch_sim) && length(batch_sim) > length(keep_sim)) {
          batch_sim <- batch_sim[keep_sim]
        }
      }
    }
  }
  
  threads <- max(1L, as.integer(threads %||% 1L))
  if (requireNamespace("RcppParallel", quietly = TRUE)) {
    try(RcppParallel::setThreadOptions(numThreads = threads), silent = TRUE)
  }
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    try(RhpcBLASctl::blas_set_num_threads(threads), silent = TRUE)
    try(RhpcBLASctl::omp_set_num_threads(threads), silent = TRUE)
  }
  if (requireNamespace("data.table", quietly = TRUE)) {
    try(data.table::setDTthreads(threads), silent = TRUE)
  }

  # Resolve modality early so ATAC uses TF-IDF transformation (simPIC/Signac standard)
  is_atac <- .is_atac_modality(val_ref$object, modality) || .is_atac_modality(val_sim$object, modality)
  modality_resolved <- if (is_atac) "scATAC-seq" else "scRNA-seq"

  # Step 1: Extract properties
  if (verbose) message("[1/5] Extracting cell-level properties...")
  ref_cell <- extract_cell_properties(ref_data, modality = modality_resolved, verbose = FALSE)
  sim_cell <- extract_cell_properties(sim_data, modality = modality_resolved, verbose = FALSE)
  
  if (verbose) message("[2/5] Extracting feature-level properties...")
  ref_feat <- extract_feature_properties(ref_data, modality = modality_resolved, verbose = FALSE)
  sim_feat <- extract_feature_properties(sim_data, modality = modality_resolved, verbose = FALSE)
  
  # Step 2: Compute Univariate Accuracy Metrics (Cell & Feature Level) in Parallel
  if (verbose) message("[3/5] Computing univariate accuracy metrics...")
  
  # For scRNA-seq, include tmm_factor if available (skip for ATAC or if tmm is empty)
  has_tmm <- !is_atac && !is.null(ref_cell$tmm_factor) && length(ref_cell$tmm_factor) > 0 &&
             !all(is.na(ref_cell$tmm_factor))
  
  cell_props_to_eval <- if (has_tmm) {
    c("library_size", "log_library_size", "tmm_factor", "zero_fraction_cell", "detection_freq_cell", "cell_cor", "pca_dist")
  } else {
    c("library_size", "log_library_size", "zero_fraction_cell", "detection_freq_cell", "cell_cor", "pca_dist")
  }
  feat_props_to_eval <- c("mean_expression", "variance", "scaled_variance", "cv", "zero_fraction_feature", "dispersion", "gene_cor")
  cell_results <- list()
  feat_results <- list()
  
  all_props <- c(
    lapply(cell_props_to_eval, function(p) list(type = "cell", name = p, r = ref_cell[[p]], s = sim_cell[[p]])),
    lapply(feat_props_to_eval, function(p) list(type = "feat", name = p, r = ref_feat[[p]], s = sim_feat[[p]]))
  )
  all_props <- all_props[vapply(all_props, function(x) !is.null(x$r) && !is.null(x$s) && length(x$r) > 0, logical(1))]

  if (threads > 1 && length(all_props) > 1 && requireNamespace("BiocParallel", quietly = TRUE)) {
    bparam_u <- if (.Platform$OS.type == "windows") {
      BiocParallel::SnowParam(workers = min(threads, length(all_props)))
    } else {
      BiocParallel::MulticoreParam(workers = min(threads, length(all_props)))
    }
    prop_results <- tryCatch({
      suppressWarnings(BiocParallel::bplapply(all_props, function(item) {
        calc_all_univariate_metrics(item$r, item$s, item$name)
      }, BPPARAM = bparam_u))
    }, error = function(e) NULL)
    
    if (!is.null(prop_results)) {
      for (idx in seq_along(all_props)) {
        item <- all_props[[idx]]
        if (item$type == "cell") cell_results[[item$name]] <- prop_results[[idx]]
        else feat_results[[item$name]] <- prop_results[[idx]]
      }
    } else {
      for (item in all_props) {
        if (item$type == "cell") cell_results[[item$name]] <- calc_all_univariate_metrics(item$r, item$s, item$name)
        else feat_results[[item$name]] <- calc_all_univariate_metrics(item$r, item$s, item$name)
      }
    }
  } else {
    for (item in all_props) {
      if (item$type == "cell") cell_results[[item$name]] <- calc_all_univariate_metrics(item$r, item$s, item$name)
      else feat_results[[item$name]] <- calc_all_univariate_metrics(item$r, item$s, item$name)
    }
  }
  
  # Step 4: Compute Bivariate Accuracy Metrics in Parallel
  bivariate_results <- list()
  if (compute_bivariate) {
    if (verbose) message("[4/5] Computing bivariate metrics (Fasano-Franceschini, KDE zstat, 2D EMD)...")
    
    biv_configs <- list(
      lib_vs_cellzero = list(ref = cbind(log10(ref_cell$library_size + 1), ref_cell$zero_fraction_cell),
                             sim = cbind(log10(sim_cell$library_size + 1), sim_cell$zero_fraction_cell),
                             prefix = "lib_vs_cellzero"),
      mean_vs_sd = list(ref = cbind(ref_feat$mean_expression, ref_feat$sd),
                        sim = cbind(sim_feat$mean_expression, sim_feat$sd),
                        prefix = "mean_vs_sd"),
      mean_vs_featzero = list(ref = cbind(ref_feat$mean_expression, ref_feat$zero_fraction_feature),
                              sim = cbind(sim_feat$mean_expression, sim_feat$zero_fraction_feature),
                              prefix = "mean_vs_featzero"),
      mean_vs_dispersion = list(ref = cbind(ref_feat$mean_expression, ref_feat$dispersion),
                                sim = cbind(sim_feat$mean_expression, sim_feat$dispersion),
                                prefix = "mean_vs_dispersion")
    )

    if (threads > 1 && requireNamespace("BiocParallel", quietly = TRUE)) {
      bparam <- if (.Platform$OS.type == "windows") {
        BiocParallel::SnowParam(workers = min(threads, length(biv_configs)))
      } else {
        BiocParallel::MulticoreParam(workers = min(threads, length(biv_configs)))
      }
      biv_list <- tryCatch({
        suppressWarnings(BiocParallel::bplapply(biv_configs, function(cfg) {
          calc_all_bivariate_metrics(cfg$ref, cfg$sim, cfg$prefix, threads = 1)
        }, BPPARAM = bparam))
      }, error = function(e) NULL)
      if (!is.null(biv_list)) {
        bivariate_results <- stats::setNames(biv_list, names(biv_configs))
      } else {
        for (nm in names(biv_configs)) {
          cfg <- biv_configs[[nm]]
          bivariate_results[[nm]] <- calc_all_bivariate_metrics(cfg$ref, cfg$sim, cfg$prefix, threads = threads)
        }
      }
    } else {
      for (nm in names(biv_configs)) {
        cfg <- biv_configs[[nm]]
        bivariate_results[[nm]] <- calc_all_bivariate_metrics(cfg$ref, cfg$sim, cfg$prefix, threads = threads)
      }
    }
  }
  
  # Step 5: Dropout, Zero-Inflation, & Generative Manifold
  if (verbose) message("[5/5] Computing zero-probability and manifold distances...")
  zero_prob_eval <- if (!is_atac) evaluate_zero_probability_curve(ref_data, sim_data) else NULL
  acc_sparsity_eval <- if (is_atac) evaluate_accessibility_sparsity_curve(ref_data, sim_data) else NULL
  ref_excess_w <- calc_excess_zero_weights(ref_data)
  sim_excess_w <- calc_excess_zero_weights(sim_data)
  
  gen_pr <- calc_generative_precision_recall(ref_data, sim_data)
  mmd_eval <- calc_mmd(ref_data, sim_data)
  
  # Build tidy summary table
  summary_rows <- list()
  all_univ_metric_names <- c(
    "Wasserstein", "KS", "MAD", "Bhattacharyya", "Runs_Statistic",
    "NN_Mismatch", "Between_Dataset_Silh_Global",
    "MAE", "RMSE", "KDE_ISE"
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
      "Fasano_Franceschini_2D_KS",
      "EMD_2D", "NN_Mismatch_2D",
      "Between_Dataset_Silh_2D_Global"
    )
    for (b_name in names(bivariate_results)) {
      for (m in all_biv_metric_names) {
        key <- paste0(b_name, "_", m)
        if (!is.null(bivariate_results[[b_name]][[key]])) {
          summary_rows[[length(summary_rows) + 1]] <- data.frame(
            Category = "Distributional Properties",
            Property = b_name,
            Metric = m,
            Value = as.numeric(bivariate_results[[b_name]][[key]]),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  if (!is.null(zero_prob_eval)) {
    for (m in names(zero_prob_eval)) {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        Category = "Sparsity & Zero-Inflation Dynamics",
        Property = "dropout_curve",
        Metric = m,
        Value = as.numeric(zero_prob_eval[[m]]),
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (!is.null(acc_sparsity_eval)) {
    atac_cat2_metrics <- c(
      curve_rmse     = acc_sparsity_eval$curve_rmse,
      r2_discrepancy = acc_sparsity_eval$r2_discrepancy
    )
    for (m in names(atac_cat2_metrics)) {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        Category = "Sparsity & Zero-Inflation Dynamics",
        Property = "accessibility_sparsity_curve",
        Metric = m,
        Value = as.numeric(atac_cat2_metrics[[m]]),
        stringsAsFactors = FALSE
      )
    }
  }
  
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    Category = "Sparsity & Zero-Inflation Dynamics",
    Property = "excess_zero_weights",
    Metric = "excess_zero_weight_discrepancy",
    Value = abs(ref_excess_w$mean_excess_zero_weight - sim_excess_w$mean_excess_zero_weight),
    stringsAsFactors = FALSE
  )
  
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    Category = "Cellular Manifold & Clustering Concordance",
    Property = "manifold",
    Metric = "MMD",
    Value = as.numeric(mmd_eval),
    stringsAsFactors = FALSE
  )
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    Category = "Cellular Manifold & Clustering Concordance",
    Property = "manifold",
    Metric = "Generative_Precision",
    Value = as.numeric(gen_pr$generative_precision),
    stringsAsFactors = FALSE
  )
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    Category = "Cellular Manifold & Clustering Concordance",
    Property = "manifold",
    Metric = "Generative_Recall",
    Value = as.numeric(gen_pr$generative_recall),
    stringsAsFactors = FALSE
  )
  
  clu_eval_uni <- tryCatch({
    evaluate_clustering_metrics(
      data = sim_data,
      cell_types = cell_types_sim,
      ref_data = ref_data,
      ref_cell_types = cell_types_ref,
      batch_info = batch_sim,
      threads = threads
    )
  }, error = function(e) NULL)
  
  if (!is.null(clu_eval_uni)) {
    clu_uni_metrics <- c(
      "silhouette_sim", "cdi",
      "neighborhood_purity", "ari", "nmi"
    )
    for (cm in clu_uni_metrics) {
      if (!is.null(clu_eval_uni[[cm]]) && !is.na(clu_eval_uni[[cm]])) {
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          Category = "Cellular Manifold & Clustering Concordance",
          Property = "clustering",
          Metric = cm,
          Value = as.numeric(clu_eval_uni[[cm]]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # Step 5a: Technical Variation & Batch Mixing Dynamics (Category IV - if batch_sim / batch_info provided)
  b_eval_uni <- NULL
  if (!is.null(batch_sim)) {
    if (verbose) message("Evaluating batch integration & confounder mixing...")
    b_eval_uni <- tryCatch({
      evaluate_batch_metrics(
        data = sim_data,
        batch_info = batch_sim,
        cell_types = cell_types_sim,
        verbose = FALSE
      )
    }, error = function(e) NULL)
    
    if (!is.null(b_eval_uni)) {
      batch_metrics_to_record <- c(
        "batch_silhouette", "pcr_r2", "cms", "isi", "cross_batch_F1"
      )
      for (bm in batch_metrics_to_record) {
        if (!is.null(b_eval_uni[[bm]])) {
          summary_rows[[length(summary_rows) + 1]] <- data.frame(
            Category = "Technical Variation & Batch Mixing Dynamics",
            Property = "batch_mixing",
            Metric = bm,
            Value = as.numeric(b_eval_uni[[bm]]),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  # Step 5b: Biological Signal Fidelity (Category V - Libra Pseudobulk edgeR-QLF)
  # Strictly requires biological replicates (batch_ref and batch_sim / batch_info)
  # for BOTH scRNA-seq and scATAC-seq (without batch_info, Category V cannot run)
  signal_eval_uni <- NULL
  if (!is.null(cell_types_ref) && !is.null(cell_types_sim)) {
    if (!is.null(batch_ref) && !is.null(batch_sim)) {
      if (is_atac) {
        if (verbose) message("Evaluating scATAC-seq biological signal fidelity (Libra Pseudobulk edgeR-QLF)...")
        signal_eval_uni <- tryCatch({
          evaluate_atac_da_fidelity(
            ref_data = ref_data,
            sim_data = sim_data,
            ref_celltypes = cell_types_ref,
            sim_celltypes = cell_types_sim,
            batch_ref = batch_ref,
            batch_sim = batch_sim
          )
        }, error = function(e) {
          if (verbose) message("scATAC-seq DA evaluation note: ", conditionMessage(e))
          NULL
        })
        
        if (!is.null(signal_eval_uni) && nrow(signal_eval_uni$da_summary_table) > 0) {
          for (row_i in seq_len(nrow(signal_eval_uni$da_summary_table))) {
            summary_rows[[length(summary_rows) + 1]] <- data.frame(
              Category = "Biological Signal & Differential Fidelity",
              Property = "scATAC_differential_accessibility",
              Metric = signal_eval_uni$da_summary_table$Metric[row_i],
              Value = as.numeric(signal_eval_uni$da_summary_table$Value[row_i]),
              stringsAsFactors = FALSE
            )
          }
        }
      } else {
        if (verbose) message("Evaluating scRNA-seq biological signal fidelity (Libra Pseudobulk edgeR-QLF)...")
        signal_eval_uni <- tryCatch({
          evaluate_deg_fidelity(
            ref_data = ref_data,
            sim_data = sim_data,
            ref_celltypes = cell_types_ref,
            sim_celltypes = cell_types_sim,
            batch_ref = batch_ref,
            batch_sim = batch_sim
          )
        }, error = function(e) {
          if (verbose) message("scRNA-seq pseudobulk DE evaluation note: ", conditionMessage(e))
          NULL
        })
        
        if (!is.null(signal_eval_uni) && nrow(signal_eval_uni$deg_summary_table) > 0) {
          for (row_i in seq_len(nrow(signal_eval_uni$deg_summary_table))) {
            summary_rows[[length(summary_rows) + 1]] <- data.frame(
              Category = "Biological Signal & Differential Fidelity",
              Property = "scRNA_differential_expression",
              Metric = signal_eval_uni$deg_summary_table$Metric[row_i],
              Value = as.numeric(signal_eval_uni$deg_summary_table$Value[row_i]),
              stringsAsFactors = FALSE
            )
          }
        }
      }
    } else {
      if (verbose) {
        mod_name_str <- if (is_atac) "scATAC-seq" else "scRNA-seq"
        message(sprintf("Biological replicate annotations ('batch_info') not provided: skipping Category V pseudobulk analysis for %s.", mod_name_str))
      }
    }
  }
  
  # Step 6: Trajectory & Lineage Topology (Category VII - dynverse workflow, scRNA-seq only)
  traj_res <- NULL
  if (!is_atac && isTRUE(include_trajectory)) {
    if (verbose) message("[6/7] Evaluating trajectory and lineage topology (HVG + dynverse, scRNA-seq only)...")
    traj_res <- tryCatch({
      evaluate_trajectory_metrics(
        ref_data = ref_data,
        sim_data = sim_data,
        cell_types_ref = cell_types_ref,
        cell_types_sim = cell_types_sim,
        num_trees = trajectory_num_trees,
        threads = threads,
        verbose = FALSE
      )
    }, error = function(e) NULL)

    summary_rows <- c(summary_rows, .append_trajectory_rows(traj_res, modality = NULL))
  } else if (is_atac && verbose) {
    message("[6/7] Category VII (Trajectory & Lineage Topology) requires scRNA-seq data: skipping for scATAC-seq.")
  }
  
  # Step 7: Computational Scalability & Efficiency (Category VIII)
  time_unit_norm <- .normalize_time_unit(time_unit)
  memory_unit_norm <- .normalize_memory_unit(memory_unit)
  elapsed_seconds <- .convert_time_to_seconds(elapsed_time, time_unit_norm)
  memory_mib <- .convert_memory_to_mib(memory_mb, memory_unit_norm)
  sim_cores_val <- .normalize_sim_cores(sim_cores)
  cores_mult <- if (!is.null(sim_cores_val)) sim_cores_val else 1L
  core_sec_val <- if (!is.null(elapsed_seconds) && !is.na(elapsed_seconds)) as.numeric(elapsed_seconds * cores_mult) else NULL
  mem_core_val <- if (!is.null(memory_mib) && !is.na(memory_mib)) as.numeric(memory_mib / cores_mult) else NULL

  if (!is.null(core_sec_val) || !is.null(mem_core_val)) {
    if (verbose) message("[7/7] Logging computational scalability & efficiency metrics...")
    if (!is.null(core_sec_val)) {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        Category = "(VIII) Computational Scalability & Efficiency",
        Property = "resource_usage",
        Metric = "core_seconds",
        Value = as.numeric(core_sec_val),
        stringsAsFactors = FALSE
      )
    }
    if (!is.null(mem_core_val)) {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        Category = "(VIII) Computational Scalability & Efficiency",
        Property = "resource_usage",
        Metric = "memory_per_core_mib",
        Value = as.numeric(mem_core_val),
        stringsAsFactors = FALSE
      )
    }
  }

  summary_table <- do.call(rbind, summary_rows)
  if (!is.null(sim_cores_val) && !is.null(summary_table) && nrow(summary_table) > 0) {
    summary_table$Sim_Cores <- sim_cores_val
  }
  if (verbose) message("Unimodal accuracy evaluation complete.")
  
  list(
    metrics_summary_table = summary_table,
    cell_metrics = cell_results,
    feature_metrics = feat_results,
    bivariate_metrics = bivariate_results,
    zero_probability_curve = zero_prob_eval,
    accessibility_sparsity_curve = acc_sparsity_eval,
    excess_zero_weights = list(ref = ref_excess_w, sim = sim_excess_w),
    generative_manifold = list(
      mmd = mmd_eval,
      precision_recall = gen_pr
    ),
    ref_properties = list(cell = ref_cell, feature = ref_feat),
    sim_properties = list(cell = sim_cell, feature = sim_feat),
    batch_effects = b_eval_uni,
    biological_signal = signal_eval_uni,
    trajectory = traj_res,
    resource_usage = list(
      core_seconds = core_sec_val,
      memory_per_core_mib = mem_core_val,
      memory_per_core_gib = if (!is.null(mem_core_val)) mem_core_val / 1024 else NULL,
      memory_mb = memory_mib,
      memory_input = memory_mb,
      memory_unit = memory_unit_norm,
      elapsed_time = elapsed_seconds,
      elapsed_time_input = elapsed_time,
      time_unit = time_unit_norm,
      sim_cores = sim_cores_val
    )
  )
}

#' Master Multiomics Simulation Benchmarking Pipeline
#'
#' Orchestrates comprehensive benchmarking across all 8 evaluation categories
#' using strictly the user's provided inputs:
#' \itemize{
#'   \item \code{ref_multi}: List of reference omics objects/matrices (e.g. \code{list(RNA = ..., ATAC = ...)}).
#'   \item \code{sim_multi}: List of simulated omics objects/matrices (e.g. \code{list(RNA = ..., ATAC = ...)}).
#'   \item \code{cell_types}: Vector/factor of cell-type identities.
#'   \item \code{batch_info}: Vector/factor of batch or biological replicate annotations.
#'   \item Computational resource usage: \code{elapsed_time}, \code{memory_mb}.
#'   \item Trajectory features inferred directly from scRNA-seq counts.
#' }
#'
#' \if{html}{\figure{workflow_diagram.png}{options: width="800" alt="scSimEval Unified Benchmarking Workflow Diagram"}}
#'
#' @details
#' \strong{Multiomics Metric Measurement Capacity Across the 56 Evaluation Measures:}
#' \itemize{
#'   \item \strong{Paired or Unpaired Count Data Alone:} Evaluates unimodal Category I (10 univariate
#'     + 4 bivariate metrics when \code{compute_bivariate = TRUE}), Category II (5 sparsity/zero-inflation
#'     metrics across RNA & ATAC), Category III manifold metrics (MMD, Generative Precision/Recall),
#'     and Category VII trajectory metrics inferred from scRNA-seq via the HVG-filtered dynverse workflow
#'     (\code{him}, \code{edge_flip}, \code{geodesic_correlation}, \code{F1_branches}, \code{wcor_features}).
#'   \item \strong{+ Cell Types:} Unlocks Category III supervised clustering concordance (Silhouette ASW,
#'     CDI, Neighborhood Purity, ARI, NMI) and Category VI (Multiomics Integration & Alignment:
#'     Batch ASW, Cell-Type ASW, Graph Connectivity, Integration ARI, Integration NMI).
#'   \item \strong{+ Batch / Biological Replicate Annotations:} Unlocks Category IV batch mixing across
#'     both layers and Category V (12 differential fidelity metrics: null P-value uniformity, accuracy,
#'     precision, recall, F1-score, ROC-AUC across both scRNA-seq DE and scATAC-seq DA).
#'   \item \strong{+ Runtime & Memory:} Logs Category VIII (\code{core_seconds}, \code{memory_per_core_mib})
#'     to complete the full 56-measure evaluation suite across all 8 canonical categories.
#' }
#'
#' @param ref_multi Named list of reference matrices for Modality 1 and 2.
#' @param sim_multi Named list of simulated matrices for Modality 1 and 2.
#' @param cell_types Optional factor or vector of cell type labels for cells.
#' @param batch_info Optional factor or vector of batch labels for cells.
#' @param memory_mb Optional numeric value of peak memory allocation (in \code{memory_unit}).
#' @param elapsed_time Optional numeric value of wall-clock elapsed time (in \code{time_unit}).
#' @param feature_pairs Optional 2-column data.frame of linked feature pairs.
#' @param pairing Character string specifying whether the multiomics dataset is
#'   \code{"paired"} (simultaneous co-assay from the same individual cells, e.g.,
#'   10x Chromium Multiome, SHARE-seq, SNARE-seq) or \code{"unpaired"} (independent
#'   profiling of modalities from separate cells of the same biological tissue).
#'   Default is \code{"paired"}. In \code{"unpaired"} mode, cell-level pairing
#'   metrics (FOSCTTM, Match@1, cross-modal generation, and direct peak-to-gene
#'   coupling) are skipped, while population-level cross-modal metrics (cross-modal
#'   label transfer, co-expression module fidelity, peak co-accessibility, and
#'   accessibility profile concordance) are computed alongside all unimodal metrics.
#' @param is_paired Optional logical alias for \code{pairing}. If provided,
#'   \code{TRUE} corresponds to \code{"paired"} and \code{FALSE} to \code{"unpaired"}.
#' @param compute_bivariate Logical, whether to compute 2D bivariate tests. Default FALSE for speed.
#' @param threads Number of CPU threads. Default is 1.
#' @param verbose Logical, whether to print execution progress. Default is TRUE.
#' @param time_unit Character string specifying the unit of \code{elapsed_time}:
#'   \code{"seconds"} (default, or \code{"s"}) or \code{"minutes"} (or \code{"min"}).
#' @param memory_unit Character string specifying the unit of \code{memory_mb}:
#'   \code{"MiB"} (default, binary mebibytes) or \code{"MB"} (decimal megabytes, converted to MiB).
#' @param sim_cores Optional positive integer specifying the number of CPU cores used
#'   by the simulator to generate \code{sim_multi}. Default is NULL.
#' @param elapsed_time_rna,elapsed_time_atac Optional separate wall-clock runtimes for RNA and ATAC layers
#'   when modalities are simulated independently (used if combined \code{elapsed_time} is NULL).
#' @param memory_mb_rna,memory_mb_atac Optional separate peak RAM allocations for RNA and ATAC layers
#'   when modalities are simulated independently (used if combined \code{memory_mb} is NULL).
#' @param time_unit_rna,time_unit_atac Optional time units (\code{"seconds"} or \code{"minutes"}) for
#'   \code{elapsed_time_rna} and \code{elapsed_time_atac}; defaults to \code{time_unit}.
#' @param memory_unit_rna,memory_unit_atac Optional memory units (\code{"MiB"} or \code{"MB"}) for
#'   \code{memory_mb_rna} and \code{memory_mb_atac}; defaults to \code{memory_unit}.
#' @param sim_cores_rna,sim_cores_atac Optional positive integers specifying simulation CPU cores for
#'   separate RNA and ATAC runs (used if combined \code{sim_cores} is NULL).
#' @param trajectory_num_trees Number of Random Forest trees for Category VII \code{wcor_features}
#'   (default \code{10000L}, per Saelens et al. 2019).
#'
#' @return A list containing detailed results and a unified tidy master summary table.
#' @examples
#' data(example_multiomics, package = "scSimEval")
#' evaluate_multiomics_accuracy(example_multiomics$ref_multi, example_multiomics$sim_multi)
#' @export
evaluate_multiomics_accuracy <- function(
  ref_multi,
  sim_multi,
  cell_types = NULL,
  batch_info = NULL,
  pairing = c("paired", "unpaired"),
  is_paired = NULL,
  memory_mb = NULL,
  elapsed_time = NULL,
  feature_pairs = NULL,
  compute_bivariate = FALSE,
  threads = 1,
  verbose = TRUE,
  time_unit = c("seconds", "minutes"),
  memory_unit = c("MiB", "MB"),
  sim_cores = NULL,
  elapsed_time_rna = NULL,
  elapsed_time_atac = NULL,
  memory_mb_rna = NULL,
  memory_mb_atac = NULL,
  time_unit_rna = NULL,
  time_unit_atac = NULL,
  memory_unit_rna = NULL,
  memory_unit_atac = NULL,
  sim_cores_rna = NULL,
  sim_cores_atac = NULL,
  trajectory_num_trees = 10000L
) {
  # Strict standardization check (Seurat / SingleCellExperiment)
  val_ref <- validate_multiomics_input(ref_multi, "ref_multi")
  val_sim <- validate_multiomics_input(sim_multi, "sim_multi")
  
  if (is.null(cell_types)) {
    if (!is.null(val_ref$cell_types)) cell_types <- val_ref$cell_types
    else if (!is.null(val_sim$cell_types)) cell_types <- val_sim$cell_types
  }
  if (is.null(batch_info)) {
    if (!is.null(val_ref$batch_info)) batch_info <- val_ref$batch_info
    else if (!is.null(val_sim$batch_info)) batch_info <- val_sim$batch_info
  }
  
  ref_multi_orig <- ref_multi
  sim_multi_orig <- sim_multi

  mod_names <- names(ref_multi_orig)
  if (is.null(mod_names) || length(mod_names) < 2 || any(mod_names == "")) {
    if (!is.null(val_ref$m1_name) && !is.null(val_ref$m2_name)) {
      mod_names <- c(val_ref$m1_name, val_ref$m2_name)
    } else {
      mod_names <- c("Modality_1", "Modality_2")
    }
  }
  m1_name <- mod_names[1]
  m2_name <- mod_names[2]

  ref_multi <- stats::setNames(list(val_ref$m1_obj, val_ref$m2_obj), c(m1_name, m2_name))
  sim_multi <- stats::setNames(list(val_sim$m1_obj, val_sim$m2_obj), c(m1_name, m2_name))
  
  if (!is.null(is_paired)) {
    pairing <- if (isTRUE(is_paired)) "paired" else "unpaired"
  } else {
    pairing <- match.arg(pairing)
  }
  
  threads <- max(1L, as.integer(threads %||% 1L))
  if (requireNamespace("RcppParallel", quietly = TRUE)) {
    try(RcppParallel::setThreadOptions(numThreads = threads), silent = TRUE)
  }
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    try(RhpcBLASctl::blas_set_num_threads(threads), silent = TRUE)
    try(RhpcBLASctl::omp_set_num_threads(threads), silent = TRUE)
  }
  if (requireNamespace("data.table", quietly = TRUE)) {
    try(data.table::setDTthreads(threads), silent = TRUE)
  }
  
  # Modality names already resolved and assigned to ref_multi / sim_multi
  
  master_rows <- list()
  
  # ============================================================================
  # 1. Unimodal Evaluation for Modality 1 and Modality 2 (Categories 1 & 2)
  # ============================================================================
  if (verbose) message(sprintf("=== [Step 1] Evaluating %s Layer Accuracy ===", m1_name))
  m1_eval <- evaluate_simulation_accuracy(
    ref_data = ref_multi[[1]],
    sim_data = sim_multi[[1]],
    modality = m1_name,
    compute_bivariate = compute_bivariate,
    threads = threads,
    verbose = verbose,
    include_trajectory = FALSE
  )
  t1 <- m1_eval$metrics_summary_table
  t1 <- t1[t1$Property != "clustering" & !grepl("Trajectory", t1$Category, ignore.case = TRUE), , drop = FALSE]
  t1$Modality <- m1_name
  master_rows[[length(master_rows) + 1]] <- t1
  
  if (verbose) message(sprintf("=== [Step 2] Evaluating %s Layer Accuracy ===", m2_name))
  m2_eval <- evaluate_simulation_accuracy(
    ref_data = ref_multi[[2]],
    sim_data = sim_multi[[2]],
    modality = m2_name,
    compute_bivariate = compute_bivariate,
    threads = threads,
    verbose = verbose,
    include_trajectory = FALSE
  )
  t2 <- m2_eval$metrics_summary_table
  t2 <- t2[t2$Property != "clustering" & !grepl("Trajectory", t2$Category, ignore.case = TRUE), , drop = FALSE]
  t2$Modality <- m2_name
  master_rows[[length(master_rows) + 1]] <- t2
  
  # ============================================================================
  # 2. Cellular Structure & Clustering Accuracy (Category 3)
  # ============================================================================
  # Extract modality specific cell types and batches (handles vector, factor, data.frame, and named list inputs)
  extract_label_vector <- function(x, n_expected, mod_name, fallback_idx = 1) {
    if (is.null(x)) return(NULL)
    if (is.data.frame(x)) {
      v <- if (!is.null(x$CellType)) x$CellType else if (!is.null(x$Batch)) x$Batch else if (ncol(x) >= 2) x[[2]] else x[[1]]
      if (length(v) == n_expected) return(as.character(v))
    }
    if (is.list(x) && !is.data.frame(x) && !is.factor(x)) {
      if (!is.null(x[[mod_name]])) return(extract_label_vector(x[[mod_name]], n_expected, mod_name, fallback_idx))
      m_idx <- match(tolower(mod_name), tolower(names(x)))
      if (!is.na(m_idx)) return(extract_label_vector(x[[m_idx]], n_expected, mod_name, fallback_idx))
      if (length(x) >= fallback_idx) return(extract_label_vector(x[[fallback_idx]], n_expected, mod_name, fallback_idx))
      return(extract_label_vector(x[[1]], n_expected, mod_name, fallback_idx))
    }
    if (length(x) == n_expected) return(as.character(x))
    return(NULL)
  }
  ct_m1 <- extract_label_vector(cell_types, ncol(ref_multi[[1]]), m1_name, fallback_idx = 1)
  ct_m2 <- extract_label_vector(cell_types, ncol(ref_multi[[2]]), m2_name, fallback_idx = 2)
  bt_m1 <- extract_label_vector(batch_info, ncol(ref_multi[[1]]), m1_name, fallback_idx = 1)
  bt_m2 <- extract_label_vector(batch_info, ncol(ref_multi[[2]]), m2_name, fallback_idx = 2)

  if (verbose) message(sprintf("=== [Step 3] Evaluating Cellular Structure & Clustering Preservation (%s & %s) ===", m1_name, m2_name))
  clu_eval_m1 <- evaluate_clustering_metrics(
    data = sim_multi[[1]],
    cell_types = ct_m1,
    ref_data = ref_multi[[1]],
    ref_cell_types = ct_m1,
    batch_info = bt_m1,
    threads = threads
  )
  clu_eval_m2 <- evaluate_clustering_metrics(
    data = sim_multi[[2]],
    cell_types = ct_m2,
    ref_data = ref_multi[[2]],
    ref_cell_types = ct_m2,
    batch_info = bt_m2,
    threads = threads
  )
  clustering_res <- list(m1 = clu_eval_m1, m2 = clu_eval_m2)
  
  clu_metrics_to_record <- c(
    "silhouette_sim", "cdi",
    "neighborhood_purity", "ari", "nmi"
  )
  for (cm in clu_metrics_to_record) {
    if (!is.null(clu_eval_m1[[cm]]) && !is.na(clu_eval_m1[[cm]])) {
      master_rows[[length(master_rows) + 1]] <- data.frame(
        Category = "Cellular Manifold & Clustering Concordance",
        Property = "clustering",
        Metric = cm,
        Value = as.numeric(clu_eval_m1[[cm]]),
        Modality = m1_name,
        stringsAsFactors = FALSE
      )
    }
    if (!is.null(clu_eval_m2[[cm]]) && !is.na(clu_eval_m2[[cm]])) {
      master_rows[[length(master_rows) + 1]] <- data.frame(
        Category = "Cellular Manifold & Clustering Concordance",
        Property = "clustering",
        Metric = cm,
        Value = as.numeric(clu_eval_m2[[cm]]),
        Modality = m2_name,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # ============================================================================
  # 3. Batch Effect & Technical Confounder Mixing (Category 4)
  # ============================================================================
  batch_res <- list()
  if (!is.null(bt_m1)) {
    if (verbose) message("=== [Step 4] Evaluating Batch Integration & Technical Confounder Mixing ===")
    b_eval <- evaluate_batch_metrics(
      data = sim_multi[[1]],
      batch_info = bt_m1,
      cell_types = ct_m1,
      verbose = FALSE
    )
    batch_res[[m1_name]] <- b_eval
    
    batch_metrics_to_record <- c(
      "batch_silhouette", "pcr_r2", "cms", "isi",
      "cross_batch_F1"
    )
    for (bm in batch_metrics_to_record) {
      if (!is.null(b_eval[[bm]])) {
        master_rows[[length(master_rows) + 1]] <- data.frame(
          Category = "Technical Variation & Batch Mixing Dynamics",
          Property = "batch_mixing",
          Metric = bm,
          Value = as.numeric(b_eval[[bm]]),
          Modality = m1_name,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!is.null(bt_m2)) {
    b_eval_m2 <- tryCatch({
      evaluate_batch_metrics(
        data = sim_multi[[2]],
        batch_info = bt_m2,
        cell_types = ct_m2,
        verbose = FALSE
      )
    }, error = function(e) NULL)
    if (!is.null(b_eval_m2)) {
      batch_res[[m2_name]] <- b_eval_m2
      batch_metrics_to_record <- c(
        "batch_silhouette", "pcr_r2", "cms", "isi",
        "cross_batch_F1"
      )
      for (bm in batch_metrics_to_record) {
        if (!is.null(b_eval_m2[[bm]])) {
          master_rows[[length(master_rows) + 1]] <- data.frame(
            Category = "Technical Variation & Batch Mixing Dynamics",
            Property = "batch_mixing",
            Metric = bm,
            Value = as.numeric(b_eval_m2[[bm]]),
            Modality = m2_name,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  # ============================================================================
  # 4. Biological Signal & Downstream Performance (Category 5)
  #    - scRNA-seq (Modality 1): Libra Pseudobulk edgeR-QLF (Requires batch_info)
  #    - scATAC-seq (Modality 2): Libra Pseudobulk edgeR-QLF (Requires batch_info)
  # ============================================================================
  signal_res <- list()
  if (!is.null(ct_m1)) {
    if (!is.null(bt_m1)) {
      if (verbose) message(sprintf("=== [Step 5] Evaluating %s Biological Signal Fidelity (Libra Pseudobulk edgeR-QLF) ===", m1_name))
      sim_ct_m1 <- if (length(ct_m1) == ncol(sim_multi[[1]])) ct_m1 else rep(ct_m1, length.out = ncol(sim_multi[[1]]))
      sim_bt_m1 <- if (length(bt_m1) == ncol(sim_multi[[1]])) bt_m1 else rep(bt_m1, length.out = ncol(sim_multi[[1]]))
      
      deg_eval <- tryCatch({
        evaluate_deg_fidelity(
          ref_data = ref_multi[[1]],
          sim_data = sim_multi[[1]],
          ref_celltypes = ct_m1,
          sim_celltypes = sim_ct_m1,
          batch_ref = bt_m1,
          batch_sim = sim_bt_m1
        )
      }, error = function(e) {
        if (verbose) message("scRNA-seq pseudobulk DE evaluation note: ", conditionMessage(e))
        NULL
      })
      
      if (!is.null(deg_eval) && nrow(deg_eval$deg_summary_table) > 0) {
        signal_res[[m1_name]] <- deg_eval
        for (row_i in seq_len(nrow(deg_eval$deg_summary_table))) {
          master_rows[[length(master_rows) + 1]] <- data.frame(
            Category = "Biological Signal & Differential Fidelity",
            Property = "scRNA_differential_expression",
            Metric = deg_eval$deg_summary_table$Metric[row_i],
            Value = as.numeric(deg_eval$deg_summary_table$Value[row_i]),
            Modality = m1_name,
            stringsAsFactors = FALSE
          )
        }
      }
    } else {
      if (verbose) message(sprintf("Biological replicate annotations ('batch_info') not provided: skipping Category V pseudobulk analysis for %s.", m1_name))
    }
  }

  # Modality 2 (scATAC-seq) Biological Signal Fidelity via Libra Pseudobulk edgeR-QLF
  # Strictly requires biological replicates (bt_m2 or batch_info)
  is_atac_m2 <- grepl("atac|peak|chromatin", tolower(m2_name)) ||
    (!is.null(rownames(ref_multi[[2]])) && length(rownames(ref_multi[[2]])) > 0 &&
     mean(grepl("^(peak|chr[0-9xyXYM]+[:-])", rownames(ref_multi[[2]]), ignore.case = TRUE)) > 0.5)

  if (is_atac_m2) {
    if (!is.null(ct_m2) && !is.null(bt_m2)) {
      if (verbose) message(sprintf("=== [Step 5] Evaluating %s Biological Signal Fidelity (Libra Pseudobulk edgeR-QLF) ===", m2_name))
      sim_ct_m2 <- if (length(ct_m2) == ncol(sim_multi[[2]])) ct_m2 else rep(ct_m2, length.out = ncol(sim_multi[[2]]))
      sim_bt_m2 <- if (length(bt_m2) == ncol(sim_multi[[2]])) bt_m2 else rep(bt_m2, length.out = ncol(sim_multi[[2]]))
      
      da_eval_m2 <- tryCatch({
        evaluate_atac_da_fidelity(
          ref_data = ref_multi[[2]],
          sim_data = sim_multi[[2]],
          ref_celltypes = ct_m2,
          sim_celltypes = sim_ct_m2,
          batch_ref = bt_m2,
          batch_sim = sim_bt_m2
        )
      }, error = function(e) {
        if (verbose) message("scATAC-seq DA evaluation note in multiomics: ", conditionMessage(e))
        NULL
      })
      
      if (!is.null(da_eval_m2)) {
        signal_res[[m2_name]] <- da_eval_m2
        if (nrow(da_eval_m2$da_summary_table) > 0) {
          for (row_i in seq_len(nrow(da_eval_m2$da_summary_table))) {
            master_rows[[length(master_rows) + 1]] <- data.frame(
              Category = "Biological Signal & Differential Fidelity",
              Property = "scATAC_differential_accessibility",
              Metric = da_eval_m2$da_summary_table$Metric[row_i],
              Value = as.numeric(da_eval_m2$da_summary_table$Value[row_i]),
              Modality = m2_name,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    } else {
      if (verbose) message("Biological replicate annotations ('batch_info') not provided: skipping Category V pseudobulk DA analysis for scATAC-seq modality.")
    }
  }
  
  # ============================================================================
  # 5. Trajectory & Lineage Topology (Category VII - dynverse workflow; scRNA-seq only)
  # ============================================================================
  m1_is_atac <- .is_atac_modality(ref_multi[[1]], m1_name) || .is_atac_modality(sim_multi[[1]], m1_name)
  m2_is_atac <- .is_atac_modality(ref_multi[[2]], m2_name) || .is_atac_modality(sim_multi[[2]], m2_name)
  rna_idx <- if (!m1_is_atac) 1L else if (!m2_is_atac) 2L else NA_integer_

  traj_eval <- NULL
  if (!is.na(rna_idx)) {
    rna_mod_name <- if (rna_idx == 1L) m1_name else m2_name
    rna_ct <- if (rna_idx == 1L) ct_m1 else ct_m2
    if (verbose) message(sprintf("=== [Step 6] Evaluating Trajectory & Lineage Topology Inferred from %s (HVG + dynverse) ===", rna_mod_name))
    sim_ct_traj <- if (!is.null(rna_ct)) (if (length(rna_ct) == ncol(sim_multi[[rna_idx]])) rna_ct else rep(rna_ct, length.out = ncol(sim_multi[[rna_idx]]))) else NULL
    traj_eval <- tryCatch({
      evaluate_trajectory_metrics(
        ref_data = ref_multi[[rna_idx]],
        sim_data = sim_multi[[rna_idx]],
        cell_types_ref = rna_ct,
        cell_types_sim = sim_ct_traj,
        num_trees = trajectory_num_trees,
        threads = threads,
        verbose = FALSE
      )
    }, error = function(e) NULL)
    master_rows <- c(master_rows, .append_trajectory_rows(traj_eval, modality = rna_mod_name))
  } else if (verbose) {
    message("=== [Step 6] Skipping Category VII Trajectory & Lineage Topology (requires scRNA-seq modality, not scATAC-seq) ===")
  }
  
  # ============================================================================
  # 6. Multiomics Integration & Alignment (Category 6)
  #    Multiomics Integration & Cross-Modal Alignment Fidelity
  #    (scib, scMultiSim, scCross, Matilda, scMoMtF)
  # ============================================================================
  if (verbose) message("=== [Step 7] Evaluating Multiomics Integration & Alignment (Category VI) ===")
  int_eval <- tryCatch({
    evaluate_multiomics_integration(
      ref_multi = ref_multi_orig,
      sim_multi = sim_multi_orig,
      cell_types = ct_m1,
      batch_info = bt_m1,
      n_pcs = 15,
      threads = threads,
      verbose = verbose
    )
  }, error = function(e) {
    if (verbose) message("Cross-modal integration evaluation note: ", conditionMessage(e))
    NULL
  })
  
  if (!is.null(int_eval) && !is.null(int_eval$integration_summary_table) && nrow(int_eval$integration_summary_table) > 0) {
    for (row_i in seq_len(nrow(int_eval$integration_summary_table))) {
      master_rows[[length(master_rows) + 1]] <- data.frame(
        Category = "(VI) Multiomics Integration & Alignment",
        Property = "multiomics_integration",
        Metric = int_eval$integration_summary_table$Metric[row_i],
        Value = as.numeric(int_eval$integration_summary_table$Value[row_i]),
        Modality = "Joint (Cross-Modal)",
        stringsAsFactors = FALSE
      )
    }
  }
  cross_cor_res <- int_eval
  fos_res <- NULL
  gen_fid <- NULL
  coupling_res <- NULL
  peak_coacc <- NULL
  coexpr_fid <- NULL
  acc_prof <- NULL
  
  # ============================================================================
  # 7. Computational Scalability & Efficiency (Category VIII)
  # ============================================================================
  time_unit_norm <- .normalize_time_unit(time_unit)
  memory_unit_norm <- .normalize_memory_unit(memory_unit)
  t_u_rna <- if (!is.null(time_unit_rna)) .normalize_time_unit(time_unit_rna) else time_unit_norm
  t_u_atac <- if (!is.null(time_unit_atac)) .normalize_time_unit(time_unit_atac) else time_unit_norm
  m_u_rna <- if (!is.null(memory_unit_rna)) .normalize_memory_unit(memory_unit_rna) else memory_unit_norm
  m_u_atac <- if (!is.null(memory_unit_atac)) .normalize_memory_unit(memory_unit_atac) else memory_unit_norm

  elapsed_rna_sec  <- .convert_time_to_seconds(elapsed_time_rna, t_u_rna)
  elapsed_atac_sec <- .convert_time_to_seconds(elapsed_time_atac, t_u_atac)
  mem_rna_mib      <- .convert_memory_to_mib(memory_mb_rna, m_u_rna)
  mem_atac_mib     <- .convert_memory_to_mib(memory_mb_atac, m_u_atac)
  cores_rna_val    <- .normalize_sim_cores(sim_cores_rna)
  cores_atac_val   <- .normalize_sim_cores(sim_cores_atac)

  has_sep_time <- (!is.null(elapsed_rna_sec) && !is.na(elapsed_rna_sec)) ||
                  (!is.null(elapsed_atac_sec) && !is.na(elapsed_atac_sec))
  has_sep_mem  <- (!is.null(mem_rna_mib) && !is.na(mem_rna_mib)) ||
                  (!is.null(mem_atac_mib) && !is.na(mem_atac_mib))

  elapsed_seconds <- .convert_time_to_seconds(elapsed_time, time_unit_norm)
  if ((is.null(elapsed_seconds) || is.na(elapsed_seconds)) && has_sep_time) {
    t_parts <- c(if (!is.null(elapsed_rna_sec)) elapsed_rna_sec else 0,
                 if (!is.null(elapsed_atac_sec)) elapsed_atac_sec else 0)
    t_parts <- t_parts[!is.na(t_parts)]
    if (length(t_parts) > 0) elapsed_seconds <- sum(t_parts)
  }

  memory_mib <- .convert_memory_to_mib(memory_mb, memory_unit_norm)
  if ((is.null(memory_mib) || is.na(memory_mib)) && has_sep_mem) {
    m_parts <- c(if (!is.null(mem_rna_mib)) mem_rna_mib else NA_real_,
                 if (!is.null(mem_atac_mib)) mem_atac_mib else NA_real_)
    m_parts <- m_parts[!is.na(m_parts)]
    if (length(m_parts) > 0) memory_mib <- max(m_parts)
  }

  sim_cores_val <- .normalize_sim_cores(sim_cores)
  if (is.null(sim_cores_val) && (!is.null(cores_rna_val) || !is.null(cores_atac_val))) {
    c_parts <- c(cores_rna_val, cores_atac_val)
    c_parts <- c_parts[!is.na(c_parts)]
    if (length(c_parts) > 0) sim_cores_val <- max(c_parts)
  }
  scal_input_mode <- if (has_sep_time || has_sep_mem) "separate" else "combined"

  # Core-Seconds and Peak RAM per Core calculations
  cores_comb_mult <- if (!is.null(sim_cores_val)) sim_cores_val else 1L
  c_rna <- if (!is.null(cores_rna_val)) cores_rna_val else cores_comb_mult
  c_atac <- if (!is.null(cores_atac_val)) cores_atac_val else cores_comb_mult

  if (has_sep_time) {
    cs_rna <- if (!is.null(elapsed_rna_sec) && !is.na(elapsed_rna_sec)) elapsed_rna_sec * c_rna else 0
    cs_atac <- if (!is.null(elapsed_atac_sec) && !is.na(elapsed_atac_sec)) elapsed_atac_sec * c_atac else 0
    core_sec_total <- cs_rna + cs_atac
  } else if (!is.null(elapsed_seconds) && !is.na(elapsed_seconds)) {
    core_sec_total <- elapsed_seconds * cores_comb_mult
    cs_rna <- if (!is.null(elapsed_rna_sec)) elapsed_rna_sec * c_rna else NULL
    cs_atac <- if (!is.null(elapsed_atac_sec)) elapsed_atac_sec * c_atac else NULL
  } else {
    core_sec_total <- NULL
    cs_rna <- NULL
    cs_atac <- NULL
  }

  if (has_sep_mem) {
    mpc_rna <- if (!is.null(mem_rna_mib) && !is.na(mem_rna_mib)) mem_rna_mib / c_rna else NA_real_
    mpc_atac <- if (!is.null(mem_atac_mib) && !is.na(mem_atac_mib)) mem_atac_mib / c_atac else NA_real_
    mpc_parts <- c(mpc_rna, mpc_atac)
    mpc_parts <- mpc_parts[!is.na(mpc_parts)]
    mem_core_total <- if (length(mpc_parts) > 0) max(mpc_parts) else NA_real_
  } else if (!is.null(memory_mib) && !is.na(memory_mib)) {
    mem_core_total <- memory_mib / cores_comb_mult
    mpc_rna <- if (!is.null(mem_rna_mib)) mem_rna_mib / c_rna else NULL
    mpc_atac <- if (!is.null(mem_atac_mib)) mem_atac_mib / c_atac else NULL
  } else {
    mem_core_total <- NULL
    mpc_rna <- NULL
    mpc_atac <- NULL
  }

  if (!is.null(core_sec_total) || !is.null(mem_core_total)) {
    if (verbose) message("=== [Step 8] Logging Computational Scalability & Efficiency ===")
    if (!is.null(core_sec_total)) {
      master_rows[[length(master_rows) + 1]] <- data.frame(
        Category = "(VIII) Computational Scalability & Efficiency",
        Property = "resource_usage",
        Metric = "core_seconds",
        Value = as.numeric(core_sec_total),
        Modality = "Computation",
        stringsAsFactors = FALSE
      )
    }
    if (!is.null(mem_core_total)) {
      master_rows[[length(master_rows) + 1]] <- data.frame(
        Category = "(VIII) Computational Scalability & Efficiency",
        Property = "resource_usage",
        Metric = "memory_per_core_mib",
        Value = as.numeric(mem_core_total),
        Modality = "Computation",
        stringsAsFactors = FALSE
      )
    }
  }
  
  master_table <- do.call(rbind, master_rows)
  if (!is.null(sim_cores_val) && !is.null(master_table) && nrow(master_table) > 0) {
    master_table$Sim_Cores <- sim_cores_val
  }
  if (verbose) message("Master multiomics benchmarking complete.")
  
  list(
    benchmark_summary_table = master_table,
    metrics_summary_table   = master_table,
    summary_table           = master_table,
    pairing = pairing,
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
    ),
    resource_usage = list(
      mode = scal_input_mode,
      core_seconds = core_sec_total,
      memory_per_core_mib = mem_core_total,
      memory_per_core_gib = if (!is.null(mem_core_total)) mem_core_total / 1024 else NULL,
      memory_mb = memory_mib,
      memory_input = memory_mb,
      memory_unit = memory_unit_norm,
      elapsed_time = elapsed_seconds,
      elapsed_time_input = elapsed_time,
      time_unit = time_unit_norm,
      sim_cores = sim_cores_val,
      rna = list(
        core_seconds = cs_rna,
        memory_per_core_mib = mpc_rna,
        elapsed_time = elapsed_rna_sec,
        elapsed_time_input = elapsed_time_rna,
        time_unit = t_u_rna,
        memory_mb = mem_rna_mib,
        memory_input = memory_mb_rna,
        memory_unit = m_u_rna,
        sim_cores = cores_rna_val
      ),
      atac = list(
        core_seconds = cs_atac,
        memory_per_core_mib = mpc_atac,
        elapsed_time = elapsed_atac_sec,
        elapsed_time_input = elapsed_time_atac,
        time_unit = t_u_atac,
        memory_mb = mem_atac_mib,
        memory_input = memory_mb_atac,
        memory_unit = m_u_atac,
        sim_cores = cores_atac_val
      )
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
#' @examples
#' data(example_scrna, package = "scSimEval")
#' res <- evaluate_simulation_accuracy(example_scrna$ref, example_scrna$sim)
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
            memory_mb = item$memory_mb,
            elapsed_time = item$elapsed_time,
            time_unit = if (!is.null(item$time_unit)) item$time_unit else "seconds",
            memory_unit = if (!is.null(item$memory_unit)) item$memory_unit else "MiB",
            sim_cores = if (!is.null(item$sim_cores)) item$sim_cores else if (!is.null(item$cores)) item$cores else item$n_cores,
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
            batch_info = item$batch_info,
            memory_mb = item$memory_mb,
            elapsed_time = item$elapsed_time,
            time_unit = if (!is.null(item$time_unit)) item$time_unit else "seconds",
            memory_unit = if (!is.null(item$memory_unit)) item$memory_unit else "MiB",
            sim_cores = if (!is.null(item$sim_cores)) item$sim_cores else if (!is.null(item$cores)) item$cores else item$n_cores
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
            pairing = if (!is.null(sub_items[[rna_idx]]$pairing)) sub_items[[rna_idx]]$pairing else if (!is.null(sub_items[[rna_idx]]$is_paired)) (if (isTRUE(sub_items[[rna_idx]]$is_paired)) "paired" else "unpaired") else "paired",
            memory_mb_rna = sub_items[[rna_idx]]$memory_mb,
            elapsed_time_rna = sub_items[[rna_idx]]$elapsed_time,
            time_unit_rna = if (!is.null(sub_items[[rna_idx]]$time_unit)) sub_items[[rna_idx]]$time_unit else "seconds",
            memory_unit_rna = if (!is.null(sub_items[[rna_idx]]$memory_unit)) sub_items[[rna_idx]]$memory_unit else "MiB",
            sim_cores_rna = if (!is.null(sub_items[[rna_idx]]$sim_cores)) sub_items[[rna_idx]]$sim_cores else if (!is.null(sub_items[[rna_idx]]$cores)) sub_items[[rna_idx]]$cores else sub_items[[rna_idx]]$n_cores,
            memory_mb_atac = sub_items[[atac_idx]]$memory_mb,
            elapsed_time_atac = sub_items[[atac_idx]]$elapsed_time,
            time_unit_atac = if (!is.null(sub_items[[atac_idx]]$time_unit)) sub_items[[atac_idx]]$time_unit else "seconds",
            memory_unit_atac = if (!is.null(sub_items[[atac_idx]]$memory_unit)) sub_items[[atac_idx]]$memory_unit else "MiB",
            sim_cores_atac = if (!is.null(sub_items[[atac_idx]]$sim_cores)) sub_items[[atac_idx]]$sim_cores else if (!is.null(sub_items[[atac_idx]]$cores)) sub_items[[atac_idx]]$cores else sub_items[[atac_idx]]$n_cores,
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
              batch_info = it$batch_info,
              memory_mb = it$memory_mb,
              elapsed_time = it$elapsed_time,
              time_unit = if (!is.null(it$time_unit)) it$time_unit else "seconds",
              memory_unit = if (!is.null(it$memory_unit)) it$memory_unit else "MiB",
              sim_cores = if (!is.null(it$sim_cores)) it$sim_cores else if (!is.null(it$cores)) it$cores else it$n_cores
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
          pairing = if (!is.null(item$pairing)) item$pairing else if (!is.null(item$is_paired)) (if (isTRUE(item$is_paired)) "paired" else "unpaired") else "paired",
          memory_mb = item$memory_mb,
          elapsed_time = item$elapsed_time,
          time_unit = if (!is.null(item$time_unit)) item$time_unit else "seconds",
          memory_unit = if (!is.null(item$memory_unit)) item$memory_unit else "MiB",
          sim_cores = if (!is.null(item$sim_cores)) item$sim_cores else if (!is.null(item$cores)) item$cores else item$n_cores,
          elapsed_time_rna = item$elapsed_time_rna,
          elapsed_time_atac = item$elapsed_time_atac,
          memory_mb_rna = item$memory_mb_rna,
          memory_mb_atac = item$memory_mb_atac,
          time_unit_rna = item$time_unit_rna,
          time_unit_atac = item$time_unit_atac,
          memory_unit_rna = item$memory_unit_rna,
          memory_unit_atac = item$memory_unit_atac,
          sim_cores_rna = item$sim_cores_rna,
          sim_cores_atac = item$sim_cores_atac,
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
          batch_info = item$batch_info,
          memory_mb = item$memory_mb,
          elapsed_time = item$elapsed_time,
          time_unit = if (!is.null(item$time_unit)) item$time_unit else "seconds",
          memory_unit = if (!is.null(item$memory_unit)) item$memory_unit else "MiB",
          sim_cores = if (!is.null(item$sim_cores)) item$sim_cores else if (!is.null(item$cores)) item$cores else item$n_cores
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
        pairing = if (!is.null(exp$pairing)) exp$pairing else "paired",
        memory_mb = exp$memory_mb,
        elapsed_time = exp$elapsed_time,
        feature_pairs = exp$feature_pairs,
        compute_bivariate = compute_bivariate,
        threads = threads,
        verbose = verbose,
        time_unit = if (!is.null(exp$time_unit)) exp$time_unit else "seconds",
        memory_unit = if (!is.null(exp$memory_unit)) exp$memory_unit else "MiB",
        sim_cores = exp$sim_cores,
        elapsed_time_rna = exp$elapsed_time_rna,
        elapsed_time_atac = exp$elapsed_time_atac,
        memory_mb_rna = exp$memory_mb_rna,
        memory_mb_atac = exp$memory_mb_atac,
        time_unit_rna = exp$time_unit_rna,
        time_unit_atac = exp$time_unit_atac,
        memory_unit_rna = exp$memory_unit_rna,
        memory_unit_atac = exp$memory_unit_atac,
        sim_cores_rna = exp$sim_cores_rna,
        sim_cores_atac = exp$sim_cores_atac
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
      keep_cols <- c("Dataset", "Data_Name", "Modality", "Category", "Property", "Metric", "Value")
      if ("Sim_Cores" %in% colnames(df)) keep_cols <- c(keep_cols, "Sim_Cores")
      df <- df[, keep_cols, drop = FALSE]
      all_master_rows[[length(all_master_rows) + 1]] <- df
      dataset_results[[exp_name]] <- res

    } else {
      res <- evaluate_simulation_accuracy(
        ref_data = exp$ref,
        sim_data = exp$sim,
        cell_types = exp$cell_types,
        compute_bivariate = compute_bivariate,
        threads = threads,
        memory_mb = exp$memory_mb,
        elapsed_time = exp$elapsed_time,
        verbose = verbose,
        time_unit = if (!is.null(exp$time_unit)) exp$time_unit else "seconds",
        memory_unit = if (!is.null(exp$memory_unit)) exp$memory_unit else "MiB",
        sim_cores = exp$sim_cores
      )

      df <- res$metrics_summary_table
      df$Dataset <- exp_name
      df$Data_Name <- exp$data_name
      df$Modality <- exp$modality
      keep_cols <- c("Dataset", "Data_Name", "Modality", "Category", "Property", "Metric", "Value")
      if ("Sim_Cores" %in% colnames(df)) keep_cols <- c(keep_cols, "Sim_Cores")
      df <- df[, keep_cols, drop = FALSE]
      all_master_rows[[length(all_master_rows) + 1]] <- df
      dataset_results[[exp$data_name]] <- res
    }
  }

  has_any_cores <- any(vapply(all_master_rows, function(d) "Sim_Cores" %in% colnames(d), logical(1)))
  if (has_any_cores) {
    all_master_rows <- lapply(all_master_rows, function(d) {
      if (!"Sim_Cores" %in% colnames(d)) d$Sim_Cores <- NA_integer_
      d
    })
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

#' Benchmark Simulation Methods Across Multiple Single-Cell Datasets
#'
#' A unified, high-performance, and beginner-friendly benchmarking runner designed
#' for multi-dataset and multi-simulator studies (e.g., 16 datasets x 6 simulators).
#' Automatically manages memory allocation with sparse matrix coercions,
#' parallel execution across datasets, progress tracking, and produces
#' both granular and cross-dataset aggregated benchmark tables ready for visualization.
#'
#' @param datasets A named list of empirical reference datasets. Each item can be a
#'   single-cell count matrix (features x cells), a path to an RDS/CSV matrix file,
#'   or a list with \code{ref}, optional \code{cell_types}, and optional \code{batch_info}.
#' @param simulators A named list of simulation methods. Each element must be a named
#'   list of simulated count matrices (or matrix file paths) corresponding to the
#'   datasets in \code{datasets}.
#'   Example: \code{list(Splatter = list(D1 = sim_splat_d1, D2 = sim_splat_d2), scDesign3 = ...)}.
#' @param compute_bivariate Logical. Whether to compute 2D bivariate tests. Default \code{FALSE}
#'   for fast high-throughput execution.
#' @param threads Number of threads / parallel workers. Default is \code{1}.
#' @param auto_sparse Logical. Automatically coerce dense count matrices to \code{dgCMatrix}
#'   to minimize RAM usage. Default \code{TRUE}.
#' @param run_gc Logical. Run garbage collection (\code{gc()}) after each dataset evaluation
#'   to prevent memory accumulation on HPC nodes. Default \code{TRUE}.
#' @param save_path Optional file path (ending in \code{.rds}) to directly persist the
#'   consolidated benchmark bundle for loading into Shiny Studio.
#' @param verbose Logical. Print informative progress logging. Default \code{TRUE}.
#'
#' @return An object of class \code{scSimEval_multidataset_benchmark} containing:
#'   \itemize{
#'     \item \code{granular_benchmark_table}: Tidy data.frame containing all evaluated metrics across every dataset and simulator.
#'     \item \code{benchmark_summary_table}: Tidy data.frame containing cross-dataset averaged fidelity scores for each simulator (ready for \code{plot_benchmark_bubble_matrix()}).
#'     \item \code{dataset_summary_table}: Summary of reference dataset dimensions, sparsity, and library properties.
#'     \item \code{methods}: Vector of evaluated simulator names.
#'     \item \code{datasets}: Vector of evaluated dataset names.
#'     \item \code{cross_dataset_stats}: Summary statistics (mean, SD, SE) for each metric across datasets.
#'   }
#' @export
#' @examples
#' \donttest{
#' data(example_scrna, package = "scSimEval")
#' refs <- list("CohortA" = example_scrna$ref)
#' sims <- list("Splatter" = list("CohortA" = example_scrna$sim))
#' res <- benchmark_simulation_suite(refs, sims, verbose = FALSE)
#' print(res)
#' }
benchmark_simulation_suite <- function(
  datasets,
  simulators,
  compute_bivariate = FALSE,
  threads           = 1,
  auto_sparse       = TRUE,
  run_gc            = TRUE,
  save_path         = NULL,
  verbose           = TRUE
) {
  if (!is.list(datasets) || length(datasets) == 0) {
    stop("`datasets` must be a non-empty named list of empirical reference datasets.")
  }
  if (!is.list(simulators) || length(simulators) == 0) {
    stop("`simulators` must be a non-empty named list of simulator outputs.")
  }
  if (is.null(names(datasets)) || any(names(datasets) == "")) {
    names(datasets) <- paste0("Dataset_", seq_along(datasets))
  }
  if (is.null(names(simulators)) || any(names(simulators) == "")) {
    names(simulators) <- paste0("Simulator_", seq_along(simulators))
  }

  d_names <- names(datasets)
  s_names <- names(simulators)
  n_datasets <- length(d_names)
  n_sims <- length(s_names)
  total_evals <- n_datasets * n_sims

  if (verbose) {
    message("=====================================================================")
    message(sprintf("  scSimEval Multi-Dataset Benchmarking Suite (HPC Optimized)"))
    message(sprintf("  Datasets (%d)   : %s", n_datasets, paste(head(d_names, 4), collapse = ", ")))
    message(sprintf("  Simulators (%d) : %s", n_sims, paste(s_names, collapse = ", ")))
    message(sprintf("  Total Runs     : %d evaluations (%d datasets x %d simulators)", total_evals, n_datasets, n_sims))
    message(sprintf("  Parallel Workers: %d", threads))
    message("=====================================================================")
  }

  # Helper for sparse matrix coercion
  .ensure_mat <- function(m) {
    if (is.character(m) && length(m) == 1 && file.exists(m)) {
      ext <- tolower(tools::file_ext(m))
      if (ext == "rds") m <- readRDS(m)
      else if (ext == "csv") m <- as.matrix(utils::read.csv(m, row.names = 1, check.names = FALSE))
      else m <- as.matrix(utils::read.table(m, header = TRUE, row.names = 1, check.names = FALSE))
    }
    if (inherits(m, "SingleCellExperiment") || inherits(m, "Seurat")) {
      return(m)
    }
    if (is.data.frame(m)) m <- as.matrix(m)
    if (auto_sparse && !inherits(m, "dgCMatrix") && (is.matrix(m) || inherits(m, "Matrix"))) {
      m <- Matrix::Matrix(m, sparse = TRUE)
    }
    SingleCellExperiment::SingleCellExperiment(assays = list(counts = m))
  }

  all_granular_rows <- list()
  dataset_summaries <- list()
  eval_count <- 0

  for (i in seq_along(d_names)) {
    d_name <- d_names[i]
    d_item <- datasets[[d_name]]

    ref_m <- if (is.list(d_item) && !inherits(d_item, "Matrix") && !is.data.frame(d_item) && !is.matrix(d_item)) {
      if (!is.null(d_item$ref)) d_item$ref else d_item[[1]]
    } else {
      d_item
    }
    ref_m <- .ensure_mat(ref_m)

    ct_labels <- if (is.list(d_item) && !is.null(d_item$cell_types)) d_item$cell_types else NULL
    batch_info <- if (is.list(d_item) && !is.null(d_item$batch_info)) d_item$batch_info else NULL

    n_c <- ncol(ref_m)
    n_f <- nrow(ref_m)
    cnts_ref <- if (inherits(ref_m, "SingleCellExperiment")) {
      SingleCellExperiment::counts(ref_m)
    } else if (inherits(ref_m, "Seurat")) {
      tryCatch(ref_m[["RNA"]]$counts, error = function(e) Seurat::GetAssayData(ref_m, slot = "counts"))
    } else {
      ref_m
    }
    sparsity_val <- round(1 - (as.numeric(Matrix::nnzero(cnts_ref)) / (as.numeric(n_c) * as.numeric(n_f))), 4)

    dataset_summaries[[d_name]] <- data.frame(
      "Dataset / Simulator" = d_name,
      "Role"                = "Biological Reference",
      "Modality"            = "Single-Cell (Counts)",
      "Cells (N)"           = n_c,
      "Features (p)"        = n_f,
      "Sparsity (%)"        = sprintf("%.1f%%", sparsity_val * 100),
      "Cell Types"          = if (!is.null(ct_labels)) length(unique(ct_labels)) else NA_integer_,
      check.names           = FALSE,
      stringsAsFactors      = FALSE
    )

    if (verbose) {
      message(sprintf("\n>>> [%d/%d] Ingesting Reference Dataset: '%s' (%d cells, %d features, %.1f%% sparsity)",
                      i, n_datasets, d_name, n_c, n_f, sparsity_val * 100))
    }

    for (j in seq_along(s_names)) {
      s_name <- s_names[j]
      eval_count <- eval_count + 1
      s_container <- simulators[[s_name]]

      sim_m <- if (is.list(s_container) && d_name %in% names(s_container)) {
        s_container[[d_name]]
      } else if (is.list(s_container) && length(s_container) >= i) {
        s_container[[i]]
      } else if (is.character(s_container) && length(s_container) >= i) {
        s_container[i]
      } else {
        s_container
      }

      if (is.null(sim_m)) {
        warning(sprintf("Simulator '%s' has no simulation output for dataset '%s'. Skipping.", s_name, d_name))
        next
      }

      sim_m <- .ensure_mat(sim_m)

      if (verbose) {
        message(sprintf("  --> [%d/%d] Evaluating Simulator '%s' on '%s'...",
                        eval_count, total_evals, s_name, d_name))
      }

      t0 <- Sys.time()
      res_ij <- tryCatch({
        evaluate_simulation_accuracy(
          ref_data          = ref_m,
          sim_data          = sim_m,
          cell_types        = ct_labels,
          compute_bivariate = compute_bivariate,
          threads           = threads,
          verbose           = FALSE
        )
      }, error = function(e) {
        warning(sprintf("Failed evaluating '%s' on '%s': %s", s_name, d_name, e$message))
        NULL
      })
      t1 <- Sys.time()
      elapsed_s <- round(as.numeric(difftime(t1, t0, units = "secs")), 2)

      if (!is.null(res_ij)) {
        tbl_ij <- res_ij$metrics_summary_table
        tbl_ij$Dataset <- d_name
        tbl_ij$Method <- s_name
        tbl_ij$Simulator <- s_name

        if (!"Score_Raw" %in% colnames(tbl_ij) && "Value" %in% colnames(tbl_ij)) {
          tbl_ij$Score_Raw <- tbl_ij$Value
        }
        if (!"Score" %in% colnames(tbl_ij)) {
          tbl_ij$Score <- tbl_ij$Score_Raw
        }

        # Normalize score
        if (!"Score_Norm" %in% colnames(tbl_ij)) {
          tbl_ij <- .normalize_bubble_scores(tbl_ij)
        }

        all_granular_rows[[length(all_granular_rows) + 1]] <- tbl_ij
        if (verbose) {
          message(sprintf("      Completed: %d metrics evaluated in %.1fs.", nrow(tbl_ij), elapsed_s))
        }
      }
    }

    if (run_gc) {
      gc(verbose = FALSE)
    }
  }

  if (length(all_granular_rows) == 0) {
    stop("No evaluations succeeded across the multi-dataset benchmark suite.")
  }

  granular_tbl <- do.call(rbind, all_granular_rows)
  rownames(granular_tbl) <- NULL

  # Standardize categories
  if ("Category" %in% colnames(granular_tbl)) {
    mc_override <- .METRIC_CATEGORY_MAP[granular_tbl$Metric]
    idx <- !is.na(mc_override)
    granular_tbl$Category[idx] <- mc_override[idx]
  }

  # Re-normalize direction-aware scores across all evaluated simulators
  if (nrow(granular_tbl) > 0 && "Score_Raw" %in% colnames(granular_tbl)) {
    granular_tbl <- .normalize_bubble_scores(granular_tbl)
    granular_tbl$Score <- granular_tbl$Score_Norm
  }

  # Build Cross-Dataset Summary Table (mean across datasets for each simulator)
  agg_norm <- stats::aggregate(
    Score_Norm ~ Method + Category + Metric,
    data = granular_tbl,
    FUN  = mean,
    na.rm = TRUE
  )
  agg_raw <- stats::aggregate(
    Score_Raw ~ Method + Category + Metric,
    data = granular_tbl,
    FUN  = mean,
    na.rm = TRUE
  )
  summary_tbl <- merge(agg_norm, agg_raw, by = c("Method", "Category", "Metric"))
  summary_tbl$Score <- round(summary_tbl$Score_Norm, 4)
  summary_tbl$Score_Norm <- round(summary_tbl$Score_Norm, 4)
  summary_tbl$Score_Raw  <- round(summary_tbl$Score_Raw, 4)
  summary_tbl$Simulator  <- summary_tbl$Method

  # Cross-dataset statistics (mean, SD, SE, min, max across datasets)
  stats_agg <- stats::aggregate(
    Score_Norm ~ Method + Category + Metric,
    data = granular_tbl,
    FUN = function(x) c(
      Mean = mean(x, na.rm = TRUE),
      SD = stats::sd(x, na.rm = TRUE),
      SE = stats::sd(x, na.rm = TRUE) / sqrt(length(stats::na.omit(x))),
      Min = min(x, na.rm = TRUE),
      Max = max(x, na.rm = TRUE),
      N = length(stats::na.omit(x))
    )
  )
  stats_df <- data.frame(
    Method     = stats_agg$Method,
    Category   = stats_agg$Category,
    Metric     = stats_agg$Metric,
    Mean_Score = round(stats_agg$Score_Norm[, "Mean"], 4),
    SD_Score   = round(stats_agg$Score_Norm[, "SD"], 4),
    SE_Score   = round(stats_agg$Score_Norm[, "SE"], 4),
    Min_Score  = round(stats_agg$Score_Norm[, "Min"], 4),
    Max_Score  = round(stats_agg$Score_Norm[, "Max"], 4),
    N_Datasets = as.integer(stats_agg$Score_Norm[, "N"]),
    stringsAsFactors = FALSE
  )

  d_summary_df <- do.call(rbind, dataset_summaries)

  out <- list(
    granular_benchmark_table = granular_tbl,
    benchmark_summary_table  = summary_tbl,
    cross_dataset_stats      = stats_df,
    dataset_summary          = d_summary_df,
    dataset_summary_df       = d_summary_df,
    methods                  = unique(summary_tbl$Method),
    datasets                 = d_names,
    is_multidataset          = TRUE,
    n_datasets               = n_datasets,
    n_simulators             = n_sims
  )
  class(out) <- c("scSimEval_multidataset_benchmark", "scSimEval_consolidated", "list")

  if (!is.null(save_path)) {
    saveRDS(out, save_path)
    if (verbose) {
      message(sprintf("\n>>> Benchmark results persisted successfully to: '%s'", save_path))
    }
  }

  if (verbose) {
    message("\n=====================================================================")
    message(sprintf("  Benchmark Suite Completed! %d evaluations across %d datasets & %d simulators.",
                    nrow(granular_tbl), n_datasets, n_sims))
    message("  Cross-dataset averages computed and ready for plot_benchmark_bubble_matrix().")
    message("=====================================================================")
  }

  out
}

#' @export
print.scSimEval_multidataset_benchmark <- function(x, ...) {
  cat("=====================================================================\n")
  cat("  scSimEval Multi-Dataset Benchmark Results\n")
  cat("=====================================================================\n")
  cat(sprintf("Datasets Evaluated     : %d (%s)\n", length(x$datasets), paste(head(x$datasets, 4), collapse = ", ")))
  cat(sprintf("Simulators Evaluated   : %d (%s)\n", length(x$methods), paste(x$methods, collapse = ", ")))
  cat(sprintf("Total Evaluated Scores : %d across all datasets & metrics\n", nrow(x$granular_benchmark_table)))
  cat(sprintf("Cross-Dataset Average  : %d averaged metric entries in summary table\n", nrow(x$benchmark_summary_table)))
  cat("=====================================================================\n")
  invisible(x)
}

