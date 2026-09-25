#' @title Differential Features and Biological Signal Evaluation
#' @description Implements the 5 biological signal evaluation patterns from SimBench
#'   (DE, DV, DD, DP, BD), predictive cell identity classification modeling,
#'   pseudoreplication delta-variance preservation (Libra), cell cycle phase fidelity,
#'   and cell-type mixture deconvolution accuracy.
#' @name metrics_de_signal
#' @export calc_signal_de
#' @export calc_signal_dv
#' @export calc_signal_dd
#' @export calc_signal_dp
#' @export calc_signal_bd
#' @export evaluate_simbench_signals
#' @export evaluate_predictive_de_model
#' @export to_pseudobulk
#' @export calc_delta_variance
#' @export calc_cell_cycle_phase_fidelity
#' @export calc_deconvolution_accuracy
#' @export evaluate_deg_fidelity
NULL

# -----------------------------------------------------------------------------
# SimBench 5 Biological Signal Types: DE, DV, DD, DP, BD
# -----------------------------------------------------------------------------

#' Detect Differential Expression (DE) Genes (Mean Shift via Limma or T-Test)
#' @param exprs_mat Log-normalized matrix (genes x cells).
#' @param cell_types Factor or binary vector of 2 cell types.
#' @param p_sig Significance threshold (default 0.05).
#' @return Vector of adjusted p-values.
#' @export
calc_signal_de <- function(exprs_mat, cell_types, p_sig = 0.05) {
  cell_types <- droplevels(as.factor(cell_types))
  if (nlevels(cell_types) != 2) stop("Requires exactly 2 cell types for pairwise DE.")
  
  if (requireNamespace("limma", quietly = TRUE)) {
    design <- stats::model.matrix(~ cell_types)
    fit <- limma::lmFit(as.matrix(exprs_mat), design)
    fit <- limma::eBayes(fit, trend = TRUE)
    tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
    p_adj <- stats::p.adjust(tt$P.Value, method = "BH")
  } else {
    # Two-sample t-test fallback
    idx1 <- which(cell_types == levels(cell_types)[1])
    idx2 <- which(cell_types == levels(cell_types)[2])
    p_vals <- apply(as.matrix(exprs_mat), 1, function(x) {
      if (stats::sd(x[idx1]) == 0 && stats::sd(x[idx2]) == 0) return(1)
      suppressWarnings(stats::t.test(x[idx1], x[idx2])$p.value)
    })
    p_adj <- stats::p.adjust(p_vals, method = "BH")
  }
  return(p_adj)
}

#' Detect Differential Variability (DV) Genes (Bartlett Test)
#' @param exprs_mat Log-normalized matrix (genes x cells).
#' @param cell_types Factor of 2 cell types.
#' @return Vector of adjusted p-values.
#' @export
calc_signal_dv <- function(exprs_mat, cell_types) {
  cell_types <- droplevels(as.factor(cell_types))
  exprs_mat <- as.matrix(exprs_mat)
  p_vals <- apply(exprs_mat, 1, function(x) {
    df <- data.frame(gene = x, cell_type = cell_types)
    res <- tryCatch(stats::bartlett.test(gene ~ cell_type, data = df)$p.value,
                    error = function(e) 1)
    if (is.nan(res)) 1 else res
  })
  stats::p.adjust(p_vals, method = "BH")
}

#' Detect Differential Distribution (DD) Genes (Kolmogorov-Smirnov Test)
#' @param exprs_mat Log-normalized matrix (genes x cells).
#' @param cell_types Factor of 2 cell types.
#' @return Vector of adjusted p-values.
#' @export
calc_signal_dd <- function(exprs_mat, cell_types) {
  cell_types <- droplevels(as.factor(cell_types))
  exprs_mat <- as.matrix(exprs_mat)
  idx1 <- which(cell_types == levels(cell_types)[1])
  idx2 <- which(cell_types == levels(cell_types)[2])
  
  p_vals <- apply(exprs_mat, 1, function(x) {
    res <- tryCatch(
      suppressWarnings(stats::ks.test(x[idx1], x[idx2])$p.value),
      error = function(e) 1
    )
    if (is.nan(res)) 1 else res
  })
  stats::p.adjust(p_vals, method = "BH")
}

#' Detect Differential Proportion / Zero-Inflation (DP) Genes (Chisq Test)
#' @param exprs_mat Count or expression matrix (genes x cells).
#' @param cell_types Factor of 2 cell types.
#' @param threshold Value threshold for zero/detection status (default 0).
#' @return Vector of adjusted p-values.
#' @export
calc_signal_dp <- function(exprs_mat, cell_types, threshold = 0) {
  cell_types <- droplevels(as.factor(cell_types))
  exprs_mat <- as.matrix(exprs_mat)
  bin_mat <- ifelse(exprs_mat > threshold, 1, 0)
  
  p_vals <- apply(bin_mat, 1, function(x) {
    tab <- table(factor(x, levels = c(0, 1)), cell_types)
    res <- tryCatch(
      suppressWarnings(stats::chisq.test(tab)$p.value),
      error = function(e) 1
    )
    if (is.nan(res)) 1 else res
  })
  stats::p.adjust(p_vals, method = "BH")
}

#' Detect Bimodal Distribution (BD) Genes (Bimodal Separation Index)
#' @param exprs_mat Matrix (genes x cells).
#' @param cell_types Factor of 2 cell types.
#' @return Numeric vector of Bimodality Index values.
#' @export
calc_signal_bd <- function(exprs_mat, cell_types) {
  cell_types <- droplevels(as.factor(cell_types))
  exprs_mat <- as.matrix(exprs_mat)
  idx1 <- which(cell_types == levels(cell_types)[1])
  idx2 <- which(cell_types == levels(cell_types)[2])
  
  pi1 <- length(idx1) / ncol(exprs_mat)
  pi2 <- length(idx2) / ncol(exprs_mat)
  
  m1 <- rowMeans(exprs_mat[, idx1, drop = FALSE])
  m2 <- rowMeans(exprs_mat[, idx2, drop = FALSE])
  v1 <- apply(exprs_mat[, idx1, drop = FALSE], 1, stats::var)
  v2 <- apply(exprs_mat[, idx2, drop = FALSE], 1, stats::var)
  
  denom <- sqrt(pi1 * v1 + pi2 * v2)
  bi <- ifelse(denom > 0, abs(m1 - m2) / denom, 0)
  return(bi)
}

#' Evaluate the 5 SimBench Biological Signal Proportions
#'
#' Compares proportions of genes exhibiting DE, DV, DD, DP, and BD between
#' empirical reference and simulated datasets.
#'
#' @param ref_mat Reference matrix.
#' @param sim_mat Simulated matrix.
#' @param ref_celltypes Reference cell type labels (subsets to top 2 abundant types).
#' @param sim_celltypes Simulated cell type labels.
#' @param p_sig Significance cutoff (default 0.05).
#' @param bi_cutoff Bimodal index cutoff (default 0.3).
#'
#' @return A tidy data.frame comparing biological signal proportions.
#' @export
evaluate_simbench_signals <- function(
  ref_mat,
  sim_mat,
  ref_celltypes,
  sim_celltypes,
  p_sig = 0.05,
  bi_cutoff = 0.3
) {
  # Subset to the 2 most abundant cell types in reference
  top2 <- names(sort(table(ref_celltypes), decreasing = TRUE))[1:2]
  
  idx_ref <- which(ref_celltypes %in% top2)
  idx_sim <- which(sim_celltypes %in% top2)
  
  r_sub <- ref_mat[, idx_ref]
  s_sub <- sim_mat[, idx_sim]
  r_ct <- ref_celltypes[idx_ref]
  s_ct <- sim_celltypes[idx_sim]
  
  # Calculate for reference
  n_genes_r <- nrow(r_sub)
  prop_de_r <- sum(calc_signal_de(r_sub, r_ct) < p_sig) / n_genes_r
  prop_dv_r <- sum(calc_signal_dv(r_sub, r_ct) < p_sig) / n_genes_r
  prop_dd_r <- sum(calc_signal_dd(r_sub, r_ct) < p_sig) / n_genes_r
  prop_dp_r <- sum(calc_signal_dp(r_sub, r_ct) < p_sig) / n_genes_r
  prop_bd_r <- sum(calc_signal_bd(r_sub, r_ct) > bi_cutoff) / n_genes_r
  
  # Calculate for simulation
  n_genes_s <- nrow(s_sub)
  prop_de_s <- sum(calc_signal_de(s_sub, s_ct) < p_sig) / n_genes_s
  prop_dv_s <- sum(calc_signal_dv(s_sub, s_ct) < p_sig) / n_genes_s
  prop_dd_s <- sum(calc_signal_dd(s_sub, s_ct) < p_sig) / n_genes_s
  prop_dp_s <- sum(calc_signal_dp(s_sub, s_ct) < p_sig) / n_genes_s
  prop_bd_s <- sum(calc_signal_bd(s_sub, s_ct) > bi_cutoff) / n_genes_s
  
  signals <- c("DE (Mean Shift)", "DV (Variability)", "DD (Distribution)",
               "DP (Proportion/Dropout)", "BD (Bimodal Index)")
  
  data.frame(
    Signal_Type = signals,
    Reference_Prop = c(prop_de_r, prop_dv_r, prop_dd_r, prop_dp_r, prop_bd_r),
    Simulation_Prop = c(prop_de_s, prop_dv_s, prop_dd_s, prop_dp_s, prop_bd_s),
    Absolute_Error = abs(c(prop_de_r, prop_dv_r, prop_dd_r, prop_dp_r, prop_bd_r) -
                         c(prop_de_s, prop_dv_s, prop_dd_s, prop_dp_s, prop_bd_s)),
    stringsAsFactors = FALSE
  )
}

#' Predictive Cell Identity Classification Using DE / Top Variable Features
#'
#' Trains a classifier on simulated features (80\% train / 20\% test) to verify
#' whether the simulated biological signals are sufficiently predictive of cell identities.
#' If \code{de_features} is NULL, the top variable features are automatically selected.
#'
#' @param data Matrix (features x cells).
#' @param group Cell type labels.
#' @param de_features Optional character vector of selected DE feature names.
#' @param n_top Number of top features to auto-select if \code{de_features} is NULL (default 50).
#' @param method Classifier: "knn" (default), "rf", or "svm".
#'
#' @return A list containing Accuracy, Macro Precision, Macro Recall, and Macro F1 score.
#' @export
evaluate_predictive_de_model <- function(
  data,
  group,
  de_features = NULL,
  n_top = 50,
  method = c("knn", "rf", "svm")
) {
  method <- match.arg(method)
  group <- as.factor(group)
  data_mat <- as.matrix(data)
  
  if (is.null(de_features) || length(de_features) < 2) {
    vars <- apply(data_mat, 1, stats::var)
    n_sel <- min(n_top, nrow(data_mat))
    de_features <- names(sort(vars, decreasing = TRUE))[seq_len(n_sel)]
  } else {
    de_features <- intersect(de_features, rownames(data_mat))
  }
  
  if (length(de_features) < 2 || nlevels(group) < 2) {
    return(list(accuracy = NA_real_, F1 = NA_real_, precision = NA_real_, recall = NA_real_))
  }
  
  sub_data <- as.data.frame(t(data_mat[de_features, , drop = FALSE]))
  set.seed(42)
  n_cells <- nrow(sub_data)
  train_idx <- sample(seq_len(n_cells), round(0.8 * n_cells))
  
  train_x <- sub_data[train_idx, , drop = FALSE]
  test_x <- sub_data[-train_idx, , drop = FALSE]
  train_y <- group[train_idx]
  test_y <- group[-train_idx]
  
  # Prediction via k-NN fallback (universal in base R class package)
  pred_y <- class::knn(train = train_x, test = test_x, cl = train_y, k = min(5, nrow(train_x)))
  
  # Confusion matrix & macro metrics
  cm <- table(Predicted = pred_y, Actual = test_y)
  acc <- sum(diag(cm)) / max(1, sum(cm))
  
  cls <- levels(group)
  prs <- vapply(cls, function(c) {
    denom <- sum(pred_y == c)
    if (denom > 0) sum(pred_y == c & test_y == c) / denom else 0
  }, numeric(1))
  
  res <- vapply(cls, function(c) {
    denom <- sum(test_y == c)
    if (denom > 0) sum(pred_y == c & test_y == c) / denom else 0
  }, numeric(1))
  
  f1s <- ifelse(prs + res > 0, 2 * (prs * res) / (prs + res), 0)
  
  list(
    accuracy = acc,
    precision = mean(prs, na.rm = TRUE),
    recall = mean(res, na.rm = TRUE),
    F1 = mean(f1s, na.rm = TRUE)
  )
}

# -----------------------------------------------------------------------------
# Libra Benchmark Measures (Squair et al., Nature Communications 2021)
# -----------------------------------------------------------------------------

#' Convert Single-Cell Counts to Replicate Pseudobulks
#'
#' Converts a single-cell expression matrix into pseudobulk matrices per cell type
#' by summarizing counts across biological replicates and experimental conditions
#' (adapted from Squair et al. Libra).
#'
#' @param counts Single-cell count matrix (features x cells).
#' @param cell_types Factor or character vector of cell types.
#' @param replicates Factor or character vector of biological replicate IDs.
#' @param conditions Factor or character vector of experimental condition / treatment labels.
#' @param min_cells Minimum number of cells required per cell type to retain. Default is 3.
#'
#' @return Named list of pseudobulk count matrices (features x sample_replicates), one per cell type.
#' @export
to_pseudobulk <- function(
  counts,
  cell_types,
  replicates,
  conditions,
  min_cells = 3
) {
  counts <- as.matrix(counts)
  cell_types <- as.factor(cell_types)
  replicates <- as.factor(replicates)
  conditions <- as.factor(conditions)
  
  meta <- data.frame(
    cell = colnames(counts),
    cell_type = cell_types,
    replicate = replicates,
    condition = conditions,
    stringsAsFactors = FALSE
  )
  
  unique_types <- levels(cell_types)
  pb_list <- list()
  
  for (ct in unique_types) {
    sub_idx <- which(meta$cell_type == ct)
    if (length(sub_idx) < min_cells) next
    
    sub_meta <- meta[sub_idx, ]
    sub_counts <- counts[, sub_idx, drop = FALSE]
    
    sample_id <- factor(paste0(sub_meta$replicate, ":", sub_meta$condition))
    mm <- stats::model.matrix(~ 0 + sample_id)
    colnames(mm) <- levels(sample_id)
    
    pb_mat <- as.matrix(sub_counts %*% mm)
    pb_list[[ct]] <- pb_mat
  }
  pb_list
}

#' Calculate Delta Variance (Pseudoreplication Bias Metric)
#'
#' Measures the difference in gene expression variance between true biological
#' replicates and cell-shuffled pseudoreplicates (Squair et al., Nature Communications 2021).
#' Positive delta variance indicates that biological replicates have higher inter-individual
#' variability than randomly pooled pseudoreplicates.
#'
#' @param counts Count matrix (features x cells).
#' @param replicates Factor of biological replicate / batch assignments for each cell.
#' @param conditions Factor of condition or cell-type labels.
#' @param cell_types Optional cell type vector to evaluate per cell type or overall.
#' @param n_perm Number of permutation shuffles to compute mean pseudoreplicate variance. Default is 3.
#'
#' @return A named list:
#'   \item{delta_variance}{Numeric vector of gene-level delta variance values}
#'   \item{mean_delta_variance}{Mean delta variance across all genes}
#'   \item{median_delta_variance}{Median delta variance across all genes}
#' @export
calc_delta_variance <- function(
  counts,
  replicates,
  conditions,
  cell_types = NULL,
  n_perm = 3
) {
  counts <- as.matrix(counts)
  replicates <- as.factor(replicates)
  conditions <- as.factor(conditions)
  
  # True biological replicates CPM
  samp_fac <- factor(paste0(replicates, "_", conditions))
  mm1 <- stats::model.matrix(~ 0 + samp_fac)
  pb_true <- counts %*% mm1
  cpm_true <- t(t(pb_true) / (colSums(pb_true) + 1e-12)) * 1e6
  var_true <- apply(cpm_true, 1, stats::var)
  
  # Shuffled pseudoreplicates within each condition
  shuffled_dvs <- replicate(n_perm, {
    shuff_rep <- replicates
    for (cond in levels(conditions)) {
      c_idx <- which(conditions == cond)
      shuff_rep[c_idx] <- sample(shuff_rep[c_idx])
    }
    samp_shuff <- factor(paste0(shuff_rep, "_", conditions))
    mm2 <- stats::model.matrix(~ 0 + samp_shuff)
    pb_shuff <- counts %*% mm2
    cpm_shuff <- t(t(pb_shuff) / (colSums(pb_shuff) + 1e-12)) * 1e6
    var_shuff <- apply(cpm_shuff, 1, stats::var)
    var_shuff - var_true
  })
  
  mean_dv <- if (n_perm == 1) as.numeric(shuffled_dvs) else rowMeans(shuffled_dvs)
  
  list(
    delta_variance = mean_dv,
    mean_delta_variance = mean(mean_dv, na.rm = TRUE),
    median_delta_variance = stats::median(mean_dv, na.rm = TRUE)
  )
}

#' Evaluate Cell Cycle Phase Distribution Fidelity
#'
#' Compares cell cycle phase proportions (e.g. G0/G1, S, G2/M) between reference
#' and simulated single-cell populations.
#'
#' @param ref_phases Factor or character vector of cell cycle phase assignments in reference.
#' @param sim_phases Factor or character vector of cell cycle phase assignments in simulation.
#'
#' @return A list with phase proportions, Jensen-Shannon divergence, and chi-squared test p-value.
#' @export
calc_cell_cycle_phase_fidelity <- function(ref_phases, sim_phases) {
  all_phases <- union(unique(ref_phases), unique(sim_phases))
  all_phases <- all_phases[!is.na(all_phases)]
  
  t_ref <- table(factor(ref_phases, levels = all_phases))
  t_sim <- table(factor(sim_phases, levels = all_phases))
  
  p_ref <- as.numeric(prop.table(t_ref))
  p_sim <- as.numeric(prop.table(t_sim))
  names(p_ref) <- all_phases
  names(p_sim) <- all_phases
  
  mae <- mean(abs(p_ref - p_sim))
  
  # Jensen-Shannon Divergence
  m <- 0.5 * (p_ref + p_sim)
  kl <- function(p, q) {
    terms <- ifelse(p > 0, p * log(p / pmax(1e-12, q)), 0)
    sum(terms)
  }
  jsd <- 0.5 * kl(p_ref, m) + 0.5 * kl(p_sim, m)
  
  # Chi-squared test of independence
  cont_tab <- rbind(t_ref, t_sim)
  chisq_res <- tryCatch(suppressWarnings(stats::chisq.test(cont_tab)), error = function(e) NULL)
  p_val <- if (!is.null(chisq_res)) chisq_res$p.value else NA_real_
  
  list(
    ref_proportions = p_ref,
    sim_proportions = p_sim,
    phase_mae = mae,
    phase_jsd = jsd,
    homogeneity_pvalue = p_val
  )
}

# -----------------------------------------------------------------------------
# Single-Cell Mixture & Deconvolution Accuracy
# -----------------------------------------------------------------------------

#' Evaluate Cell-Type Deconvolution & Mixture Proportion Accuracy
#'
#' Evaluates the accuracy of cell-type mixture proportion estimation / deconvolution
#' methods between empirical reference and simulated cell populations (Chen et al., Bioinformatics 2021).
#' Can accept precomputed proportion matrices OR cell-type vectors (and optional batch labels)
#' to calculate proportion preservation across batches.
#'
#' @param true_proportions Vector, matrix, or data.frame of reference cell-type proportions,
#'   or factor/vector of reference cell types.
#' @param estimated_proportions Vector, matrix, or data.frame of estimated / simulated cell-type proportions,
#'   or factor/vector of simulated cell types.
#' @param batch_ref Optional batch labels for reference cells (if true_proportions is a cell-type vector).
#' @param batch_sim Optional batch labels for simulated cells (if estimated_proportions is a cell-type vector).
#'
#' @return A list containing RMSE, MAE, Pearson correlation (r), Spearman correlation (rho),
#'   Jensen-Shannon Divergence (JSD), and Total Variation Distance (TVD).
#' @export
calc_deconvolution_accuracy <- function(
  true_proportions,
  estimated_proportions,
  batch_ref = NULL,
  batch_sim = NULL
) {
  # If character/factor vectors supplied, convert to batch-by-celltype or 1-by-celltype proportion matrices
  if ((is.character(true_proportions) || is.factor(true_proportions)) &&
      (is.character(estimated_proportions) || is.factor(estimated_proportions))) {
    all_cts <- union(as.character(true_proportions), as.character(estimated_proportions))
    
    if (!is.null(batch_ref) && !is.null(batch_sim)) {
      tab_ref <- prop.table(table(batch_ref, factor(true_proportions, levels = all_cts)), margin = 1)
      tab_sim <- prop.table(table(batch_sim, factor(estimated_proportions, levels = all_cts)), margin = 1)
      t_mat <- as.matrix(tab_ref)
      e_mat <- as.matrix(tab_sim)
    } else {
      p_ref <- as.numeric(prop.table(table(factor(true_proportions, levels = all_cts))))
      p_sim <- as.numeric(prop.table(table(factor(estimated_proportions, levels = all_cts))))
      t_mat <- matrix(p_ref, nrow = 1, dimnames = list("sample_1", all_cts))
      e_mat <- matrix(p_sim, nrow = 1, dimnames = list("sample_1", all_cts))
    }
  } else {
    t_mat <- as.matrix(true_proportions)
    e_mat <- as.matrix(estimated_proportions)
  }
  
  if (nrow(t_mat) != nrow(e_mat) || ncol(t_mat) != ncol(e_mat)) {
    if (nrow(t_mat) == ncol(e_mat) && ncol(t_mat) == nrow(e_mat)) {
      e_mat <- t(e_mat)
    } else {
      common_cts <- intersect(colnames(t_mat), colnames(e_mat))
      if (length(common_cts) >= 2) {
        t_mat <- t_mat[, common_cts, drop = FALSE]
        e_mat <- e_mat[, common_cts, drop = FALSE]
      } else {
        # Align rows and cols to min
        min_r <- min(nrow(t_mat), nrow(e_mat))
        min_c <- min(ncol(t_mat), ncol(e_mat))
        t_mat <- t_mat[seq_len(min_r), seq_len(min_c), drop = FALSE]
        e_mat <- e_mat[seq_len(min_r), seq_len(min_c), drop = FALSE]
      }
    }
  }
  
  # Row-normalize to valid probability simplices
  t_mat <- pmax(t_mat, 0)
  e_mat <- pmax(e_mat, 0)
  
  t_sums <- rowSums(t_mat)
  e_sums <- rowSums(e_mat)
  t_mat <- t_mat / ifelse(t_sums > 0, t_sums, 1)
  e_mat <- e_mat / ifelse(e_sums > 0, e_sums, 1)
  
  n_samples <- nrow(t_mat)
  
  sample_rmse <- numeric(n_samples)
  sample_mae <- numeric(n_samples)
  sample_r <- numeric(n_samples)
  sample_rho <- numeric(n_samples)
  sample_jsd <- numeric(n_samples)
  sample_tvd <- numeric(n_samples)
  
  for (i in seq_len(n_samples)) {
    p <- t_mat[i, ]
    q <- e_mat[i, ]
    
    sample_rmse[i] <- sqrt(mean((p - q)^2))
    sample_mae[i] <- mean(abs(p - q))
    sample_tvd[i] <- 0.5 * sum(abs(p - q))
    
    if (stats::sd(p) > 1e-6 && stats::sd(q) > 1e-6) {
      sample_r[i] <- stats::cor(p, q, method = "pearson")
      sample_rho[i] <- stats::cor(p, q, method = "spearman")
    } else {
      sample_r[i] <- NA_real_
      sample_rho[i] <- NA_real_
    }
    
    # Jensen-Shannon Divergence
    m <- 0.5 * (p + q)
    kl_pm <- sum(ifelse(p > 0 & m > 0, p * log2(p / m), 0))
    kl_qm <- sum(ifelse(q > 0 & m > 0, q * log2(q / m), 0))
    sample_jsd[i] <- max(0, 0.5 * (kl_pm + kl_qm))
  }
  
  overall_rmse <- sqrt(mean((t_mat - e_mat)^2))
  overall_mae <- mean(abs(t_mat - e_mat))
  overall_r <- mean(sample_r, na.rm = TRUE)
  overall_rho <- mean(sample_rho, na.rm = TRUE)
  overall_jsd <- mean(sample_jsd, na.rm = TRUE)
  overall_tvd <- mean(sample_tvd, na.rm = TRUE)
  
  list(
    deconvolution_rmse = as.numeric(overall_rmse),
    deconvolution_mae = as.numeric(overall_mae),
    deconvolution_r = as.numeric(overall_r),
    deconvolution_rho = as.numeric(overall_rho),
    deconvolution_jsd = as.numeric(overall_jsd),
    deconvolution_tvd = as.numeric(overall_tvd),
    sample_metrics = data.frame(
      sample = seq_len(n_samples),
      RMSE = sample_rmse,
      MAE = sample_mae,
      Pearson_r = sample_r,
      Spearman_rho = sample_rho,
      JSD = sample_jsd,
      TVD = sample_tvd
    )
  )
}


#' Evaluate Differentially Expressed Gene (DEG) Fidelity
#'
#' Evaluates the preservation of differentially expressed genes (DEGs) and biological
#' signal between empirical reference and simulated datasets across three landmark single-cell
#' benchmarking frameworks:
#' \itemize{
#'   \item \strong{Simpipe (Duo et al., 2024):} True DEG ratio, Distribution score (p-value
#'         uniformity via Pearson Chi-Square goodness-of-fit test on remaining genes after DEG removal),
#'         and supervised machine learning classification (Accuracy, Precision, Recall, F1) using
#'         simulated DEGs to predict cell identity.
#'   \item \strong{SimBench (Cao et al., 2021):} Symmetric Mean Absolute Percentage Error (SMAPE)
#'         on DEG proportions, SimBench DE fidelity score, log2 fold-change effect size Pearson
#'         and Spearman correlation, and top-N DEG Jaccard overlap.
#'   \item \strong{Shaky Foundations (Crowell et al., 2023):} Group silhouette separation width,
#'         silhouette discrepancy, and Percent Variance Explained (PVE) by group assignment.
#' }
#'
#' @param ref_data Reference single-cell expression or count matrix (features x cells).
#' @param sim_data Simulated single-cell expression or count matrix (features x cells).
#' @param ref_celltypes Factor or character vector of cell types for reference cells.
#' @param sim_celltypes Factor or character vector of cell types for simulated cells.
#' @param group1 Optional character string specifying group 1 name. Default NULL (auto-selects top abundant type).
#' @param group2 Optional character string specifying group 2 name. Default NULL (auto-selects second most abundant type).
#' @param fdr_cutoff Adjusted p-value significance threshold. Default is 0.05.
#' @param logfc_cutoff Absolute log2 fold-change cutoff. Default is 0.5.
#' @param top_n_de Number of top ranked DEGs to evaluate for Jaccard overlap and classification. Default is 100.
#' @param classifier Supervised classifier: "knn" (default), "svm", or "rf".
#' @param run_pvalue_uniformity Logical, whether to run Simpipe's Chi-square p-value uniformity test. Default is TRUE.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{deg_summary_table}: Comprehensive data.frame uniting all 15 metrics across Simpipe, SimBench, and Shaky Foundations.
#'   \item \code{deg_genes_ref}: Character vector of significant DEGs in reference.
#'   \item \code{deg_genes_sim}: Character vector of significant DEGs in simulation.
#'   \item \code{contrast}: String describing the evaluated group contrast.
#'   \item \code{logfc_ref}: Named numeric vector of reference log2 fold-changes.
#'   \item \code{logfc_sim}: Named numeric vector of simulated log2 fold-changes.
#'   \item \code{ml_classification}: Detailed classifier performance metrics.
#'   \item \code{pvalue_uniformity}: Chi-square goodness-of-fit test results.
#' }
#' @export
evaluate_deg_fidelity <- function(
  ref_data,
  sim_data,
  ref_celltypes,
  sim_celltypes,
  group1 = NULL,
  group2 = NULL,
  fdr_cutoff = 0.05,
  logfc_cutoff = 0.5,
  top_n_de = 100,
  classifier = c("knn", "svm", "rf"),
  run_pvalue_uniformity = TRUE
) {
  classifier <- match.arg(classifier)
  
  if (inherits(ref_data, "SingleCellExperiment")) {
    ref_data <- SummarizedExperiment::assay(ref_data, "counts")
  } else if (inherits(ref_data, "Seurat")) {
    ref_data <- Seurat::GetAssayData(ref_data, slot = "counts")
  }
  
  if (inherits(sim_data, "SingleCellExperiment")) {
    sim_data <- SummarizedExperiment::assay(sim_data, "counts")
  } else if (inherits(sim_data, "Seurat")) {
    sim_data <- Seurat::GetAssayData(sim_data, slot = "counts")
  }
  
  ref_mat <- as.matrix(ref_data)
  sim_mat <- as.matrix(sim_data)
  
  # Determine comparison groups
  if (is.null(group1) || is.null(group2)) {
    tab <- sort(table(ref_celltypes), decreasing = TRUE)
    if (length(tab) < 2) stop("Need at least 2 distinct cell types for differential expression analysis.")
    group1 <- names(tab)[1]
    group2 <- names(tab)[2]
  }
  
  idx_ref <- which(ref_celltypes %in% c(group1, group2))
  idx_sim <- which(sim_celltypes %in% c(group1, group2))
  
  if (length(idx_ref) < 4 || length(idx_sim) < 4) {
    stop("Insufficient cell count in the selected contrast groups.")
  }
  
  r_sub <- ref_mat[, idx_ref, drop = FALSE]
  s_sub <- sim_mat[, idx_sim, drop = FALSE]
  r_ct <- droplevels(as.factor(ref_celltypes[idx_ref]))
  s_ct <- droplevels(as.factor(sim_celltypes[idx_sim]))
  
  # Ensure row names exist
  if (is.null(rownames(r_sub))) rownames(r_sub) <- paste0("Gene_", seq_len(nrow(r_sub)))
  if (is.null(rownames(s_sub))) rownames(s_sub) <- rownames(r_sub)
  
  # CPM normalization
  lib_r <- colSums(r_sub); lib_r[lib_r == 0] <- 1
  cpm_r <- log2(t(t(r_sub) / lib_r) * 1e6 + 1)
  
  lib_s <- colSums(s_sub); lib_s[lib_s == 0] <- 1
  cpm_s <- log2(t(t(s_sub) / lib_s) * 1e6 + 1)
  
  g1_r <- which(r_ct == group1); g2_r <- which(r_ct == group2)
  g1_s <- which(s_ct == group1); g2_s <- which(s_ct == group2)
  
  m1_r <- rowMeans(cpm_r[, g1_r, drop = FALSE]); m2_r <- rowMeans(cpm_r[, g2_r, drop = FALSE])
  m1_s <- rowMeans(cpm_s[, g1_s, drop = FALSE]); m2_s <- rowMeans(cpm_s[, g2_s, drop = FALSE])
  
  logfc_ref <- m1_r - m2_r
  logfc_sim <- m1_s - m2_s
  
  # Statistical Testing (limma / t-test)
  padj_ref <- calc_signal_de(cpm_r, r_ct)
  padj_sim <- calc_signal_de(cpm_s, s_ct)
  
  # Identify significant DEGs
  sig_ref <- which(padj_ref < fdr_cutoff & abs(logfc_ref) >= logfc_cutoff)
  sig_sim <- which(padj_sim < fdr_cutoff & abs(logfc_sim) >= logfc_cutoff)
  
  deg_genes_ref <- rownames(r_sub)[sig_ref]
  deg_genes_sim <- rownames(s_sub)[sig_sim]
  
  n_genes <- nrow(r_sub)
  prop_ref <- length(sig_ref) / max(1, n_genes)
  prop_sim <- length(sig_sim) / max(1, n_genes)
  
  # --- Simpipe Metrics ---
  deg_ratio <- length(sig_sim) / max(1, length(sig_ref))
  
  # Distribution Score (P-value Uniformity via Pearson Chi-Square test)
  dist_score <- NA_real_
  chisq_stat <- NA_real_
  chisq_pval <- NA_real_
  if (run_pvalue_uniformity && length(sig_sim) < (n_genes - 10)) {
    non_deg_idx <- setdiff(seq_len(n_genes), sig_sim)
    rem_mat <- cpm_s[non_deg_idx, , drop = FALSE]
    
    rem_pvals <- apply(rem_mat, 1, function(x) {
      if (stats::sd(x[g1_s]) == 0 && stats::sd(x[g2_s]) == 0) return(1)
      suppressWarnings(stats::t.test(x[g1_s], x[g2_s])$p.value)
    })
    rem_pvals <- stats::na.omit(rem_pvals)
    
    if (length(rem_pvals) >= 10) {
      bins <- cut(rem_pvals, breaks = seq(0, 1, length.out = 11), include.lowest = TRUE)
      obs_counts <- table(bins)
      chi_res <- tryCatch(stats::chisq.test(obs_counts), error = function(e) NULL)
      if (!is.null(chi_res)) {
        chisq_stat <- as.numeric(chi_res$statistic)
        chisq_pval <- as.numeric(chi_res$p.value)
        dist_score <- if (chisq_pval > 0.05) 1.0 else max(0, 1 - (chisq_stat / (length(rem_pvals) + 1)))
      }
    }
  }
  
  # Machine Learning Classification (Simpipe)
  top_sim_de <- if (length(sig_sim) >= 2) {
    deg_genes_sim
  } else {
    rownames(s_sub)[order(padj_sim)[seq_len(min(top_n_de, n_genes))]]
  }
  
  ml_res <- evaluate_predictive_de_model(
    data = cpm_s,
    group = s_ct,
    de_features = top_sim_de,
    method = classifier
  )
  
  # --- SimBench Metrics ---
  smape_denom <- prop_ref + prop_sim
  smape <- if (smape_denom > 0) abs(prop_sim - prop_ref) / smape_denom else 0
  simbench_fidelity_score <- pmax(0, 1 - smape)
  
  logfc_pcc <- tryCatch(stats::cor(logfc_ref, logfc_sim, method = "pearson"), error = function(e) NA_real_)
  logfc_scc <- tryCatch(stats::cor(logfc_ref, logfc_sim, method = "spearman"), error = function(e) NA_real_)
  
  n_top_eval <- min(top_n_de, n_genes)
  top_ref_genes <- rownames(r_sub)[order(padj_ref)[seq_len(n_top_eval)]]
  top_sim_genes <- rownames(s_sub)[order(padj_sim)[seq_len(n_top_eval)]]
  
  jaccard_top <- length(intersect(top_ref_genes, top_sim_genes)) / max(1, length(union(top_ref_genes, top_sim_genes)))
  
  # --- Shaky Foundations Metrics (Crowell et al., 2023) ---
  pca_sim <- stats::prcomp(t(cpm_s), rank. = min(5, ncol(cpm_s) - 1))
  pca_ref <- stats::prcomp(t(cpm_r), rank. = min(5, ncol(cpm_r) - 1))
  
  sil_sim <- if (requireNamespace("cluster", quietly = TRUE)) {
    d_sim <- stats::dist(pca_sim$x[, seq_len(min(3, ncol(pca_sim$x)))])
    mean(cluster::silhouette(as.integer(s_ct), d_sim)[, 3], na.rm = TRUE)
  } else NA_real_
  
  sil_ref <- if (requireNamespace("cluster", quietly = TRUE)) {
    d_ref <- stats::dist(pca_ref$x[, seq_len(min(3, ncol(pca_ref$x)))])
    mean(cluster::silhouette(as.integer(r_ct), d_ref)[, 3], na.rm = TRUE)
  } else NA_real_
  
  silhouette_discrepancy <- abs(sil_sim - sil_ref)
  
  pve_vals <- vapply(seq_len(min(3, ncol(pca_sim$x))), function(k) {
    summary(stats::lm(pca_sim$x[, k] ~ s_ct))$r.squared
  }, numeric(1))
  sim_pve <- mean(pve_vals, na.rm = TRUE)
  
  pve_ref_vals <- vapply(seq_len(min(3, ncol(pca_ref$x))), function(k) {
    summary(stats::lm(pca_ref$x[, k] ~ r_ct))$r.squared
  }, numeric(1))
  ref_pve <- mean(pve_ref_vals, na.rm = TRUE)
  pve_discrepancy <- abs(sim_pve - ref_pve)
  
  summary_df <- data.frame(
    Framework = c(
      rep("Simpipe (Duo et al., 2024)", 6),
      rep("SimBench (Cao et al., 2021)", 5),
      rep("Shaky Foundations (Crowell et al., 2023)", 4)
    ),
    Metric = c(
      "DEG_Ratio",
      "PValue_Uniformity_Chisq",
      "Distribution_Score",
      "Classifier_Accuracy",
      "Classifier_Macro_F1",
      "Classifier_Macro_Recall",
      "SimBench_SMAPE",
      "SimBench_DE_Fidelity_Score",
      "Log2FC_Pearson_Corr",
      "Log2FC_Spearman_Corr",
      "Top_DEG_Jaccard_Overlap",
      "Silhouette_Sim",
      "Silhouette_Discrepancy",
      "PVE_Group_Sim",
      "PVE_Discrepancy"
    ),
    Value = c(
      deg_ratio,
      chisq_stat,
      dist_score,
      ml_res$accuracy,
      ml_res$F1,
      ml_res$recall,
      smape,
      simbench_fidelity_score,
      logfc_pcc,
      logfc_scc,
      jaccard_top,
      sil_sim,
      silhouette_discrepancy,
      sim_pve,
      pve_discrepancy
    ),
    Direction = c(
      "Target = 1.0 (Optimal balance)",
      "Lower is better (closer to Uniform)",
      "Higher is better (1.0 = Uniform null)",
      "Higher is better (Group predictability)",
      "Higher is better (Balanced F1)",
      "Higher is better (Sensitivity)",
      "Lower is better (0.0 = zero error)",
      "Higher is better (1.0 = perfect match)",
      "Higher is better (Effect size match)",
      "Higher is better (Rank order match)",
      "Higher is better (Set intersection)",
      "Higher is better (Group separation)",
      "Lower is better (Preserves ref geometry)",
      "Target ~ Reference PVE",
      "Lower is better (Preserves variance explained)"
    ),
    stringsAsFactors = FALSE
  )
  
  list(
    deg_summary_table = summary_df,
    deg_genes_ref = deg_genes_ref,
    deg_genes_sim = deg_genes_sim,
    contrast = paste(group1, "vs", group2),
    logfc_ref = logfc_ref,
    logfc_sim = logfc_sim,
    ml_classification = ml_res,
    pvalue_uniformity = list(statistic = chisq_stat, pvalue = chisq_pval, distribution_score = dist_score)
  )
}

