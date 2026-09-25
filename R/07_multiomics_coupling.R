#' @title Single-Cell Multiomics Cross-Modality Coupling
#' @description Evaluates inter-modality coupling, cross-modal label transfer,
#'   FOSCTTM alignment, cross-modal generation fidelity, joint embedding alignment,
#'   peak-to-gene ATAC-RNA regulatory linkage, chromatin peak co-accessibility,
#'   and co-regulation/co-expression modularity in multiomics simulations.
#' @name multiomics_coupling
#' @export calc_cross_modality_correlation
#' @export evaluate_cross_modal_prediction
#' @export calc_foscttm
#' @export calc_cross_modal_generation
#' @export calc_modality_alignment
#' @export calc_network_jaccard
#' @export calc_atac_rna_coupling
#' @export calc_coregulation_fidelity
#' @export calc_peak_coaccessibility_fidelity
#' @export calc_coexpression_module_fidelity
#' @export calc_accessibility_profile_concordance
NULL

#' Calculate Cross-Modality Correlation Fidelity
#'
#' Evaluates whether cross-modality feature relationships (e.g. peak-to-gene links,
#' promoter accessibility vs. gene expression, or mRNA vs. surface protein) are accurately
#' captured by the simulation method compared to the empirical reference.
#'
#' @param ref_mod1 Reference matrix for Modality 1 (features x cells).
#' @param ref_mod2 Reference matrix for Modality 2 (features x cells).
#' @param sim_mod1 Simulated matrix for Modality 1.
#' @param sim_mod2 Simulated matrix for Modality 2.
#' @param feature_pairs Optional 2-column data.frame/matrix of paired feature names/indices.
#' @param method Correlation method: "spearman" (default) or "pearson".
#'
#' @return A named list of the 7 univariate accuracy metrics comparing the reference vs.
#'   simulated cross-modality correlation distributions.
#' @export
calc_cross_modality_correlation <- function(
  ref_mod1,
  ref_mod2,
  sim_mod1,
  sim_mod2,
  feature_pairs = NULL,
  method = c("spearman", "pearson")
) {
  method <- match.arg(method)
  
  ref_mod1 <- as.matrix(ref_mod1)
  ref_mod2 <- as.matrix(ref_mod2)
  sim_mod1 <- as.matrix(sim_mod1)
  sim_mod2 <- as.matrix(sim_mod2)
  
  if (is.null(feature_pairs)) {
    n_pairs <- min(nrow(ref_mod1), nrow(ref_mod2), nrow(sim_mod1), nrow(sim_mod2))
    idx1 <- seq_len(n_pairs)
    idx2 <- seq_len(n_pairs)
  } else {
    feature_pairs <- as.data.frame(feature_pairs)
    idx1 <- feature_pairs[[1]]
    idx2 <- feature_pairs[[2]]
  }
  
  # Compute pairwise correlations in reference
  ref_cors <- sapply(seq_along(idx1), function(i) {
    x <- ref_mod1[idx1[i], ]
    y <- ref_mod2[idx2[i], ]
    if (stats::sd(x) > 0 && stats::sd(y) > 0) {
      stats::cor(x, y, method = method)
    } else {
      0
    }
  })
  
  # Compute pairwise correlations in simulation
  sim_cors <- sapply(seq_along(idx1), function(i) {
    x <- sim_mod1[idx1[i], ]
    y <- sim_mod2[idx2[i], ]
    if (stats::sd(x) > 0 && stats::sd(y) > 0) {
      stats::cor(x, y, method = method)
    } else {
      0
    }
  })
  
  calc_all_univariate_metrics(ref_cors, sim_cors, metric_prefix = "cross_modality_cor")
}

#' Cross-Modal Cell-Type Label Transfer Accuracy
#'
#' Trains a model using ONLY Modality 1 (e.g., scATAC-seq) and predicts cell-type
#' identities defined in Modality 2 (e.g., scRNA-seq) to test cross-modal biological coherence.
#'
#' @param mod1_data Matrix for Modality 1 (features x cells).
#' @param cell_types Ground truth cell-type vector from Modality 2.
#' @return A list with cross-modal classification accuracy and macro F1 score.
#' @export
evaluate_cross_modal_prediction <- function(mod1_data, cell_types) {
  cell_types <- as.factor(cell_types)
  data_t <- t(as.matrix(mod1_data))
  
  set.seed(42)
  n <- nrow(data_t)
  train_idx <- sample(seq_len(n), round(0.8 * n))
  
  train_x <- data_t[train_idx, , drop = FALSE]
  test_x <- data_t[-train_idx, , drop = FALSE]
  train_y <- cell_types[train_idx]
  test_y <- cell_types[-train_idx]
  
  # k-NN prediction
  pred_y <- class::knn(train = train_x, test = test_x, cl = train_y, k = min(5, nrow(train_x)))
  
  cm <- table(Predicted = pred_y, Actual = test_y)
  acc <- sum(diag(cm)) / max(1, sum(cm))
  
  cls <- levels(cell_types)
  prs <- vapply(cls, function(c) {
    d <- sum(pred_y == c)
    if (d > 0) sum(pred_y == c & test_y == c) / d else 0
  }, numeric(1))
  
  res <- vapply(cls, function(c) {
    d <- sum(test_y == c)
    if (d > 0) sum(pred_y == c & test_y == c) / d else 0
  }, numeric(1))
  
  f1s <- ifelse(prs + res > 0, 2 * (prs * res) / (prs + res), 0)
  
  list(
    cross_modal_accuracy = acc,
    cross_modal_F1 = mean(f1s, na.rm = TRUE)
  )
}

# ==============================================================================
# Cross-Modality Integration & In Silico Generation Benchmarks (scCross Metrics)
# Reference: Zhai et al., Genome Biology 2024, 25:178 (scCross framework)
# ==============================================================================

#' Fraction of Samples Closer Than The True Match (FOSCTTM)
#'
#' Computes the Fraction of Samples Closer Than The True Match (FOSCTTM) to quantify
#' cell alignment error across single-cell multi-omics modalities in a shared latent
#' or integrated space (Zhai et al., Genome Biology 2024; Liu et al., Nat Biotechnol 2023).
#' A score of 0 represents perfect alignment (where the true paired cell is the nearest neighbor),
#' while 0.5 corresponds to random chance.
#'
#' @param x Coordinates matrix for Modality 1 (cells x dimensions).
#' @param y Coordinates matrix for Modality 2 (cells x dimensions), where row i of y
#'   corresponds to the true match of row i of x.
#' @param metric Distance metric to evaluate: "euclidean" (default) or "cosine".
#'
#' @return A list containing:
#'   \item{foscttm}{Bidirectional mean FOSCTTM score across all cells (lower is better, 0 to 0.5).}
#'   \item{foscttm_xy}{Directional FOSCTTM from Modality 1 to Modality 2.}
#'   \item{foscttm_yx}{Directional FOSCTTM from Modality 2 to Modality 1.}
#'   \item{match_at_1}{Top-1 match rate (proportion of cells where the true match is rank 1).}
#'   \item{match_at_5}{Top-5 match rate (proportion of cells where the true match is in top 5).}
#'   \item{cell_foscttm}{Vector of bidirectional FOSCTTM scores for individual cells.}
#' @export
calc_foscttm <- function(x, y, metric = c("euclidean", "cosine")) {
  metric <- match.arg(metric)
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  if (nrow(x) != nrow(y)) {
    min_r <- min(nrow(x), nrow(y))
    x <- x[seq_len(min_r), , drop = FALSE]
    y <- y[seq_len(min_r), , drop = FALSE]
  }
  if (ncol(x) != ncol(y)) {
    min_c <- min(ncol(x), ncol(y))
    x <- x[, seq_len(min_c), drop = FALSE]
    y <- y[, seq_len(min_c), drop = FALSE]
  }
  
  N <- nrow(x)
  if (N < 2) {
    return(list(foscttm = 0, foscttm_xy = 0, foscttm_yx = 0, match_at_1 = 1, match_at_5 = 1, cell_foscttm = 0))
  }
  
  if (metric == "cosine") {
    norm_x <- sqrt(rowSums(x^2))
    norm_y <- sqrt(rowSums(y^2))
    norm_x[norm_x == 0] <- 1e-10
    norm_y[norm_y == 0] <- 1e-10
    x_norm <- sweep(x, 1, norm_x, "/")
    y_norm <- sweep(y, 1, norm_y, "/")
    d_mat <- 1 - (x_norm %*% t(y_norm))
  } else {
    x_sq <- rowSums(x^2)
    y_sq <- rowSums(y^2)
    d_sq <- outer(x_sq, y_sq, "+") - 2 * (x %*% t(y))
    d_sq[d_sq < 0] <- 0
    d_mat <- sqrt(d_sq)
  }
  
  true_dist <- diag(d_mat)
  
  comp_xy <- d_mat < true_dist
  foscttm_xy <- rowMeans(comp_xy)
  
  comp_yx <- t(d_mat) < true_dist
  foscttm_yx <- rowMeans(comp_yx)
  
  cell_foscttm <- (foscttm_xy + foscttm_yx) / 2
  mean_foscttm <- mean(cell_foscttm)
  
  ranks_xy <- rowSums(comp_xy)
  match_at_1 <- mean(ranks_xy == 0)
  match_at_5 <- mean(ranks_xy < 5)
  
  list(
    foscttm = as.numeric(mean_foscttm),
    foscttm_xy = as.numeric(mean(foscttm_xy)),
    foscttm_yx = as.numeric(mean(foscttm_yx)),
    match_at_1 = as.numeric(match_at_1),
    match_at_5 = as.numeric(match_at_5),
    cell_foscttm = cell_foscttm
  )
}

#' Cross-Modality In Silico Generation & Translation Fidelity
#'
#' Evaluates the biological accuracy and reconstruction fidelity of cross-modality
#' generation (e.g. predicting scATAC from scRNA or vice versa; Zhai et al., Genome Biology 2024).
#' Computes cell-wise and feature-wise Pearson and Spearman correlations, cosine similarity,
#' RMSE, and MAE between predicted and measured multi-omics profiles.
#'
#' @param true_data Matrix or data.frame of measured features across cells (features x cells).
#' @param pred_data Matrix or data.frame of in-silico generated / predicted features (features x cells).
#'
#' @return A list containing cell-wise and feature-wise correlation, cosine similarity, RMSE, and MAE.
#' @export
calc_cross_modal_generation <- function(true_data, pred_data) {
  true_data <- as.matrix(true_data)
  pred_data <- as.matrix(pred_data)
  
  if (all(dim(true_data) != dim(pred_data))) {
    if (nrow(true_data) == ncol(pred_data) && ncol(true_data) == nrow(pred_data)) {
      pred_data <- t(pred_data)
    } else {
      min_r <- min(nrow(true_data), nrow(pred_data))
      min_c <- min(ncol(true_data), ncol(pred_data))
      true_data <- true_data[seq_len(min_r), seq_len(min_c), drop = FALSE]
      pred_data <- pred_data[seq_len(min_r), seq_len(min_c), drop = FALSE]
    }
  }
  
  G <- nrow(true_data)
  N <- ncol(true_data)
  
  pcc_cell <- vapply(seq_len(N), function(i) {
    u <- true_data[, i]
    v <- pred_data[, i]
    if (stats::sd(u) > 0 && stats::sd(v) > 0) {
      stats::cor(u, v, method = "pearson")
    } else 0
  }, numeric(1))
  
  scc_cell <- vapply(seq_len(N), function(i) {
    u <- true_data[, i]
    v <- pred_data[, i]
    if (stats::sd(u) > 0 && stats::sd(v) > 0) {
      stats::cor(u, v, method = "spearman")
    } else 0
  }, numeric(1))
  
  pcc_feat <- vapply(seq_len(G), function(g) {
    u <- true_data[g, ]
    v <- pred_data[g, ]
    if (stats::sd(u) > 0 && stats::sd(v) > 0) {
      stats::cor(u, v, method = "pearson")
    } else 0
  }, numeric(1))
  
  scc_feat <- vapply(seq_len(G), function(g) {
    u <- true_data[g, ]
    v <- pred_data[g, ]
    if (stats::sd(u) > 0 && stats::sd(v) > 0) {
      stats::cor(u, v, method = "spearman")
    } else 0
  }, numeric(1))
  
  diff_mat <- true_data - pred_data
  rmse <- sqrt(mean(diff_mat^2))
  mae <- mean(abs(diff_mat))
  
  cos_sim <- vapply(seq_len(N), function(i) {
    u <- true_data[, i]
    v <- pred_data[, i]
    norm_u <- sqrt(sum(u^2))
    norm_v <- sqrt(sum(v^2))
    if (norm_u > 0 && norm_v > 0) {
      sum(u * v) / (norm_u * norm_v)
    } else 0
  }, numeric(1))
  
  list(
    mean_cell_pcc = as.numeric(mean(pcc_cell, na.rm = TRUE)),
    median_cell_pcc = as.numeric(stats::median(pcc_cell, na.rm = TRUE)),
    mean_cell_scc = as.numeric(mean(scc_cell, na.rm = TRUE)),
    median_cell_scc = as.numeric(stats::median(scc_cell, na.rm = TRUE)),
    mean_feat_pcc = as.numeric(mean(pcc_feat, na.rm = TRUE)),
    median_feat_pcc = as.numeric(stats::median(pcc_feat, na.rm = TRUE)),
    mean_feat_scc = as.numeric(mean(scc_feat, na.rm = TRUE)),
    median_feat_scc = as.numeric(stats::median(scc_feat, na.rm = TRUE)),
    mean_cell_cosine = as.numeric(mean(cos_sim, na.rm = TRUE)),
    rmse = as.numeric(rmse),
    mae = as.numeric(mae)
  )
}

#' Modality Alignment & Omics Layer Mixing in Joint Latent Space
#'
#' Quantifies how well different single-cell omics layers mix in a unified latent embedding space
#' without sacrificing biological clustering (Zhai et al., Genome Biology 2024; Luecken et al., Nat Methods 2022).
#'
#' @param embedding Joint coordinates matrix (cells x latent_dimensions).
#' @param modalities Factor or character vector indicating the modality of each cell (e.g. "RNA" vs "ATAC").
#' @param cell_types Optional factor or character vector of cell type labels.
#' @param k Number of nearest neighbors for neighborhood connectivity calculation (default 15).
#'
#' @return A list containing modality_asw, modality_mixing_score, and mean_cross_modality_neighbor_frac.
#' @export
calc_modality_alignment <- function(embedding, modalities, cell_types = NULL, k = 15) {
  embedding <- as.matrix(embedding)
  modalities <- as.factor(modalities)
  N <- nrow(embedding)
  
  dist_mat <- stats::dist(embedding)
  asw_modality <- if (requireNamespace("cluster", quietly = TRUE)) {
    tryCatch({
      sil <- cluster::silhouette(as.integer(modalities), dist_mat)
      mean(sil[, 3], na.rm = TRUE)
    }, error = function(e) NA_real_)
  } else {
    calc_silhouette(dist_mat, modalities)
  }
  
  modality_mixing_score <- if (!is.na(asw_modality)) 1 - abs(asw_modality) else NA_real_
  
  dist_full <- as.matrix(dist_mat)
  diag(dist_full) <- Inf
  
  k <- min(k, N - 1)
  knn_idx <- t(apply(dist_full, 1, order)[seq_len(k), , drop = FALSE])
  
  cross_mod_frac <- vapply(seq_len(N), function(i) {
    mean(modalities[knn_idx[i, ]] != modalities[i])
  }, numeric(1))
  
  list(
    modality_asw = as.numeric(asw_modality),
    modality_mixing_score = as.numeric(modality_mixing_score),
    mean_cross_modality_neighbor_frac = as.numeric(mean(cross_mod_frac))
  )
}

# Standardize edge list helper
standardize_edges <- function(edges, directed = TRUE) {
  if (is.data.frame(edges)) {
    cols <- tolower(colnames(edges))
    g1_col <- which(cols %in% c("gene1", "from", "source", "regulator", "tf"))[1]
    g2_col <- which(cols %in% c("gene2", "to", "target"))[1]
    w_col  <- which(cols %in% c("edgeweight", "weight", "score", "importance", "abs_weight"))[1]
    
    if (is.na(g1_col) || is.na(g2_col)) {
      if (ncol(edges) >= 2) {
        g1_col <- 1
        g2_col <- 2
      } else {
        stop("edges data.frame must contain at least 2 columns")
      }
    }
    
    g1 <- as.character(edges[[g1_col]])
    g2 <- as.character(edges[[g2_col]])
    w <- if (!is.na(w_col)) as.numeric(edges[[w_col]]) else rep(1, length(g1))
    
    valid <- g1 != g2
    g1 <- g1[valid]
    g2 <- g2[valid]
    w <- w[valid]
    
    if (!directed) {
      swap <- g1 > g2
      tmp <- g1[swap]
      g1[swap] <- g2[swap]
      g2[swap] <- tmp
    }
    
    sep <- if (directed) "->" else "--"
    edge_id <- paste0(g1, sep, g2)
    df <- data.frame(
      Gene1 = g1,
      Gene2 = g2,
      EdgeWeight = w,
      EdgeID = edge_id,
      stringsAsFactors = FALSE
    )
    df$abs_weight <- abs(df$EdgeWeight)
    df <- df[order(df$abs_weight, decreasing = TRUE), ]
    df <- df[!duplicated(df$EdgeID), ]
    return(df)
  } else if (is.character(edges)) {
    return(unique(edges))
  } else {
    stop("Unsupported edge format: must be data.frame or character vector")
  }
}

#' Jaccard Similarity of Predicted Regulatory Networks / Links
#'
#' Computes the Jaccard similarity index of top-k edge predictions between two networks
#' or top correlated feature pairs between reference and simulated datasets.
#'
#' @param edges1 Data frame or character vector representing the first network.
#' @param edges2 Data frame or character vector representing the second network.
#' @param k Integer cutoff for top edge selection (optional).
#' @param directed Logical, whether edges are directed (default TRUE) or undirected (FALSE).
#'
#' @return Jaccard similarity index in [0, 1] (|E1 cap E2| / |E1 cup E2|).
#' @export
calc_network_jaccard <- function(edges1, edges2, k = NULL, directed = TRUE) {
  std1 <- standardize_edges(edges1, directed = directed)
  std2 <- standardize_edges(edges2, directed = directed)
  
  ids1 <- if (is.data.frame(std1)) std1$EdgeID else std1
  ids2 <- if (is.data.frame(std2)) std2$EdgeID else std2
  
  if (!is.null(k)) {
    ids1 <- ids1[seq_len(min(k, length(ids1)))]
    ids2 <- ids2[seq_len(min(k, length(ids2)))]
  }
  
  set1 <- unique(ids1)
  set2 <- unique(ids2)
  
  inter <- length(intersect(set1, set2))
  uni <- length(union(set1, set2))
  
  if (uni == 0) 1.0 else inter / uni
}

#' Chromatin Accessibility to RNA Regulatory Effect Coupling
#'
#' Evaluates whether simulated paired scATAC-seq and scRNA-seq profiles capture the true
#' regulatory coupling where regional chromatin opening gates target gene transcription
#' (scMultiSim; Li et al., Nat Methods 2023).
#'
#' @param atac_data Matrix of chromatin peak accessibilities (peaks x cells).
#' @param rna_data Matrix of gene expression counts (genes x cells).
#' @param linked_pairs Optional 2-column data.frame of known linked peak-gene pairs.
#'   If NULL, assumes 1-to-1 matching by row index.
#'
#' @return A list containing mean coupling correlation, positive coupling ratio, and mean R-squared.
#' @export
calc_atac_rna_coupling <- function(atac_data, rna_data, linked_pairs = NULL) {
  atac_data <- as.matrix(atac_data)
  rna_data <- as.matrix(rna_data)
  
  if (is.null(linked_pairs)) {
    n_pairs <- min(nrow(atac_data), nrow(rna_data))
    idx_atac <- seq_len(n_pairs)
    idx_rna <- seq_len(n_pairs)
  } else {
    linked_pairs <- as.data.frame(linked_pairs)
    idx_atac <- linked_pairs[[1]]
    idx_rna <- linked_pairs[[2]]
  }
  
  n_eval <- length(idx_atac)
  
  cors <- vapply(seq_len(n_eval), function(i) {
    a <- atac_data[idx_atac[i], ]
    r <- rna_data[idx_rna[i], ]
    if (stats::sd(a) > 0 && stats::sd(r) > 0) {
      stats::cor(a, r, method = "pearson")
    } else 0
  }, numeric(1))
  
  r2s <- vapply(seq_len(n_eval), function(i) {
    a <- atac_data[idx_atac[i], ]
    r <- rna_data[idx_rna[i], ]
    if (stats::sd(a) > 0 && stats::sd(r) > 0) {
      summary(stats::lm(r ~ a))$r.squared
    } else 0
  }, numeric(1))
  
  list(
    mean_coupling_cor = as.numeric(mean(cors)),
    median_coupling_cor = as.numeric(stats::median(cors)),
    positive_coupling_ratio = as.numeric(mean(cors > 0)),
    mean_r_squared = as.numeric(mean(r2s)),
    coupling_correlations = cors
  )
}

#' Multi-Omics Co-Regulation & Modularity Fidelity
#'
#' Evaluates the preservation of co-regulated gene/feature modules between reference
#' and simulated multi-omics datasets (Monzo et al., 2025; Arzalluz-Luque et al., 2022).
#'
#' @param ref_data Reference feature-by-cell matrix or data frame.
#' @param sim_data Simulated feature-by-cell matrix or data frame.
#' @param modules A list of character vectors representing feature clusters/modules,
#'   or a named vector/factor of module assignments per feature.
#' @param method Correlation method: "pearson" (default) or "spearman".
#'
#' @return A list containing module correlation r, RMSE, MAE, and modularity fidelity.
#' @export
calc_coregulation_fidelity <- function(
  ref_data,
  sim_data,
  modules,
  method = c("pearson", "spearman")
) {
  method <- match.arg(method)
  ref_mat <- as.matrix(ref_data)
  sim_mat <- as.matrix(sim_data)
  
  if (!is.list(modules)) {
    if (!is.null(names(modules))) {
      modules <- split(names(modules), as.character(modules))
    } else {
      stop("modules must be a list of feature vectors or a named factor/vector.")
    }
  }
  
  all_mod_genes <- unlist(modules)
  common_genes <- intersect(all_mod_genes, intersect(rownames(ref_mat), rownames(sim_mat)))
  if (length(common_genes) < 4) {
    # If gene names don't match, use top rows
    n_use <- min(length(all_mod_genes), nrow(ref_mat), nrow(sim_mat))
    common_genes <- rownames(ref_mat)[seq_len(n_use)]
    modules <- list(Module_1 = common_genes[seq_len(n_use %/% 2)],
                    Module_2 = common_genes[(n_use %/% 2 + 1):n_use])
  }
  
  modules_filt <- lapply(modules, function(m) intersect(m, common_genes))
  modules_filt <- modules_filt[vapply(modules_filt, length, integer(1)) >= 2]
  
  if (length(modules_filt) == 0) {
    stop("No modules have at least 2 common features between ref and sim data.")
  }
  
  genes_eval <- unlist(modules_filt)
  c_ref <- stats::cor(t(ref_mat[genes_eval, ]), method = method)
  c_sim <- stats::cor(t(sim_mat[genes_eval, ]), method = method)
  c_ref[is.na(c_ref)] <- 0
  c_sim[is.na(c_sim)] <- 0
  
  gene_to_mod <- rep(names(modules_filt), times = vapply(modules_filt, length, integer(1)))
  names(gene_to_mod) <- genes_eval
  
  is_intra <- outer(gene_to_mod, gene_to_mod, "==")
  diag(is_intra) <- FALSE
  
  is_inter <- outer(gene_to_mod, gene_to_mod, "!=")
  
  ref_intra_vals <- c_ref[is_intra]
  sim_intra_vals <- c_sim[is_intra]
  
  r_val <- stats::cor(ref_intra_vals, sim_intra_vals, method = "pearson")
  rmse_val <- sqrt(mean((ref_intra_vals - sim_intra_vals)^2))
  mae_val <- mean(abs(ref_intra_vals - sim_intra_vals))
  
  ref_mean_intra <- mean(ref_intra_vals, na.rm = TRUE)
  sim_mean_intra <- mean(sim_intra_vals, na.rm = TRUE)
  
  ref_mean_inter <- if (any(is_inter)) mean(abs(c_ref[is_inter]), na.rm = TRUE) else 1e-4
  sim_mean_inter <- if (any(is_inter)) mean(abs(c_sim[is_inter]), na.rm = TRUE) else 1e-4
  
  ref_modularity <- ref_mean_intra / max(ref_mean_inter, 1e-4)
  sim_modularity <- sim_mean_intra / max(sim_mean_inter, 1e-4)
  
  modularity_fid <- if (ref_modularity > 0 && sim_modularity > 0) {
    min(ref_modularity, sim_modularity) / max(ref_modularity, sim_modularity)
  } else 0
  
  list(
    module_correlation_r = as.numeric(r_val),
    module_correlation_rmse = as.numeric(rmse_val),
    module_correlation_mae = as.numeric(mae_val),
    ref_modularity_ratio = as.numeric(ref_modularity),
    sim_modularity_ratio = as.numeric(sim_modularity),
    modularity_fidelity = as.numeric(modularity_fid)
  )
}

#' Evaluate Single-Cell ATAC Peak Co-Accessibility Fidelity (SCRIP)
#'
#' Evaluates the preservation of chromatin peak-to-peak co-accessibility correlation
#' matrices (cis-regulatory interactions) between reference and simulated scATAC-seq datasets.
#'
#' @param ref_atac Reference scATAC-seq matrix (peaks x cells).
#' @param sim_atac Simulated scATAC-seq matrix (peaks x cells).
#' @param top_n_peaks Number of highest variance peaks to evaluate (default 300).
#' @param method Correlation method: "spearman" (default) or "pearson".
#'
#' @return A list containing RV coefficient, matrix correlation, Frobenius distance, and MAE.
#' @export
calc_peak_coaccessibility_fidelity <- function(
  ref_atac,
  sim_atac,
  top_n_peaks = 300,
  method = c("spearman", "pearson")
) {
  method <- match.arg(method)
  ref_mat <- as.matrix(ref_atac)
  sim_mat <- as.matrix(sim_atac)
  
  common_peaks <- intersect(rownames(ref_mat), rownames(sim_mat))
  if (length(common_peaks) < 10) {
    n_use <- min(top_n_peaks, nrow(ref_mat), nrow(sim_mat))
    ref_vars <- apply(ref_mat, 1, stats::var)
    sim_vars <- apply(sim_mat, 1, stats::var)
    ref_sub <- ref_mat[order(ref_vars, decreasing = TRUE)[seq_len(n_use)], ]
    sim_sub <- sim_mat[order(sim_vars, decreasing = TRUE)[seq_len(n_use)], ]
  } else {
    ref_vars <- apply(ref_mat[common_peaks, , drop = FALSE], 1, stats::var)
    n_use <- min(top_n_peaks, length(common_peaks))
    top_p <- names(sort(ref_vars, decreasing = TRUE))[seq_len(n_use)]
    ref_sub <- ref_mat[top_p, , drop = FALSE]
    sim_sub <- sim_mat[top_p, , drop = FALSE]
  }
  
  cor_ref <- stats::cor(t(ref_sub), method = method)
  cor_sim <- stats::cor(t(sim_sub), method = method)
  cor_ref[is.na(cor_ref)] <- 0
  cor_sim[is.na(cor_sim)] <- 0
  
  tr <- function(m) sum(diag(m))
  rv <- tr(cor_ref %*% cor_sim) / sqrt(tr(cor_ref %*% cor_ref) * tr(cor_sim %*% cor_sim) + 1e-12)
  
  mask <- upper.tri(cor_ref)
  vals_ref <- cor_ref[mask]
  vals_sim <- cor_sim[mask]
  
  r_val <- stats::cor(vals_ref, vals_sim, method = "pearson")
  rho_val <- stats::cor(vals_ref, vals_sim, method = "spearman")
  frobenius_dist <- sqrt(sum((cor_ref - cor_sim)^2)) / max(1, length(vals_ref))
  mae_val <- mean(abs(vals_ref - vals_sim))
  
  list(
    rv_coefficient = as.numeric(rv),
    coaccessibility_pearson = as.numeric(r_val),
    coaccessibility_spearman = as.numeric(rho_val),
    frobenius_distance = as.numeric(frobenius_dist),
    coaccessibility_mae = as.numeric(mae_val)
  )
}

#' Evaluate Gene Co-Expression Module Fidelity (ESCO)
#'
#' Automatically discovers or evaluates co-expression modules and measures
#' the preservation of intra-module vs. inter-module correlation structure.
#'
#' @param ref_mat Reference gene expression matrix (genes x cells).
#' @param sim_mat Simulated gene expression matrix (genes x cells).
#' @param modules Optional named list of gene modules. If NULL, auto-detected via hierarchical clustering.
#' @param n_modules Number of modules to detect if modules is NULL (default 5).
#' @param top_genes Number of top variable genes to consider for clustering (default 200).
#'
#' @return A list summarizing module correlation preservation, RMSE, and modularity ratio fidelity.
#' @export
calc_coexpression_module_fidelity <- function(
  ref_mat,
  sim_mat,
  modules = NULL,
  n_modules = 5,
  top_genes = 200
) {
  ref_m <- as.matrix(ref_mat)
  sim_m <- as.matrix(sim_mat)
  
  common_genes <- intersect(rownames(ref_m), rownames(sim_m))
  if (length(common_genes) < 20) {
    n_use <- min(top_genes, nrow(ref_m), nrow(sim_m))
    ref_m <- ref_m[seq_len(n_use), ]
    sim_m <- sim_m[seq_len(n_use), ]
    common_genes <- rownames(ref_m)
  }
  
  if (is.null(modules)) {
    gene_vars <- apply(ref_m[common_genes, , drop = FALSE], 1, stats::var)
    n_sel <- min(top_genes, length(common_genes))
    top_g <- names(sort(gene_vars, decreasing = TRUE))[seq_len(n_sel)]
    
    sub_ref <- ref_m[top_g, , drop = FALSE]
    cor_ref <- stats::cor(t(sub_ref), method = "spearman")
    cor_ref[is.na(cor_ref)] <- 0
    dist_mat <- stats::as.dist(1 - cor_ref)
    hc <- stats::hclust(dist_mat, method = "average")
    clusters <- stats::cutree(hc, k = min(n_modules, max(2, length(top_g) %/% 4)))
    modules <- split(names(clusters), paste0("Module_", clusters))
  }
  
  calc_coregulation_fidelity(ref_m, sim_m, modules = modules, method = "spearman")
}

#' Evaluate Chromatin Accessibility Profile Concordance (DiTSim)
#'
#' Computes global and cell-type-stratified Pearson and Spearman correlations
#' and Kullback-Leibler (KL) divergence between empirical reference and simulated
#' single-cell chromatin accessibility profiles.
#'
#' @param ref_data Reference count/accessibility matrix (peaks x cells).
#' @param sim_data Simulated count/accessibility matrix (peaks x cells).
#' @param cell_types Optional factor or vector of cell-type annotations for cells.
#' @param use_tfidf Logical, whether to apply TF-IDF transformation prior to mean profile calculation (default TRUE).
#'
#' @return A list containing global PCC, global SCC, per-cell-type mean PCC/SCC,
#'   and mean accessibility KL divergence.
#' @export
calc_accessibility_profile_concordance <- function(
  ref_data,
  sim_data,
  cell_types = NULL,
  use_tfidf = TRUE
) {
  r_mat <- as.matrix(ref_data)
  s_mat <- as.matrix(sim_data)
  
  common_feats <- intersect(rownames(r_mat), rownames(s_mat))
  if (length(common_feats) >= 10) {
    r_mat <- r_mat[common_feats, , drop = FALSE]
    s_mat <- s_mat[common_feats, , drop = FALSE]
  } else {
    n_use <- min(nrow(r_mat), nrow(s_mat))
    r_mat <- r_mat[seq_len(n_use), , drop = FALSE]
    s_mat <- s_mat[seq_len(n_use), , drop = FALSE]
  }
  
  apply_tfidf <- function(mat) {
    tf <- sweep(mat, 2, pmax(colSums(mat), 1), "/")
    n_cells <- ncol(mat)
    peak_counts <- rowSums(mat > 0)
    idf <- log(1 + n_cells / pmax(peak_counts, 1))
    tfidf <- sweep(tf, 1, idf, "*")
    log(1 + 1e4 * tfidf)
  }
  
  if (use_tfidf) {
    r_proc <- apply_tfidf(r_mat)
    s_proc <- apply_tfidf(s_mat)
  } else {
    r_proc <- r_mat
    s_proc <- s_mat
  }
  
  mean_ref <- rowMeans(r_proc)
  mean_sim <- rowMeans(s_proc)
  
  global_pcc <- stats::cor(mean_ref, mean_sim, method = "pearson")
  global_scc <- stats::cor(mean_ref, mean_sim, method = "spearman")
  
  p <- mean_ref / max(sum(mean_ref), 1e-12)
  q <- mean_sim / max(sum(mean_sim), 1e-12)
  eps <- 1e-10
  p <- (p + eps) / sum(p + eps)
  q <- (q + eps) / sum(q + eps)
  kl_div <- sum(p * log(p / q))
  
  ct_pcc <- NA_real_
  ct_scc <- NA_real_
  ct_summary <- NULL
  
  if (!is.null(cell_types)) {
    cts <- as.factor(cell_types)
    if (length(cts) == ncol(r_proc) && length(cts) == ncol(s_proc)) {
      u_cts <- levels(cts)
      pccs <- numeric(length(u_cts))
      sccs <- numeric(length(u_cts))
      names(pccs) <- names(sccs) <- u_cts
      
      for (cl in u_cts) {
        idx <- which(cts == cl)
        if (length(idx) >= 2) {
          m_r <- rowMeans(r_proc[, idx, drop = FALSE])
          m_s <- rowMeans(s_proc[, idx, drop = FALSE])
          pccs[cl] <- stats::cor(m_r, m_s, method = "pearson")
          sccs[cl] <- stats::cor(m_r, m_s, method = "spearman")
        } else {
          pccs[cl] <- NA_real_
          sccs[cl] <- NA_real_
        }
      }
      ct_pcc <- mean(pccs, na.rm = TRUE)
      ct_scc <- mean(sccs, na.rm = TRUE)
      ct_summary <- data.frame(
        CellType = u_cts,
        PCC = pccs,
        SCC = sccs,
        stringsAsFactors = FALSE
      )
    }
  }
  
  list(
    global_pcc = as.numeric(global_pcc),
    global_scc = as.numeric(global_scc),
    kl_divergence = as.numeric(kl_div),
    mean_celltype_pcc = as.numeric(ct_pcc),
    mean_celltype_scc = as.numeric(ct_scc),
    celltype_summary = ct_summary
  )
}
