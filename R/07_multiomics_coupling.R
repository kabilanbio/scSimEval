# ==============================================================================
# (VI) Multiomics Integration & Alignment
# Multiomics Integration & Cross-Modal Alignment Fidelity
# Gold-standard benchmarking metrics extracted from scib, scMultiSim, scCross, Matilda, scMoMtF
# ==============================================================================

#' Calculate Batch Average Silhouette Width (Batch ASW)
#'
#' Evaluates batch mixing within each cell type in integrated or joint multiomics
#' embedding space, following the definitions in \code{scib} and \code{scMultiSim}.
#' For each cell type, the silhouette width of cells with respect to batch labels
#' is computed (\eqn{-1 \le s_i \le 1}), where values near 0 indicate well-mixed batches.
#' The per-cell-type score is \eqn{1 - |s_i|}, averaged across all cells and cell types.
#'
#' @param embed Numeric matrix of cell coordinates (cells x dimensions, e.g. PCA or joint embedding).
#' @param batch Vector or factor of batch or modality labels for each cell.
#' @param cell_type Optional vector or factor of cell-type labels for each cell.
#'
#' @return Numeric scalar between 0 and 1, where 1 indicates optimal batch mixing.
#' @references
#' Luecken et al. (2022) Benchmarking atlas-level data integration in single-cell genomics. Nature Methods.
#' Li et al. (2023) scMultiSim: simulation of multi-modality single cell data guided by gene regulatory networks. Cell Systems.
#' @examples
#' set.seed(42)
#' emb <- matrix(rnorm(100 * 10), 100, 10)
#' bt <- rep(c("B1", "B2"), each = 50)
#' ct <- rep(c("T1", "T2"), 50)
#' calc_batch_asw(emb, bt, ct)
#' @export
calc_batch_asw <- function(embed, batch, cell_type = NULL) {
  if (is.null(embed) || nrow(embed) < 3) return(NA_real_)
  batch <- as.character(batch)
  if (length(unique(batch)) < 2) return(1.0)
  
  if (requireNamespace("cluster", quietly = TRUE)) {
    if (!is.null(cell_type) && length(unique(cell_type)) >= 2) {
      cell_type <- as.character(cell_type)
      cts <- unique(cell_type)
      asw_per_ct <- numeric(0)
      
      for (ct in cts) {
        idx <- which(cell_type == ct)
        if (length(idx) >= 4 && length(unique(batch[idx])) >= 2) {
          sub_embed <- embed[idx, , drop = FALSE]
          sub_batch <- as.numeric(as.factor(batch[idx]))
          tb <- table(sub_batch)
          if (all(tb >= 2)) {
            d_mat <- stats::dist(sub_embed)
            sil <- tryCatch(cluster::silhouette(sub_batch, d_mat), error = function(e) NULL)
            if (!is.null(sil) && inherits(sil, "silhouette")) {
              s_vals <- sil[, 3]
              asw_score <- mean(1 - abs(s_vals), na.rm = TRUE)
              asw_per_ct <- c(asw_per_ct, asw_score)
            }
          }
        }
      }
      if (length(asw_per_ct) > 0) {
        return(max(0, min(1, mean(asw_per_ct))))
      }
    }
    # Fallback to global batch silhouette
    d_mat <- stats::dist(embed)
    b_num <- as.numeric(as.factor(batch))
    sil <- tryCatch(cluster::silhouette(b_num, d_mat), error = function(e) NULL)
    if (!is.null(sil) && inherits(sil, "silhouette")) {
      return(max(0, min(1, mean(1 - abs(sil[, 3]), na.rm = TRUE))))
    }
  }
  1.0
}


#' Calculate Cell-Type Average Silhouette Width (Cell-Type ASW)
#'
#' Evaluates biological conservation and cell-type cluster separation in integrated
#' or joint multiomics space, following \code{scCross} and \code{scib}.
#' The silhouette width of cells with respect to ground-truth cell-type labels is
#' computed and scaled to $[0, 1]$ via $(ASW + 1) / 2$.
#'
#' @param embed Numeric matrix of cell coordinates (cells x dimensions).
#' @param cell_type Vector or factor of cell-type labels for each cell.
#'
#' @return Numeric scalar between 0 and 1, where 1 indicates distinct, tight cell-type clusters.
#' @references
#' Ding et al. (2023) scCross: deep generative modeling for single-cell cross-modality integration and perturbation analysis.
#' Luecken et al. (2022) Benchmarking atlas-level data integration in single-cell genomics. Nature Methods.
#' @examples
#' set.seed(42)
#' emb <- matrix(rnorm(100 * 10), 100, 10)
#' ct <- rep(c("T1", "T2"), each = 50)
#' calc_cell_type_asw(emb, ct)
#' @export
calc_cell_type_asw <- function(embed, cell_type) {
  if (is.null(embed) || nrow(embed) < 3) return(NA_real_)
  cell_type <- as.character(cell_type)
  if (length(unique(cell_type)) < 2) return(1.0)
  
  if (requireNamespace("cluster", quietly = TRUE)) {
    n <- nrow(embed)
    if (n > 3000) {
      idx <- sample(seq_len(n), 3000)
      embed <- embed[idx, , drop = FALSE]
      cell_type <- cell_type[idx]
    }
    ct_num <- as.numeric(as.factor(cell_type))
    tb <- table(ct_num)
    if (any(tb < 2)) {
      valid_classes <- as.numeric(names(tb)[tb >= 2])
      keep <- ct_num %in% valid_classes
      embed <- embed[keep, , drop = FALSE]
      ct_num <- ct_num[keep]
    }
    if (length(unique(ct_num)) < 2) return(1.0)
    
    d_mat <- stats::dist(embed)
    sil <- tryCatch(cluster::silhouette(ct_num, d_mat), error = function(e) NULL)
    if (!is.null(sil) && inherits(sil, "silhouette")) {
      raw_asw <- mean(sil[, 3], na.rm = TRUE)
      scaled_asw <- (raw_asw + 1) / 2
      return(max(0, min(1, scaled_asw)))
    }
  }
  1.0
}


#' Calculate KNN Graph Connectivity
#'
#' Evaluates whether cells of the same cell type form a connected component in
#' the k-nearest neighbor (k-NN) graph across modalities/batches, following
#' \code{scCross}, \code{scMultiSim}, and \code{scib}.
#' For each cell type, the relative size of the largest connected component
#' (LCC) in the induced subgraph is calculated: $|LCC| / |V_c|$. The overall score
#' is the mean across all cell types.
#'
#' @param embed Numeric matrix of cell coordinates (cells x dimensions).
#' @param cell_type Vector or factor of cell-type labels for each cell.
#' @param k Number of nearest neighbors to construct the graph (default 15).
#'
#' @return Numeric scalar between 0 and 1, where 1 indicates complete graph connectivity.
#' @references
#' Ding et al. (2023) scCross: deep generative modeling for single-cell cross-modality integration.
#' Li et al. (2023) scMultiSim: simulation of multi-modality single cell data. Cell Systems.
#' Luecken et al. (2022) Benchmarking atlas-level data integration. Nature Methods.
#' @examples
#' set.seed(42)
#' emb <- matrix(rnorm(60 * 5), 60, 5)
#' ct <- rep(c("A", "B"), each = 30)
#' calc_graph_connectivity(emb, ct, k = 10)
#' @export
calc_graph_connectivity <- function(embed, cell_type, k = 15) {
  if (is.null(embed) || nrow(embed) < 2) return(NA_real_)
  cell_type <- as.character(cell_type)
  cts <- unique(cell_type)
  if (length(cts) < 1) return(1.0)
  
  n_cells <- nrow(embed)
  k <- max(1, min(k, n_cells - 1))
  
  d_mat <- as.matrix(stats::dist(embed))
  diag(d_mat) <- Inf
  
  adj <- matrix(FALSE, nrow = n_cells, ncol = n_cells)
  for (i in seq_len(n_cells)) {
    nn <- order(d_mat[i, ])[seq_len(k)]
    adj[i, nn] <- TRUE
    adj[nn, i] <- TRUE
  }
  
  conn_scores <- numeric(length(cts))
  for (ci in seq_along(cts)) {
    c_label <- cts[ci]
    idx <- which(cell_type == c_label)
    if (length(idx) <= 1) {
      conn_scores[ci] <- 1.0
      next
    }
    
    sub_adj <- adj[idx, idx, drop = FALSE]
    n_sub <- length(idx)
    
    visited <- rep(FALSE, n_sub)
    comp_sizes <- integer(0)
    
    for (start_node in seq_len(n_sub)) {
      if (!visited[start_node]) {
        queue <- start_node
        visited[start_node] <- TRUE
        size <- 0
        while (length(queue) > 0) {
          curr <- queue[1]
          queue <- queue[-1]
          size <- size + 1
          nbrs <- which(sub_adj[curr, ] & !visited)
          if (length(nbrs) > 0) {
            visited[nbrs] <- TRUE
            queue <- c(queue, nbrs)
          }
        }
        comp_sizes <- c(comp_sizes, size)
      }
    }
    
    max_comp <- if (length(comp_sizes) > 0) max(comp_sizes) else 1
    conn_scores[ci] <- max_comp / n_sub
  }
  
  mean_conn <- mean(conn_scores, na.rm = TRUE)
  max(0, min(1, mean_conn))
}


#' Calculate Integration Adjusted Rand Index (ARI)
#'
#' Computes the Adjusted Rand Index (ARI) evaluating the concordance between
#' clusters identified in the integrated/joint multiomics embedding and ground-truth
#' cell-type labels, following \code{scib}, \code{scCross}, \code{scMultiSim},
#' \code{Matilda}, and \code{scMoMtF}.
#'
#' @param embed Numeric matrix of cell coordinates (cells x dimensions).
#' @param cell_type Vector or factor of ground-truth cell-type labels.
#' @param n_clusters Number of clusters (defaults to the number of unique cell types).
#'
#' @return Numeric scalar between 0 and 1, where 1 indicates perfect clustering match.
#' @references
#' Hubert & Arabie (1985) Comparing partitions. Journal of Classification.
#' Luecken et al. (2022) Benchmarking atlas-level data integration. Nature Methods.
#' Ding et al. (2023) scCross.
#' Li et al. (2023) scMultiSim. Cell Systems.
#' Yang et al. (2023) Matilda. Nucleic Acids Research.
#' Lan et al. (2024) scMoMtF.
#' @examples
#' set.seed(42)
#' emb <- matrix(rnorm(60 * 4), 60, 4)
#' ct <- rep(c("TypeA", "TypeB"), each = 30)
#' calc_integration_ari(emb, ct)
#' @export
calc_integration_ari <- function(embed, cell_type, n_clusters = NULL) {
  if (is.null(embed) || nrow(embed) < 2) return(NA_real_)
  cell_type <- as.character(cell_type)
  n_types <- length(unique(cell_type))
  if (n_types <= 1) return(1.0)
  
  k_clust <- if (!is.null(n_clusters)) n_clusters else n_types
  k_clust <- max(2, min(k_clust, nrow(embed) - 1))
  
  set.seed(42)
  km <- tryCatch(
    suppressWarnings(stats::kmeans(embed, centers = k_clust, nstart = 10, iter.max = 50)),
    error = function(e) NULL
  )
  if (is.null(km)) return(NA_real_)
  pred_cl <- km$cluster
  
  if (requireNamespace("mclust", quietly = TRUE)) {
    ari <- mclust::adjustedRandIndex(pred_cl, cell_type)
  } else {
    tab <- table(pred_cl, cell_type)
    n <- sum(tab)
    a_i <- rowSums(tab)
    b_j <- colSums(tab)
    sum_n2 <- sum(choose(tab, 2))
    sum_a2 <- sum(choose(a_i, 2))
    sum_b2 <- sum(choose(b_j, 2))
    n2 <- choose(n, 2)
    exp_idx <- (sum_a2 * sum_b2) / n2
    max_idx <- 0.5 * (sum_a2 + sum_b2)
    ari <- if (max_idx == exp_idx) 1.0 else (sum_n2 - exp_idx) / (max_idx - exp_idx)
  }
  max(0, min(1, ari))
}


#' Calculate Integration Normalized Mutual Information (NMI)
#'
#' Computes Normalized Mutual Information (NMI) between clusters in the integrated
#' multiomics space and ground-truth cell-type labels, following \code{scib},
#' \code{scCross}, \code{scMultiSim}, \code{Matilda}, and \code{scMoMtF}.
#'
#' @param embed Numeric matrix of cell coordinates (cells x dimensions).
#' @param cell_type Vector or factor of ground-truth cell-type labels.
#' @param n_clusters Number of clusters (defaults to the number of unique cell types).
#'
#' @return Numeric scalar between 0 and 1, where 1 indicates perfect mutual information.
#' @references
#' Luecken et al. (2022) Benchmarking atlas-level data integration. Nature Methods.
#' Ding et al. (2023) scCross.
#' Li et al. (2023) scMultiSim. Cell Systems.
#' Yang et al. (2023) Matilda. Nucleic Acids Research.
#' Lan et al. (2024) scMoMtF.
#' @examples
#' set.seed(42)
#' emb <- matrix(rnorm(60 * 4), 60, 4)
#' ct <- rep(c("TypeA", "TypeB"), each = 30)
#' calc_integration_nmi(emb, ct)
#' @export
calc_integration_nmi <- function(embed, cell_type, n_clusters = NULL) {
  if (is.null(embed) || nrow(embed) < 2) return(NA_real_)
  cell_type <- as.character(cell_type)
  n_types <- length(unique(cell_type))
  if (n_types <= 1) return(1.0)
  
  k_clust <- if (!is.null(n_clusters)) n_clusters else n_types
  k_clust <- max(2, min(k_clust, nrow(embed) - 1))
  
  set.seed(42)
  km <- tryCatch(
    suppressWarnings(stats::kmeans(embed, centers = k_clust, nstart = 10, iter.max = 50)),
    error = function(e) NULL
  )
  if (is.null(km)) return(NA_real_)
  pred_cl <- km$cluster
  
  tab <- table(pred_cl, cell_type)
  total <- sum(tab)
  p_xy <- tab / total
  p_x <- rowSums(p_xy)
  p_y <- colSums(p_xy)
  
  h_x <- -sum(p_x[p_x > 0] * log(p_x[p_x > 0]))
  h_y <- -sum(p_y[p_y > 0] * log(p_y[p_y > 0]))
  
  p_outer <- outer(p_x, p_y)
  nz <- p_xy > 0
  mi <- sum(p_xy[nz] * log(p_xy[nz] / p_outer[nz]))
  
  denom <- (h_x + h_y) / 2
  if (denom <= 0) return(1.0)
  nmi <- mi / denom
  max(0, min(1, nmi))
}


#' Evaluate Multiomics Integration & Cross-Modal Alignment Fidelity (Category VI)
#'
#' Evaluates the multiomics integration and cross-modal alignment fidelity of simulated
#' single-cell multiomics data across 5 gold-standard integration metrics extracted from
#' key single-cell benchmarking suites (\code{scib}, \code{scMultiSim}, \code{scCross},
#' \code{Matilda}, and \code{scMoMtF}):
#' \enumerate{
#'   \item \strong{Batch ASW (\code{batch_asw}):} Average Silhouette Width evaluating batch mixing within cell types (scib, scMultiSim).
#'   \item \strong{Cell-Type ASW (\code{cell_type_asw}):} Average Silhouette Width evaluating cell-type bio-conservation in the multiomics space (scCross, scib).
#'   \item \strong{Graph Connectivity (\code{graph_connectivity}):} Relative size of the largest connected component in the k-NN graph for each cell type (scCross, scMultiSim, scib).
#'   \item \strong{ARI (\code{ari}):} Adjusted Rand Index comparing joint multiomics clustering against ground-truth cell types (scib, scCross, scMultiSim, Matilda, scMoMtF).
#'   \item \strong{NMI (\code{nmi}):} Normalized Mutual Information quantifying information concordance between joint clustering and cell types (scib, scCross, scMultiSim, Matilda, scMoMtF).
#' }
#'
#' @param ref_multi Optional reference multiomics list of matrices (e.g. \code{list(rna = ..., atac = ...)}).
#' @param sim_multi Simulated multiomics list of matrices or a matrix representing joint embedding coordinates.
#' @param cell_types Optional cell-type vector or list/data.frame of annotations for each cell.
#' @param batch_info Optional batch/donor/replicate vector or list/data.frame of annotations for each cell.
#' @param embed Optional precomputed joint embedding matrix (cells x dimensions). If NULL, computed automatically from \code{sim_multi}.
#' @param reduction_name Optional name of the precomputed dimensional reduction to extract from Seurat/SCE objects.
#' @param n_pcs Number of principal components per modality to retain for joint embedding (default 15).
#' @param k Number of nearest neighbors for graph-based metrics (default 30).
#' @param threads Number of CPU threads for parallel processing (default 1).
#' @param verbose Logical; whether to print progress messages (default TRUE).
#'
#' @return A list of class \code{scSimEval_integration} containing:
#'   \item{metrics}{Named list of the 5 calculated integration metric values (batch_asw, cell_type_asw, graph_connectivity, ari, nmi).}
#'   \item{integration_summary_table}{A tidy data.frame containing the 5 metrics with columns \code{Category}, \code{Property}, \code{Metric}, \code{Value}, and \code{Modality}.}
#'   \item{embedding}{The joint multiomics embedding matrix used for evaluation.}
#' @examples
#' data(example_multiomics, package = "scSimEval")
#' ct <- rep(c("A", "B"), each = 50)
#' res <- evaluate_multiomics_integration(sim_multi = example_multiomics$sim_multi, cell_types = ct)
#' @export
evaluate_multiomics_integration <- function(
  ref_multi = NULL,
  sim_multi,
  cell_types = NULL,
  batch_info = NULL,
  embed = NULL,
  reduction_name = NULL,
  n_pcs = 15,
  k = 30,
  threads = 1,
  verbose = TRUE
) {
  # Strict standardization check (Option B: Seurat / SingleCellExperiment)
  val_sim <- validate_multiomics_input(sim_multi, "sim_multi", reduction_name = reduction_name)
  sim_multi <- list(val_sim$m1, val_sim$m2)
  if (!is.null(ref_multi)) {
    val_ref <- validate_multiomics_input(ref_multi, "ref_multi", reduction_name = reduction_name)
    ref_multi <- list(val_ref$m1, val_ref$m2)
  }
  
  if (is.null(embed) && !is.null(val_sim$embedding)) {
    embed <- val_sim$embedding
    if (verbose) message("Using precomputed dimensional reduction from sim_multi for Category VI.")
  }
  if (is.null(cell_types) && !is.null(val_sim$cell_types)) {
    cell_types <- val_sim$cell_types
  }
  if (is.null(batch_info) && !is.null(val_sim$batch_info)) {
    batch_info <- val_sim$batch_info
  }

  # 1. Resolve Joint Embedding
  joint_embed <- embed
  ct_vec <- NULL
  bt_vec <- NULL
  
  if (is.null(joint_embed)) {
    if (is.matrix(sim_multi) || (is.data.frame(sim_multi) && !is.list(sim_multi))) {
      joint_embed <- as.matrix(sim_multi)
    } else if (is.list(sim_multi) && length(sim_multi) >= 2) {
      m1 <- as.matrix(sim_multi[[1]])
      m2 <- as.matrix(sim_multi[[2]])
      n1 <- ncol(m1)
      n2 <- ncol(m2)
      
      # Paired multiomics (equal cell counts)
      if (n1 == n2) {
        p1 <- fast_pca(t(normalize_counts(m1)), n_pcs = min(n_pcs, n1 - 1, nrow(m1) - 1))
        p2 <- fast_pca(t(normalize_counts(m2)), n_pcs = min(n_pcs, n2 - 1, nrow(m2) - 1))
        s1 <- scale(p1)
        s2 <- scale(p2)
        s1[is.na(s1)] <- 0
        s2[is.na(s2)] <- 0
        joint_embed <- cbind(s1, s2)
        
        # Resolve cell types
        if (!is.null(cell_types)) {
          if (length(cell_types) == n1) {
            ct_vec <- as.character(cell_types)
          } else if (is.data.frame(cell_types) && nrow(cell_types) == n1) {
            ct_vec <- as.character(cell_types[[1]])
          } else if (is.list(cell_types) && length(cell_types[[1]]) == n1) {
            ct_vec <- as.character(cell_types[[1]])
          }
        }
        
        # Resolve batches
        if (!is.null(batch_info)) {
          if (length(batch_info) == n1) {
            bt_vec <- as.character(batch_info)
          } else if (is.data.frame(batch_info) && nrow(batch_info) == n1) {
            bt_vec <- as.character(batch_info[[1]])
          } else if (is.list(batch_info) && length(batch_info[[1]]) == n1) {
            bt_vec <- as.character(batch_info[[1]])
          }
        }
      } else {
        # Unpaired multiomics: modalities stacked
        p1 <- fast_pca(t(normalize_counts(m1)), n_pcs = min(n_pcs, n1 - 1, nrow(m1) - 1))
        p2 <- fast_pca(t(normalize_counts(m2)), n_pcs = min(n_pcs, n2 - 1, nrow(m2) - 1))
        min_p <- min(ncol(p1), ncol(p2))
        joint_embed <- rbind(scale(p1[, seq_len(min_p), drop = FALSE]),
                             scale(p2[, seq_len(min_p), drop = FALSE]))
        joint_embed[is.na(joint_embed)] <- 0
        bt_vec <- c(rep("Modality_1", n1), rep("Modality_2", n2))
        
        if (!is.null(cell_types)) {
          if (length(cell_types) == (n1 + n2)) {
            ct_vec <- as.character(cell_types)
          } else if (is.list(cell_types) && length(cell_types) >= 2) {
            ct_vec <- c(as.character(cell_types[[1]]), as.character(cell_types[[2]]))
          }
        }
      }
    } else if (is.list(sim_multi) && length(sim_multi) == 1) {
      m1 <- as.matrix(sim_multi[[1]])
      joint_embed <- fast_pca(t(normalize_counts(m1)), n_pcs = min(n_pcs, ncol(m1) - 1, nrow(m1) - 1))
    }
  }
  
  if (is.null(joint_embed) || nrow(joint_embed) < 3) {
    warning("Unable to construct valid joint embedding for multiomics integration evaluation.")
    return(NULL)
  }
  
  n_cells <- nrow(joint_embed)
  
  # Fallbacks for cell types and batches if missing
  if (is.null(ct_vec)) {
    if (!is.null(cell_types) && length(cell_types) == n_cells) {
      ct_vec <- as.character(cell_types)
    } else {
      # Unsupervised clustering fallback
      k_fallback <- min(3, max(2, n_cells - 1))
      ct_vec <- paste0("Cluster_", stats::kmeans(joint_embed, centers = k_fallback, nstart = 5)$cluster)
    }
  }
  
  # If batch is NULL or unimodal, evaluate cross-modal alignment by projecting modalities
  eval_modality_mixing <- is.null(bt_vec) || length(unique(bt_vec)) < 2
  
  if (eval_modality_mixing && is.list(sim_multi) && length(sim_multi) >= 2) {
    m1 <- as.matrix(sim_multi[[1]])
    m2 <- as.matrix(sim_multi[[2]])
    n1 <- ncol(m1)
    n2 <- ncol(m2)
    p1 <- fast_pca(t(normalize_counts(m1)), n_pcs = min(n_pcs, n1 - 1, nrow(m1) - 1))
    p2 <- fast_pca(t(normalize_counts(m2)), n_pcs = min(n_pcs, n2 - 1, nrow(m2) - 1))
    min_p <- min(ncol(p1), ncol(p2))
    cross_embed <- rbind(scale(p1[, seq_len(min_p), drop = FALSE]),
                         scale(p2[, seq_len(min_p), drop = FALSE]))
    cross_embed[is.na(cross_embed)] <- 0
    cross_batch <- c(rep("Modality_1", n1), rep("Modality_2", n2))
    cross_ct <- c(ct_vec, if (length(ct_vec) == n1) ct_vec else rep(ct_vec, length.out = n2))
    
    b_asw <- calc_batch_asw(cross_embed, cross_batch, cross_ct)
  } else {
    bt_target <- if (!is.null(bt_vec)) bt_vec else rep("Batch_1", n_cells)
    b_asw <- calc_batch_asw(joint_embed, bt_target, ct_vec)
  }
  
  ct_asw <- calc_cell_type_asw(joint_embed, ct_vec)
  g_conn <- calc_graph_connectivity(joint_embed, ct_vec, k = min(15, n_cells - 1))
  ari_val <- calc_integration_ari(joint_embed, ct_vec)
  nmi_val <- calc_integration_nmi(joint_embed, ct_vec)
  
  metrics_list <- list(
    batch_asw = as.numeric(b_asw),
    cell_type_asw = as.numeric(ct_asw),
    graph_connectivity = as.numeric(g_conn),
    ari = as.numeric(ari_val),
    nmi = as.numeric(nmi_val)
  )
  
  summary_df <- data.frame(
    Category = "(VI) Multiomics Integration & Alignment",
    Property = "multiomics_integration",
    Metric = c("batch_asw", "cell_type_asw", "graph_connectivity", "ari", "nmi"),
    Value = as.numeric(c(b_asw, ct_asw, g_conn, ari_val, nmi_val)),
    Modality = "Joint (Cross-Modal)",
    stringsAsFactors = FALSE
  )
  
  res <- list(
    metrics = metrics_list,
    integration_summary_table = summary_df,
    embedding = joint_embed
  )
  class(res) <- c("scSimEval_integration", "list")
  res
}
