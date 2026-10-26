
#' Automatically Infer Pseudotime Trajectory from scRNA-seq Counts
#'
#' Estimates cell differentiation pseudotime directly from single-cell expression
#' counts using diffusion/principal curve projection along the first principal component
#' of top highly variable genes.
#'
#' @param data Expression or count matrix (genes x cells), or SingleCellExperiment object.
#' @param n_top Number of top variable genes to use (default 500).
#' @return Numeric vector of pseudotime values in [0, 1] for each cell.
#' @examples
#' data(example_scrna, package = "scSimEval")
#' infer_scrna_pseudotime(example_scrna$ref)
#' @export
infer_scrna_pseudotime <- function(data, n_top = 500) {
  if (inherits(data, "SingleCellExperiment") && requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    data <- SummarizedExperiment::assay(data, "counts")
  }
  n_cells <- ncol(data)
  if (is.null(n_cells) || n_cells < 3) return(seq(0, 1, length.out = max(1, as.integer(n_cells %||% length(data)))))
  
  data_norm <- normalize_counts(data)
  gene_vars <- fast_row_vars(data_norm)
  n_top_use <- min(n_top, length(gene_vars))
  top_idx <- order(gene_vars, decreasing = TRUE)[seq_len(n_top_use)]
  
  # Run PCA on top variable genes
  sub_mat <- as.matrix(data_norm[top_idx, , drop = FALSE])
  pca_pts <- fast_pca(t(sub_mat), n_pcs = min(2, ncol(sub_mat) - 1), n_hvg = n_top, scale = FALSE)
  pc1 <- pca_pts[, 1]
  
  # Normalize to [0, 1]
  rng <- max(pc1) - min(pc1)
  if (rng > 0) {
    pt <- (pc1 - min(pc1)) / rng
  } else {
    pt <- seq(0, 1, length.out = n_cells)
  }
  names(pt) <- colnames(data)
  return(pt)
}

#' Automatically Infer Lineage Tree (hclust) from scRNA-seq Counts & Cell Types
#'
#' Builds a hierarchical clustering lineage tree between cell-type centroids
#' from scRNA-seq expression data. If cell types are omitted or have fewer than 2
#' classes, cells are clustered automatically along the principal trajectory.
#'
#' @param data Expression or count matrix (genes x cells).
#' @param cell_types Factor or vector of cell-type annotations for each cell (optional).
#' @param method Linkage method for hierarchical clustering (default "ward.D2").
#' @return An \code{hclust} tree object representing the cell-type lineage tree.
#' @examples
#' data(example_scrna, package = "scSimEval")
#' infer_scrna_lineage_tree(example_scrna$ref, example_scrna$cell_types)
#' @export
infer_scrna_lineage_tree <- function(data, cell_types = NULL, method = "ward.D2") {
  if (inherits(data, "SingleCellExperiment") && requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    data <- SummarizedExperiment::assay(data, "counts")
  }
  mat <- as.matrix(normalize_counts(data))
  n_cells <- ncol(mat)
  
  if (is.data.frame(cell_types)) {
    cell_types <- if (!is.null(cell_types$CellType)) {
      cell_types$CellType
    } else if (ncol(cell_types) >= 2) {
      cell_types[[2]]
    } else {
      cell_types[[1]]
    }
  } else if (is.list(cell_types) && !is.factor(cell_types)) {
    cell_types <- cell_types[[1]]
  }
  
  # If cell types not provided, mismatched in length, or < 2 types, auto-cluster
  if (is.null(cell_types) || length(cell_types) != n_cells || length(unique(stats::na.omit(cell_types))) < 2) {
    k_use <- min(5, max(2, round(n_cells / 20)))
    pca_clu <- fast_pca(t(mat[seq_len(min(100, nrow(mat))), , drop = FALSE]), n_pcs = min(2, ncol(mat) - 1), scale = FALSE)
    set.seed(42)
    km <- stats::kmeans(pca_clu, centers = k_use, nstart = 3)
    cell_types <- as.factor(paste0("Cluster_", km$cluster))
  } else {
    cell_types <- as.factor(cell_types)
  }
  
  unique_cts <- levels(cell_types)
  if (length(unique_cts) < 2) {
    unique_cts <- unique(as.character(cell_types))
  }
  
  # Compute centroid mean expression per cell type
  centroids <- vapply(unique_cts, function(ct) {
    idx <- which(cell_types == ct)
    if (length(idx) > 1) {
      rowMeans(mat[, idx, drop = FALSE])
    } else if (length(idx) == 1) {
      mat[, idx]
    } else {
      rep(0, nrow(mat))
    }
  }, numeric(nrow(mat)))
  
  # Distance matrix between centroids
  dist_mat <- stats::dist(t(centroids))
  tree <- stats::hclust(dist_mat, method = method)
  return(tree)
}

#' Calculate Correlation of Geodesic Pseudotime Distances
#'
#' Evaluates whether the relative cell-to-cell ordering along a differentiation
#' trajectory is preserved between reference and simulated data.
#' Accepts either precomputed pseudotime vectors OR scRNA-seq expression matrices
#' (which will be automatically inferred via PCA trajectory).
#'
#' @param ref_pseudotime Numeric vector of reference pseudotimes OR reference count matrix (genes x cells).
#' @param sim_pseudotime Numeric vector of simulated pseudotimes OR simulated count matrix (genes x cells).
#' @param method Correlation method: "spearman" (default) or "pearson".
#'
#' @return Correlation coefficient between reference and simulated pseudotime trajectories.
#' @examples
#' data(example_scrna, package = "scSimEval")
#' calc_pseudotime_correlation(example_scrna$ref, example_scrna$sim)
#' @export
calc_pseudotime_correlation <- function(ref_pseudotime, sim_pseudotime, method = c("spearman", "pearson")) {
  method <- match.arg(method)
  
  # Auto-infer if 2D expression matrices provided
  if (!is.null(dim(ref_pseudotime)) && length(dim(ref_pseudotime)) == 2) {
    ref_pseudotime <- infer_scrna_pseudotime(ref_pseudotime)
  }
  if (!is.null(dim(sim_pseudotime)) && length(dim(sim_pseudotime)) == 2) {
    sim_pseudotime <- infer_scrna_pseudotime(sim_pseudotime)
  }
  
  aligned <- align_distributions(ref_pseudotime, sim_pseudotime)
  r <- stats::cor(aligned$ref, aligned$sim, method = method)
  if (is.na(r) || !is.finite(r)) 0 else r
}

#' Calculate Lineage Tree Branch Height Discrepancy (RMSE)
#'
#' Computes the Root Mean Squared Error (RMSE) between sorted branch heights of
#' hierarchical lineage trees. Accepts either precomputed \code{hclust} objects
#' OR scRNA-seq count matrices + \code{cell_types} (which will be automatically inferred).
#'
#' @param ref_tree Hierarchical clustering (\code{hclust}) object from reference data,
#'   OR reference count matrix.
#' @param sim_tree Hierarchical clustering (\code{hclust}) object from simulated data,
#'   OR simulated count matrix.
#' @param cell_types_ref Optional cell type labels for reference (if count matrix provided).
#' @param cell_types_sim Optional cell type labels for simulation (if count matrix provided).
#'
#' @return Root mean squared error between branch heights.
#' @examples
#' data(example_scrna, package = "scSimEval")
#' calc_tree_height_discrepancy(example_scrna$ref, example_scrna$sim,
#'                              cell_types_ref = example_scrna$cell_types,
#'                              cell_types_sim = example_scrna$cell_types)
#' @export
calc_tree_height_discrepancy <- function(
  ref_tree,
  sim_tree,
  cell_types_ref = NULL,
  cell_types_sim = NULL
) {
  # Auto-infer lineage tree if expression matrices provided
  if (!inherits(ref_tree, "hclust")) {
    ref_tree <- infer_scrna_lineage_tree(ref_tree, cell_types_ref)
  }
  if (!inherits(sim_tree, "hclust")) {
    sim_tree <- infer_scrna_lineage_tree(sim_tree, cell_types_sim)
  }
  
  h_ref <- sort(ref_tree$height)
  h_sim <- sort(sim_tree$height)
  
  if (length(h_ref) == 0 || length(h_sim) == 0) return(NA_real_)
  
  # Align tree branch lengths if different number of clusters
  if (length(h_ref) != length(h_sim)) {
    n_target <- max(length(h_ref), length(h_sim))
    h_ref <- stats::approx(seq_along(h_ref), h_ref, n = n_target)$y
    h_sim <- stats::approx(seq_along(h_sim), h_sim, n = n_target)$y
  }
  
  calc_rmse(h_ref, h_sim)
}

#' Calculate Trajectory Graph Alignment Score
#'
#' Evaluates trajectory topology preservation by measuring the concordance of
#' pairwise graph distances between landmark cells along the inferred trajectory.
#'
#' @param ref_data Reference expression/count matrix OR pseudotime vector.
#' @param sim_data Simulated expression/count matrix OR pseudotime vector.
#' @param ref_pseudotime Optional precomputed reference pseudotime.
#' @param sim_pseudotime Optional precomputed simulated pseudotime.
#' @param n_landmarks Number of trajectory landmark points (default 50).
#' @return Numeric trajectory graph alignment score in [0, 1] (higher is better).
#' @examples
#' data(example_scrna, package = "scSimEval")
#' calc_trajectory_graph_alignment(example_scrna$ref, example_scrna$sim)
#' @export
calc_trajectory_graph_alignment <- function(
  ref_data,
  sim_data,
  ref_pseudotime = NULL,
  sim_pseudotime = NULL,
  n_landmarks = 50
) {
  if (is.null(ref_pseudotime)) {
    ref_pseudotime <- if (!is.null(dim(ref_data)) && length(dim(ref_data)) == 2) {
      infer_scrna_pseudotime(ref_data)
    } else {
      as.numeric(ref_data)
    }
  }
  if (is.null(sim_pseudotime)) {
    sim_pseudotime <- if (!is.null(dim(sim_data)) && length(dim(sim_data)) == 2) {
      infer_scrna_pseudotime(sim_data)
    } else {
      as.numeric(sim_data)
    }
  }
  
  # Extract or compute PCA coordinates
  get_pca_coords <- function(mat) {
    if (!is.null(dim(mat)) && length(dim(mat)) == 2) {
      norm_m <- normalize_counts(mat)
      n_genes_use <- min(100, nrow(norm_m))
      g_vars <- fast_row_vars(norm_m)
      top_g <- order(g_vars, decreasing = TRUE)[seq_len(n_genes_use)]
      sub_m <- as.matrix(norm_m[top_g, , drop = FALSE])
      fast_pca(t(sub_m), n_pcs = min(3, ncol(sub_m) - 1), scale = FALSE)
    } else {
      matrix(as.numeric(mat), ncol = 1)
    }
  }
  
  ref_pca <- tryCatch(get_pca_coords(ref_data), error = function(e) matrix(as.numeric(ref_pseudotime), ncol = 1))
  sim_pca <- tryCatch(get_pca_coords(sim_data), error = function(e) matrix(as.numeric(sim_pseudotime), ncol = 1))
  
  n_pts <- min(n_landmarks, length(ref_pseudotime), length(sim_pseudotime))
  if (n_pts < 5) return(1.0)
  
  probs <- seq(0.02, 0.98, length.out = n_pts)
  ref_q <- stats::quantile(ref_pseudotime, probs = probs, na.rm = TRUE)
  sim_q <- stats::quantile(sim_pseudotime, probs = probs, na.rm = TRUE)
  
  ref_landmarks <- matrix(0, nrow = n_pts, ncol = ncol(ref_pca))
  sim_landmarks <- matrix(0, nrow = n_pts, ncol = ncol(sim_pca))
  
  window <- 0.15
  for (i in seq_len(n_pts)) {
    ref_idx <- which(abs(ref_pseudotime - ref_q[i]) <= window)
    if (length(ref_idx) == 0) ref_idx <- order(abs(ref_pseudotime - ref_q[i]))[seq_len(min(5, length(ref_pseudotime)))]
    ref_landmarks[i, ] <- colMeans(ref_pca[ref_idx, , drop = FALSE])
    
    sim_idx <- which(abs(sim_pseudotime - sim_q[i]) <= window)
    if (length(sim_idx) == 0) sim_idx <- order(abs(sim_pseudotime - sim_q[i]))[seq_len(min(5, length(sim_pseudotime)))]
    sim_landmarks[i, ] <- colMeans(sim_pca[sim_idx, , drop = FALSE])
  }
  
  d_ref <- as.matrix(stats::dist(ref_landmarks))
  d_sim <- as.matrix(stats::dist(sim_landmarks))
  
  v_ref <- d_ref[upper.tri(d_ref)]
  v_sim <- d_sim[upper.tri(d_sim)]
  
  r <- stats::cor(v_ref, v_sim, method = "spearman")
  if (is.na(r) || !is.finite(r)) r <- 0
  max(0, min(1, r))
}

#' Calculate Lineage Gene Expression Pattern Fidelity
#'
#' Evaluates the preservation of continuous gene expression trends along developmental
#' trajectories by calculating the mean correlation of expression profiles across pseudotime
#' for top dynamic lineage genes.
#'
#' @param ref_data Reference expression/count matrix.
#' @param sim_data Simulated expression/count matrix.
#' @param ref_pseudotime Optional precomputed reference pseudotime.
#' @param sim_pseudotime Optional precomputed simulated pseudotime.
#' @param n_genes Number of dynamic genes to evaluate (default 50).
#' @param n_bins Number of pseudotime trajectory bins (default 20).
#' @return Mean correlation of gene expression trends along lineages in [0, 1].
#' @examples
#' data(example_scrna, package = "scSimEval")
#' calc_lineage_expression_fidelity(example_scrna$ref, example_scrna$sim)
#' @export
calc_lineage_expression_fidelity <- function(
  ref_data,
  sim_data,
  ref_pseudotime = NULL,
  sim_pseudotime = NULL,
  n_genes = 50,
  n_bins = 20
) {
  if (is.null(ref_pseudotime)) {
    ref_pseudotime <- if (!is.null(dim(ref_data)) && length(dim(ref_data)) == 2) {
      infer_scrna_pseudotime(ref_data)
    } else {
      as.numeric(ref_data)
    }
  }
  if (is.null(sim_pseudotime)) {
    sim_pseudotime <- if (!is.null(dim(sim_data)) && length(dim(sim_data)) == 2) {
      infer_scrna_pseudotime(sim_data)
    } else {
      as.numeric(sim_data)
    }
  }
  
  if (is.null(dim(ref_data)) || is.null(dim(sim_data))) return(1.0)
  
  ref_norm <- normalize_counts(ref_data)
  sim_norm <- normalize_counts(sim_data)
  
  common_genes <- intersect(rownames(ref_norm), rownames(sim_norm))
  if (length(common_genes) < 5) {
    n_use <- min(nrow(ref_norm), nrow(sim_norm), n_genes)
    common_genes <- rownames(ref_norm)[seq_len(n_use)]
  }
  
  g_vars <- fast_row_vars(ref_norm[common_genes, , drop = FALSE])
  top_genes <- common_genes[order(g_vars, decreasing = TRUE)[seq_len(min(n_genes, length(g_vars)))]]
  
  ref_bins <- cut(ref_pseudotime, breaks = seq(0, 1, length.out = n_bins + 1), include.lowest = TRUE, labels = FALSE)
  sim_bins <- cut(sim_pseudotime, breaks = seq(0, 1, length.out = n_bins + 1), include.lowest = TRUE, labels = FALSE)
  
  ref_mat <- as.matrix(ref_norm[top_genes, , drop = FALSE])
  sim_mat <- as.matrix(sim_norm[top_genes, , drop = FALSE])
  
  gene_cors <- vapply(seq_along(top_genes), function(i) {
    ref_profile <- vapply(seq_len(n_bins), function(b) {
      idx <- which(ref_bins == b)
      if (length(idx) > 0) mean(ref_mat[i, idx]) else NA_real_
    }, numeric(1))
    
    sim_profile <- vapply(seq_len(n_bins), function(b) {
      idx <- which(sim_bins == b)
      if (length(idx) > 0) mean(sim_mat[i, idx]) else NA_real_
    }, numeric(1))
    
    valid <- !is.na(ref_profile) & !is.na(sim_profile)
    if (sum(valid) < 3) return(0)
    
    sd_r <- stats::sd(ref_profile[valid])
    sd_s <- stats::sd(sim_profile[valid])
    if (is.na(sd_r) || sd_r < 1e-8 || is.na(sd_s) || sd_s < 1e-8) return(0)
    
    r <- stats::cor(ref_profile[valid], sim_profile[valid], method = "spearman")
    if (is.na(r) || !is.finite(r)) 0 else max(0, min(1, r))
  }, numeric(1))
  
  val <- mean(gene_cors, na.rm = TRUE)
  if (is.na(val) || !is.finite(val)) val <- 0
  max(0, min(1, val))
}

#' Full Trajectory Accuracy Evaluation
#'
#' Evaluates trajectory preservation using pseudotime correlation, branch height
#' discrepancy (RMSE), trajectory graph alignment, and gene expression patterns along lineages.
#' If raw count matrices are provided, pseudotime, graphs, and lineage trees are automatically
#' inferred directly from the scRNA-seq data.
#'
#' @param ref_data Reference count matrix (genes x cells) OR numeric vector of pseudotimes.
#' @param sim_data Simulated count matrix (genes x cells) OR numeric vector of pseudotimes.
#' @param cell_types_ref Optional vector of reference cell types (for lineage tree inference).
#' @param cell_types_sim Optional vector of simulated cell types (for lineage tree inference).
#'
#' @return A named list of trajectory accuracy metrics:
#'   \describe{
#'     \item{pseudotime_correlation}{Spearman correlation between reference and simulated pseudotime.}
#'     \item{tree_height_rmse}{Root mean squared error of lineage tree branch heights.}
#'     \item{graph_alignment}{Concordance of pairwise graph distances along trajectory landmarks.}
#'     \item{lineage_expression_fidelity}{Mean correlation of gene expression trends along lineages.}
#'   }
#' @examples
#' data(example_scrna, package = "scSimEval")
#' evaluate_trajectory_metrics(example_scrna$ref, example_scrna$sim)
#' @export
evaluate_trajectory_metrics <- function(
  ref_data,
  sim_data,
  cell_types_ref = NULL,
  cell_types_sim = NULL
) {
  is_ref_mat <- !is.null(dim(ref_data)) && length(dim(ref_data)) == 2
  is_sim_mat <- !is.null(dim(sim_data)) && length(dim(sim_data)) == 2
  
  # Extract or infer pseudotime
  if (is_ref_mat) {
    ref_pseudotime <- infer_scrna_pseudotime(ref_data)
  } else {
    ref_pseudotime <- as.numeric(ref_data)
  }
  
  if (is_sim_mat) {
    sim_pseudotime <- infer_scrna_pseudotime(sim_data)
  } else {
    sim_pseudotime <- as.numeric(sim_data)
  }
  
  # 1. Pseudotime correlation
  cor_val <- calc_pseudotime_correlation(ref_pseudotime, sim_pseudotime)
  if (is.na(cor_val) || !is.finite(cor_val)) cor_val <- 0
  
  # 2. 1D distance metrics on pseudotime distributions
  dist_metrics <- calc_all_univariate_metrics(ref_pseudotime, sim_pseudotime, "pseudotime")
  
  # 3. Tree height discrepancy (RMSE)
  tree_rmse <- tryCatch(
    calc_tree_height_discrepancy(ref_data, sim_data, cell_types_ref, cell_types_sim),
    error = function(e) NA_real_
  )
  
  # 4. Trajectory graph alignment
  graph_align <- tryCatch(
    calc_trajectory_graph_alignment(ref_data, sim_data, ref_pseudotime, sim_pseudotime),
    error = function(e) NA_real_
  )
  
  # 5. Lineage gene expression pattern fidelity
  lineage_expr <- tryCatch(
    calc_lineage_expression_fidelity(ref_data, sim_data, ref_pseudotime, sim_pseudotime),
    error = function(e) NA_real_
  )
  
  c(list(
    pseudotime_correlation      = cor_val,
    tree_height_rmse            = tree_rmse,
    graph_alignment             = graph_align,
    lineage_expression_fidelity = lineage_expr
  ), dist_metrics)
}
