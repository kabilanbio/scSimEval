
#' Automatically Infer Pseudotime Trajectory from scRNA-seq Counts
#'
#' Estimates cell differentiation pseudotime directly from single-cell expression
#' counts using diffusion/principal curve projection along the first principal component
#' of top highly variable genes.
#'
#' @param data Expression or count matrix (genes x cells).
#' @param n_top Number of top variable genes to use (default 500).
#' @return Numeric vector of pseudotime values in [0, 1] for each cell.
#' @export
infer_scrna_pseudotime <- function(data, n_top = 500) {
  mat <- as.matrix(data)
  n_cells <- ncol(mat)
  if (n_cells < 3) return(seq(0, 1, length.out = n_cells))
  
  # Select top variable genes
  gene_vars <- apply(mat, 1, stats::var)
  n_top_use <- min(n_top, length(gene_vars))
  top_idx <- order(gene_vars, decreasing = TRUE)[seq_len(n_top_use)]
  
  # Run PCA
  sub_mat <- mat[top_idx, , drop = FALSE]
  pca_res <- stats::prcomp(t(sub_mat), center = TRUE, scale. = FALSE, rank. = 2)
  pc1 <- pca_res$x[, 1]
  
  # Normalize to [0, 1]
  rng <- max(pc1) - min(pc1)
  if (rng > 0) {
    pt <- (pc1 - min(pc1)) / rng
  } else {
    pt <- seq(0, 1, length.out = n_cells)
  }
  names(pt) <- colnames(mat)
  return(pt)
}

#' Automatically Infer Lineage Tree (hclust) from scRNA-seq Counts & Cell Types
#'
#' Builds a hierarchical clustering lineage tree between cell-type centroids
#' from scRNA-seq expression data.
#'
#' @param data Expression or count matrix (genes x cells).
#' @param cell_types Factor or vector of cell-type annotations for each cell.
#' @param method Linkage method for hierarchical clustering (default "ward.D2").
#' @return An \code{hclust} tree object representing the cell-type lineage tree.
#' @export
infer_scrna_lineage_tree <- function(data, cell_types, method = "ward.D2") {
  mat <- as.matrix(data)
  cell_types <- as.factor(cell_types)
  unique_cts <- levels(cell_types)
  
  if (length(unique_cts) < 2) {
    stop("Need at least 2 cell types to build a lineage tree.")
  }
  
  # Compute centroid mean expression per cell type
  centroids <- vapply(unique_cts, function(ct) {
    idx <- which(cell_types == ct)
    if (length(idx) > 1) {
      rowMeans(mat[, idx, drop = FALSE])
    } else {
      mat[, idx]
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
#' @export
calc_pseudotime_correlation <- function(ref_pseudotime, sim_pseudotime, method = c("spearman", "pearson")) {
  method <- match.arg(method)
  
  # Auto-infer if count matrices provided
  if (is.matrix(ref_pseudotime) || is.data.frame(ref_pseudotime)) {
    ref_pseudotime <- infer_scrna_pseudotime(ref_pseudotime)
  }
  if (is.matrix(sim_pseudotime) || is.data.frame(sim_pseudotime)) {
    sim_pseudotime <- infer_scrna_pseudotime(sim_pseudotime)
  }
  
  aligned <- align_distributions(ref_pseudotime, sim_pseudotime)
  stats::cor(aligned$ref, aligned$sim, method = method)
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
#' @export
calc_tree_height_discrepancy <- function(
  ref_tree,
  sim_tree,
  cell_types_ref = NULL,
  cell_types_sim = NULL
) {
  # Auto-infer lineage tree if expression matrices provided
  if (!inherits(ref_tree, "hclust")) {
    if (is.null(cell_types_ref)) stop("cell_types_ref must be provided when ref_tree is a matrix.")
    ref_tree <- infer_scrna_lineage_tree(ref_tree, cell_types_ref)
  }
  if (!inherits(sim_tree, "hclust")) {
    if (is.null(cell_types_sim)) stop("cell_types_sim must be provided when sim_tree is a matrix.")
    sim_tree <- infer_scrna_lineage_tree(sim_tree, cell_types_sim)
  }
  
  h_ref <- sort(ref_tree$height)
  h_sim <- sort(sim_tree$height)
  calc_rmse(h_ref, h_sim)
}

#' Full Trajectory Accuracy Evaluation
#'
#' Evaluates trajectory preservation using pseudotime correlation, branch height
#' discrepancy, and univariate distribution distance metrics on pseudotime.
#' If raw count matrices are provided, pseudotime and lineage trees are automatically
#' inferred directly from the scRNA-seq data.
#'
#' @param ref_data Reference count matrix (genes x cells) OR numeric vector of pseudotimes.
#' @param sim_data Simulated count matrix (genes x cells) OR numeric vector of pseudotimes.
#' @param cell_types_ref Optional vector of reference cell types (for lineage tree inference).
#' @param cell_types_sim Optional vector of simulated cell types (for lineage tree inference).
#'
#' @return A named list of trajectory accuracy metrics.
#' @export
evaluate_trajectory_metrics <- function(
  ref_data,
  sim_data,
  cell_types_ref = NULL,
  cell_types_sim = NULL
) {
  # Extract or use pseudotime
  if (is.matrix(ref_data) || is.data.frame(ref_data)) {
    ref_pseudotime <- infer_scrna_pseudotime(ref_data)
  } else {
    ref_pseudotime <- as.numeric(ref_data)
  }
  
  if (is.matrix(sim_data) || is.data.frame(sim_data)) {
    sim_pseudotime <- infer_scrna_pseudotime(sim_data)
  } else {
    sim_pseudotime <- as.numeric(sim_data)
  }
  
  # 1. Pseudotime correlation
  cor_val <- calc_pseudotime_correlation(ref_pseudotime, sim_pseudotime)
  
  # 2. 1D distance metrics on pseudotime distributions
  dist_metrics <- calc_all_univariate_metrics(ref_pseudotime, sim_pseudotime, "pseudotime")
  
  # 3. Tree height discrepancy if cell types are provided or if inputs were trees
  tree_rmse <- NA_real_
  if (!is.null(cell_types_ref) && !is.null(cell_types_sim) &&
      (is.matrix(ref_data) || is.data.frame(ref_data))) {
    tree_rmse <- tryCatch(
      calc_tree_height_discrepancy(ref_data, sim_data, cell_types_ref, cell_types_sim),
      error = function(e) NA_real_
    )
  }
  
  c(list(
    pseudotime_correlation = cor_val,
    tree_height_rmse = tree_rmse
  ), dist_metrics)
}
