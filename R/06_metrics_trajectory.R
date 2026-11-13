# ==============================================================================
# Saelens et al. (2019) / dynverse Category VII Workflow (scRNA-seq Only)
# ==============================================================================

#' Extract Sparse Count Matrix and Optional Cell Types for Trajectory Analysis
#' @keywords internal
#' @noRd
.traj_extract_input <- function(x, cell_types = NULL, arg_name = "data") {
  if (inherits(x, "SingleCellExperiment") || inherits(x, "Seurat")) {
    v <- validate_sc_object(x, arg_name = arg_name)
    mat <- v$counts
    if (is.null(cell_types) && !is.null(v$cell_types)) {
      cell_types <- v$cell_types
    }
  } else if (is.list(x) && !is.null(x$counts)) {
    mat <- if (inherits(x$counts, "dgCMatrix")) x$counts else methods::as(as.matrix(x$counts), "dgCMatrix")
    if (is.null(cell_types) && !is.null(x$cell_types)) cell_types <- x$cell_types
  } else if (!is.null(dim(x)) && length(dim(x)) == 2) {
    mat <- if (inherits(x, "dgCMatrix")) x else methods::as(as.matrix(x), "dgCMatrix")
  } else {
    stop(sprintf("Trajectory input '%s' must be a Seurat, SingleCellExperiment, or 2D count matrix.", arg_name),
         call. = FALSE)
  }

  if (is.null(rownames(mat)) || any(!nzchar(rownames(mat)))) {
    rownames(mat) <- paste0("Gene_", seq_len(nrow(mat)))
  } else {
    rownames(mat) <- make.unique(as.character(rownames(mat)))
  }
  if (is.null(colnames(mat)) || any(!nzchar(colnames(mat)))) {
    colnames(mat) <- paste0("Cell_", seq_len(ncol(mat)))
  } else {
    colnames(mat) <- make.unique(as.character(colnames(mat)))
  }

  if (!is.null(cell_types)) {
    cell_types <- as.character(cell_types)
    if (length(cell_types) != ncol(mat)) {
      cell_types <- rep(cell_types, length.out = ncol(mat))
    }
  }
  list(counts = mat, cell_types = cell_types)
}

#' Filter Outlier Cells via 3-MAD Quality Control (Saelens et al., 2019)
#' @keywords internal
#' @noRd
.traj_filter_outlier_cells <- function(counts, cell_types = NULL, nmads = 3, min_cells = 15L) {
  n_cells <- ncol(counts)
  if (n_cells <= min_cells) {
    keep <- Matrix::colSums(counts) > 0
    if (sum(keep) < 3L) keep <- rep(TRUE, n_cells)
    return(list(
      counts = counts[, keep, drop = FALSE],
      cell_types = if (!is.null(cell_types)) cell_types[keep] else NULL,
      keep = keep
    ))
  }

  lib_size <- as.numeric(Matrix::colSums(counts))
  n_detected <- as.numeric(Matrix::colSums(counts > 0))
  mito_idx <- grep("^MT-|^mt-", rownames(counts))
  pct_mito <- if (length(mito_idx) > 0) {
    mito_sum <- as.numeric(Matrix::colSums(counts[mito_idx, , drop = FALSE]))
    ifelse(lib_size > 0, 100 * mito_sum / lib_size, 0)
  } else {
    rep(0, n_cells)
  }

  mad_outlier <- function(val, type = "lower", use_log = TRUE) {
    if (requireNamespace("scater", quietly = TRUE)) {
      res <- tryCatch(
        suppressWarnings(as.logical(scater::isOutlier(val, nmads = nmads, type = type, log = use_log))),
        error = function(e) NULL
      )
      if (!is.null(res) && length(res) == length(val)) {
        res[is.na(res)] <- FALSE
        return(res)
      }
    }
    x <- if (use_log) log1p(pmax(val, 0)) else val
    med <- stats::median(x, na.rm = TRUE)
    m <- stats::mad(x, center = med, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(m) || m <= 0) return(rep(FALSE, length(val)))
    if (type == "lower") x < (med - nmads * m) else x > (med + nmads * m)
  }

  out_lib <- mad_outlier(lib_size, type = "lower", use_log = TRUE)
  out_det <- mad_outlier(n_detected, type = "lower", use_log = TRUE)
  out_mito <- if (length(mito_idx) > 0 && any(pct_mito > 0)) {
    mad_outlier(pct_mito, type = "higher", use_log = FALSE)
  } else {
    rep(FALSE, n_cells)
  }

  keep <- !(out_lib | out_det | out_mito) & (lib_size > 0) & (n_detected > 0)
  if (sum(keep) < max(min_cells, floor(0.5 * n_cells))) {
    keep <- (lib_size > 0) & (n_detected > 0)
  }
  if (sum(keep) < 3L) {
    keep <- rep(TRUE, n_cells)
  }

  list(
    counts = counts[, keep, drop = FALSE],
    cell_types = if (!is.null(cell_types)) cell_types[keep] else NULL,
    keep = keep
  )
}

#' Library-Size Log-Normalization for Trajectory Analysis
#' @keywords internal
#' @noRd
.traj_lognorm <- function(counts) {
  lib <- as.numeric(Matrix::colSums(counts))
  pos_lib <- lib[lib > 0]
  med_lib <- if (length(pos_lib) > 0) stats::median(pos_lib) else 1
  if (!is.finite(med_lib) || med_lib <= 0) med_lib <- 1
  sf <- ifelse(lib > 0, lib / med_lib, 1)
  scaled <- Matrix::t(Matrix::t(counts) / sf)
  log1p(scaled)
}

#' Empirical Mean-Variance Decomposition Fallback for HVG Selection
#' @keywords internal
#' @noRd
.traj_empirical_hvg_stats <- function(log_mat) {
  mu <- Matrix::rowMeans(log_mat)
  v <- as.numeric(fast_row_vars(log_mat))
  v[!is.finite(v) | v < 0] <- 0

  ok <- is.finite(mu) & is.finite(v) & (mu > 0)
  tech <- rep(0, nrow(log_mat))
  if (sum(ok) >= 10) {
    fit <- tryCatch(
      suppressWarnings(stats::loess(v[ok] ~ mu[ok], span = 0.75, degree = 1)),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      pred <- tryCatch(stats::predict(fit, newdata = mu), error = function(e) NULL)
      if (!is.null(pred)) {
        pred[is.na(pred)] <- stats::median(v[ok], na.rm = TRUE)
        tech <- pmax(pred, 0)
      }
    }
  }
  if (all(tech == 0)) {
    tech <- stats::median(v, na.rm = TRUE)
  }
  bio <- v - tech
  pos_sd <- stats::mad(bio, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(pos_sd) || pos_sd <= 1e-8) pos_sd <- stats::sd(bio, na.rm = TRUE)
  if (!is.finite(pos_sd) || pos_sd <= 1e-8) pos_sd <- 1
  pvals <- stats::pnorm(bio / pos_sd, lower.tail = FALSE)
  fdr <- stats::p.adjust(pvals, method = "BH")
  data.frame(
    gene = rownames(log_mat),
    mean = as.numeric(mu),
    total = v,
    tech = tech,
    bio = bio,
    FDR = fdr,
    stringsAsFactors = FALSE
  )
}

#' Preprocess and Select Highly Variable Genes for Trajectory Evaluation (Saelens et al., 2019)
#'
#' Implements the standardized scRNA-seq quality control and HVG selection workflow
#' from the \emph{dynverse} benchmarking study (\emph{Saelens et al., Nature Biotechnology, 2019}):
#' \enumerate{
#'   \item Harsh cell outlier filtering (3 MADs from median on library size, detected genes, and mitochondrial percentage).
#'   \item Low-expression gene filtering (expressed in \eqn{\ge 5\%} of cells and average expression \eqn{> 0.02}).
#'   \item HVG selection via mean-variance decomposition (\code{scran::modelGeneVar} with \code{FDR < 0.05}
#'         and biological variance \code{> 0.5}, with adaptive top-HVG fallback for sparse/small benchmark sets).
#'   \item Iterative removal of all-zero cells and uninformative zero features until convergence.
#' }
#'
#' @param ref_data Reference count matrix, \code{SingleCellExperiment}, or \code{Seurat} object.
#' @param sim_data Optional simulated count matrix, \code{SingleCellExperiment}, or \code{Seurat} object.
#' @param cell_types_ref Optional vector of reference cell-type labels.
#' @param cell_types_sim Optional vector of simulated cell-type labels.
#' @param nmads MAD threshold for outlier cell removal (default \code{3}).
#' @param min_cell_frac Minimum fraction of cells expressing a retained gene (default \code{0.05}).
#' @param min_mean_expr Minimum average expression count for a retained gene (default \code{0.02}).
#' @param fdr_threshold FDR cutoff for HVG biological variance (default \code{0.05}).
#' @param bio_var_threshold Minimum biological variance component for strict HVG selection (default \code{0.5}).
#' @param min_hvgs Minimum HVG floor with adaptive fallback if fewer pass strict thresholds (default \code{30L}).
#' @param max_hvgs Maximum HVG cap for tractable Random Forest importance fitting (default \code{1500L}).
#'
#' @return A list containing HVG-filtered log-normalized expression matrices (\code{ref_expr}, \code{sim_expr}),
#'   selected HVG names (\code{hvgs}), synced cell-type vectors (\code{cell_types_ref}, \code{cell_types_sim}),
#'   HVG statistics table (\code{hvg_stats}), and filtering counts (\code{summary}).
#' @examples
#' data(example_scrna, package = "scSimEval")
#' prep <- preprocess_trajectory_hvgs(example_scrna$ref, example_scrna$sim,
#'   cell_types_ref = example_scrna$cell_types, cell_types_sim = example_scrna$cell_types)
#' length(prep$hvgs)
#' @export
preprocess_trajectory_hvgs <- function(
  ref_data,
  sim_data = NULL,
  cell_types_ref = NULL,
  cell_types_sim = NULL,
  nmads = 3,
  min_cell_frac = 0.05,
  min_mean_expr = 0.02,
  fdr_threshold = 0.05,
  bio_var_threshold = 0.5,
  min_hvgs = 30L,
  max_hvgs = 1500L
) {
  ref_in <- .traj_extract_input(ref_data, cell_types_ref, arg_name = "ref_data")
  ref_counts <- ref_in$counts
  ct_ref <- ref_in$cell_types

  sim_counts <- NULL
  ct_sim <- NULL
  if (!is.null(sim_data)) {
    sim_in <- .traj_extract_input(sim_data, cell_types_sim, arg_name = "sim_data")
    sim_counts <- sim_in$counts
    ct_sim <- sim_in$cell_types

    common_genes <- intersect(rownames(ref_counts), rownames(sim_counts))
    if (length(common_genes) >= 10L) {
      ref_counts <- ref_counts[common_genes, , drop = FALSE]
      sim_counts <- sim_counts[common_genes, , drop = FALSE]
    } else if (nrow(ref_counts) == nrow(sim_counts)) {
      rownames(sim_counts) <- rownames(ref_counts)
    } else {
      stop("Reference and simulated datasets share fewer than 10 genes for trajectory evaluation.", call. = FALSE)
    }
  }

  # 1. Harsh cell QC (3 MADs)
  ref_qc <- .traj_filter_outlier_cells(ref_counts, ct_ref, nmads = nmads)
  ref_counts <- ref_qc$counts
  ct_ref <- ref_qc$cell_types

  if (!is.null(sim_counts)) {
    sim_qc <- .traj_filter_outlier_cells(sim_counts, ct_sim, nmads = nmads)
    sim_counts <- sim_qc$counts
    ct_sim <- sim_qc$cell_types
  }

  # 2. Low-expression gene filtering (>= 5% cells, mean > 0.02)
  det_frac <- as.numeric(Matrix::rowMeans(ref_counts > 0))
  mean_cnt <- as.numeric(Matrix::rowMeans(ref_counts))
  gene_keep <- (det_frac >= min_cell_frac) & (mean_cnt > min_mean_expr)
  target_floor <- min(as.integer(min_hvgs), nrow(ref_counts))
  if (sum(gene_keep) < target_floor) {
    gene_keep <- det_frac > 0
  }
  if (sum(gene_keep) < 2L) {
    gene_keep <- rep(TRUE, nrow(ref_counts))
  }
  ref_counts <- ref_counts[gene_keep, , drop = FALSE]
  if (!is.null(sim_counts)) {
    sim_counts <- sim_counts[gene_keep, , drop = FALSE]
  }

  # 3. HVG selection on reference log-normalized expression (scran + fallback)
  ref_log <- .traj_lognorm(ref_counts)
  hvg_stats <- NULL
  hvg_method <- "empirical_mean_var"
  if (requireNamespace("scran", quietly = TRUE) && nrow(ref_log) >= 20L && ncol(ref_log) >= 10L) {
    hvg_df <- tryCatch({
      suppressWarnings(scran::modelGeneVar(ref_log))
    }, error = function(e) NULL)
    if (!is.null(hvg_df) && all(c("bio", "FDR") %in% colnames(hvg_df))) {
      hvg_stats <- as.data.frame(hvg_df)
      hvg_stats$gene <- rownames(ref_log)
      hvg_method <- "scran_modelGeneVar"
    }
  }
  if (is.null(hvg_stats)) {
    hvg_stats <- .traj_empirical_hvg_stats(ref_log)
  }

  strict_idx <- which(is.finite(hvg_stats$FDR) & is.finite(hvg_stats$bio) &
                        hvg_stats$FDR < fdr_threshold & hvg_stats$bio > bio_var_threshold)
  sel_idx <- strict_idx
  if (length(sel_idx) < target_floor) {
    relaxed_idx <- which(is.finite(hvg_stats$FDR) & is.finite(hvg_stats$bio) &
                           hvg_stats$FDR < fdr_threshold & hvg_stats$bio > 0)
    sel_idx <- if (length(relaxed_idx) >= target_floor) {
      relaxed_idx
    } else {
      ord <- order(hvg_stats$bio, hvg_stats$total, decreasing = TRUE, na.last = TRUE)
      utils::head(ord, min(max(target_floor, min(300L, nrow(hvg_stats))), nrow(hvg_stats)))
    }
  }
  if (length(sel_idx) > max_hvgs) {
    sub_ord <- sel_idx[order(hvg_stats$bio[sel_idx], decreasing = TRUE, na.last = TRUE)]
    sel_idx <- utils::head(sub_ord, max_hvgs)
  }
  hvgs <- hvg_stats$gene[sel_idx]

  # 4. Subset to HVGs & iterative zero convergence filtering
  ref_sub <- ref_counts[hvgs, , drop = FALSE]
  sim_sub <- if (!is.null(sim_counts)) sim_counts[hvgs, , drop = FALSE] else NULL

  repeat {
    n_g0 <- nrow(ref_sub)
    n_c_ref0 <- ncol(ref_sub)
    n_c_sim0 <- if (!is.null(sim_sub)) ncol(sim_sub) else 0L

    # Remove all-zero reference genes (genes silent in sim remain so sim receives 0 RF importance)
    g_nz <- Matrix::rowSums(ref_sub) > 0
    if (sum(g_nz) >= 2L && any(!g_nz)) {
      ref_sub <- ref_sub[g_nz, , drop = FALSE]
      if (!is.null(sim_sub)) sim_sub <- sim_sub[g_nz, , drop = FALSE]
    }

    c_ref_nz <- Matrix::colSums(ref_sub) > 0
    if (sum(c_ref_nz) >= 3L && any(!c_ref_nz)) {
      ref_sub <- ref_sub[, c_ref_nz, drop = FALSE]
      if (!is.null(ct_ref)) ct_ref <- ct_ref[c_ref_nz]
    }

    if (!is.null(sim_sub)) {
      c_sim_nz <- Matrix::colSums(sim_sub) > 0
      if (sum(c_sim_nz) >= 3L && any(!c_sim_nz)) {
        sim_sub <- sim_sub[, c_sim_nz, drop = FALSE]
        if (!is.null(ct_sim)) ct_sim <- ct_sim[c_sim_nz]
      }
    }

    changed <- (nrow(ref_sub) != n_g0) ||
      (ncol(ref_sub) != n_c_ref0) ||
      (!is.null(sim_sub) && ncol(sim_sub) != n_c_sim0)
    if (!changed) break
  }

  hvgs_final <- rownames(ref_sub)
  ref_expr <- .traj_lognorm(ref_sub)
  sim_expr <- if (!is.null(sim_sub)) .traj_lognorm(sim_sub) else NULL

  list(
    ref_expr = ref_expr,
    sim_expr = sim_expr,
    hvgs = hvgs_final,
    cell_types_ref = ct_ref,
    cell_types_sim = ct_sim,
    hvg_stats = hvg_stats,
    summary = list(
      n_hvgs = length(hvgs_final),
      n_strict_hvgs = length(strict_idx),
      n_ref_cells = ncol(ref_expr),
      n_sim_cells = if (!is.null(sim_expr)) ncol(sim_expr) else NA_integer_,
      hvg_method = hvg_method
    )
  )
}

#' Project Cells onto Continuous Minimum Spanning Tree Branches
#' @keywords internal
#' @noRd
.traj_build_mst_and_progressions <- function(pca_coords, cell_ids, groups, max_milestones = 25L) {
  groups <- as.character(groups)
  tab <- table(groups)
  if (length(tab) > max_milestones) {
    keep_top <- names(sort(tab, decreasing = TRUE))[seq_len(max_milestones)]
    other <- !(groups %in% keep_top)
    if (any(other)) groups[other] <- keep_top[1]
  }
  # Merge singleton or tiny (< 2 cell) groups into nearest centroid group
  tab <- table(groups)
  u_groups <- names(tab)
  if (length(u_groups) > 2L && any(tab < 2L)) {
    big_groups <- names(tab[tab >= 2L])
    if (length(big_groups) >= 2L) {
      big_cents <- t(vapply(big_groups, function(g) colMeans(pca_coords[groups == g, , drop = FALSE]),
                            numeric(ncol(pca_coords))))
      for (sg in names(tab[tab < 2L])) {
        idx <- which(groups == sg)
        for (i in idx) {
          dists <- rowSums(sweep(big_cents, 2L, pca_coords[i, ], `-`)^2)
          groups[i] <- big_groups[which.min(dists)]
        }
      }
    }
  }

  milestone_ids <- sort(unique(groups))
  M <- length(milestone_ids)
  if (M < 2L) {
    # Fallback 2-state split along PC1 so every trajectory has a valid edge
    ord_pc1 <- order(pca_coords[, 1])
    groups <- ifelse(seq_len(nrow(pca_coords)) %in% ord_pc1[seq_len(floor(nrow(pca_coords) / 2))],
                     "M1", "M2")
    milestone_ids <- c("M1", "M2")
    M <- 2L
  }

  centroids <- t(vapply(milestone_ids, function(m) {
    colMeans(pca_coords[groups == m, , drop = FALSE])
  }, numeric(ncol(pca_coords))))
  rownames(centroids) <- milestone_ids

  cdist <- as.matrix(stats::dist(centroids, method = "euclidean"))
  cdist[cdist < 1e-6] <- 1e-6
  diag(cdist) <- 0

  g_full <- igraph::graph_from_adjacency_matrix(cdist, mode = "undirected", weighted = TRUE)
  mst_g <- igraph::mst(g_full, weights = igraph::E(g_full)$weight)
  edf <- igraph::as_data_frame(mst_g, what = "edges")

  milestone_network <- data.frame(
    from = as.character(edf$from),
    to = as.character(edf$to),
    length = pmax(as.numeric(edf$weight), 1e-4),
    directed = FALSE,
    stringsAsFactors = FALSE
  )

  # Map each cell continuously onto its best incident MST segment [from, to]
  n_cells <- nrow(pca_coords)
  prog_from <- character(n_cells)
  prog_to <- character(n_cells)
  prog_pct <- numeric(n_cells)

  for (i in seq_len(n_cells)) {
    g_i <- groups[i]
    inc_idx <- which(milestone_network$from == g_i | milestone_network$to == g_i)
    if (length(inc_idx) == 0L) inc_idx <- seq_len(nrow(milestone_network))

    x_i <- pca_coords[i, ]
    best_dist2 <- Inf
    best_from <- milestone_network$from[inc_idx[1]]
    best_to <- milestone_network$to[inc_idx[1]]
    best_t <- 0.5

    for (e_id in inc_idx) {
      u <- milestone_network$from[e_id]
      v <- milestone_network$to[e_id]
      cu <- centroids[u, ]
      cv <- centroids[v, ]
      vec <- cv - cu
      denom <- sum(vec^2)
      t_val <- if (denom > 1e-12) sum((x_i - cu) * vec) / denom else 0.5
      t_val <- max(0, min(1, t_val))
      proj <- cu + t_val * vec
      d2 <- sum((x_i - proj)^2)
      if (d2 < best_dist2) {
        best_dist2 <- d2
        best_from <- u
        best_to <- v
        best_t <- t_val
      }
    }
    prog_from[i] <- best_from
    prog_to[i] <- best_to
    prog_pct[i] <- best_t
  }

  progressions <- data.frame(
    cell_id = cell_ids,
    from = prog_from,
    to = prog_to,
    percentage = prog_pct,
    stringsAsFactors = FALSE
  )

  # Compute exact cell x milestone geodesic distances along the MST
  m_geodesic <- igraph::distances(mst_g, weights = igraph::E(mst_g)$weight)
  m_geodesic <- m_geodesic[milestone_ids, milestone_ids, drop = FALSE]
  edge_len_lookup <- stats::setNames(milestone_network$length,
                                     paste0(milestone_network$from, "||", milestone_network$to))
  cell_edge_len <- edge_len_lookup[paste0(prog_from, "||", prog_to)]

  d_to_from <- outer(prog_pct * cell_edge_len, rep(1, M))
  d_to_to <- outer((1 - prog_pct) * cell_edge_len, rep(1, M))
  cell_milestone_geodesic <- pmin(
    d_to_from + m_geodesic[prog_from, , drop = FALSE],
    d_to_to + m_geodesic[prog_to, , drop = FALSE]
  )
  rownames(cell_milestone_geodesic) <- cell_ids
  colnames(cell_milestone_geodesic) <- milestone_ids

  list(
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    progressions = progressions,
    centroids = centroids,
    cell_milestone_geodesic = cell_milestone_geodesic
  )
}

#' Infer Docker-Free dynwrap Milestone Trajectory from Normalized Expression
#'
#' Constructs a continuous Minimum Spanning Tree (MST) milestone trajectory directly from
#' log-normalized expression counts (typically HVG-filtered via \code{\link{preprocess_trajectory_hvgs}}),
#' projects cells continuously onto MST edges, and computes pairwise cell-to-milestone geodesic
#' distances along the milestone network.
#'
#' @param expr_mat Log-normalized expression matrix (\code{genes x cells}).
#' @param cell_types Optional character vector or factor of cell-type identities.
#' @param n_pcs Number of PCA components for MST construction (default \code{15L}).
#' @param pca_coords Optional precomputed cell PCA matrix (\code{cells x dims}).
#' @param milestone_prefix Prefix for unsupervised cluster milestones when \code{cell_types} is \code{NULL}.
#'
#' @return A list (augmented \code{dynwrap} trajectory object when \pkg{dynwrap} is available)
#'   containing \code{milestone_ids}, \code{milestone_network}, \code{progressions},
#'   \code{milestone_percentages}, \code{pca_coords}, and \code{cell_milestone_geodesic}.
#' @examples
#' data(example_scrna, package = "scSimEval")
#' prep <- preprocess_trajectory_hvgs(example_scrna$ref, cell_types_ref = example_scrna$cell_types)
#' traj <- infer_dynverse_trajectory(prep$ref_expr, prep$cell_types_ref)
#' @export
infer_dynverse_trajectory <- function(
  expr_mat,
  cell_types = NULL,
  n_pcs = 15L,
  pca_coords = NULL,
  milestone_prefix = "M"
) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for MST trajectory inference.", call. = FALSE)
  }

  cell_ids <- colnames(expr_mat)
  if (is.null(cell_ids)) {
    cell_ids <- paste0("Cell_", seq_len(ncol(expr_mat)))
    colnames(expr_mat) <- cell_ids
  }

  if (is.null(pca_coords)) {
    k <- min(as.integer(n_pcs), nrow(expr_mat) - 1L, ncol(expr_mat) - 1L)
    k <- max(2L, k)
    pca_coords <- suppressWarnings(fast_pca(Matrix::t(expr_mat), n_pcs = k))
    rownames(pca_coords) <- cell_ids
  } else {
    pca_coords <- as.matrix(pca_coords)
    rownames(pca_coords) <- cell_ids
  }

  groups <- if (!is.null(cell_types) && length(unique(stats::na.omit(cell_types))) >= 2L) {
    ct <- as.character(cell_types)
    ct[is.na(ct) | !nzchar(ct)] <- "Unknown"
    ct
  } else {
    k_cl <- min(max(3L, floor(sqrt(nrow(pca_coords) / 6))), 8L, nrow(pca_coords) - 1L)
    ord_pc1 <- order(pca_coords[, 1])
    init_idx <- ord_pc1[unique(round(seq(1, nrow(pca_coords), length.out = k_cl)))]
    init_centers <- pca_coords[init_idx, , drop = FALSE]
    if (nrow(init_centers) < 2L || any(duplicated(init_centers))) init_centers <- k_cl
    km <- stats::kmeans(pca_coords, centers = init_centers, iter.max = 50)
    paste0(milestone_prefix, km$cluster)
  }

  mst_res <- .traj_build_mst_and_progressions(pca_coords, cell_ids, groups)

  # Wrap cleanly via dynwrap when available (Docker-free)
  traj_obj <- NULL
  if (requireNamespace("dynwrap", quietly = TRUE)) {
    traj_obj <- tryCatch({
      ds <- dynwrap::wrap_data(cell_ids = cell_ids)
      dynwrap::add_trajectory(
        dataset = ds,
        milestone_ids = mst_res$milestone_ids,
        milestone_network = mst_res$milestone_network,
        progressions = mst_res$progressions
      )
    }, error = function(e) NULL)
  }

  if (is.null(traj_obj)) {
    traj_obj <- list(
      cell_ids = cell_ids,
      milestone_ids = mst_res$milestone_ids,
      milestone_network = mst_res$milestone_network,
      progressions = mst_res$progressions
    )
  }

  traj_obj$pca_coords <- pca_coords
  traj_obj$centroids <- mst_res$centroids
  traj_obj$cell_milestone_geodesic <- mst_res$cell_milestone_geodesic
  traj_obj$branch_assignment <- paste(
    pmin(mst_res$progressions$from, mst_res$progressions$to),
    pmax(mst_res$progressions$from, mst_res$progressions$to),
    sep = "->"
  )
  traj_obj
}

#' Align Unsupervised Trajectory Milestones Between Reference and Simulated Graphs
#' @keywords internal
#' @noRd
.traj_align_unsupervised_milestones <- function(ref_traj, sim_traj) {
  ref_c <- ref_traj$centroids
  sim_c <- sim_traj$centroids
  k_match <- min(ncol(ref_c), ncol(sim_c))
  d_mat <- fast_cross_distance_sq(ref_c[, seq_len(k_match), drop = FALSE],
                                  sim_c[, seq_len(k_match), drop = FALSE])
  n_r <- nrow(ref_c)
  n_s <- nrow(sim_c)

  mapping <- stats::setNames(rownames(sim_c), rownames(sim_c))
  if (requireNamespace("clue", quietly = TRUE) && n_r >= 2L && n_s >= 2L) {
    if (n_s <= n_r) {
      perm <- as.integer(clue::solve_LSAP(t(d_mat), maximum = FALSE))
      mapping <- stats::setNames(rownames(ref_c)[perm], rownames(sim_c))
    } else {
      perm <- as.integer(clue::solve_LSAP(d_mat, maximum = FALSE))
      used_r <- stats::setNames(rownames(ref_c), rownames(sim_c)[perm])
      for (nm in names(used_r)) mapping[[nm]] <- used_r[[nm]]
    }
  } else {
    used_ref <- character(0)
    for (s_idx in seq_len(n_s)) {
      ord <- order(d_mat[, s_idx])
      pick <- rownames(ref_c)[ord[!(rownames(ref_c)[ord] %in% used_ref)][1]]
      if (is.na(pick)) pick <- rownames(ref_c)[ord[1]]
      used_ref <- c(used_ref, pick)
      mapping[rownames(sim_c)[s_idx]] <- pick
    }
  }

  # Rename simulated milestone IDs in network, progressions, and geodesic columns
  sim_traj$milestone_ids <- unname(mapping[sim_traj$milestone_ids])
  sim_traj$milestone_network$from <- unname(mapping[sim_traj$milestone_network$from])
  sim_traj$milestone_network$to <- unname(mapping[sim_traj$milestone_network$to])
  sim_traj$progressions$from <- unname(mapping[sim_traj$progressions$from])
  sim_traj$progressions$to <- unname(mapping[sim_traj$progressions$to])
  rownames(sim_traj$centroids) <- unname(mapping[rownames(sim_traj$centroids)])
  colnames(sim_traj$cell_milestone_geodesic) <- unname(mapping[colnames(sim_traj$cell_milestone_geodesic)])
  sim_traj$branch_assignment <- paste(
    pmin(sim_traj$progressions$from, sim_traj$progressions$to),
    pmax(sim_traj$progressions$from, sim_traj$progressions$to),
    sep = "->"
  )
  sim_traj
}

#' Pure-R Hamming-Ipsen-Mikhailov (HIM) Topology Similarity (Jurman et al., 2015)
#' @keywords internal
#' @noRd
.traj_him_pure_r <- function(net1, net2, ga = 0.1) {
  nodes <- sort(unique(c(net1$from, net1$to, net2$from, net2$to)))
  n <- length(nodes)
  if (n < 2L) return(1)

  build_adj <- function(net) {
    A <- matrix(0, nrow = n, ncol = n, dimnames = list(nodes, nodes))
    for (r in seq_len(nrow(net))) {
      u <- as.character(net$from[r])
      v <- as.character(net$to[r])
      w <- if ("length" %in% colnames(net)) max(as.numeric(net$length[r]), 1e-6) else 1
      if (u != v) {
        A[u, v] <- A[u, v] + w
        A[v, u] <- A[v, u] + w
      }
    }
    A
  }

  A1 <- build_adj(net1)
  A2 <- build_adj(net2)
  if (max(A1) == 0 || max(A2) == 0) return(0)
  A1 <- A1 / sum(A1)
  A2 <- A2 / sum(A2)

  # Normalized Hamming distance over normalized weighted edges
  H <- 0.5 * sum(abs(A1 - A2))
  H <- max(0, min(1, H))

  # Ipsen-Mikhailov spectral distance on Laplacian frequencies
  lap_freqs <- function(A) {
    L <- diag(rowSums(A)) - A
    ev <- sort(eigen(L, symmetric = TRUE, only.values = TRUE)$values)
    sqrt(pmax(ev[-1L], 0))
  }

  im_norm_const <- function(N) {
    # Max spectral L2 distance: empty graph (all omega = 0) vs K_N (equal eigenvalues sum=1/(N-1) or scaled)
    K_N <- matrix(1 / (N * (N - 1)), nrow = N, ncol = N)
    diag(K_N) <- 0
    wK <- lap_freqs(K_N)
    wE <- rep(0, N - 1L)
    spec_l2(wK, wE)
  }

  lorentz_K <- function(w) {
    1 / sum(pi / 2 + atan(w / ga))
  }

  spec_l2 <- function(w1, w2) {
    K1 <- lorentz_K(w1)
    K2 <- lorentz_K(w2)
    integrand <- function(x) {
      r1 <- vapply(x, function(xi) K1 * sum(ga / ((xi - w1)^2 + ga^2)), numeric(1))
      r2 <- vapply(x, function(xi) K2 * sum(ga / ((xi - w2)^2 + ga^2)), numeric(1))
      (r1 - r2)^2
    }
    val <- tryCatch(stats::integrate(integrand, lower = 0, upper = Inf, subdivisions = 200L)$value,
                    error = function(e) 0)
    sqrt(max(val, 0))
  }

  w1 <- lap_freqs(A1)
  w2 <- lap_freqs(A2)
  IM_raw <- spec_l2(w1, w2)
  IM_den <- im_norm_const(n)
  IM <- if (IM_den > 1e-12) max(0, min(1, IM_raw / IM_den)) else 0

  HIM <- sqrt((H^2 + IM^2) / 2)
  max(0, min(1, 1 - HIM))
}

#' Compute dynverse Hamming-Ipsen-Mikhailov (HIM) Trajectory Topology Similarity
#'
#' Quantifies both local edge correspondence (Hamming) and global Laplacian spectral similarity
#' (Ipsen-Mikhailov) between reference and simulated trajectory milestone networks.
#' Uses \code{dyneval::calculate_him} when \pkg{netdist} is available and falls back to an
#' exact pure-R implementation otherwise.
#'
#' @param ref_traj Reference trajectory list returned by \code{\link{infer_dynverse_trajectory}}.
#' @param sim_traj Simulated trajectory list returned by \code{\link{infer_dynverse_trajectory}}.
#' @return Numeric HIM similarity in \code{[0, 1]} (higher is better).
#' @export
calc_dynverse_him <- function(ref_traj, sim_traj) {
  net1 <- ref_traj$milestone_network
  net2 <- sim_traj$milestone_network
  if (nrow(net1) == 0L || nrow(net2) == 0L) return(0)

  if (requireNamespace("dyneval", quietly = TRUE)) {
    val <- tryCatch(
      suppressWarnings(as.numeric(dyneval::calculate_him(net2, net1, simplify = TRUE))),
      error = function(e) NA_real_
    )
    if (is.finite(val)) return(max(0, min(1, val)))
  }
  .traj_him_pure_r(net1, net2)
}

#' Compute dynverse Edge Flip Trajectory Topology Similarity
#'
#' Quantifies topological similarity between reference and simulated milestone graphs via
#' minimum edge flip edit distance (\code{dyneval::calculate_edge_flip}).
#'
#' @param ref_traj Reference trajectory list returned by \code{\link{infer_dynverse_trajectory}}.
#' @param sim_traj Simulated trajectory list returned by \code{\link{infer_dynverse_trajectory}}.
#' @return Numeric edge-flip similarity in \code{[0, 1]} (higher is better).
#' @export
calc_dynverse_edge_flip <- function(ref_traj, sim_traj) {
  net1 <- ref_traj$milestone_network
  net2 <- sim_traj$milestone_network
  if (nrow(net1) == 0L || nrow(net2) == 0L) return(0)

  if (requireNamespace("dyneval", quietly = TRUE)) {
    val <- tryCatch(
      suppressWarnings(as.numeric(dyneval::calculate_edge_flip(net2, net1, simplify = TRUE))),
      error = function(e) NA_real_
    )
    if (is.finite(val)) return(max(0, min(1, val)))
  }

  # Normalized binary edge disagreement fallback
  e1 <- unique(paste(pmin(net1$from, net1$to), pmax(net1$from, net1$to), sep = "||"))
  e2 <- unique(paste(pmin(net2$from, net2$to), pmax(net2$from, net2$to), sep = "||"))
  u_len <- length(union(e1, e2))
  if (u_len == 0L) return(1)
  length(intersect(e1, e2)) / u_len
}

#' Nearest-Neighbor Cell Matching in Shared PCA Space
#' @keywords internal
#' @noRd
.traj_match_cells_pca <- function(ref_pca, sim_pca) {
  k <- min(ncol(ref_pca), ncol(sim_pca))
  r_mat <- ref_pca[, seq_len(k), drop = FALSE]
  s_mat <- sim_pca[, seq_len(k), drop = FALSE]

  if (requireNamespace("RANN", quietly = TRUE)) {
    r_to_s <- RANN::nn2(data = s_mat, query = r_mat, k = 1L)$nn.idx[, 1L]
    s_to_r <- RANN::nn2(data = r_mat, query = s_mat, k = 1L)$nn.idx[, 1L]
  } else {
    d2 <- fast_cross_distance_sq(r_mat, s_mat)
    r_to_s <- max.col(-d2, ties.method = "first")
    s_to_r <- max.col(-t(d2), ties.method = "first")
  }
  list(r_to_s = r_to_s, s_to_r = s_to_r)
}

#' Compute dynverse Geodesic Distance Ordering Correlation
#'
#' Computes pairwise cell-to-milestone geodesic distances along each trajectory network,
#' aligns cells across the shared HVG PCA manifold, and computes the Spearman rank
#' correlation clipped to \code{[0, 1]} (matching \code{dyneval}'s \code{correlation} metric).
#'
#' @param ref_traj Reference trajectory list.
#' @param sim_traj Simulated trajectory list.
#' @param nn_map Precomputed nearest-neighbor index map between reference and simulated cells.
#' @return Numeric geodesic distance correlation in \code{[0, 1]} (higher is better).
#' @export
calc_dynverse_geodesic_correlation <- function(ref_traj, sim_traj, nn_map = NULL) {
  D_ref <- ref_traj$cell_milestone_geodesic
  D_sim <- sim_traj$cell_milestone_geodesic
  common_m <- intersect(colnames(D_ref), colnames(D_sim))
  if (length(common_m) < 2L) return(0)

  D_ref <- D_ref[, common_m, drop = FALSE]
  D_sim <- D_sim[, common_m, drop = FALSE]

  if (is.null(nn_map)) {
    nn_map <- .traj_match_cells_pca(ref_traj$pca_coords, sim_traj$pca_coords)
  }

  v_ref1 <- as.numeric(D_ref)
  v_sim1 <- as.numeric(D_sim[nn_map$r_to_s, , drop = FALSE])
  v_ref2 <- as.numeric(D_ref[nn_map$s_to_r, , drop = FALSE])
  v_sim2 <- as.numeric(D_sim)

  cor1 <- if (stats::sd(v_ref1) > 0 && stats::sd(v_sim1) > 0) {
    max(0, stats::cor(v_ref1, v_sim1, method = "spearman"))
  } else 0
  cor2 <- if (stats::sd(v_ref2) > 0 && stats::sd(v_sim2) > 0) {
    max(0, stats::cor(v_ref2, v_sim2, method = "spearman"))
  } else 0

  max(0, min(1, 0.5 * (cor1 + cor2)))
}

#' Compute dynverse Branch Assignment Overlap (F1_branches)
#'
#' Evaluates recovery and relevance of cellular MST branch assignments between reference
#' and simulated trajectories via bidirectional Jaccard overlap and harmonic mean F1
#' (matching \code{dyneval::calculate_mapping} \code{F1_branches}).
#'
#' @param ref_traj Reference trajectory list.
#' @param sim_traj Simulated trajectory list.
#' @param nn_map Precomputed nearest-neighbor index map between reference and simulated cells.
#' @return Numeric \code{F1_branches} score in \code{[0, 1]} (higher is better).
#' @export
calc_dynverse_f1_branches <- function(ref_traj, sim_traj, nn_map = NULL) {
  b_ref <- as.character(ref_traj$branch_assignment)
  b_sim <- as.character(sim_traj$branch_assignment)
  if (length(b_ref) == 0L || length(b_sim) == 0L) return(0)

  if (is.null(nn_map)) {
    nn_map <- .traj_match_cells_pca(ref_traj$pca_coords, sim_traj$pca_coords)
  }

  calc_pair_f1 <- function(g1, g2) {
    tab <- table(g1, g2)
    r_tot <- rowSums(tab)
    c_tot <- colSums(tab)
    union_mat <- outer(r_tot, c_tot, `+`) - tab
    jacc <- ifelse(union_mat > 0, tab / union_mat, 0)
    recovery <- mean(apply(jacc, 1L, max))
    relevance <- mean(apply(jacc, 2L, max))
    if ((recovery + relevance) <= 0) return(0)
    2 * recovery * relevance / (recovery + relevance)
  }

  f1_r <- calc_pair_f1(b_ref, b_sim[nn_map$r_to_s])
  f1_s <- calc_pair_f1(b_ref[nn_map$s_to_r], b_sim)
  max(0, min(1, 0.5 * (f1_r + f1_s)))
}

#' Compute Weighted Feature Importance Correlation (wcor_features, Saelens et al., 2019)
#'
#' Trains \pkg{ranger} Random Forest regression models to predict cell-to-milestone geodesic
#' distances from HVG expression values (\code{num.trees = 10000}, \code{mtry_prop = 0.01} (1 percent of features),
#' \code{importance = "impurity"}), averages Gini impurity importance across trajectory milestones
#' for both reference and simulated trajectories, and computes the reference-weighted Pearson
#' correlation (\code{wcor_features}) clipped to \code{[0, 1]}.
#'
#' @param ref_traj Reference trajectory list containing \code{cell_milestone_geodesic}.
#' @param sim_traj Simulated trajectory list containing \code{cell_milestone_geodesic}.
#' @param ref_expr Reference log-normalized HVG matrix (\code{hvgs x cells_ref}).
#' @param sim_expr Simulated log-normalized HVG matrix (\code{hvgs x cells_sim}).
#' @param num_trees Number of Random Forest trees per milestone (default \code{10000L}, per Saelens et al. 2019).
#' @param mtry_prop Fraction of features sampled per split (default \code{0.01}, i.e., 1 percent).
#' @param threads Number of threads for \pkg{ranger} (default \code{1L}).
#' @param seed Random seed for reproducible forest fitting (default \code{42L}).
#'
#' @return A list containing \code{wcor_features} (weighted correlation clipped to \code{[0, 1]}),
#'   \code{cor_features} (unweighted correlation clipped to \code{[0, 1]}),
#'   \code{ref_importance}, and \code{sim_importance}.
#' @export
calc_wcor_features <- function(
  ref_traj,
  sim_traj,
  ref_expr,
  sim_expr,
  num_trees = 10000L,
  mtry_prop = 0.01,
  threads = 1L,
  seed = 42L
) {
  genes <- intersect(rownames(ref_expr), rownames(sim_expr))
  if (length(genes) < 3L) {
    return(list(wcor_features = 0, cor_features = 0, ref_importance = numeric(0), sim_importance = numeric(0)))
  }

  X_ref <- as.matrix(Matrix::t(ref_expr[genes, colnames(ref_expr), drop = FALSE]))
  X_sim <- as.matrix(Matrix::t(sim_expr[genes, colnames(sim_expr), drop = FALSE]))
  D_ref <- ref_traj$cell_milestone_geodesic[rownames(X_ref), , drop = FALSE]
  D_sim <- sim_traj$cell_milestone_geodesic[rownames(X_sim), , drop = FALSE]

  P <- length(genes)
  mtry_val <- max(1L, as.integer(floor(mtry_prop * P)))
  use_ranger <- requireNamespace("ranger", quietly = TRUE)

  compute_imp <- function(X, D, seed_offset = 0L) {
    n_milestones <- ncol(D)
    imp_mat <- matrix(0, nrow = P, ncol = n_milestones, dimnames = list(genes, colnames(D)))
    valid_m <- 0L
    for (m in seq_len(n_milestones)) {
      y <- as.numeric(D[, m])
      if (length(y) < 5L || !is.finite(stats::sd(y)) || stats::sd(y) <= 1e-8) next
      valid_m <- valid_m + 1L
      if (use_ranger) {
        fit <- ranger::ranger(
          x = X,
          y = y,
          num.trees = as.integer(num_trees),
          mtry = mtry_val,
          importance = "impurity",
          write.forest = FALSE,
          num.threads = as.integer(max(1L, threads)),
          seed = as.integer(seed + seed_offset + m),
          verbose = FALSE
        )
        vi <- as.numeric(fit$variable.importance[genes])
        vi[!is.finite(vi)] <- 0
        imp_mat[, m] <- vi
      } else {
        cs <- apply(X, 2L, function(col) {
          if (stats::sd(col) <= 1e-8) return(0)
          abs(stats::cor(col, y))
        })
        cs[!is.finite(cs)] <- 0
        imp_mat[, m] <- cs
      }
    }
    if (valid_m == 0L) return(stats::setNames(rep(0, P), genes))
    rowMeans(imp_mat[, seq_len(n_milestones), drop = FALSE], na.rm = TRUE)
  }

  imp_ref <- compute_imp(X_ref, D_ref, seed_offset = 0L)
  imp_sim <- compute_imp(X_sim, D_sim, seed_offset = 1000L)

  if (!is.finite(stats::sd(imp_ref)) || !is.finite(stats::sd(imp_sim)) ||
      stats::sd(imp_ref) <= 1e-12 || stats::sd(imp_sim) <= 1e-12 ||
      sum(imp_ref > 0) < 2L) {
    return(list(wcor_features = 0, cor_features = 0, ref_importance = imp_ref, sim_importance = imp_sim))
  }

  cor_unw <- max(0, min(1, stats::cor(imp_ref, imp_sim)))
  wt <- pmax(imp_ref, 0)
  wcor <- tryCatch({
    cw <- stats::cov.wt(cbind(imp_ref, imp_sim), wt = wt, cor = TRUE)
    val <- cw$cor[1, 2]
    if (!is.finite(val)) 0 else max(0, min(1, val))
  }, error = function(e) 0)

  list(
    wcor_features = wcor,
    cor_features = cor_unw,
    ref_importance = imp_ref,
    sim_importance = imp_sim
  )
}

#' Full Category VII Trajectory & Lineage Topology Evaluation (dynverse Workflow)
#'
#' Evaluates trajectory and lineage preservation using the standardized \emph{dynverse}
#' benchmarking framework (\emph{Saelens et al., Nature Biotechnology, 2019}):
#' \enumerate{
#'   \item Standardized cell outlier QC (3 MADs), low-expression gene filtering
#'         (\eqn{\ge 5\%} cells, mean \eqn{> 0.02}), and HVG selection (\code{scran::modelGeneVar}).
#'   \item Shared HVG PCA projection and Docker-free Minimum Spanning Tree (MST) milestone
#'         trajectory inference with continuous cell-to-milestone geodesic distances.
#'   \item Five complementary trajectory metrics: \code{him} (Hamming-Ipsen-Mikhailov similarity),
#'         \code{edge_flip} (minimum edge flip similarity), \code{geodesic_correlation},
#'         \code{F1_branches}, and \code{wcor_features} (\pkg{ranger} Random Forest feature
#'         importance weighted correlation with \code{num_trees = 10000} and 1 percent \code{mtry}).
#' }
#'
#' @param ref_data Reference count matrix, \code{SingleCellExperiment}, or \code{Seurat} object.
#' @param sim_data Simulated count matrix, \code{SingleCellExperiment}, or \code{Seurat} object.
#' @param cell_types_ref Optional vector of reference cell types.
#' @param cell_types_sim Optional vector of simulated cell types.
#' @param num_trees Number of \pkg{ranger} trees per milestone for \code{wcor_features} (default \code{10000L}).
#' @param nmads MAD cutoff for outlier cell QC (default \code{3}).
#' @param min_cell_frac Minimum cell fraction per gene (default \code{0.05}).
#' @param min_mean_expr Minimum average gene count (default \code{0.02}).
#' @param fdr_threshold FDR threshold for HVG selection (default \code{0.05}).
#' @param bio_var_threshold Biological variance threshold for HVG selection (default \code{0.5}).
#' @param threads Number of CPU threads for \pkg{ranger} (default \code{1L}).
#' @param seed Random seed (default \code{42L}).
#' @param verbose Logical; print progress messages.
#'
#' @return A named list of Category VII \emph{dynverse} trajectory accuracy metrics:
#'   \describe{
#'     \item{him}{Hamming-Ipsen-Mikhailov topology similarity (\code{[0, 1]}).}
#'     \item{edge_flip}{Edge flip topology similarity (\code{[0, 1]}).}
#'     \item{geodesic_correlation}{Spearman correlation of aligned geodesic distances (\code{[0, 1]}).}
#'     \item{F1_branches}{Macro F1 overlap of MST branch assignments (\code{[0, 1]}).}
#'     \item{wcor_features}{Reference-weighted Pearson correlation of Random Forest gene importances (\code{[0, 1]}).}
#'     \item{cor_features}{Unweighted Pearson correlation of Random Forest gene importances.}
#'     \item{n_hvgs}{Number of HVGs retained for trajectory evaluation.}
#'     \item{ref_trajectory}{Inferred reference milestone trajectory object.}
#'     \item{sim_trajectory}{Inferred simulated milestone trajectory object.}
#'     \item{preprocessing_summary}{Cell and HVG filtering summary list.}
#'   }
#' @examples
#' data(example_scrna, package = "scSimEval")
#' evaluate_trajectory_metrics(
#'   example_scrna$ref, example_scrna$sim,
#'   cell_types_ref = example_scrna$cell_types,
#'   cell_types_sim = example_scrna$cell_types,
#'   num_trees = 500L
#' )
#' @export
evaluate_trajectory_metrics <- function(
  ref_data,
  sim_data,
  cell_types_ref = NULL,
  cell_types_sim = NULL,
  num_trees = 10000L,
  nmads = 3,
  min_cell_frac = 0.05,
  min_mean_expr = 0.02,
  fdr_threshold = 0.05,
  bio_var_threshold = 0.5,
  threads = 1L,
  seed = 42L,
  verbose = FALSE
) {
  # Numeric pseudotime vector guard
  if (is.numeric(ref_data) && is.numeric(sim_data) &&
      is.null(dim(ref_data)) && is.null(dim(sim_data))) {
    aligned <- align_distributions(as.numeric(ref_data), as.numeric(sim_data))
    pt_cor <- suppressWarnings(stats::cor(aligned$ref, aligned$sim, method = "spearman"))
    if (is.na(pt_cor) || !is.finite(pt_cor)) pt_cor <- 0
    return(list(
      him = NA_real_,
      edge_flip = NA_real_,
      geodesic_correlation = max(0, pt_cor),
      F1_branches = NA_real_,
      wcor_features = NA_real_,
      cor_features = NA_real_,
      n_hvgs = 0L
    ))
  }

  # 1. Saelens et al. (2019) QC & HVG selection
  if (verbose) message("  Running Saelens et al. (2019) QC & HVG preprocessing...")
  prep <- preprocess_trajectory_hvgs(
    ref_data = ref_data,
    sim_data = sim_data,
    cell_types_ref = cell_types_ref,
    cell_types_sim = cell_types_sim,
    nmads = nmads,
    min_cell_frac = min_cell_frac,
    min_mean_expr = min_mean_expr,
    fdr_threshold = fdr_threshold,
    bio_var_threshold = bio_var_threshold
  )

  # 2. Shared HVG PCA projection (fit on reference, project both ref & sim)
  n_pcs <- min(15L, nrow(prep$ref_expr) - 1L, ncol(prep$ref_expr) - 1L, ncol(prep$sim_expr) - 1L)
  n_pcs <- max(2L, n_pcs)
  ref_dense_t <- as.matrix(Matrix::t(prep$ref_expr))
  sim_dense_t <- as.matrix(Matrix::t(prep$sim_expr))
  keep_cols <- which(apply(ref_dense_t, 2L, stats::sd) > 1e-8)
  if (length(keep_cols) < 2L) keep_cols <- seq_len(ncol(ref_dense_t))
  n_pcs <- min(n_pcs, length(keep_cols), nrow(ref_dense_t) - 1L, nrow(sim_dense_t) - 1L)
  n_pcs <- max(2L, n_pcs)

  pca_fit <- stats::prcomp(ref_dense_t[, keep_cols, drop = FALSE], center = TRUE, scale. = FALSE, rank. = n_pcs)
  ref_pca <- pca_fit$x[, seq_len(n_pcs), drop = FALSE]
  sim_pca <- stats::predict(pca_fit, newdata = sim_dense_t[, keep_cols, drop = FALSE])[, seq_len(n_pcs), drop = FALSE]
  rownames(ref_pca) <- colnames(prep$ref_expr)
  rownames(sim_pca) <- colnames(prep$sim_expr)

  # 3. Infer Docker-free dynwrap MST milestone trajectories
  has_shared_ct <- !is.null(prep$cell_types_ref) && !is.null(prep$cell_types_sim) &&
    length(intersect(unique(prep$cell_types_ref), unique(prep$cell_types_sim))) >= 2L

  ref_traj <- infer_dynverse_trajectory(
    expr_mat = prep$ref_expr,
    cell_types = if (has_shared_ct) prep$cell_types_ref else NULL,
    pca_coords = ref_pca
  )
  sim_traj <- infer_dynverse_trajectory(
    expr_mat = prep$sim_expr,
    cell_types = if (has_shared_ct) prep$cell_types_sim else NULL,
    pca_coords = sim_pca
  )

  if (!has_shared_ct) {
    sim_traj <- .traj_align_unsupervised_milestones(ref_traj, sim_traj)
  }

  # 4. Shared PCA nearest-neighbor matching between reference & simulated cells
  nn_map <- .traj_match_cells_pca(ref_pca, sim_pca)

  # 5. Compute the 5 Category VII dynverse metrics
  him_val <- tryCatch(calc_dynverse_him(ref_traj, sim_traj), error = function(e) NA_real_)
  ef_val <- tryCatch(calc_dynverse_edge_flip(ref_traj, sim_traj), error = function(e) NA_real_)
  geo_val <- tryCatch(calc_dynverse_geodesic_correlation(ref_traj, sim_traj, nn_map = nn_map),
                      error = function(e) NA_real_)
  f1_val <- tryCatch(calc_dynverse_f1_branches(ref_traj, sim_traj, nn_map = nn_map),
                     error = function(e) NA_real_)
  wcor_res <- tryCatch(
    calc_wcor_features(
      ref_traj = ref_traj,
      sim_traj = sim_traj,
      ref_expr = prep$ref_expr,
      sim_expr = prep$sim_expr,
      num_trees = num_trees,
      threads = threads,
      seed = seed
    ),
    error = function(e) list(wcor_features = NA_real_, cor_features = NA_real_)
  )

  list(
    him = him_val,
    edge_flip = ef_val,
    geodesic_correlation = geo_val,
    F1_branches = f1_val,
    wcor_features = wcor_res$wcor_features,
    cor_features = wcor_res$cor_features,
    n_hvgs = prep$summary$n_hvgs,
    ref_trajectory = ref_traj,
    sim_trajectory = sim_traj,
    preprocessing_summary = prep$summary
  )
}
