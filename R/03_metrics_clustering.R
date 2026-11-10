#' Calculate Average Silhouette Width (ASW)
#' @param dist_mat Distance matrix among cells or numeric matrix (features x cells).
#' @param cluster_labels Vector of cluster or cell-type labels.
#' @return Mean silhouette width.
#' @examples
#' data <- matrix(stats::rnorm(200), 20, 10)
#' cl <- factor(rep(c("A", "B"), each = 10))
#' dist_mat <- stats::dist(data)
#' calc_silhouette(dist_mat, cl)
#' @export
calc_silhouette <- function(dist_mat, cluster_labels) {
  if (!inherits(dist_mat, "dist")) {
    dist_mat <- stats::dist(t(as.matrix(dist_mat)))
  }
  cluster_labels <- as.numeric(as.factor(cluster_labels))
  if (requireNamespace("cluster", quietly = TRUE)) {
    sil <- tryCatch(cluster::silhouette(cluster_labels, dist_mat), error = function(e) NULL)
    if (is.null(sil) || is.na(sil)[1] || (!is.matrix(sil) && !inherits(sil, "silhouette"))) {
      NA_real_
    } else {
      mean(sil[, 3], na.rm = TRUE)
    }
  } else {
    NA_real_
  }
}

#' Calculate Davies-Bouldin Index (DB)
#' @param data Matrix with features in rows, cells in columns.
#' @param cluster_labels Cluster assignments.
#' @return Davies-Bouldin index (lower is better).
#' @examples
#' data <- matrix(stats::rnorm(200), 20, 10)
#' cl <- factor(rep(c("A", "B"), each = 5))
#' dist_mat <- stats::dist(data)
#' calc_davies_bouldin(data, cl)
#' @export
calc_davies_bouldin <- function(data, cluster_labels) {
  cluster_labels <- as.numeric(as.factor(cluster_labels))
  if (length(unique(cluster_labels)) < 2 || length(cluster_labels) < 3) return(NA_real_)
  val <- tryCatch({
    if (requireNamespace("clusterSim", quietly = TRUE)) {
      res <- clusterSim::index.DB(t(as.matrix(data)), cluster_labels)
      if (is.list(res) && "DB" %in% names(res)) as.numeric(res$DB) else as.numeric(res)
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)
  if (is.na(val)) {
    val <- tryCatch({
      X <- t(as.matrix(data))
      u_cl <- unique(cluster_labels)
      k <- length(u_cl)
      if (k < 2) return(NA_real_)
      centroids <- list()
      s_i <- numeric(k)
      for (idx_c in seq_len(k)) {
        cl <- u_cl[idx_c]
        idx <- which(cluster_labels == cl)
        cl_mat <- X[idx, , drop = FALSE]
        cent <- colMeans(cl_mat)
        centroids[[idx_c]] <- cent
        diffs <- sweep(cl_mat, 2, cent, "-")
        s_i[idx_c] <- mean(sqrt(rowSums(diffs^2)))
      }
      r_ij <- matrix(0, k, k)
      for (i in seq_len(k)) {
        for (j in seq_len(k)) {
          if (i != j) {
            m_ij <- sqrt(sum((centroids[[i]] - centroids[[j]])^2))
            r_ij[i, j] <- if (m_ij > 0) (s_i[i] + s_i[j]) / m_ij else 0
          }
        }
      }
      mean(apply(r_ij, 1, max))
    }, error = function(e) NA_real_)
  }
  val
}

#' Calculate Adjusted Rand Index (ARI)
#' @param pred Predicted cluster labels.
#' @param truth Ground truth cell type labels.
#' @return ARI value between -1 and 1.
#' @examples
#' pred <- factor(rep(c("A", "B"), each = 20))
#' truth <- factor(rep(c("A", "B"), each = 20))
#' calc_ari(pred, truth)
#' @export
calc_ari <- function(pred, truth) {
  if (requireNamespace("mclust", quietly = TRUE)) {
    mclust::adjustedRandIndex(pred, truth)
  } else {
    # Fallback ARI calculation via contingency table
    tab <- table(pred, truth)
    a <- sum(choose(tab, 2))
    b <- sum(choose(rowSums(tab), 2))
    c <- sum(choose(colSums(tab), 2))
    d <- choose(sum(tab), 2)
    expected <- (b * c) / d
    max_idx <- (b + c) / 2
    (a - expected) / (max_idx - expected)
  }
}

#' Calculate Normalized Mutual Information (NMI)
#' @param pred Predicted cluster labels.
#' @param truth Ground truth cell type labels.
#' @return NMI value between 0 and 1.
#' @examples
#' pred <- factor(rep(c("A", "B"), each = 20))
#' truth <- factor(rep(c("A", "B"), each = 20))
#' calc_nmi(pred, truth)
#' @export
calc_nmi <- function(pred, truth) {
  tab <- table(pred, truth)
  n <- sum(tab)
  p_ij <- tab / n
  p_i <- rowSums(tab) / n
  p_j <- colSums(tab) / n
  
  # Entropies
  h_i <- -sum(p_i * log2(p_i + 1e-12))
  h_j <- -sum(p_j * log2(p_j + 1e-12))
  
  # Mutual info
  mi <- 0
  for (i in seq_len(nrow(tab))) {
    for (j in seq_len(ncol(tab))) {
      if (tab[i, j] > 0) {
        mi <- mi + p_ij[i, j] * log2(p_ij[i, j] / (p_i[i] * p_j[j]))
      }
    }
  }
  if (h_i + h_j == 0) return(1)
  2 * mi / (h_i + h_j)
}

#' Calculate Neighborhood Purity
#'
#' Integrated from simPIC (Chugh et al., 2024; bluster::neighborPurity).
#' For each cell, evaluates the proportion of its k-nearest neighbors that
#' share the same cluster or cell-type identity.
#'
#' @param data Matrix of coordinates (cells x dimensions) or expression matrix (features x cells) or dist matrix.
#' @param cluster_labels Vector of cluster or cell-type labels.
#' @param k Number of nearest neighbors (default 5\% of cells).
#' @param is_distance Logical, whether data is already a distance matrix. Default FALSE.
#'
#' @return A numeric vector of neighborhood purity scores (values from 0 to 1),
#'   with an attribute "mean_purity" containing the overall average.
#' @examples
#' data <- matrix(stats::rnorm(200), 20, 10)
#' cl <- factor(rep(c("A", "B"), each = 10))
#' dist_mat <- stats::dist(data)
#' calc_neighborhood_purity(data, cl)
#' @export
calc_neighborhood_purity <- function(data, cluster_labels, k = NULL, is_distance = FALSE) {
  cluster_labels <- as.character(cluster_labels)
  n_cells <- length(cluster_labels)
  if (n_cells < 2) return(numeric(0))
  
  if (is.null(k)) {
    k <- max(5, min(50, round(0.05 * n_cells)))
  }
  k <- min(k, n_cells - 1)
  
  if (requireNamespace("bluster", quietly = TRUE) && !is_distance) {
    mat <- if (nrow(data) == n_cells) as.matrix(data) else t(as.matrix(data))
    pur <- tryCatch({
      res <- bluster::neighborPurity(mat, clusters = cluster_labels, k = k)
      as.numeric(res$purity)
    }, error = function(e) NULL)
    if (!is.null(pur)) {
      attr(pur, "mean_purity") <- mean(pur, na.rm = TRUE)
      return(pur)
    }
  }
  
  if (is_distance) {
    dist_mat <- as.matrix(data)
  } else {
    mat <- if (nrow(data) == n_cells) as.matrix(data) else t(as.matrix(data))
    if (requireNamespace("RANN", quietly = TRUE)) {
      nn_res <- RANN::nn2(mat, k = k + 1)
      nn_idx <- nn_res$nn.idx[, 2:(k + 1), drop = FALSE]
      purities <- vapply(seq_len(n_cells), function(i) {
        mean(cluster_labels[nn_idx[i, ]] == cluster_labels[i], na.rm = TRUE)
      }, numeric(1))
      attr(purities, "mean_purity") <- mean(purities, na.rm = TRUE)
      return(purities)
    } else {
      dist_mat <- as.matrix(stats::dist(mat))
    }
  }
  
  diag(dist_mat) <- Inf
  purities <- vapply(seq_len(n_cells), function(i) {
    neighbors <- order(dist_mat[i, ])[seq_len(k)]
    mean(cluster_labels[neighbors] == cluster_labels[i], na.rm = TRUE)
  }, numeric(1))
  attr(purities, "mean_purity") <- mean(purities, na.rm = TRUE)
  purities
}

#' Compute Clustering Deviation Index (CDI) for Single-Cell Clustering Evaluation
#'
#' Evaluates the unsupervised goodness-of-fit of single-cell clustering labels directly
#' on count matrices using the official Bioconductor \pkg{CDI} package
#' (Fang et al., \emph{Genome Biology}, 2022; DOI: 10.1186/s13059-022-02825-5).
#' CDI models UMI counts under cluster- and gene-specific Negative Binomial distributions
#' adjusted by cell size factors, returning AIC- and BIC-penalized deviation scores.
#'
#' To minimize time and space complexity from \eqn{O(G \cdot N)} dense allocations down to
#' \eqn{O(\mathrm{nnz}(X) + G_0 \cdot N_{\mathrm{eval}})}, sparse \code{dgCMatrix} inputs are
#' stratified-subsampled across cells and pre-filtered by zero-proportion in sparse space prior
#' to Working Dispersion Score (WDS) feature selection (\code{CDI::feature_gene_selection}),
#' median-of-ratios size factor estimation (\code{CDI::size_factor}), and Negative Binomial MLE
#' evaluation (\code{CDI::calculate_CDI}).
#'
#' @param counts Count matrix (\code{matrix}, sparse \code{dgCMatrix}, or \code{SingleCellExperiment};
#'   genes x cells or cells x genes).
#' @param cluster_labels Vector or factor of cluster / cell-type assignments for each cell.
#' @param batch_label Optional vector or factor of batch assignments for each cell (enables
#'   \pkg{CDI}'s batch-common vs. batch-specific NB likelihood ratio test).
#' @param size_factors Optional numeric vector of positive cell-specific size factors. If \code{NULL},
#'   computed via \code{CDI::size_factor()}.
#' @param nfeature Integer number of informative feature genes to select via
#'   \code{CDI::feature_gene_selection()} (default \code{500L}, per Fang et al., 2022).
#' @param method Feature selection method passed to \code{CDI::feature_gene_selection()}:
#'   \code{"wds"} (Working Dispersion Score, default) or \code{"vst"}.
#' @param zp_threshold Numeric zero-proportion cutoff in \code{(0, 1)} used to pre-filter
#'   ultra-sparse features prior to WDS selection (default \code{0.95}). Automatically relaxed
#'   for ultra-sparse modalities such as scATAC-seq.
#' @param lrt_pval_threshold Numeric p-value threshold in \code{(0, 1)} for multi-batch LRT
#'   in \code{CDI::calculate_CDI()} (default \code{0.01}).
#' @param max_cells Maximum number of cells evaluated via proportional stratified subsampling
#'   across clusters (default \code{1000L}; set \code{Inf} or \code{NULL} to evaluate all cells).
#' @param normalize_per_obs Logical; if \code{TRUE} (default), divides raw \code{CDI_AIC} and
#'   \code{CDI_BIC} by \eqn{N_{\mathrm{eval}} \times G_{\mathrm{selected}}} so scores are
#'   scale-invariant across simulated datasets of differing dimensions.
#' @param threads Integer number of parallel worker threads for \code{BiocParallel} (default \code{1L}).
#' @param seed Integer random seed used only during stratified cell subsampling (default \code{42L}).
#' @param top_features Optional deprecated alias for \code{nfeature}.
#'
#' @return A list containing:
#'   \item{cdi_aic}{Per-observation normalized (or raw if \code{normalize_per_obs = FALSE}) CDI-AIC score (lower indicates better clustering, sensitive to subtypes).}
#'   \item{cdi_bic}{Per-observation normalized (or raw if \code{normalize_per_obs = FALSE}) CDI-BIC score (lower indicates better clustering, favors main cell types).}
#'   \item{cdi_aic_raw}{Unnormalized \code{CDI_AIC} returned directly by \code{CDI::calculate_CDI()}.}
#'   \item{cdi_bic_raw}{Unnormalized \code{CDI_BIC} returned directly by \code{CDI::calculate_CDI()}.}
#'   \item{loglik}{Total maximized NB log-likelihood (\code{-neg_llk_val}) across selected features.}
#'   \item{n_clusters}{Number of clusters evaluated.}
#'   \item{n_features}{Number of selected feature genes evaluated.}
#'   \item{n_cells}{Number of cells evaluated.}
#' @importFrom CDI calculate_CDI feature_gene_selection size_factor
#' @examples
#' counts <- matrix(stats::rpois(1200, lambda = 3), nrow = 30, ncol = 40)
#' cl <- factor(rep(c("TypeA", "TypeB"), each = 20))
#' calc_cdi(counts, cl)
#' @export
calc_cdi <- function(
  counts,
  cluster_labels,
  batch_label = NULL,
  size_factors = NULL,
  nfeature = 500L,
  method = c("wds", "vst"),
  zp_threshold = 0.95,
  lrt_pval_threshold = 0.01,
  max_cells = 1000L,
  normalize_per_obs = TRUE,
  threads = 1L,
  seed = 42L,
  top_features = NULL
) {
  method <- match.arg(method)
  if (!is.null(top_features)) {
    nfeature <- as.integer(top_features)
  }
  nfeature <- max(2L, as.integer(nfeature))

  na_out <- function(k_val = NA_integer_) {
    list(
      cdi_aic = NA_real_,
      cdi_bic = NA_real_,
      cdi_aic_raw = NA_real_,
      cdi_bic_raw = NA_real_,
      loglik = NA_real_,
      n_clusters = k_val,
      n_features = 0L,
      n_cells = 0L
    )
  }

  if (inherits(counts, "SingleCellExperiment") || inherits(counts, "SummarizedExperiment")) {
    assay_nms <- SummarizedExperiment::assayNames(counts)
    slot_use <- if ("counts" %in% assay_nms) "counts" else assay_nms[1]
    counts <- SummarizedExperiment::assay(counts, slot_use)
  } else if (is.data.frame(counts)) {
    counts <- as.matrix(counts)
  }

  n_lbl <- length(cluster_labels)
  if (ncol(counts) != n_lbl && nrow(counts) == n_lbl) {
    counts <- Matrix::t(counts)
  }
  if (ncol(counts) != n_lbl) {
    return(na_out(0L))
  }

  # 1. Clean NA labels & enforce CDI's >= 3 cells per cluster requirement in O(N)
  valid_cell <- !is.na(cluster_labels)
  if (!is.null(batch_label)) {
    if (length(batch_label) != n_lbl) {
      batch_label <- NULL
    } else {
      valid_cell <- valid_cell & !is.na(batch_label)
    }
  }
  if (sum(valid_cell) < 6L) {
    return(na_out(length(unique(stats::na.omit(cluster_labels)))))
  }
  if (!all(valid_cell)) {
    counts <- counts[, valid_cell, drop = FALSE]
    cluster_labels <- cluster_labels[valid_cell]
    if (!is.null(batch_label)) batch_label <- batch_label[valid_cell]
    if (!is.null(size_factors) && length(size_factors) == n_lbl) {
      size_factors <- size_factors[valid_cell]
    }
  }

  cl_chr <- as.character(cluster_labels)
  cl_tab <- table(cl_chr)
  keep_clusters <- names(cl_tab)[cl_tab >= 3L]
  if (length(keep_clusters) < 2L) {
    return(na_out(length(cl_tab)))
  }
  if (length(keep_clusters) < length(cl_tab)) {
    keep_idx <- which(cl_chr %in% keep_clusters)
    prev_ncol <- ncol(counts)
    counts <- counts[, keep_idx, drop = FALSE]
    cl_chr <- cl_chr[keep_idx]
    if (!is.null(batch_label)) batch_label <- batch_label[keep_idx]
    if (!is.null(size_factors) && length(size_factors) == prev_ncol) {
      size_factors <- size_factors[keep_idx]
    }
  }

  # Validate multi-batch cell counts; fall back to single-batch CDI if batch splits are too small
  if (!is.null(batch_label)) {
    bt_chr <- as.character(batch_label)
    if (length(unique(bt_chr)) < 2L) {
      batch_label <- NULL
    } else {
      comb_tab <- table(paste0(cl_chr, "||", bt_chr))
      valid_comb <- comb_tab >= 3L
      surviving_cl <- unique(sub("\\|\\|.*$", "", names(comb_tab)[valid_comb]))
      if (length(surviving_cl) < 2L) {
        batch_label <- NULL
      } else {
        batch_label <- bt_chr
      }
    }
  }

  # 2. Proportional stratified cell subsampling in sparse space: O(N)
  N_cur <- ncol(counts)
  if (!is.null(max_cells) && is.finite(max_cells) && N_cur > max_cells && max_cells >= 6L) {
    has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit({
      if (has_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(seed))

    idx_by_cl <- split(seq_len(N_cur), cl_chr)
    sub_idx <- unlist(lapply(idx_by_cl, function(ix) {
      n_take <- max(3L, as.integer(round(length(ix) * (max_cells / N_cur))))
      n_take <- min(length(ix), n_take)
      if (n_take == length(ix)) ix else sample(ix, size = n_take, replace = FALSE)
    }), use.names = FALSE)
    sub_idx <- sort(sub_idx)

    counts <- counts[, sub_idx, drop = FALSE]
    cl_chr <- cl_chr[sub_idx]
    if (!is.null(batch_label)) batch_label <- batch_label[sub_idx]
    if (!is.null(size_factors) && length(size_factors) == N_cur) {
      size_factors <- size_factors[sub_idx]
    }
  }

  # 3. Sparse zero-proportion pre-filtering in O(nnz(X)) before any dense coercion
  #    Matches CDI::feature_gene_selection's zp < zp_threshold rule identically while avoiding
  #    densifying thousands of ultra-sparse/all-zero rows.
  zp_vec <- fast_sparse_zero_fraction_feature(counts)
  zp_eff <- zp_threshold
  active_genes <- which(zp_vec < zp_eff)
  if (length(active_genes) < min(20L, nrow(counts))) {
    zp_eff <- 0.995
    active_genes <- which(zp_vec < zp_eff)
  }
  if (length(active_genes) < 2L) {
    active_genes <- which(zp_vec < 1)
  }
  if (length(active_genes) < 2L) {
    return(na_out(length(unique(cl_chr))))
  }
  counts <- counts[active_genes, , drop = FALSE]

  # Drop empty cells (all-zero library across active genes)
  cell_sums <- Matrix::colSums(counts)
  pos_cells <- which(cell_sums > 0)
  if (length(pos_cells) < 6L) {
    return(na_out(length(unique(cl_chr))))
  }
  if (length(pos_cells) < ncol(counts)) {
    prev_n <- ncol(counts)
    counts <- counts[, pos_cells, drop = FALSE]
    cl_chr <- cl_chr[pos_cells]
    if (!is.null(batch_label)) batch_label <- batch_label[pos_cells]
    if (!is.null(size_factors) && length(size_factors) == prev_n) {
      size_factors <- size_factors[pos_cells]
    }
    cl_tab2 <- table(cl_chr)
    keep_cl2 <- names(cl_tab2)[cl_tab2 >= 3L]
    if (length(keep_cl2) < 2L) return(na_out(length(cl_tab2)))
    if (length(keep_cl2) < length(cl_tab2)) {
      k2_idx <- which(cl_chr %in% keep_cl2)
      counts <- counts[, k2_idx, drop = FALSE]
      cl_chr <- cl_chr[k2_idx]
      if (!is.null(batch_label)) batch_label <- batch_label[k2_idx]
      if (!is.null(size_factors) && length(size_factors) == length(pos_cells)) {
        size_factors <- size_factors[k2_idx]
      }
    }
  }

  # 4. Materialize compact active matrix & ensure non-negative integer counts for NB PMF
  counts_mat <- as.matrix(counts)
  rm(counts)
  if (anyNA(counts_mat)) counts_mat[is.na(counts_mat)] <- 0
  if (any(counts_mat < 0)) counts_mat[counts_mat < 0] <- 0
  if (is.double(counts_mat) && any(counts_mat != floor(counts_mat))) {
    counts_mat <- round(counts_mat)
  }
  storage.mode(counts_mat) <- "integer"

  # Ensure unique gene/cell names (prevents rowname collision inside CDI::calculate_CDI)
  rownames(counts_mat) <- paste0("g", seq_len(nrow(counts_mat)))
  colnames(counts_mat) <- paste0("c", seq_len(ncol(counts_mat)))

  # 5. Official CDI::size_factor() on active count matrix
  if (is.null(size_factors) || length(size_factors) != ncol(counts_mat)) {
    size_factors <- CDI::size_factor(X = counts_mat)
  }
  size_factors <- as.numeric(size_factors)
  bad_sf <- is.na(size_factors) | !is.finite(size_factors) | size_factors <= 0
  if (any(bad_sf)) {
    min_pos <- if (any(!bad_sf)) min(size_factors[!bad_sf]) else 1
    size_factors[bad_sf] <- max(min_pos, 1e-4)
  }

  # 6. Official CDI::feature_gene_selection() (WDS default, Fang et al., 2022)
  n_sel <- min(nfeature, nrow(counts_mat))
  if (nrow(counts_mat) > n_sel) {
    fg_idx <- tryCatch(
      CDI::feature_gene_selection(
        X = counts_mat,
        batch_label = batch_label,
        method = method,
        nfeature = n_sel,
        zp_threshold = zp_eff
      ),
      error = function(e) integer(0)
    )
    fg_idx <- unique(as.integer(fg_idx))
    fg_idx <- fg_idx[!is.na(fg_idx) & fg_idx >= 1L & fg_idx <= nrow(counts_mat)]
    if (length(fg_idx) >= 2L) {
      sub_mat <- counts_mat[fg_idx, , drop = FALSE]
    } else {
      sub_mat <- counts_mat[seq_len(n_sel), , drop = FALSE]
    }
  } else {
    sub_mat <- counts_mat
  }
  rm(counts_mat)
  rownames(sub_mat) <- paste0("g", seq_len(nrow(sub_mat)))

  # 7. Adaptive BiocParallel backend (avoids socket startup latency on small tasks)
  K_eval <- length(unique(cl_chr))
  G_eval <- nrow(sub_mat)
  N_eval <- ncol(sub_mat)
  threads <- max(1L, as.integer(threads))

  bpparam <- if (threads > 1L && (G_eval * K_eval >= 1500L || .Platform$OS.type != "windows")) {
    if (.Platform$OS.type == "windows") {
      BiocParallel::SnowParam(workers = threads)
    } else {
      BiocParallel::MulticoreParam(workers = threads)
    }
  } else {
    BiocParallel::SerialParam()
  }

  # 8. Official CDI::calculate_CDI()
  cand_df <- data.frame(cand = cl_chr, stringsAsFactors = FALSE)
  cdi_df <- CDI::calculate_CDI(
    X = sub_mat,
    feature_gene_index = NULL,
    cand_lab_df = cand_df,
    cell_size_factor = size_factors,
    batch_label = batch_label,
    lrt_pval_threshold = lrt_pval_threshold,
    BPPARAM = bpparam
  )

  aic_raw <- as.numeric(cdi_df$CDI_AIC[1])
  bic_raw <- as.numeric(cdi_df$CDI_BIC[1])
  neg_llk <- if ("neg_llk_val" %in% colnames(cdi_df)) as.numeric(cdi_df$neg_llk_val[1]) else NA_real_
  n_clust_out <- if ("N_cluster" %in% colnames(cdi_df)) as.integer(cdi_df$N_cluster[1]) else K_eval

  denom <- if (isTRUE(normalize_per_obs)) as.numeric(N_eval) * as.numeric(G_eval) else 1.0

  list(
    cdi_aic = aic_raw / denom,
    cdi_bic = bic_raw / denom,
    cdi_aic_raw = aic_raw,
    cdi_bic_raw = bic_raw,
    loglik = -neg_llk,
    n_clusters = n_clust_out,
    n_features = G_eval,
    n_cells = N_eval
  )
}

#' Full Clustering Performance Evaluation
#'
#' Evaluates both unsupervised cluster separation (ASW, DB, Neighborhood Purity, CDI)
#' and supervised concordance against ground truth labels (ARI, NMI).
#'
#' @param data Count or normalized matrix (features x cells).
#' @param cluster_info Cluster or cell type labels.
#' @param pred_clusters Optional predicted cluster labels (if different from ground truth). If NULL, k-means is automatically run.
#' @param dist_mat Optional precomputed distance matrix.
#' @param cell_types Optional alias for cluster_info.
#' @param ref_data Optional reference data to compute reference clustering quality.
#' @param ref_cell_types Optional reference cell type labels (if ref_data has different cells than data).
#' @param batch_info Optional batch labels per cell passed to \code{\link{calc_cdi}}.
#' @param threads Integer number of parallel CPU threads passed to \code{\link{calc_cdi}} (default \code{1L}).
#'
#' @return A named list of all clustering metrics.
#' @examples
#' data <- matrix(stats::rpois(200, lambda = 3), 20, 10)
#' cl <- factor(rep(c("A", "B"), each = 5))
#' dist_mat <- stats::dist(t(data))
#' evaluate_clustering_metrics(data, cl)
#' @export
evaluate_clustering_metrics <- function(
  data,
  cluster_info = NULL,
  pred_clusters = NULL,
  dist_mat = NULL,
  cell_types = NULL,
  ref_data = NULL,
  ref_cell_types = NULL,
  batch_info = NULL,
  threads = 1L
) {
  n_cells_sim <- ncol(data)
  if (is.null(cluster_info) && !is.null(cell_types)) {
    if (is.list(cell_types) && !is.data.frame(cell_types) && !is.factor(cell_types)) {
      lens <- vapply(cell_types, length, integer(1))
      match_idx <- which(lens == n_cells_sim)
      if (length(match_idx) > 0) cluster_info <- cell_types[[match_idx[1]]] else cluster_info <- cell_types[[1]]
    } else {
      cluster_info <- cell_types
    }
  }
  
  if (is.null(cluster_info) || length(cluster_info) != n_cells_sim || length(unique(stats::na.omit(cluster_info))) < 2) {
    k_fallback <- min(5, max(2, n_cells_sim - 1))
    d_norm_temp <- normalize_counts(data)
    if (inherits(d_norm_temp, "Matrix")) d_norm_temp <- as.matrix(d_norm_temp)
    pca_temp <- fast_pca(t(d_norm_temp), n_pcs = min(15, n_cells_sim - 1, nrow(d_norm_temp) - 1), n_hvg = 2000, scale = FALSE)
    km_fallback <- tryCatch(
      suppressWarnings(stats::kmeans(pca_temp, centers = k_fallback, iter.max = 50, nstart = 5)),
      error = function(e) list(cluster = rep(seq_len(k_fallback), length.out = n_cells_sim))
    )
    cluster_info <- factor(paste0("Cluster_", km_fallback$cluster))
    cluster_info <- as.factor(cluster_info)
    # SimBench standard: skip cell types with fewer than 10 cells
    ct_tab <- table(cluster_info)
    rare_cts <- names(ct_tab)[ct_tab < 10]
    if (length(rare_cts) > 0 && length(ct_tab) > length(rare_cts)) {
      keep_idx <- which(!cluster_info %in% rare_cts)
      if (length(keep_idx) >= 10 && length(unique(cluster_info[keep_idx])) >= 2) {
        if (!is.null(batch_info) && length(batch_info) == ncol(data)) {
          batch_info <- batch_info[keep_idx]
        }
        data <- data[, keep_idx, drop = FALSE]
        cluster_info <- droplevels(factor(cluster_info[keep_idx]))
        n_cells_sim <- ncol(data)
        if (!is.null(dist_mat)) dist_mat <- NULL
      }
    }
  }
  
  d_use <- normalize_counts(data)
  if (inherits(d_use, "Matrix")) d_use <- as.matrix(d_use)

  # Fast low-dimensional PCA projection for high-dimensional single-cell clustering (scDesign3, scVI standard)
  pca_coords <- NULL
  if (nrow(d_use) > 30) {
    pca_coords <- fast_pca(t(d_use), n_pcs = min(30, ncol(d_use) - 1, nrow(d_use) - 1), n_hvg = 2000, scale = FALSE)
  }
  
  if (is.null(dist_mat)) {
    dist_mat <- if (!is.null(pca_coords)) stats::dist(pca_coords) else stats::dist(t(as.matrix(d_use)))
  }
  
  db_sim <- if (length(unique(stats::na.omit(cluster_info))) > 1) {
    tryCatch(calc_davies_bouldin(if (!is.null(pca_coords)) t(pca_coords) else d_use, cluster_info), error = function(e) NA_real_)
  } else {
    NA_real_
  }

  res <- list(
    silhouette_sim = calc_silhouette(dist_mat, cluster_info),
    silhouette = calc_silhouette(dist_mat, cluster_info),
    davies_bouldin = db_sim,
    davies_bouldin_sim = db_sim,
    neighborhood_purity = mean(calc_neighborhood_purity(dist_mat, cluster_info, is_distance = TRUE), na.rm = TRUE)
  )
  
  if (!is.null(ref_data)) {
    r_use <- normalize_counts(ref_data)
    if (inherits(r_use, "Matrix")) r_use <- as.matrix(r_use)
    ref_pca <- if (nrow(r_use) > 30) fast_pca(t(r_use), n_pcs = min(30, ncol(r_use) - 1, nrow(r_use) - 1), n_hvg = 2000, scale = FALSE) else NULL
    ref_dist <- if (!is.null(ref_pca)) stats::dist(ref_pca) else stats::dist(t(as.matrix(r_use)))
    ref_cl <- if (!is.null(ref_cell_types) && length(ref_cell_types) == ncol(r_use)) {
      as.factor(ref_cell_types)
    } else if (length(cluster_info) == ncol(r_use)) {
      cluster_info
    } else {
      NULL
    }
    if (!is.null(ref_cl) && length(ref_cl) == ncol(r_use)) {
      ref_tab <- table(ref_cl)
      rare_ref <- names(ref_tab)[ref_tab < 10]
      if (length(rare_ref) > 0 && length(ref_tab) > length(rare_ref)) {
        keep_ref <- which(!ref_cl %in% rare_ref)
        if (length(keep_ref) >= 10 && length(unique(ref_cl[keep_ref])) >= 2) {
          r_use <- r_use[, keep_ref, drop = FALSE]
          ref_cl <- droplevels(factor(ref_cl[keep_ref]))
          ref_dist <- stats::dist(t(as.matrix(r_use)))
        }
      }
    }
    if (!is.null(ref_cl) && length(unique(stats::na.omit(ref_cl))) > 1) {
      res$silhouette_ref <- calc_silhouette(ref_dist, ref_cl)
      res$davies_bouldin_ref <- calc_davies_bouldin(if (!is.null(ref_pca)) t(ref_pca) else r_use, ref_cl)
    }
  }
  
  # Unsupervised Clustering Deviation Index (CDI; Fang et al., Genome Biology 2022)
  cdi_res <- tryCatch(
    calc_cdi(data, cluster_info, batch_label = batch_info, threads = threads),
    error = function(e) NULL
  )
  if (!is.null(cdi_res) && !is.na(cdi_res$cdi_aic)) {
    res$cdi <- cdi_res$cdi_aic
    res$CDI_AIC <- cdi_res$cdi_aic
    res$CDI_BIC <- cdi_res$cdi_bic
  }
  
  # If pred_clusters is not provided, run k-means with k = nlevels(cluster_info)
  if (is.null(pred_clusters) && nlevels(cluster_info) >= 2) {
    k_centers <- nlevels(cluster_info)
    km_input <- if (!is.null(pca_coords)) pca_coords else t(as.matrix(d_use))
    km_fit <- tryCatch(suppressWarnings(stats::kmeans(km_input, centers = k_centers, iter.max = 50, nstart = 5)), error = function(e) NULL)
    if (!is.null(km_fit)) {
      pred_clusters <- factor(km_fit$cluster)
    }
  }
  
  # Supervised concordance if predicted clusters are supplied or inferred
  if (!is.null(pred_clusters)) {
    res$ari <- calc_ari(pred_clusters, cluster_info)
    res$ARI <- res$ari
    res$nmi <- calc_nmi(pred_clusters, cluster_info)
    res$NMI <- res$nmi
  }
  
  return(res)
}

# -----------------------------------------------------------------------------
# Deep Generative Manifold Quality: Precision, Recall & Diversity (ACTIVA, scDiffusion)
# -----------------------------------------------------------------------------

#' Generative Precision and Recall for Single-Cell Manifolds
#'
#' Evaluates the realism (Precision) and coverage/diversity (Recall) of generated
#' single-cell distributions relative to empirical reference data using k-NN hyperspheres.
#'
#' @param ref_mat Matrix of reference cells (cells x features or features x cells).
#' @param sim_mat Matrix of simulated cells (cells x features or features x cells).
#' @param cells_as_cols Logical, whether cells are columns (default TRUE).
#' @param k Number of nearest neighbors to define the local manifold radius (default 5).
#' @param max_cells Maximum number of cells to subsample (default 500).
#' @param max_features Maximum number of features to retain for speed (default 1000).
#'
#' @return A list containing Generative Precision, Generative Recall, and Generative F1.
#' @examples
#' ref <- matrix(stats::rnorm(200), 20, 10)
#' sim <- matrix(stats::rnorm(200), 20, 10)
#' calc_generative_precision_recall(ref, sim)
#' @export
calc_generative_precision_recall <- function(
  ref_mat,
  sim_mat,
  cells_as_cols = TRUE,
  k = 5,
  max_cells = 500,
  max_features = 1000
) {
  if (inherits(ref_mat, "SingleCellExperiment")) {
    ref_mat <- SummarizedExperiment::assay(ref_mat, "counts")
  }
  if (inherits(sim_mat, "SingleCellExperiment")) {
    sim_mat <- SummarizedExperiment::assay(sim_mat, "counts")
  }
  
  if (!is.matrix(ref_mat) && !inherits(ref_mat, "Matrix")) ref_mat <- as.matrix(ref_mat)
  if (!is.matrix(sim_mat) && !inherits(sim_mat, "Matrix")) sim_mat <- as.matrix(sim_mat)
  
  if (!cells_as_cols) {
    ref_mat <- t(ref_mat)
    sim_mat <- t(sim_mat)
    cells_as_cols <- TRUE
  }
  
  # Align features between reference and simulation
  r_names <- rownames(ref_mat)
  s_names <- rownames(sim_mat)
  if (!is.null(r_names) && !is.null(s_names) && length(intersect(r_names, s_names)) >= 5) {
    common_feats <- intersect(r_names, s_names)
    ref_mat <- ref_mat[common_feats, , drop = FALSE]
    sim_mat <- sim_mat[common_feats, , drop = FALSE]
  } else {
    n_common <- min(nrow(ref_mat), nrow(sim_mat))
    if (n_common < 2) return(list(generative_precision = NA_real_, generative_recall = NA_real_, generative_f1 = NA_real_))
    ref_mat <- ref_mat[seq_len(n_common), , drop = FALSE]
    sim_mat <- sim_mat[seq_len(n_common), , drop = FALSE]
  }
  
  # Standard protocol: library-size normalize & log-transform counts
  ref_norm <- normalize_counts(ref_mat)
  sim_norm <- normalize_counts(sim_mat)
  
  # Subsample features to top HVGs if large
  if (nrow(ref_norm) > max_features) {
    rv <- fast_row_vars(ref_norm)
    top_f <- order(rv, decreasing = TRUE)[seq_len(max_features)]
    ref_norm <- ref_norm[top_f, , drop = FALSE]
    sim_norm <- sim_norm[top_f, , drop = FALSE]
  }
  
  x_raw <- t(as.matrix(ref_norm))
  y_raw <- t(as.matrix(sim_norm))
  
  if (nrow(x_raw) > max_cells) x_raw <- x_raw[sample(nrow(x_raw), max_cells), , drop = FALSE]
  if (nrow(y_raw) > max_cells) y_raw <- y_raw[sample(nrow(y_raw), max_cells), , drop = FALSE]
  
  # Joint PCA projection (30 PCs) for density sphere calculation
  combined <- rbind(x_raw, y_raw)
  n_pcs_use <- min(30, ncol(combined) - 1, nrow(combined) - 1)
  if (n_pcs_use >= 2 && ncol(combined) > n_pcs_use) {
    pca_comb <- fast_pca(combined, n_pcs = n_pcs_use, scale = FALSE)
    x <- pca_comb[seq_len(nrow(x_raw)), , drop = FALSE]
    y <- pca_comb[(nrow(x_raw) + 1):nrow(combined), , drop = FALSE]
  } else {
    x <- x_raw
    y <- y_raw
  }
  
  m <- nrow(x)
  n <- nrow(y)
  
  k_x <- min(k, m - 1)
  k_y <- min(k, n - 1)
  
  if (k_x < 1 || k_y < 1) {
    return(list(generative_precision = NA_real_, generative_recall = NA_real_, generative_f1 = NA_real_))
  }
  
  # Pairwise distance within reference: dist_xx
  dist_xx <- as.matrix(stats::dist(x))
  diag(dist_xx) <- Inf
  radius_x <- apply(dist_xx, 1, function(row) sort(row)[k_x])
  
  # Pairwise distance within simulation: dist_yy
  dist_yy <- as.matrix(stats::dist(y))
  diag(dist_yy) <- Inf
  radius_y <- apply(dist_yy, 1, function(row) sort(row)[k_y])
  
  # Fast BLAS cross distance between reference and simulation: dist_xy (m x n)
  dist_xy <- sqrt(fast_cross_distance_sq(x, y))
  
  # Precision: For each simulated cell y_j, is it within at least one reference sphere?
  precision_hits <- apply(dist_xy, 2, function(col_dists) {
    any(col_dists <= radius_x)
  })
  precision_val <- mean(precision_hits)
  
  # Recall: For each reference cell x_i, is it within at least one simulated sphere?
  recall_hits <- apply(dist_xy, 1, function(row_dists) {
    any(row_dists <= radius_y)
  })
  recall_val <- mean(recall_hits)
  
  f1_val <- if ((precision_val + recall_val) > 0) {
    2 * (precision_val * recall_val) / (precision_val + recall_val)
  } else 0
  
  list(
    generative_precision = as.numeric(precision_val),
    generative_recall = as.numeric(recall_val),
    generative_f1 = as.numeric(f1_val)
  )
}

#' Calculate Maximum Mean Discrepancy (MMD) on Single-Cell PCA Manifolds
#'
#' Integrated from ACTIVA (Heydari et al., Bioinformatics 2022; SindiLab/ACTIVA)
#' and scGAN (Marouf et al., Nature Communications 2020). Evaluates how closely
#' the multivariate distribution of simulated single-cell profiles matches the
#' reference biological manifold in a shared Principal Component Analysis (PCA)
#' latent space using a multi-scale Radial Basis Function (RBF) Gaussian kernel
#' in a Reproducing Kernel Hilbert Space (RKHS).
#'
#' @param ref_mat Matrix of reference cells (features x cells by default, or SingleCellExperiment).
#' @param sim_mat Matrix of simulated cells (features x cells by default, or SingleCellExperiment).
#' @param cells_as_cols Logical, whether cells are columns (default TRUE).
#' @param n_pcs Number of principal components for shared latent projection (default 50, as in ACTIVA).
#' @param scales Numeric vector of bandwidth multipliers for the multi-scale RBF kernel
#'   (default \code{c(0.01, 0.1, 1, 10, 100)}).
#' @param max_cells Maximum number of cells per dataset to subsample for speed (default 1000).
#' @param max_features Maximum number of highly variable features to retain prior to PCA (default 2000).
#'
#' @return Numeric scalar representing the multi-scale RBF Maximum Mean Discrepancy
#'   (\code{MMD >= 0}; lower values indicate higher manifold fidelity).
#' @examples
#' ref <- matrix(stats::rnorm(400), 20, 20)
#' sim <- matrix(stats::rnorm(400), 20, 20)
#' calc_mmd(ref, sim)
#' @export
calc_mmd <- function(
  ref_mat,
  sim_mat,
  cells_as_cols = TRUE,
  n_pcs = 50,
  scales = c(0.01, 0.1, 1, 10, 100),
  max_cells = 1000,
  max_features = 2000
) {
  if (inherits(ref_mat, "SingleCellExperiment")) {
    ref_mat <- SummarizedExperiment::assay(ref_mat, "counts")
  }
  if (inherits(sim_mat, "SingleCellExperiment")) {
    sim_mat <- SummarizedExperiment::assay(sim_mat, "counts")
  }
  
  if (!is.matrix(ref_mat) && !inherits(ref_mat, "Matrix")) ref_mat <- as.matrix(ref_mat)
  if (!is.matrix(sim_mat) && !inherits(sim_mat, "Matrix")) sim_mat <- as.matrix(sim_mat)
  
  if (!cells_as_cols) {
    ref_mat <- t(ref_mat)
    sim_mat <- t(sim_mat)
  }
  
  # Align features between reference and simulation
  r_names <- rownames(ref_mat)
  s_names <- rownames(sim_mat)
  if (!is.null(r_names) && !is.null(s_names) && length(intersect(r_names, s_names)) >= 5) {
    common_feats <- intersect(r_names, s_names)
    ref_mat <- ref_mat[common_feats, , drop = FALSE]
    sim_mat <- sim_mat[common_feats, , drop = FALSE]
  } else {
    n_common <- min(nrow(ref_mat), nrow(sim_mat))
    if (n_common < 2) return(NA_real_)
    ref_mat <- ref_mat[seq_len(n_common), , drop = FALSE]
    sim_mat <- sim_mat[seq_len(n_common), , drop = FALSE]
  }
  
  if (ncol(ref_mat) < 2 || ncol(sim_mat) < 2) return(NA_real_)
  
  # Standard library-size log-normalization
  ref_norm <- normalize_counts(ref_mat)
  sim_norm <- normalize_counts(sim_mat)
  
  # Retain top HVGs if feature space is large
  if (nrow(ref_norm) > max_features) {
    rv <- fast_row_vars(ref_norm)
    top_f <- order(rv, decreasing = TRUE)[seq_len(max_features)]
    ref_norm <- ref_norm[top_f, , drop = FALSE]
    sim_norm <- sim_norm[top_f, , drop = FALSE]
  }
  
  x_raw <- t(as.matrix(ref_norm))
  y_raw <- t(as.matrix(sim_norm))
  
  if (nrow(x_raw) > max_cells) {
    idx_x <- round(seq(1, nrow(x_raw), length.out = max_cells))
    x_raw <- x_raw[idx_x, , drop = FALSE]
  }
  if (nrow(y_raw) > max_cells) {
    idx_y <- round(seq(1, nrow(y_raw), length.out = max_cells))
    y_raw <- y_raw[idx_y, , drop = FALSE]
  }
  
  # Project into shared PCA space (ACTIVA default: 50 PCs)
  combined <- rbind(x_raw, y_raw)
  n_pcs_use <- min(as.integer(n_pcs), ncol(combined) - 1L, nrow(combined) - 1L)
  if (n_pcs_use >= 2L && ncol(combined) > n_pcs_use) {
    pca_comb <- fast_pca(combined, n_pcs = n_pcs_use, scale = FALSE)
    x <- pca_comb[seq_len(nrow(x_raw)), , drop = FALSE]
    y <- pca_comb[(nrow(x_raw) + 1L):nrow(combined), , drop = FALSE]
  } else {
    x <- x_raw
    y <- y_raw
  }
  
  # Pairwise squared Euclidean distances in PCA space
  d_xx_sq <- fast_cross_distance_sq(x, x)
  d_yy_sq <- fast_cross_distance_sq(y, y)
  d_xy_sq <- fast_cross_distance_sq(x, y)
  
  # Median heuristic base bandwidth
  med_sq <- stats::median(d_xy_sq[d_xy_sq > 0], na.rm = TRUE)
  if (is.na(med_sq) || med_sq <= 0) med_sq <- 1.0
  
  # Multi-scale RBF kernel evaluation across bandwidth scales
  k_xx <- matrix(0, nrow = nrow(x), ncol = nrow(x))
  k_yy <- matrix(0, nrow = nrow(y), ncol = nrow(y))
  k_xy <- matrix(0, nrow = nrow(x), ncol = nrow(y))
  
  for (s in scales) {
    bw <- 2 * s * med_sq
    if (bw <= 0) bw <- 1e-8
    k_xx <- k_xx + exp(-d_xx_sq / bw)
    k_yy <- k_yy + exp(-d_yy_sq / bw)
    k_xy <- k_xy + exp(-d_xy_sq / bw)
  }
  
  n_scales <- length(scales)
  k_xx <- k_xx / n_scales
  k_yy <- k_yy / n_scales
  k_xy <- k_xy / n_scales
  
  mmd_sq <- mean(k_xx, na.rm = TRUE) + mean(k_yy, na.rm = TRUE) - 2 * mean(k_xy, na.rm = TRUE)
  as.numeric(sqrt(max(0, mmd_sq)))
}
