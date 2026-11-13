#' Align Two Empirical Distributions
#'
#' Sorts both distributions and aligns their lengths using either quantile interpolation
#' or nearest-neighbor resampling, ensuring element-wise distance metrics (MAD, MAE, RMSE)
#' are mathematically well-defined even when the number of cells or features differs.
#'
#' @param ref Numeric vector of reference values.
#' @param sim Numeric vector of simulated values.
#' @param method Alignment method: "quantile" (default, interpolated quantiles) or "resample" (resampling).
#' @param n_points Number of points to evaluate when method = "quantile". Default is max(length(ref), length(sim)).
#'
#' @return A list with aligned numeric vectors: \code{list(ref = ..., sim = ...)}.
#' @keywords internal
#' @noRd
align_distributions <- function(ref, sim, method = c("quantile", "resample"), n_points = NULL) {
  method <- match.arg(method)
  ref <- sort(stats::na.omit(as.numeric(ref)))
  sim <- sort(stats::na.omit(as.numeric(sim)))
  
  n_ref <- length(ref)
  n_sim <- length(sim)
  
  if (n_ref == 0 || n_sim == 0) {
    return(list(ref = numeric(0), sim = numeric(0)))
  }
  
  if (n_ref == n_sim) {
    return(list(ref = ref, sim = sim))
  }
  
  if (method == "quantile") {
    if (is.null(n_points)) {
      n_points <- min(max(n_ref, n_sim), 2500)
    }
    probs <- seq(0, 1, length.out = n_points)
    aligned_ref <- stats::quantile(ref, probs = probs, names = FALSE, type = 7)
    aligned_sim <- stats::quantile(sim, probs = probs, names = FALSE, type = 7)
    return(list(ref = aligned_ref, sim = aligned_sim))
  } else {
    diff_len <- abs(n_ref - n_sim)
    if (n_ref > n_sim) {
      add_vals <- sample(sim, diff_len, replace = TRUE)
      sim <- sort(c(sim, add_vals))
    } else {
      add_vals <- sample(ref, diff_len, replace = TRUE)
      ref <- sort(c(ref, add_vals))
    }
    return(list(ref = ref, sim = sim))
  }
}

#' Calculate Proportion of Outliers in a Numeric Vector
#'
#' Uses the standard Tukey 1.5 * IQR criterion on quantiles, as used in simpipe.
#'
#' @param x Numeric vector.
#' @return Outlier proportion (between 0 and 1).
#' @examples
#' ref <- stats::rnorm(50)
#' sim <- stats::rnorm(50)
#' calc_outlier_proportion(ref)
#' @export
calc_outlier_proportion <- function(x) {
  x <- stats::na.omit(as.numeric(x))
  if (length(x) < 4) return(0)
  q <- stats::quantile(x, probs = c(0.25, 0.5, 0.75))
  iqr <- stats::IQR(x)
  lower_bound <- q[1] - 1.5 * iqr
  upper_bound <- q[3] + 1.5 * iqr
  sum(x < lower_bound | x > upper_bound) / length(x)
}

#' Hungarian Maximum Weight Matching for Cluster Assignment
#'
#' Matches predicted cluster labels to ground truth labels using the Hungarian
#' algorithm (Kuhn-Munkres) on the contingency/F1 matrix.
#' Adapted from HelenaLC/simulation-comparison and simpipe.
#'
#' @param pred Vector of predicted cluster labels.
#' @param truth Vector of ground truth cell-type labels.
#'
#' @return A list with matched cluster pairs, precision, recall, and macro F1 score.
#' @keywords internal
#' @noRd
hungarian_match <- function(pred, truth) {
  valid <- !is.na(pred) & !is.na(truth)
  pred <- as.character(pred[valid])
  truth <- as.character(truth[valid])
  
  u_pred <- unique(pred)
  u_truth <- unique(truth)
  
  # Compute pairwise F1 matrix
  f1_mat <- matrix(0, nrow = length(u_pred), ncol = length(u_truth),
                   dimnames = list(u_pred, u_truth))
  pr_mat <- re_mat <- f1_mat
  
  for (i in seq_along(u_pred)) {
    p_mask <- pred == u_pred[i]
    n_p <- sum(p_mask)
    for (j in seq_along(u_truth)) {
      t_mask <- truth == u_truth[j]
      n_t <- sum(t_mask)
      n_overlap <- sum(p_mask & t_mask)
      
      pr <- if (n_p > 0) n_overlap / n_p else 0
      re <- if (n_t > 0) n_overlap / n_t else 0
      f1 <- if (pr + re > 0) 2 * (pr * re) / (pr + re) else 0
      
      pr_mat[i, j] <- pr
      re_mat[i, j] <- re
      f1_mat[i, j] <- f1
    }
  }
  
  # Optimal matching
  matched_pairs <- data.frame()
  if (requireNamespace("clue", quietly = TRUE)) {
    # Maximize total F1 score using solve_LSAP
    if (nrow(f1_mat) <= ncol(f1_mat)) {
      assignment <- clue::solve_LSAP(f1_mat, maximum = TRUE)
      for (i in seq_along(assignment)) {
        j <- as.numeric(assignment[i])
        matched_pairs <- rbind(matched_pairs, data.frame(
          pred = u_pred[i], truth = u_truth[j],
          precision = pr_mat[i, j], recall = re_mat[i, j], F1 = f1_mat[i, j],
          stringsAsFactors = FALSE
        ))
      }
    } else {
      assignment <- clue::solve_LSAP(t(f1_mat), maximum = TRUE)
      for (j in seq_along(assignment)) {
        i <- as.numeric(assignment[j])
        matched_pairs <- rbind(matched_pairs, data.frame(
          pred = u_pred[i], truth = u_truth[j],
          precision = pr_mat[i, j], recall = re_mat[i, j], F1 = f1_mat[i, j],
          stringsAsFactors = FALSE
        ))
      }
    }
  } else {
    # Greedy fallback if clue is not installed
    used_t <- c()
    for (i in seq_along(u_pred)) {
      avail_j <- setdiff(seq_along(u_truth), used_t)
      if (length(avail_j) == 0) break
      best_j <- avail_j[which.max(f1_mat[i, avail_j])]
      used_t <- c(used_t, best_j)
      matched_pairs <- rbind(matched_pairs, data.frame(
        pred = u_pred[i], truth = u_truth[best_j],
        precision = pr_mat[i, best_j], recall = re_mat[i, best_j], F1 = f1_mat[i, best_j],
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Align predicted cluster labels to ground truth using matched pairs
  mapping <- stats::setNames(matched_pairs$truth, matched_pairs$pred)
  aligned_pred <- unname(mapping[as.character(pred)])
  acc <- mean(aligned_pred == as.character(truth), na.rm = TRUE)
  
  list(
    matched_pairs = matched_pairs,
    accuracy = acc,
    mean_precision = mean(matched_pairs$precision, na.rm = TRUE),
    mean_recall = mean(matched_pairs$recall, na.rm = TRUE),
    mean_F1 = mean(matched_pairs$F1, na.rm = TRUE),
    aligned_pred = aligned_pred
  )
}

# -----------------------------------------------------------------------------
# High-Performance Computational Acceleration Helpers
# (Vectorized Sparse Matrix Arithmetic, BLAS Distance, & Truncated SVD)
# -----------------------------------------------------------------------------

#' Fast Row-Wise Variances for Dense or Sparse Single-Cell Matrices
#' @param mat Matrix or dgCMatrix (features x cells).
#' @return Numeric vector of row variances.
#' @export
fast_row_vars <- function(mat) {
  if (inherits(mat, "dgCMatrix")) {
    n_cols <- ncol(mat)
    if (n_cols <= 1) return(rep(0, nrow(mat)))
    rm <- Matrix::rowMeans(mat)
    mat_sq <- mat
    mat_sq@x <- mat_sq@x^2
    sq_rm <- Matrix::rowMeans(mat_sq)
    v <- (sq_rm - rm^2) * (n_cols / (n_cols - 1))
    v[v < 0] <- 0
    return(as.numeric(v))
  } else if (requireNamespace("matrixStats", quietly = TRUE)) {
    return(as.numeric(matrixStats::rowVars(as.matrix(mat))))
  } else {
    n_cols <- ncol(mat)
    if (n_cols <= 1) return(rep(0, nrow(mat)))
    rm <- rowMeans(mat)
    sq_rm <- rowMeans(mat^2)
    v <- (sq_rm - rm^2) * (n_cols / (n_cols - 1))
    v[v < 0] <- 0
    return(as.numeric(v))
  }
}

#' Fast Column-Wise Variances for Dense or Sparse Single-Cell Matrices
#' @param mat Matrix or dgCMatrix (features x cells).
#' @return Numeric vector of column variances.
#' @export
fast_col_vars <- function(mat) {
  if (inherits(mat, "dgCMatrix")) {
    n_rows <- nrow(mat)
    if (n_rows <= 1) return(rep(0, ncol(mat)))
    cm <- Matrix::colMeans(mat)
    mat_sq <- mat
    mat_sq@x <- mat_sq@x^2
    sq_cm <- Matrix::colMeans(mat_sq)
    v <- (sq_cm - cm^2) * (n_rows / (n_rows - 1))
    v[v < 0] <- 0
    return(as.numeric(v))
  } else if (requireNamespace("matrixStats", quietly = TRUE)) {
    return(as.numeric(matrixStats::colVars(as.matrix(mat))))
  } else {
    n_rows <- nrow(mat)
    if (n_rows <= 1) return(rep(0, ncol(mat)))
    cm <- colMeans(mat)
    sq_cm <- colMeans(mat^2)
    v <- (sq_cm - cm^2) * (n_rows / (n_rows - 1))
    v[v < 0] <- 0
    return(as.numeric(v))
  }
}

#' Fast Zero Fraction per Cell (Column-wise)
#' @param mat Matrix or dgCMatrix (features x cells).
#' @return Numeric vector of zero fraction per cell.
#' @export
fast_sparse_zero_fraction_cell <- function(mat) {
  if (inherits(mat, "dgCMatrix")) {
    nz_per_col <- diff(mat@p)
    return(as.numeric((nrow(mat) - nz_per_col) / nrow(mat)))
  }
  as.numeric(colMeans(mat == 0))
}

#' Fast Zero Fraction per Feature (Row-wise)
#' @param mat Matrix or dgCMatrix (features x cells).
#' @return Numeric vector of zero fraction per feature.
#' @export
fast_sparse_zero_fraction_feature <- function(mat) {
  if (inherits(mat, "dgCMatrix")) {
    nz_idx <- which(mat@x != 0)
    nz_per_row <- tabulate(mat@i[nz_idx] + 1, nbins = nrow(mat))
    return(as.numeric((ncol(mat) - nz_per_row) / ncol(mat)))
  }
  as.numeric(rowMeans(mat == 0))
}

#' Fast BLAS-Accelerated Cross-Distance Squared Matrix
#'
#' Computes ||X_i - Y_j||^2 using BLAS matrix product: ||x||^2 + ||y||^2 - 2 <x, y>.
#' @param X Matrix (m x p).
#' @param Y Matrix (n x p).
#' @return Matrix of squared Euclidean distances (m x n).
#' @export
fast_cross_distance_sq <- function(X, Y) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  x_norm <- rowSums(X^2)
  y_norm <- rowSums(Y^2)
  d2 <- outer(x_norm, y_norm, "+") - 2 * (X %*% t(Y))
  d2[d2 < 0] <- 0
  d2
}

#' Fast Truncated PCA on Top Highly Variable Features
#' @param mat Matrix or dgCMatrix (cells x features).
#' @param n_pcs Number of principal components. Default 20.
#' @param n_hvg Number of top highly variable genes to retain. Default 2000.
#' @param scale Logical, whether to scale features. Default FALSE.
#' @param normalize Logical, whether to apply library-size log-normalization if raw counts. Default FALSE.
#' @return Matrix of cell principal component embeddings (cells x n_pcs).
#' @export
fast_pca <- function(mat, n_pcs = 20, n_hvg = 2000, scale = FALSE, normalize = FALSE) {
  if (inherits(mat, "Matrix")) {
    mat <- as.matrix(mat)
  }
  if (normalize) {
    # mat is cells x features; normalize_counts expects features x cells
    mat <- t(as.matrix(normalize_counts(t(mat))))
  }
  n_feats <- ncol(mat)
  if (n_feats > n_hvg) {
    feat_vars <- fast_col_vars(mat)
    top_feats <- order(feat_vars, decreasing = TRUE)[seq_len(n_hvg)]
    mat <- mat[, top_feats, drop = FALSE]
  }
  rank_use <- min(n_pcs, nrow(mat) - 1, ncol(mat) - 1)
  if (rank_use < 2) return(as.matrix(mat)[, seq_len(min(1, ncol(mat))), drop = FALSE])
  
  min_dim <- min(nrow(mat), ncol(mat))
  if (requireNamespace("irlba", quietly = TRUE) && rank_use >= 3 && min_dim > 50 && rank_use < 0.5 * min_dim) {
    res <- tryCatch(suppressWarnings(irlba::prcomp_irlba(as.matrix(mat), n = rank_use, scale. = scale)), error = function(e) NULL)
    if (!is.null(res)) return(res$x)
  }
  stats::prcomp(as.matrix(mat), rank. = rank_use, scale. = scale)$x
}

#' Subsample Large Vectors for Accelerated Distributional Distance Estimation
#' @param x Numeric vector.
#' @param max_n Maximum sample size (default 10000, matching SimBench framework).
#' @return Subsampled numeric vector if length(x) > max_n, else x.
#' @export
subsample_dist_vector <- function(x, max_n = 10000) {
  x <- stats::na.omit(as.numeric(x))
  if (length(x) > max_n) {
    sample(x, max_n)
  } else {
    x
  }
}

#' Normalize Single-Cell Count Matrix (Library Size + Log Transformation)
#'
#' Converts raw single-cell sequencing count matrices to library-size normalized,
#' log-transformed values (e.g. log2(CPM + 1) or log(normalized counts + 1)),
#' mitigating library size confounding across distance, correlation, and manifold metrics
#' as established across benchmark literature (SimCH, scDesign3, SPsimSeq, GDSim, SymBench, cfDiffusion, scCross).
#'
#' @param mat Expression or count matrix (features x cells), dgCMatrix, SingleCellExperiment, or Seurat object.
#' Detect Whether Single-Cell Input Represents Chromatin Accessibility (scATAC-seq)
#'
#' Checks explicit modality strings, Seurat default assay names, SingleCellExperiment
#' experiment names, or peak coordinate feature names (e.g., \code{Peak_01},
#' \code{chr1:100-200}).
#'
#' @param data Matrix, dgCMatrix, SingleCellExperiment, or Seurat object.
#' @param modality Optional character modality hint.
#' @return Logical scalar (\code{TRUE} for scATAC-seq, \code{FALSE} otherwise).
#' @keywords internal
#' @noRd
.is_atac_modality <- function(data = NULL, modality = NULL) {
  if (!is.null(modality) && length(modality) > 0 && !is.na(modality[1]) && nzchar(as.character(modality[1]))) {
    m_low <- tolower(as.character(modality[1]))
    if (grepl("atac|peak|chromatin", m_low)) return(TRUE)
    if (grepl("rna|gene|transcript", m_low)) return(FALSE)
  }

  if (is.null(data)) return(FALSE)

  if (inherits(data, "Seurat") && requireNamespace("SeuratObject", quietly = TRUE)) {
    def_assay <- tryCatch(SeuratObject::DefaultAssay(data), error = function(e) NULL)
    if (!is.null(def_assay) && grepl("atac|peak|chromatin", tolower(def_assay))) {
      return(TRUE)
    }
  } else if (inherits(data, "SingleCellExperiment") && requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    exp_nm <- tryCatch(SingleCellExperiment::mainExpName(data), error = function(e) NULL)
    if (!is.null(exp_nm) && grepl("atac|peak|chromatin", tolower(exp_nm))) {
      return(TRUE)
    }
  }

  rn <- tryCatch(rownames(data), error = function(e) NULL)
  if (!is.null(rn) && length(rn) > 0) {
    is_peak_name <- grepl("^(peak|chr[0-9xyXYM]+[:-]|[A-Za-z0-9_.-]+[:-][0-9]+-[0-9]+)", rn, ignore.case = TRUE)
    if (mean(is_peak_name) > 0.5) return(TRUE)
  }

  FALSE
}

#' Perform TF-IDF Normalization for scATAC-seq Peak Matrices
#'
#' Computes Term Frequency–Inverse Document Frequency (TF-IDF) normalization
#' following the Signac (\code{RunTFIDF(..., method = 1)}) and simPIC
#' (Chugh et al., 2024) standard for scATAC-seq chromatin accessibility matrices:
#' \deqn{\text{TF-IDF}_{i,j} = \log\left(1 + \frac{X_{i,j}}{N_j} \cdot \frac{N_{\text{cells}}}{\sum_k X_{i,k}} \cdot \text{scale\_factor}\right)}
#'
#' @param mat A peak-by-cell count matrix, sparse \code{dgCMatrix}, \code{SingleCellExperiment}, or \code{Seurat} object.
#' @param scale_factor Numeric scaling factor applied prior to \code{log1p}. Default is \code{1e4} (10,000).
#' @return TF-IDF transformed matrix (preserving \code{dgCMatrix} sparsity when sparse input is supplied).
#' @examples
#' set.seed(42)
#' m <- matrix(rpois(100, 0.5), nrow = 10)
#' tfidf_mat <- tfidf_transform(m)
#' @export
tfidf_transform <- function(mat, scale_factor = 1e4) {
  if (is.null(mat)) return(NULL)
  if (inherits(mat, "SingleCellExperiment")) {
    mat <- SummarizedExperiment::assay(mat, "counts")
  } else if (inherits(mat, "Seurat")) {
    mat <- tryCatch(
      Seurat::GetAssayData(mat, slot = "counts"),
      error = function(e) Seurat::GetAssayData(mat, layer = "counts")
    )
  }

  libs <- as.numeric(Matrix::colSums(mat))
  libs[libs == 0] <- 1

  rs <- as.numeric(Matrix::rowSums(mat))
  idf <- ifelse(rs > 0, ncol(mat) / rs, 0)

  if (inherits(mat, "dgCMatrix")) {
    norm_mat <- mat
    col_idx <- rep(seq_len(ncol(mat)), diff(mat@p))
    row_idx <- mat@i + 1L
    norm_mat@x <- log1p((norm_mat@x / libs[col_idx]) * idf[row_idx] * scale_factor)
    return(norm_mat)
  } else {
    mat <- as.matrix(mat)
    tf <- t(t(mat) / libs)
    scaled <- sweep(tf, 1L, idf, `*`) * scale_factor
    return(log1p(scaled))
  }
}

.get_seurat_counts <- function(obj) {
  tryCatch(
    Seurat::GetAssayData(obj, layer = "counts"),
    error = function(e) Seurat::GetAssayData(obj, slot = "counts")
  )
}

#' Fast Count Normalization (Log2-CPM or TF-IDF)
#'
#' Converts raw single-cell count matrices into log-normalized CPM (\code{"lognorm"})
#' for scRNA-seq or TF-IDF (\code{"tfidf"}) for scATAC-seq, preserving
#' \code{dgCMatrix} sparsity.
#'
#' @param mat A count matrix (genes/features x cells), \code{dgCMatrix}, \code{SingleCellExperiment}, or \code{Seurat} object.
#' @param scale_factor Numeric scaling factor for counts per cell. Default is 1e6 (CPM) for \code{"lognorm"}, or 1e4 for \code{"tfidf"}.
#' @param log_base Numeric base for the logarithm in \code{"lognorm"} mode: 2 (default, for log2(CPM+1)) or natural log \code{exp(1)}.
#' @param pseudo_count Numeric pseudo-count added before logarithm in \code{"lognorm"} mode. Default is 1.
#' @param method Normalization method: \code{"lognorm"} (default) or \code{"tfidf"} (Signac/simPIC scATAC-seq standard).
#' @return Normalized matrix of transformed values.
#' @examples
#' data(example_scrna, package = "scSimEval")
#' norm_mat <- normalize_counts(example_scrna$ref)
#' @export
normalize_counts <- function(mat, scale_factor = 1e6, log_base = 2, pseudo_count = 1,
                             method = c("lognorm", "tfidf")) {
  method <- match.arg(method)
  if (is.null(mat)) return(NULL)
  if (method == "tfidf") {
    sf <- if (missing(scale_factor)) 1e4 else scale_factor
    return(tfidf_transform(mat, scale_factor = sf))
  }
  if (inherits(mat, "SingleCellExperiment")) {
    mat <- SummarizedExperiment::assay(mat, "counts")
  } else if (inherits(mat, "Seurat")) {
    mat <- .get_seurat_counts(mat)
  }
  
  # Check if matrix is already log-normalized (max < 35 and has fractional values)
  sample_vals <- if (inherits(mat, "dgCMatrix")) head(mat@x, 100) else as.numeric(mat[seq_len(min(10, nrow(mat))), seq_len(min(10, ncol(mat)))])
  sample_vals <- sample_vals[!is.na(sample_vals)]
  is_already_lognorm <- length(sample_vals) > 0 && max(sample_vals) < 35 && any(sample_vals %% 1 != 0)
  if (is_already_lognorm) {
    return(mat)
  }
  
  libs <- as.numeric(Matrix::colSums(mat))
  libs[libs == 0] <- 1
  
  if (inherits(mat, "dgCMatrix")) {
    norm_mat <- mat
    cell_factors <- scale_factor / libs
    col_idx <- rep(seq_len(ncol(mat)), diff(mat@p))
    scaled_x <- norm_mat@x * cell_factors[col_idx] + pseudo_count
    norm_mat@x <- if (log_base == 2) log2(scaled_x) else log(scaled_x)
    return(norm_mat)
  } else {
    mat <- as.matrix(mat)
    scaled <- t(t(mat) / libs) * scale_factor + pseudo_count
    norm_mat <- if (log_base == 2) log2(scaled) else log(scaled)
    return(norm_mat)
  }
}

#' Extract Dataset Dimensional Properties and Summary Statistics
#'
#' Computes key dimensional and biological properties from single-cell or multiomics count matrices,
#' including cell counts, feature counts, sparsity percentage, biological group/cell type counts,
#' batch counts, median library size, median detected features per cell, and mean expression level.
#'
#' @param mat A count matrix (genes/features x cells), dgCMatrix, SingleCellExperiment, or Seurat object.
#' @param role Character string describing the role of the dataset (e.g. \code{"Biological Reference"} or \code{"Simulated"}). Default is \code{"Reference"}.
#' @param method_name Character string indicating the method or dataset name (e.g. \code{"Empirical Reference"}, \code{"Splatter"}). Default is \code{"Dataset"}.
#' @param modality Character string describing the molecular modality (e.g. \code{"scRNA-seq"}, \code{"scATAC-seq"}). Default is \code{"scRNA-seq"}.
#' @param cell_types Optional factor or character vector of cell type annotations.
#' @param batch_info Optional factor or character vector of batch annotations.
#'
#' @return A data.frame with 1 row summarizing the dataset properties.
#' @examples
#' data(example_scrna, package = "scSimEval")
#' extract_dataset_summary(example_scrna$ref, role = "Reference", method_name = "Empirical")
#' @export
extract_dataset_summary <- function(mat, role = "Reference", method_name = "Dataset",
                                    modality = "scRNA-seq", cell_types = NULL, batch_info = NULL) {
  if (is.null(mat)) return(NULL)
  
  if (inherits(mat, "SingleCellExperiment")) {
    if (is.null(cell_types) && "cell_type" %in% colnames(SummarizedExperiment::colData(mat))) {
      cell_types <- SummarizedExperiment::colData(mat)$cell_type
    }
    if (is.null(batch_info) && "batch" %in% colnames(SummarizedExperiment::colData(mat))) {
      batch_info <- SummarizedExperiment::colData(mat)$batch
    }
    mat <- SingleCellExperiment::counts(mat)
  } else if (inherits(mat, "Seurat")) {
    if (is.null(cell_types)) cell_types <- Seurat::Idents(mat)
    if (is.null(batch_info) && "batch" %in% colnames(mat@meta.data)) {
      batch_info <- mat@meta.data$batch
    }
    mat <- mat[["RNA"]]$counts
  }
  
  n_cells <- ncol(mat)
  n_features <- nrow(mat)
  
  if (inherits(mat, "dgCMatrix")) {
    n_zeros <- (as.numeric(n_cells) * as.numeric(n_features)) - length(mat@x)
    sparsity_pct <- round((n_zeros / (as.numeric(n_cells) * as.numeric(n_features))) * 100, 2)
    lib_sizes <- Matrix::colSums(mat)
    det_feats <- Matrix::colSums(mat > 0)
    mean_expr <- mean(mat@x) * (length(mat@x) / (as.numeric(n_cells) * as.numeric(n_features)))
  } else {
    mat_num <- as.matrix(mat)
    sparsity_pct <- round(mean(mat_num == 0) * 100, 2)
    lib_sizes <- colSums(mat_num)
    det_feats <- colSums(mat_num > 0)
    mean_expr <- mean(mat_num)
  }
  
  med_lib <- stats::median(lib_sizes)
  med_det <- stats::median(det_feats)
  
  ct_str <- if (!is.null(cell_types) && length(cell_types) == n_cells) {
    u_ct <- unique(stats::na.omit(as.character(cell_types)))
    formatC(length(u_ct), format = "d", big.mark = ",")
  } else {
    "Not specified"
  }
  
  b_str <- if (!is.null(batch_info) && length(batch_info) == n_cells) {
    u_b <- unique(stats::na.omit(as.character(batch_info)))
    formatC(length(u_b), format = "d", big.mark = ",")
  } else {
    "1"
  }
  
  data.frame(
    "Dataset / Simulator" = as.character(method_name),
    "Role" = as.character(role),
    "Modality" = as.character(modality),
    "Cells (N)" = formatC(as.integer(n_cells), format = "d", big.mark = ","),
    "Features (P)" = formatC(as.integer(n_features), format = "d", big.mark = ","),
    "Sparsity" = sprintf("%.2f%%", sparsity_pct),
    "Cell Types (Groups)" = ct_str,
    "Batches" = b_str,
    "Median Lib Size" = formatC(round(med_lib, 1), format = "f", digits = 1, big.mark = ","),
    "Median Detected Features" = formatC(round(med_det, 0), format = "d", big.mark = ","),
    "Mean Expression" = formatC(round(mean_expr, 3), format = "f", digits = 3),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

.SC_CELLTYPE_COLS <- c(
  "cell_type", "celltype", "cell_types", "celltypes", "cluster", "clusters",
  "group", "groups", "label", "labels", "cell_identity", "cell_line",
  "cellstate", "state", "annotation", "subclass", "cell_ontology_class"
)

.SC_BATCH_COLS <- c(
  "batch", "batch_info", "batches", "replicate", "rep", "donor", "sample",
  "dataset", "protocol", "run", "experiment", "plate", "library", "gem_group"
)

.seurat_default_assay <- function(obj) {
  if (requireNamespace("SeuratObject", quietly = TRUE)) {
    tryCatch(as.character(SeuratObject::DefaultAssay(obj)), error = function(e) names(obj@assays)[1])
  } else {
    names(obj@assays)[1]
  }
}

.seurat_assay_counts <- function(obj, assay_name = NULL) {
  if (is.null(assay_name) || !nzchar(assay_name)) assay_name <- .seurat_default_assay(obj)
  res <- tryCatch({
    assay_obj <- obj[[assay_name]]
    if (!is.null(assay_obj)) {
      if (inherits(assay_obj, "Assay5")) {
        try(SeuratObject::JoinLayers(assay_obj, layers = "counts"), silent = TRUE)
        assay_obj$counts
      } else if (inherits(assay_obj, "Assay")) {
        assay_obj@counts
      } else {
        Seurat::GetAssayData(obj, assay = assay_name, layer = "counts")
      }
    } else {
      NULL
    }
  }, error = function(e) {
    tryCatch(
      Seurat::GetAssayData(obj, assay = assay_name, layer = "counts"),
      error = function(e2) {
        tryCatch(
          Seurat::GetAssayData(obj, assay = assay_name, slot = "counts"),
          error = function(e3) NULL
        )
      }
    )
  })
  if (is.null(res) || nrow(res) == 0) {
    res <- tryCatch(Seurat::GetAssayData(obj, slot = "counts"), error = function(e) obj@assays[[1]]@counts)
  }
  res
}

#' Validate and Standardize Single-Cell Input Objects
#'
#' Strictly enforces standard input object compliance in \code{scSimEval} (Seurat / SingleCellExperiment).
#' Validates that input datasets are either \code{SingleCellExperiment} or
#' \code{Seurat} objects, extracting underlying count matrices, cell annotations,
#' batch labels, and dimensional reductions.
#'
#' @param obj An object to validate, strictly required to be a \code{SingleCellExperiment} or \code{Seurat} object.
#' @param arg_name Name of the argument being validated (for clear error messages).
#' @param reduction_name Optional name of precomputed dimensional reduction (e.g. \code{"PCA"}, \code{"integrated"}).
#'
#' @return A list with standardized components:
#'   \item{object}{The original \code{SingleCellExperiment} or \code{Seurat} object.}
#'   \item{counts}{Extracted raw count matrix (features x cells).}
#'   \item{cell_types}{Extracted or parsed cell-type factor/vector, or NULL if absent.}
#'   \item{batch_info}{Extracted or parsed batch factor/vector, or NULL if absent.}
#'   \item{embedding}{Extracted dimensional reduction matrix (cells x dims), or NULL if absent.}
#'
#' @export
validate_sc_object <- function(obj, arg_name = "data", reduction_name = NULL) {
  if (is.null(obj)) return(NULL)
  
  is_sce <- inherits(obj, "SingleCellExperiment")
  is_seurat <- inherits(obj, "Seurat")
  
  if (!is_sce && !is_seurat) {
    stop(sprintf(
      "scSimEval strictly requires a 'Seurat' or 'SingleCellExperiment' object for '%s' (received '%s').\nTo benchmark your simulation and evaluate all categories (including Category VI Multiomics Integration & Alignment reductions), please package your counts and metadata:\n  - Seurat: Seurat::CreateSeuratObject(counts = <matrix>, meta.data = <data.frame>)\n  - SingleCellExperiment: SingleCellExperiment::SingleCellExperiment(assays = list(counts = <matrix>), colData = <data.frame>)",
      arg_name, paste(class(obj), collapse = "/")
    ), call. = FALSE)
  }
  
  counts_mat <- NULL
  ct_vec <- NULL
  bt_vec <- NULL
  embed_mat <- NULL
  
  if (is_sce) {
    # 1. Counts
    if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      assays_list <- SummarizedExperiment::assays(obj)
      if ("counts" %in% names(assays_list)) {
        counts_mat <- assays_list[["counts"]]
      } else if (length(assays_list) > 0) {
        counts_mat <- assays_list[[1]]
      }
    }
    # 2. Metadata (cell_type and batch)
    cd <- SummarizedExperiment::colData(obj)
    cd_cols <- tolower(colnames(cd))
    ct_col <- which(cd_cols %in% tolower(.SC_CELLTYPE_COLS))[1]
    if (!is.na(ct_col)) ct_vec <- as.character(cd[[ct_col]])
    
    bt_col <- which(cd_cols %in% tolower(.SC_BATCH_COLS))[1]
    if (!is.na(bt_col)) bt_vec <- as.character(cd[[bt_col]])
    
    # 3. Dimensional Reductions
    rd_names <- SingleCellExperiment::reducedDimNames(obj)
    if (!is.null(rd_names) && length(rd_names) > 0) {
      target_rd <- if (!is.null(reduction_name) && reduction_name %in% rd_names) {
        reduction_name
      } else {
        pref <- c("integrated", "joint", "harmony", "wnn", "PCA", "pca")
        found <- intersect(pref, rd_names)
        if (length(found) > 0) found[1] else rd_names[1]
      }
      embed_mat <- as.matrix(SingleCellExperiment::reducedDim(obj, target_rd))
    }
  } else if (is_seurat) {
    # 1. Counts (respect DefaultAssay first so ATAC/secondary assays work cleanly;
    #    Seurat v5 split 'counts.*' layers are joined automatically)
    counts_mat <- tryCatch(
      .seurat_assay_counts(obj, .seurat_default_assay(obj)),
      error = function(e0) tryCatch(
      Seurat::GetAssayData(obj, slot = "counts"),
      error = function(e) {
        tryCatch(
          Seurat::GetAssayData(obj, layer = "counts"),
          error = function(e2) {
            tryCatch(
              obj[["RNA"]]$counts,
              error = function(e3) obj@assays[[1]]@counts
            )
          }
        )
      }
    ))
    # 2. Metadata
    md <- obj@meta.data
    md_cols <- tolower(colnames(md))
    ct_col <- which(md_cols %in% tolower(.SC_CELLTYPE_COLS))[1]
    if (!is.na(ct_col)) {
      ct_vec <- as.character(md[[ct_col]])
    } else {
      ident_vec <- tryCatch(as.character(Seurat::Idents(obj)), error = function(e) NULL)
      if (!is.null(ident_vec) && length(ident_vec) == ncol(counts_mat)) ct_vec <- ident_vec
    }
    
    bt_col <- which(md_cols %in% tolower(.SC_BATCH_COLS))[1]
    if (!is.na(bt_col)) bt_vec <- as.character(md[[bt_col]])
    
    # 3. Dimensional Reductions
    red_names <- names(obj@reductions)
    if (!is.null(red_names) && length(red_names) > 0) {
      target_red <- if (!is.null(reduction_name) && reduction_name %in% red_names) {
        reduction_name
      } else {
        pref <- c("integrated", "joint", "harmony", "wnn", "pca", "PCA")
        found <- intersect(pref, red_names)
        if (length(found) > 0) found[1] else red_names[1]
      }
      embed_mat <- as.matrix(Seurat::Embeddings(obj, reduction = target_red))
    }
  }
  
  list(
    object = obj,
    counts = counts_mat,
    cell_types = ct_vec,
    batch_info = bt_vec,
    embedding = embed_mat
  )
}



#' Validate and Standardize Multiomics Input Objects
#'
#' Strictly enforces standard multiomics input compliance in \code{scSimEval}.
#' Validates that multiomics inputs conform to \code{Seurat} or \code{SingleCellExperiment} structures.
#' Accepts a single multi-modal object or a list of two modality objects:
#' \itemize{
#'   \item A single multi-modal \code{Seurat} object containing both RNA and ATAC assays (with metadata and reductions).
#'   \item A single multi-modal \code{SingleCellExperiment} object with an \code{altExp} for the secondary modality.
#'   \item A named list of 2 \code{Seurat} or \code{SingleCellExperiment} objects (e.g. \code{list(rna = ..., atac = ...)}).
#'   \item A file path to an \code{.rds} file containing any of the above.
#' }
#'
#' @param obj Multiomics dataset object or file path to validate.
#' @param arg_name Name of argument for diagnostic error messages.
#' @param reduction_name Optional reduction name to extract.
#' @return A list with standardized components: \code{m1}, \code{m2}, \code{m1_obj}, \code{m2_obj}, \code{m1_name}, \code{m2_name}, \code{cell_types}, \code{batch_info}, \code{embedding}.
#' @export
validate_multiomics_input <- function(obj, arg_name = "multi_data", reduction_name = NULL) {
  if (is.null(obj)) return(NULL)

  # Auto-load if file path is provided
  if (is.character(obj) && length(obj) == 1 && file.exists(obj)) {
    obj <- readRDS(obj)
  }

  # Case A: Named list of 2 objects
  if (is.list(obj) && !inherits(obj, "Seurat") && !inherits(obj, "SingleCellExperiment")) {
    if (length(obj) == 1) {
      obj <- obj[[1]]
    } else if (length(obj) >= 2) {
      o1 <- obj[[1]]
      o2 <- obj[[2]]
      nm1 <- names(obj)[1] %||% "RNA"
      nm2 <- names(obj)[2] %||% "ATAC"
      v1 <- validate_sc_object(o1, paste0(arg_name, "[[1]]"), reduction_name = reduction_name)
      v2 <- validate_sc_object(o2, paste0(arg_name, "[[2]]"), reduction_name = reduction_name)

      ct <- if (!is.null(v1$cell_types)) v1$cell_types else v2$cell_types
      bt <- if (!is.null(v1$batch_info)) v1$batch_info else v2$batch_info
      emb <- if (!is.null(v1$embedding)) v1$embedding else v2$embedding

      return(list(
        m1 = v1$counts,
        m2 = v2$counts,
        m1_obj = v1$object,
        m2_obj = v2$object,
        m1_name = nm1,
        m2_name = nm2,
        cell_types = ct,
        batch_info = bt,
        embedding = emb
      ))
    }
  }

  # Case B: Multi-modal Seurat object
  if (inherits(obj, "Seurat")) {
    assay_names <- names(obj@assays)
    if (length(assay_names) >= 2) {
      is_atac <- grepl("atac|peak|chromatin", tolower(assay_names))
      is_rna  <- grepl("rna|gene|transcr", tolower(assay_names))

      if (any(is_rna) && any(is_atac)) {
        a1_name <- assay_names[is_rna][1]
        a2_name <- assay_names[is_atac][1]
      } else {
        a1_name <- assay_names[1]
        a2_name <- assay_names[2]
      }

      m1 <- .seurat_assay_counts(obj, a1_name)
      m2 <- .seurat_assay_counts(obj, a2_name)
      m1_obj <- obj
      m2_obj <- obj
      if (requireNamespace("SeuratObject", quietly = TRUE)) {
        try({ SeuratObject::DefaultAssay(m1_obj) <- a1_name }, silent = TRUE)
        try({ SeuratObject::DefaultAssay(m2_obj) <- a2_name }, silent = TRUE)
      }
      v <- validate_sc_object(obj, arg_name, reduction_name = reduction_name)
      return(list(
        m1 = m1,
        m2 = m2,
        m1_obj = m1_obj,
        m2_obj = m2_obj,
        m1_name = a1_name,
        m2_name = a2_name,
        cell_types = v$cell_types,
        batch_info = v$batch_info,
        embedding = v$embedding
      ))
    } else {
      stop(sprintf("Multiomics Seurat object '%s' must contain at least 2 assays (found only '%s': expected e.g. RNA and ATAC).", arg_name, assay_names[1]), call. = FALSE)
    }
  }

  # Case C: Multi-modal SingleCellExperiment with altExp
  if (inherits(obj, "SingleCellExperiment")) {
    alt_names <- SingleCellExperiment::altExpNames(obj)
    if (length(alt_names) >= 1) {
      m1 <- SummarizedExperiment::assay(obj, "counts")
      a2 <- SingleCellExperiment::altExp(obj, alt_names[1])
      m2 <- SummarizedExperiment::assay(a2, "counts")
      m1_nm <- SingleCellExperiment::mainExpName(obj)
      if (is.null(m1_nm) || !nzchar(m1_nm)) m1_nm <- "RNA"
      v <- validate_sc_object(obj, arg_name, reduction_name = reduction_name)
      return(list(
        m1 = m1,
        m2 = m2,
        m1_obj = obj,
        m2_obj = a2,
        m1_name = m1_nm,
        m2_name = alt_names[1],
        cell_types = v$cell_types,
        batch_info = v$batch_info,
        embedding = v$embedding
      ))
    } else {
      stop(sprintf("Multiomics SingleCellExperiment object '%s' must contain at least one altExp for the secondary modality (e.g. ATAC).", arg_name), call. = FALSE)
    }
  }

  stop(sprintf(
    "scSimEval strictly requires multiomics inputs for '%s' to be standardized 'Seurat' or 'SingleCellExperiment' objects containing multiomics assays, cell metadata, and reductions. Received '%s'.",
    arg_name, paste(class(obj), collapse = "/")
  ), call. = FALSE)
}
