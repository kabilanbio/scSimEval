#' Cell-Specific Mixing Score (CMS)
#'
#' Tests the hypothesis that group-specific distance distributions of k-nearest
#' neighbor cells have the same underlying unspecified distribution.
#' Computes cell-specific p-values via Anderson-Darling tests (if kSamples
#' installed) or two-sample Kolmogorov-Smirnov / Kruskal-Wallis tests, with
#' optional delegation to CellMixS::cms.
#'
#' @param data Count or log-normalized expression matrix (features x cells) or NULL if coords provided.
#' @param batch_info Factor or character vector of batch assignments for each cell.
#' @param k Neighborhood size k for k-NN graph. Default is min(table(batch_info))/2 or 30.
#' @param n_pcs Number of principal components to compute if coords is NULL. Default is 30.
#' @param coords Optional precomputed embedding matrix (cells x dimensions, e.g. PCA or UMAP).
#' @param cell_min Minimum number of cells per batch in neighborhood to perform test. Default is 4.
#' @param unbalanced Logical, whether to set pure neighborhoods to NA. Default is FALSE.
#'
#' @return A named list:
#'   \item{cms_scores}{Numeric vector of CMS p-values for each cell (higher means better mixed)}
#'   \item{mean_cms}{Mean CMS score across cells}
#'   \item{median_cms}{Median CMS score across cells}
#' @examples
#' coords <- matrix(stats::rnorm(100), 50, 2)
#' batch <- factor(rep(c("B1", "B2"), length.out = 50))
#' calc_cms(coords, batch)
#' @export
calc_cms <- function(
  data = NULL,
  batch_info,
  k = NULL,
  n_pcs = 30,
  coords = NULL,
  cell_min = 4,
  unbalanced = FALSE
) {
  batch_info <- as.factor(batch_info)
  
  if (is.null(coords)) {
    d_use <- normalize_counts(data)
    d_use <- if (inherits(d_use, "Matrix")) as.matrix(d_use) else d_use
    coords <- fast_pca(t(d_use), n_pcs = n_pcs, n_hvg = 2000, scale = FALSE)
  }
  
  n_cells <- nrow(coords)
  if (is.null(k)) {
    k <- max(10, round(min(table(batch_info)) / 2))
  }
  k <- min(k, n_cells - 1)
  
  # Check if CellMixS is available
  if (requireNamespace("CellMixS", quietly = TRUE) &&
      requireNamespace("SingleCellExperiment", quietly = TRUE) &&
      !is.null(data)) {
    res_cms <- tryCatch({
      sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = as.matrix(data)),
        colData = data.frame(batch = batch_info)
      )
      SingleCellExperiment::reducedDim(sce, "PCA") <- coords
      eval_res <- CellMixS::cms(sce, k = k, group = "batch", dim_red = "PCA", n_dim = ncol(coords))
      p_vals <- eval_res$cms
      return(list(
        cms_scores = p_vals,
        mean_cms = mean(p_vals, na.rm = TRUE),
        median_cms = stats::median(p_vals, na.rm = TRUE)
      ))
    }, error = function(e) NULL)
    if (!is.null(res_cms)) return(res_cms)
  }
  
  # Native implementation
  dist_mat <- as.matrix(stats::dist(coords))
  diag(dist_mat) <- Inf
  has_ksamples <- requireNamespace("kSamples", quietly = TRUE)
  
  p_vals <- vapply(seq_len(n_cells), function(i) {
    nn_idx <- order(dist_mat[i, ])[seq_len(k)]
    nn_dists <- dist_mat[i, nn_idx]
    nn_batches <- batch_info[nn_idx]
    
    dist_by_b <- split(nn_dists, nn_batches)
    dist_by_b <- dist_by_b[vapply(dist_by_b, length, integer(1)) >= cell_min]
    
    if (length(dist_by_b) <= 1) {
      return(if (unbalanced) NA_real_ else 0)
    }
    
    if (has_ksamples) {
      ad_res <- tryCatch(kSamples::ad.test(dist_by_b), error = function(e) NULL)
      if (!is.null(ad_res)) {
        return(mean(ad_res$ad[, " asympt. P-value"], na.rm = TRUE))
      }
    }
    
    if (length(dist_by_b) == 2) {
      ks_res <- tryCatch(stats::ks.test(dist_by_b[[1]], dist_by_b[[2]]), error = function(e) NULL)
      return(if (!is.null(ks_res)) ks_res$p.value else NA_real_)
    } else {
      flat_dists <- unlist(dist_by_b, use.names = FALSE)
      flat_groups <- factor(rep(names(dist_by_b), lengths(dist_by_b)))
      kw_res <- tryCatch(stats::kruskal.test(flat_dists ~ flat_groups), error = function(e) NULL)
      return(if (!is.null(kw_res)) kw_res$p.value else NA_real_)
    }
  }, numeric(1))
  
  list(
    cms_scores = p_vals,
    mean_cms = mean(p_vals, na.rm = TRUE),
    median_cms = stats::median(p_vals, na.rm = TRUE)
  )
}

#' Inverse Simpson Index (ISI) for Batch Mixing
#'
#' Evaluates the diversity of batch labels within each cell's k-nearest neighborhood
#' using the Inverse Simpson Index (1 / sum(p_b^2)).
#' Implements a pure-R, lightweight formulation of LISI with optional distance weighting.
#'
#' @param coords Matrix of cell coordinates (cells x dimensions).
#' @param batch_info Factor or character vector of batch assignments.
#' @param k Number of nearest neighbors. Default is 30.
#' @param weighted Logical, whether to weight neighbor contributions by inverse distance. Default is TRUE.
#'
#' @return A named list:
#'   \item{isi_scores}{Numeric vector of cell-level ISI values}
#'   \item{mean_isi}{Mean ISI score (ranges from 1 to number of batches; higher indicates better mixing)}
#'   \item{median_isi}{Median ISI score}
#' @examples
#' coords <- matrix(stats::rnorm(100), 50, 2)
#' batch <- factor(rep(c("B1", "B2"), length.out = 50))
#' calc_isi(coords, batch)
#' @export
calc_isi <- function(coords, batch_info, k = 30, weighted = TRUE) {
  batch_info <- as.factor(batch_info)
  coords <- as.matrix(coords)
  n_cells <- nrow(coords)
  k <- min(k, n_cells - 1)
  
  dist_mat <- as.matrix(stats::dist(coords))
  diag(dist_mat) <- Inf
  batch_levels <- levels(batch_info)
  n_batch <- length(batch_levels)
  
  isi_vals <- vapply(seq_len(n_cells), function(i) {
    nn_idx <- order(dist_mat[i, ])[seq_len(k)]
    nn_b <- batch_info[nn_idx]
    nn_d <- dist_mat[i, nn_idx]
    
    if (weighted) {
      weights <- 1 / (nn_d + 1)
      p_b <- vapply(batch_levels, function(b) {
        sum(weights[nn_b == b]) / sum(weights)
      }, numeric(1))
    } else {
      p_b <- vapply(batch_levels, function(b) {
        sum(nn_b == b) / k
      }, numeric(1))
    }
    
    simpson <- sum(p_b^2)
    if (simpson <= 0) as.numeric(n_batch) else 1 / simpson
  }, numeric(1))
  
  list(
    isi_scores = isi_vals,
    mean_isi = mean(isi_vals, na.rm = TRUE),
    median_isi = stats::median(isi_vals, na.rm = TRUE)
  )
}

#' Calculate Batch Integration Metrics
#'
#' Evaluates how realistically simulated technical batches mix or separate,
#' computing the five canonical Category IV (Technical Variation & Batch Mixing
#' Dynamics) metrics: Batch Silhouette, PC Regression R2, CMS (Cell-Specific
#' Mixing Score), Inverse Simpson Index (ISI), and cross-batch cell type
#' classification macro F1.
#'
#' @param data Count or log-normalized expression matrix (features x cells).
#' @param batch_info Factor or character vector of batch assignments for each cell.
#' @param k Neighborhood size k for k-NN graph and mixing metrics. Default is min(table(batch_info))/2.
#' @param n_pcs Number of principal components for embedding evaluation. Default is 30.
#' @param cell_types Optional factor or character vector of cell types for cross-batch evaluation.
#' @param verbose Logical, whether to print progress.
#'
#' @return A named list of batch metrics:
#'   \item{batch_silhouette}{Average silhouette width using batch labels (lower is better mixed)}
#'   \item{pcr_r2}{Total variance explained by batch in principal components (PCR)}
#'   \item{cms}{Mean Cell-specific mixing score (CellMixS)}
#'   \item{isi}{Mean Inverse Simpson Index for batch mixing (CellMixS)}
#'   \item{cross_batch_F1}{Mean cross-batch cell type transfer macro F1 (if cell_types supplied)}
#' @examples
#' data(example_scrna, package = "scSimEval")
#' evaluate_batch_metrics(example_scrna$ref, example_scrna$batch_info,
#'                        cell_types = example_scrna$cell_types)
#' @export
evaluate_batch_metrics <- function(
  data,
  batch_info,
  k = NULL,
  n_pcs = 30,
  cell_types = NULL,
  verbose = FALSE
) {
  batch_info <- as.factor(batch_info)
  n_cells <- ncol(data)
  
  if (is.null(k)) {
    k <- max(5, round(min(table(batch_info)) / 2))
  }
  
  if (verbose) message("Running PCA for batch evaluation...")
  data_norm <- normalize_counts(data)
  sub_mat <- as.matrix(data_norm)
  gene_vars <- apply(sub_mat, 1, stats::var)
  n_hvg <- min(2000, length(gene_vars))
  hvg_idx <- order(gene_vars, decreasing = TRUE)[seq_len(n_hvg)]
  
  pca_res <- stats::prcomp(t(sub_mat[hvg_idx, , drop = FALSE]),
                           center = TRUE, scale. = FALSE,
                           rank. = min(n_pcs, n_cells - 1))
  pca_coords <- pca_res$x
  
  # 1. Batch Silhouette Width (kBET::batch_sil or cluster::silhouette)
  if (verbose) message("Calculating batch silhouette width...")
  batch_sil <- NA_real_
  pkg_kbet <- "kBET"
  if (requireNamespace(pkg_kbet, quietly = TRUE)) {
    batch_sil <- tryCatch({
      kbet_batch_sil <- getExportedValue(pkg_kbet, "batch_sil")
      kbet_batch_sil(pca.data = pca_res, batch = batch_info, nPCs = min(n_pcs, ncol(pca_coords)))
    }, error = function(e) NA_real_)
  } else if (requireNamespace("cluster", quietly = TRUE)) {
    dist_pca <- stats::dist(pca_coords)
    sil_res <- tryCatch(cluster::silhouette(as.numeric(as.factor(batch_info)), dist_pca), error = function(e) NULL)
    if (!is.null(sil_res) && (is.matrix(sil_res) || inherits(sil_res, "silhouette"))) {
      batch_sil <- mean(sil_res[, 3], na.rm = TRUE)
    }
  }
  
  # 2. Principal Component Regression (PCR) R2
  if (verbose) message("Calculating principal component regression (PCR)...")
  pcr_r2 <- NA_real_
  if (requireNamespace(pkg_kbet, quietly = TRUE)) {
    pcr_res <- tryCatch({
      kbet_pcr <- getExportedValue(pkg_kbet, "pcRegression")
      kbet_pcr(pca.data = pca_res, batch = as.character(batch_info), n_top = min(n_pcs, ncol(pca_coords)))
    }, error = function(e) NULL)
    if (!is.null(pcr_res)) pcr_r2 <- pcr_res$R2Var
  } else {
    r2_vals <- vapply(seq_len(ncol(pca_coords)), function(pc_idx) {
      fit <- stats::lm(pca_coords[, pc_idx] ~ batch_info)
      summary(fit)$r.squared
    }, numeric(1))
    pc_vars <- pca_res$sdev[seq_len(ncol(pca_coords))]^2
    pcr_r2 <- sum(r2_vals * pc_vars) / sum(pc_vars)
  }
  
  # 3. Cell-Specific Mixing Score (CMS, CellMixS)
  if (verbose) message("Calculating CMS (CellMixS)...")
  cms_res <- calc_cms(data = sub_mat, batch_info = batch_info, k = k, coords = pca_coords)
  cms_score <- cms_res$mean_cms
  
  # 4. Inverse Simpson Index (ISI, CellMixS)
  if (verbose) message("Calculating Inverse Simpson Index (ISI)...")
  isi_res <- calc_isi(coords = pca_coords, batch_info = batch_info, k = k, weighted = TRUE)
  isi_score <- isi_res$mean_isi
  
  # 5. Cross-Batch Cell Type Classification Macro F1
  cross_batch_f1 <- NA_real_
  if (!is.null(cell_types)) {
    if (verbose) message("Evaluating cross-batch cell type classification...")
    cb_res <- evaluate_cross_batch_prediction(
      coords = pca_coords,
      cell_types = cell_types,
      batch_info = batch_info,
      k = min(5, k)
    )
    cross_batch_f1 <- cb_res$mean_cross_batch_F1
  }
  
  list(
    batch_silhouette = batch_sil,
    pcr_r2 = pcr_r2,
    cms = cms_score,
    isi = isi_score,
    cross_batch_F1 = cross_batch_f1
  )
}

#' Evaluate Cross-Batch Cell Classification Transfer F1
#'
#' Evaluates whether simulated or integrated multi-batch single-cell data preserve
#' cell-type identity across distinct experimental batches or donors.
#' Trains a k-NN classifier on each individual batch and tests the macro F1 score
#' of predictions on all other remaining batches.
#'
#' @param coords Embedding coordinates (cells x dimensions).
#' @param cell_types Factor or character vector of ground-truth cell type labels.
#' @param batch_info Factor or character vector of batch assignments.
#' @param k Integer number of nearest neighbors for the k-NN classifier (default 5).
#'
#' @return A list containing:
#'   \item{mean_cross_batch_F1}{Overall mean macro-averaged F1 score across all directed batch pairs}
#' @examples
#' coords <- matrix(stats::rnorm(100), 50, 2)
#' batch <- factor(rep(c("B1", "B2"), length.out = 50))
#' ct <- factor(rep(c("TypeA", "TypeB"), length.out = 50))
#' evaluate_cross_batch_prediction(coords, ct, batch)
#' @export
evaluate_cross_batch_prediction <- function(coords, cell_types, batch_info, k = 5) {
  cell_types <- as.factor(cell_types)
  batch_info <- as.factor(batch_info)
  coords <- as.matrix(coords)
  
  batches <- levels(batch_info)
  n_b <- length(batches)
  pair_f1 <- matrix(NA_real_, nrow = n_b, ncol = n_b)
  cls <- levels(cell_types)
  
  for (i in seq_along(batches)) {
    train_idx <- which(batch_info == batches[i])
    train_x <- coords[train_idx, , drop = FALSE]
    train_y <- cell_types[train_idx]
    
    if (length(unique(train_y)) < 2) next
    
    for (j in seq_along(batches)) {
      if (i == j) next
      test_idx <- which(batch_info == batches[j])
      test_x <- coords[test_idx, , drop = FALSE]
      test_y <- cell_types[test_idx]
      
      k_use <- min(k, length(train_idx) - 1)
      pred_y <- class::knn(train = train_x, test = test_x, cl = train_y, k = k_use)
      
      f1s <- vapply(cls, function(c) {
        tp <- sum(pred_y == c & test_y == c)
        fp <- sum(pred_y == c & test_y != c)
        fn <- sum(pred_y != c & test_y == c)
        if ((2 * tp + fp + fn) > 0) (2 * tp) / (2 * tp + fp + fn) else 0
      }, numeric(1))
      pair_f1[i, j] <- mean(f1s[cls %in% test_y], na.rm = TRUE)
    }
  }
  
  list(
    mean_cross_batch_F1 = mean(pair_f1, na.rm = TRUE)
  )
}
