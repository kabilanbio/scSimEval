#' @title Batch Effect & Integration Evaluation Metrics
#' @description Comprehensive suite of batch mixing, cell-specific mixing,
#'   and local structure preservation metrics integrated from simpipe,
#'   HelenaLC/simulation-comparison, and CellMixS (Lütge et al., Life Sci Alliance 2021).
#' @name metrics_batch
NULL

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
    if (is.null(data)) stop("Either 'data' or 'coords' must be provided.")
    sub_mat <- as.matrix(data)
    n_cells <- ncol(sub_mat)
    gene_vars <- apply(sub_mat, 1, stats::var)
    n_hvg <- min(2000, length(gene_vars))
    hvg_idx <- order(gene_vars, decreasing = TRUE)[seq_len(n_hvg)]
    pca_res <- stats::prcomp(t(sub_mat[hvg_idx, , drop = FALSE]),
                             center = TRUE, scale. = FALSE,
                             rank. = min(n_pcs, n_cells - 1))
    coords <- pca_res$x
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

#' Local Density Factor (LDF)
#'
#' Computes the Local Density Estimate (LDE) and Local Density Factor (LDF)
#' using a Gaussian kernel over reachability distances in the k-NN neighborhood,
#' adapted from Latecki et al. and CellMixS.
#'
#' @param coords Matrix of coordinates / embeddings (cells x dimensions).
#' @param k Number of nearest neighbors. Default is 15.
#' @param h Bandwidth parameter for Gaussian kernel. Default is 1.
#' @param c Scaling constant for comparison of LDE to neighboring observations. Default is 1.
#'
#' @return A named list:
#'   \item{lde}{Local density estimate for each cell}
#'   \item{ldf}{Local density factor for each cell}
#' @export
calc_ldf <- function(coords, k = 15, h = 1, c = 1) {
  coords <- as.matrix(coords)
  n_cells <- nrow(coords)
  k <- min(k, n_cells - 2)
  dim_sub <- ncol(coords)
  
  dist_mat <- as.matrix(stats::dist(coords))
  diag(dist_mat) <- Inf
  
  knn_idx <- t(apply(dist_mat, 1, function(row) order(row)[seq_len(k)]))
  knn_dist <- t(apply(dist_mat, 1, function(row) sort(row)[seq_len(k)]))
  k_dist <- knn_dist[, k]
  
  # Local Density Estimate (LDE)
  lde_vals <- vapply(seq_len(n_cells), function(i) {
    neighbors <- knn_idx[i, ]
    nbr_k_dist <- k_dist[neighbors]
    dist_to_nbrs <- dist_mat[i, neighbors]
    reach_dist <- pmax(nbr_k_dist, dist_to_nbrs)
    
    denom <- ((2 * pi)^(dim_sub / 2)) * ((h * nbr_k_dist)^dim_sub)
    denom[denom == 0] <- 1e-12
    kernel <- (1 / denom) * exp(- (reach_dist^2) / (2 * ((h * nbr_k_dist)^2) + 1e-12))
    sum(kernel) / k
  }, numeric(1))
  
  # Local Density Factor (LDF)
  ldf_vals <- vapply(seq_len(n_cells), function(i) {
    neighbors <- knn_idx[i, ]
    mean_nbr_lde <- mean(lde_vals[neighbors])
    mean_nbr_lde / (lde_vals[i] + c * mean_nbr_lde + 1e-12)
  }, numeric(1))
  
  list(lde = lde_vals, ldf = ldf_vals)
}

#' Local Density Differences (ldfDiff)
#'
#' Quantifies cell-specific changes in the Local Density Factor (LDF) before
#' and after data integration or between reference and simulated datasets.
#'
#' @param coords_pre Matrix of coordinates before integration / reference (cells x dimensions).
#' @param coords_post Matrix of coordinates after integration / simulated (cells x dimensions).
#' @param k Number of nearest neighbors. Default is 15.
#' @param h Bandwidth parameter for Gaussian kernel. Default is 1.
#' @param c Scaling constant. Default is 1.
#'
#' @return A named list:
#'   \item{diff}{Numeric vector of absolute LDF differences per cell}
#'   \item{mean_ldf_diff}{Mean absolute LDF difference (lower indicates better structure preservation)}
#'   \item{median_ldf_diff}{Median absolute LDF difference}
#'   \item{ldf_pre}{LDF in pre-integration / reference space}
#'   \item{ldf_post}{LDF in post-integration / simulated space}
#' @export
calc_ldf_diff <- function(coords_pre, coords_post, k = 15, h = 1, c = 1) {
  res_pre <- calc_ldf(coords_pre, k = k, h = h, c = c)
  res_post <- calc_ldf(coords_post, k = k, h = h, c = c)
  
  diff_vals <- abs(res_post$ldf - res_pre$ldf)
  
  list(
    diff = diff_vals,
    mean_ldf_diff = mean(diff_vals, na.rm = TRUE),
    median_ldf_diff = stats::median(diff_vals, na.rm = TRUE),
    ldf_pre = res_pre$ldf,
    ldf_post = res_post$ldf
  )
}

#' Seurat Mixing Metric
#'
#' Calculates the mixing metric originally proposed in Seurat (Stuart et al., Cell 2019)
#' and adapted by CellMixS. For each cell, finds the rank of the k_pos-th neighbor
#' from each batch in the sorted neighborhood, and computes the median rank across batches.
#'
#' @param coords Matrix of cell coordinates (cells x dimensions).
#' @param batch_info Factor or character vector of batch assignments.
#' @param k Maximum neighborhood size to search. Default is 300.
#' @param k_pos The rank position to extract per batch. Default is 5.
#'
#' @return A named list:
#'   \item{mixing_metrics}{Numeric vector of median ranks per cell}
#'   \item{mean_mixing_metric}{Mean mixing rank (lower indicates better mixing)}
#'   \item{median_mixing_metric}{Median mixing rank}
#' @export
calc_seurat_mixing_metric <- function(coords, batch_info, k = 300, k_pos = 5) {
  batch_info <- as.factor(batch_info)
  coords <- as.matrix(coords)
  n_cells <- nrow(coords)
  k <- min(k, n_cells)
  
  dist_mat <- as.matrix(stats::dist(coords))
  diag(dist_mat) <- 0 # include self per Seurat convention
  batch_levels <- levels(batch_info)
  
  med_ranks <- vapply(seq_len(n_cells), function(i) {
    ordered_idx <- order(dist_mat[i, ])[seq_len(k)]
    ordered_batches <- batch_info[ordered_idx]
    
    pos_per_batch <- vapply(batch_levels, function(b) {
      match_pos <- which(ordered_batches == b)
      if (length(match_pos) < k_pos) {
        k
      } else {
        match_pos[k_pos]
      }
    }, numeric(1))
    
    stats::median(pos_per_batch)
  }, numeric(1))
  
  list(
    mixing_metrics = med_ranks,
    mean_mixing_metric = mean(med_ranks, na.rm = TRUE),
    median_mixing_metric = stats::median(med_ranks, na.rm = TRUE)
  )
}

#' Local Structure Preservation Metric
#'
#' Calculates the proportion of overlapping k-nearest neighbors between
#' an original/reference embedding and an integrated/simulated embedding
#' (adapted from Seurat LocalStruct and CellMixS locStructure).
#'
#' @param coords_pre Matrix of coordinates before integration / reference (cells x dimensions).
#' @param coords_post Matrix of coordinates after integration / simulated (cells x dimensions).
#' @param k Number of nearest neighbors. Default is 30.
#'
#' @return A named list:
#'   \item{cell_overlaps}{Numeric vector of overlap fractions per cell}
#'   \item{mean_local_structure}{Mean overlap fraction across cells (higher indicates better preservation)}
#'   \item{median_local_structure}{Median overlap fraction}
#' @export
calc_local_structure_metric <- function(coords_pre, coords_post, k = 30) {
  coords_pre <- as.matrix(coords_pre)
  coords_post <- as.matrix(coords_post)
  n_cells <- nrow(coords_pre)
  k <- min(k, n_cells - 1)
  
  dist_pre <- as.matrix(stats::dist(coords_pre))
  dist_post <- as.matrix(stats::dist(coords_post))
  diag(dist_pre) <- Inf
  diag(dist_post) <- Inf
  
  overlaps <- vapply(seq_len(n_cells), function(i) {
    nn_pre <- order(dist_pre[i, ])[seq_len(k)]
    nn_post <- order(dist_post[i, ])[seq_len(k)]
    length(intersect(nn_pre, nn_post)) / k
  }, numeric(1))
  
  list(
    cell_overlaps = overlaps,
    mean_local_structure = mean(overlaps, na.rm = TRUE),
    median_local_structure = stats::median(overlaps, na.rm = TRUE)
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
#' computing kBET, LISI, Batch Silhouette, Shannon Entropy, PC Regression,
#' CMS (Cell-Specific Mixing Score), Inverse Simpson Index (ISI), Seurat Mixing Metric,
#' and optionally Local Density Differences (ldfDiff) and Local Structure Preservation.
#'
#' @param data Count or log-normalized expression matrix (features x cells).
#' @param batch_info Factor or character vector of batch assignments for each cell.
#' @param k Neighborhood size k for k-NN graph and mixing metrics. Default is min(table(batch_info))/2.
#' @param n_pcs Number of principal components for embedding evaluation. Default is 30.
#' @param pre_data Optional matrix of pre-integration / reference expression data (for ldfDiff and local structure).
#' @param cell_types Optional factor or character vector of cell types for cross-batch evaluation.
#' @param verbose Logical, whether to print progress.
#'
#' @return A named list of batch metrics:
#'   \item{batch_silhouette}{Average silhouette width using batch labels (lower is better mixed)}
#'   \item{shannon_entropy}{Entropy of batch frequencies in local neighborhoods}
#'   \item{pcr_r2}{Total variance explained by batch in principal components (PCR)}
#'   \item{cms}{Mean Cell-specific mixing score (CellMixS)}
#'   \item{isi}{Mean Inverse Simpson Index for batch mixing (CellMixS)}
#'   \item{seurat_mixing_metric}{Mean Seurat mixing metric (Stuart et al. / CellMixS)}
#'   \item{ldf_diff}{Mean Local Density Factor difference (if pre_data supplied)}
#'   \item{local_structure}{Mean Local structure preservation overlap (if pre_data supplied)}
#'   \item{kbet_rejection}{kBET rejection rate (if kBET installed)}
#'   \item{lisi_batch}{Average batch LISI (if lisi installed)}
#'   \item{cross_batch_accuracy}{Mean cross-batch cell type transfer accuracy (if cell_types supplied)}
#'   \item{cross_batch_F1}{Mean cross-batch cell type transfer macro F1 (if cell_types supplied)}
#' @export
evaluate_batch_metrics <- function(
  data,
  batch_info,
  k = NULL,
  n_pcs = 30,
  pre_data = NULL,
  cell_types = NULL,
  verbose = FALSE
) {
  batch_info <- as.factor(batch_info)
  n_cells <- ncol(data)
  
  if (is.null(k)) {
    k <- max(5, round(min(table(batch_info)) / 2))
  }
  
  if (verbose) message("Running PCA for batch evaluation...")
  sub_mat <- as.matrix(data)
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
  if (requireNamespace("kBET", quietly = TRUE)) {
    batch_sil <- tryCatch(
      kBET::batch_sil(pca.data = pca_res, batch = batch_info, nPCs = min(n_pcs, ncol(pca_coords))),
      error = function(e) NA_real_
    )
  } else if (requireNamespace("cluster", quietly = TRUE)) {
    dist_pca <- stats::dist(pca_coords)
    sil_res <- cluster::silhouette(as.numeric(batch_info), dist_pca)
    batch_sil <- mean(sil_res[, 3], na.rm = TRUE)
  }
  
  # 2. Shannon Entropy of Batch Proportions in k-NN neighborhoods
  if (verbose) message("Calculating local neighborhood Shannon entropy...")
  shannon_entropy <- NA_real_
  tryCatch({
    dist_pca <- as.matrix(stats::dist(pca_coords))
    diag(dist_pca) <- Inf
    entropies <- vapply(seq_len(n_cells), function(i) {
      nn_idx <- order(dist_pca[i, ])[seq_len(k)]
      b_freqs <- table(batch_info[nn_idx]) / k
      -sum(b_freqs * log2(b_freqs + 1e-12))
    }, numeric(1))
    shannon_entropy <- mean(entropies, na.rm = TRUE)
  }, error = function(e) NULL)
  
  # 3. Principal Component Regression (PCR) R2
  if (verbose) message("Calculating principal component regression (PCR)...")
  pcr_r2 <- NA_real_
  if (requireNamespace("kBET", quietly = TRUE)) {
    pcr_res <- tryCatch(
      kBET::pcRegression(pca.data = pca_res, batch = as.character(batch_info), n_top = min(n_pcs, ncol(pca_coords))),
      error = function(e) NULL
    )
    if (!is.null(pcr_res)) pcr_r2 <- pcr_res$R2Var
  } else {
    r2_vals <- vapply(seq_len(ncol(pca_coords)), function(pc_idx) {
      fit <- stats::lm(pca_coords[, pc_idx] ~ batch_info)
      summary(fit)$r.squared
    }, numeric(1))
    pc_vars <- pca_res$sdev[seq_len(ncol(pca_coords))]^2
    pcr_r2 <- sum(r2_vals * pc_vars) / sum(pc_vars)
  }
  
  # 4. Cell-Specific Mixing Score (CMS, CellMixS)
  if (verbose) message("Calculating CMS (CellMixS)...")
  cms_res <- calc_cms(data = sub_mat, batch_info = batch_info, k = k, coords = pca_coords)
  cms_score <- cms_res$mean_cms
  
  # 5. Inverse Simpson Index (ISI, CellMixS)
  if (verbose) message("Calculating Inverse Simpson Index (ISI)...")
  isi_res <- calc_isi(coords = pca_coords, batch_info = batch_info, k = k, weighted = TRUE)
  isi_score <- isi_res$mean_isi
  
  # 6. Seurat Mixing Metric (CellMixS)
  if (verbose) message("Calculating Seurat Mixing Metric...")
  mm_res <- calc_seurat_mixing_metric(coords = pca_coords, batch_info = batch_info, k = max(k, 50), k_pos = 5)
  mixing_metric <- mm_res$mean_mixing_metric
  
  # 7. Local Density Differences (ldfDiff) & Local Structure Preservation (if pre_data provided)
  ldf_diff_val <- NA_real_
  loc_struct_val <- NA_real_
  if (!is.null(pre_data)) {
    if (verbose) message("Calculating LDF difference and Local Structure Preservation...")
    pre_mat <- as.matrix(pre_data)
    pre_vars <- apply(pre_mat, 1, stats::var)
    pre_hvg <- order(pre_vars, decreasing = TRUE)[seq_len(min(2000, length(pre_vars)))]
    pre_pca <- stats::prcomp(t(pre_mat[pre_hvg, , drop = FALSE]),
                             center = TRUE, scale. = FALSE,
                             rank. = min(n_pcs, ncol(pre_mat) - 1))
    
    ldf_res <- tryCatch(calc_ldf_diff(pre_pca$x, pca_coords, k = k), error = function(e) NULL)
    if (!is.null(ldf_res)) ldf_diff_val <- ldf_res$mean_ldf_diff
    
    struct_res <- tryCatch(calc_local_structure_metric(pre_pca$x, pca_coords, k = k), error = function(e) NULL)
    if (!is.null(struct_res)) loc_struct_val <- struct_res$mean_local_structure
  }
  
  # 8. kBET Rejection Rate
  kbet_rejection <- NA_real_
  if (requireNamespace("kBET", quietly = TRUE)) {
    if (verbose) message("Calculating kBET rejection rate...")
    kbet_res <- tryCatch(
      kBET::kBET(df = t(sub_mat), batch = batch_info, k0 = k, plot = FALSE),
      error = function(e) NULL
    )
    if (!is.null(kbet_res)) {
      kbet_rejection <- mean(kbet_res$results$kBET.observed, na.rm = TRUE)
    }
  }
  
  # 9. LISI (Local Inverse Simpson Index)
  lisi_batch <- NA_real_
  if (requireNamespace("lisi", quietly = TRUE)) {
    if (verbose) message("Calculating LISI...")
    meta_df <- data.frame(batch = batch_info)
    lisi_res <- tryCatch(
      lisi::compute_lisi(pca_coords, meta_data = meta_df, label_colnames = "batch", perplexity = k),
      error = function(e) NULL
    )
    if (!is.null(lisi_res)) lisi_batch <- mean(lisi_res$batch, na.rm = TRUE)
  }
  
  # 10. Cross-Batch Label Transfer Accuracy
  cross_batch_acc <- NA_real_
  cross_batch_f1 <- NA_real_
  if (!is.null(cell_types)) {
    if (verbose) message("Evaluating cross-batch cell type classification...")
    cb_res <- evaluate_cross_batch_prediction(
      coords = pca_coords,
      cell_types = cell_types,
      batch_info = batch_info,
      k = min(5, k)
    )
    cross_batch_acc <- cb_res$mean_cross_batch_accuracy
    cross_batch_f1 <- cb_res$mean_cross_batch_F1
  }
  
  list(
    batch_silhouette = batch_sil,
    shannon_entropy = shannon_entropy,
    pcr_r2 = pcr_r2,
    cms = cms_score,
    isi = isi_score,
    seurat_mixing_metric = mixing_metric,
    ldf_diff = ldf_diff_val,
    local_structure = loc_struct_val,
    kbet_rejection = kbet_rejection,
    lisi_batch = lisi_batch,
    cross_batch_accuracy = cross_batch_acc,
    cross_batch_F1 = cross_batch_f1
  )
}

#' Evaluate Cross-Batch Cell Classification Transfer Accuracy
#'
#' Evaluates whether simulated or integrated multi-batch single-cell data preserve
#' cell-type identity across distinct experimental batches or donors.
#' Trains a k-NN classifier on each individual batch and tests prediction accuracy 
#' and macro F1 score on all other remaining batches.
#'
#' @param coords Embedding coordinates (cells x dimensions).
#' @param cell_types Factor or character vector of ground-truth cell type labels.
#' @param batch_info Factor or character vector of batch assignments.
#' @param k Integer number of nearest neighbors for the k-NN classifier (default 5).
#'
#' @return A list containing:
#'   \item{mean_cross_batch_accuracy}{Overall mean classification accuracy across all directed batch pairs}
#'   \item{mean_cross_batch_F1}{Overall mean macro-averaged F1 score across all directed batch pairs}
#' @export
evaluate_cross_batch_prediction <- function(coords, cell_types, batch_info, k = 5) {
  cell_types <- as.factor(cell_types)
  batch_info <- as.factor(batch_info)
  coords <- as.matrix(coords)
  
  batches <- levels(batch_info)
  n_b <- length(batches)
  pair_acc <- matrix(NA_real_, nrow = n_b, ncol = n_b)
  pair_f1  <- matrix(NA_real_, nrow = n_b, ncol = n_b)
  
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
      
      pair_acc[i, j] <- mean(pred_y == test_y)
      
      cls <- levels(cell_types)
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
    mean_cross_batch_accuracy = mean(pair_acc, na.rm = TRUE),
    mean_cross_batch_F1 = mean(pair_f1, na.rm = TRUE)
  )
}

# -----------------------------------------------------------------------------
# BASiCS, muscat, hierarchicell & rescueSim Hierarchical Multi-Level Variance
# -----------------------------------------------------------------------------

#' Variance Component Decomposition (BASiCS & muscat)
#'
#' Decomposes gene expression variance into biological cell-type heterogeneity,
#' donor / biological replicate variation, and residual technical noise.
#'
#' @param counts Matrix or data.frame (genes x cells) or SingleCellExperiment.
#' @param cell_metadata Data.frame containing donor and cell-type columns.
#' @param donor_col Character, column name for donor/replicate (default "donor").
#' @param celltype_col Character, column name for cell identity (default "cell_type").
#' @param n_genes Number of top variable genes to analyze (default 200).
#'
#' @return A list of mean percentage of variance explained by cell-type, donor, and residual noise.
#' @export
calc_variance_decomposition <- function(
  counts,
  cell_metadata,
  donor_col = "donor",
  celltype_col = "cell_type",
  n_genes = 200
) {
  if (inherits(counts, "SingleCellExperiment")) {
    counts <- SummarizedExperiment::assay(counts, "counts")
  }
  counts <- as.matrix(counts)
  
  if (!all(c(donor_col, celltype_col) %in% colnames(cell_metadata))) {
    stop(paste("cell_metadata must contain columns:", donor_col, "and", celltype_col))
  }
  
  # Select top variable genes
  gene_vars <- apply(counts, 1, stats::var)
  n_use <- min(n_genes, length(gene_vars))
  top_genes <- order(gene_vars, decreasing = TRUE)[seq_len(n_use)]
  sub_counts <- counts[top_genes, , drop = FALSE]
  
  # Normalization: log2(CPM+1)
  libs <- colSums(sub_counts)
  libs[libs == 0] <- 1
  cpm <- t(t(sub_counts) / libs) * 1e6
  log_cpm <- log2(cpm + 1)
  
  donors <- as.factor(cell_metadata[[donor_col]])
  celltypes <- as.factor(cell_metadata[[celltype_col]])
  
  var_celltype <- numeric(n_use)
  var_donor <- numeric(n_use)
  var_residual <- numeric(n_use)
  
  for (i in seq_len(n_use)) {
    y <- log_cpm[i, ]
    fit <- tryCatch(
      stats::aov(y ~ celltypes + donors),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      ss <- summary(fit)[[1]][, "Sum Sq"]
      total_ss <- sum(ss)
      if (total_ss > 0 && length(ss) >= 3) {
        var_celltype[i] <- ss[1] / total_ss
        var_donor[i] <- ss[2] / total_ss
        var_residual[i] <- ss[3] / total_ss
      } else {
        var_residual[i] <- 1
      }
    } else {
      var_residual[i] <- 1
    }
  }
  
  list(
    mean_var_celltype = mean(var_celltype, na.rm = TRUE),
    mean_var_donor = mean(var_donor, na.rm = TRUE),
    mean_var_residual = mean(var_residual, na.rm = TRUE),
    gene_var_celltype = var_celltype,
    gene_var_donor = var_donor,
    gene_var_residual = var_residual
  )
}

#' Calculate Intraclass Correlation Coefficient (ICC) (hierarchicell & rescueSim)
#'
#' Computes the one-way random effects Intraclass Correlation Coefficient (ICC)
#' per gene across donors or biological replicates to quantify hierarchical clustering.
#'
#' @param counts Matrix (genes x cells) or SingleCellExperiment.
#' @param donor_labels Vector or factor of donor / subject IDs.
#' @param n_genes Number of top variable genes to evaluate (default 200).
#'
#' @return A list with median ICC, mean ICC, and vector of per-gene ICC values.
#' @export
calc_intraclass_correlation <- function(counts, donor_labels, n_genes = 200) {
  if (inherits(counts, "SingleCellExperiment")) {
    counts <- SummarizedExperiment::assay(counts, "counts")
  }
  counts <- as.matrix(counts)
  donor_labels <- as.factor(donor_labels)
  
  if (nlevels(donor_labels) < 2) {
    return(list(median_icc = NA_real_, mean_icc = NA_real_, gene_icc = numeric(0)))
  }
  
  # Select top variable genes
  gene_vars <- apply(counts, 1, stats::var)
  n_use <- min(n_genes, length(gene_vars))
  top_genes <- order(gene_vars, decreasing = TRUE)[seq_len(n_use)]
  sub_counts <- counts[top_genes, , drop = FALSE]
  
  # log2(CPM+1)
  libs <- colSums(sub_counts)
  libs[libs == 0] <- 1
  cpm <- t(t(sub_counts) / libs) * 1e6
  log_cpm <- log2(cpm + 1)
  
  k <- nlevels(donor_labels)
  n_total <- ncol(counts)
  n_per_group <- table(donor_labels)
  k0 <- (1 / (k - 1)) * (n_total - sum(n_per_group^2) / n_total)
  
  icc_vals <- numeric(n_use)
  for (i in seq_len(n_use)) {
    y <- log_cpm[i, ]
    fit <- tryCatch(stats::aov(y ~ donor_labels), error = function(e) NULL)
    if (!is.null(fit)) {
      ms <- summary(fit)[[1]][, "Mean Sq"]
      msb <- ms[1]
      msw <- ms[2]
      sigma2_b <- max(0, (msb - msw) / k0)
      sigma2_w <- msw
      denom <- sigma2_b + sigma2_w
      icc_vals[i] <- if (denom > 0) sigma2_b / denom else 0
    } else {
      icc_vals[i] <- 0
    }
  }
  
  list(
    median_icc = stats::median(icc_vals, na.rm = TRUE),
    mean_icc = mean(icc_vals, na.rm = TRUE),
    gene_icc = icc_vals
  )
}

