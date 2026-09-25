#' @title Extraction of Single-Cell and Multiomics Summary Properties
#' @description Comprehensive feature- and cell-level property extraction combining
#'   simpipe, HelenaLC, and SimBench benchmarking methodologies.
#' @name extract_properties
NULL

#' Extract Comprehensive Cell-Level Properties
#'
#' Computes library size, log-library size, zero fraction / detection frequency,
#' pairwise cell correlation (on top HVGs), cell k-NN hubness, pairwise PCA distance,
#' and outlier proportions.
#'
#' @param data Count matrix (features x cells).
#' @param max_cells_cor Maximum number of cells for computing pairwise correlation (default 1000).
#' @param n_top_hvgs Number of highly variable genes to use for correlation and PCA (default 500).
#' @param verbose Logical, whether to print progress messages.
#'
#' @return A named list of cell-level summary vectors.
#' @export
extract_cell_properties <- function(data, max_cells_cor = 1000, n_top_hvgs = 500, verbose = FALSE) {
  if (inherits(data, "SingleCellExperiment")) {
    data <- SummarizedExperiment::assay(data, "counts")
  } else if (inherits(data, "Seurat")) {
    data <- Seurat::GetAssayData(data, slot = "counts")
  }
  
  if (verbose) message("Extracting cell library size & log-library size...")
  lib_size <- as.numeric(Matrix::colSums(data))
  log_lib_size <- log(lib_size + 1)
  
  if (verbose) message("Extracting cell zero fraction & detection rate...")
  if (inherits(data, "dgCMatrix")) {
    n_features <- nrow(data)
    non_zeros <- diff(data@p)
    zero_fraction_cell <- (n_features - non_zeros) / n_features
  } else {
    zero_fraction_cell <- as.numeric(colMeans(data == 0))
  }
  detection_freq_cell <- 1 - zero_fraction_cell
  
  # Normalize to log2(CPM+1) for correlation and PCA
  if (verbose) message("Extracting cell-cell correlations on top HVGs...")
  lib_adj <- lib_size
  lib_adj[lib_adj == 0] <- 1
  cpm <- t(t(as.matrix(data)) / lib_adj) * 1e6
  log_cpm <- log2(cpm + 1)
  
  gene_vars <- apply(log_cpm, 1, stats::var)
  n_hvg <- min(n_top_hvgs, length(gene_vars))
  top_hvg_idx <- order(gene_vars, decreasing = TRUE)[seq_len(n_hvg)]
  sub_cpm <- log_cpm[top_hvg_idx, , drop = FALSE]
  
  n_cells <- ncol(data)
  if (n_cells > max_cells_cor) {
    sub_cell_idx <- sample(seq_len(n_cells), max_cells_cor)
    cor_cpm <- sub_cpm[, sub_cell_idx, drop = FALSE]
  } else {
    cor_cpm <- sub_cpm
  }
  
  cell_cor_mat <- stats::cor(cor_cpm, method = "spearman")
  cell_cor <- cell_cor_mat[upper.tri(cell_cor_mat)]
  
  # PCA and pairwise cell distance (cell_pcd from HelenaLC)
  if (verbose) message("Extracting cell PCA distance and k-NN hubness...")
  pca_dist <- numeric(0)
  knn_hubness <- numeric(0)
  tryCatch({
    pca_res <- stats::prcomp(t(sub_cpm), rank. = min(20, ncol(sub_cpm) - 1, nrow(sub_cpm)))
    pca_pts <- pca_res$x
    
    # Subsample for distance if large
    if (n_cells > 200) {
      d_idx <- sample(seq_len(n_cells), 200)
      pca_dist <- as.numeric(stats::dist(pca_pts[d_idx, ]))
    } else {
      pca_dist <- as.numeric(stats::dist(pca_pts))
    }
    
    # k-NN hubness (count how often a cell is in the 5% nearest neighbors)
    k_val <- max(2, round(0.05 * n_cells))
    if (requireNamespace("RANN", quietly = TRUE)) {
      nn <- RANN::nn2(pca_pts, k = k_val + 1)
      idx_nn <- nn$nn.idx[, seq(2, k_val + 1)]
      knn_hubness <- vapply(seq_len(n_cells), function(i) sum(idx_nn == i), numeric(1))
    }
  }, error = function(e) NULL)
  
  # TMM normalization factor & effective library size (countsimQC & simpipe)
  tmm_factors <- rep(1, n_cells)
  if (requireNamespace("edgeR", quietly = TRUE)) {
    tryCatch({
      dge <- suppressMessages(suppressWarnings(edgeR::DGEList(counts = data)))
      dge <- suppressMessages(suppressWarnings(edgeR::calcNormFactors(dge, method = "TMM")))
      tmm_factors <- as.numeric(dge$samples$norm.factors)
    }, error = function(e) NULL)
  }
  effective_lib_size <- lib_size * tmm_factors
  
  prop_outliers_cell <- calc_outlier_proportion(lib_size)
  
  list(
    library_size = lib_size,
    log_library_size = log_lib_size,
    tmm_factor = tmm_factors,
    effective_library_size = effective_lib_size,
    zero_fraction_cell = zero_fraction_cell,
    detection_freq_cell = detection_freq_cell,
    cell_cor = cell_cor,
    pca_dist = pca_dist,
    knn_hubness = knn_hubness,
    prop_outliers_cell = prop_outliers_cell
  )
}

#' Extract Comprehensive Feature-Level (Gene / Peak) Properties
#'
#' Computes mean abundance, standard deviation, variance, coefficient of variation (CV),
#' feature dropout (zero fraction), dispersion, Biological Coefficient of Variation (BCV),
#' gene-gene correlations, and outlier proportions.
#'
#' @param data Count matrix (features x cells).
#' @param n_top_cor Number of top HVGs to use for gene-gene correlation matrix (default 400).
#' @param verbose Logical, whether to print progress messages.
#'
#' @return A named list of feature-level summary vectors.
#' @export
extract_feature_properties <- function(data, n_top_cor = 400, verbose = FALSE) {
  if (inherits(data, "SingleCellExperiment")) {
    data <- SummarizedExperiment::assay(data, "counts")
  } else if (inherits(data, "Seurat")) {
    data <- Seurat::GetAssayData(data, slot = "counts")
  }
  
  raw_mat <- as.matrix(data)
  
  if (verbose) message("Extracting feature zero fraction...")
  zero_fraction_feature <- as.numeric(rowMeans(raw_mat == 0))
  detection_freq_feature <- 1 - zero_fraction_feature
  
  # log2(CPM+1)
  lib_sizes <- colSums(raw_mat)
  lib_sizes[lib_sizes == 0] <- 1
  cpm <- t(t(raw_mat) / lib_sizes) * 1e6
  norm_mat <- log2(cpm + 1)
  
  if (verbose) message("Extracting mean, variance, and SD...")
  mean_val <- as.numeric(rowMeans(norm_mat))
  var_val <- apply(norm_mat, 1, stats::var)
  sd_val <- sqrt(var_val)
  
  raw_mean <- as.numeric(rowMeans(raw_mat))
  raw_sd <- apply(raw_mat, 1, stats::sd)
  cv_val <- ifelse(raw_mean > 0, (raw_sd / raw_mean) * 100, 0)
  dispersion_val <- ifelse(raw_mean > 0, (raw_sd^2) / raw_mean, 0)
  
  # Biological Coefficient of Variation (BCV) from countsimQC
  if (verbose) message("Extracting Biological Coefficient of Variation (BCV)...")
  bcv_val <- numeric(0)
  if (requireNamespace("edgeR", quietly = TRUE)) {
    tryCatch({
      dge <- suppressMessages(suppressWarnings(edgeR::DGEList(counts = raw_mat)))
      dge <- suppressMessages(suppressWarnings(edgeR::calcNormFactors(dge)))
      dge <- suppressMessages(suppressWarnings(edgeR::estimateDisp(dge)))
      bcv_val <- sqrt(dge$tagwise.dispersion)
    }, error = function(e) NULL)
  }
  if (length(bcv_val) == 0) {
    # Fallback empirical BCV approximation
    bcv_val <- sqrt(pmax(0, dispersion_val - 1) / pmax(1e-6, raw_mean))
  }
  
  # Gene-gene correlation on top variable genes (SimBench, countsimQC, & HelenaLC)
  if (verbose) message("Extracting gene-gene correlation...")
  n_genes_cor <- min(n_top_cor, length(var_val))
  top_g_idx <- order(var_val, decreasing = TRUE)[seq_len(n_genes_cor)]
  sub_g_mat <- norm_mat[top_g_idx, , drop = FALSE]
  
  g_cor_mat <- stats::cor(t(sub_g_mat), method = "spearman")
  gene_cor <- g_cor_mat[upper.tri(g_cor_mat)]
  
  prop_outliers_feature <- calc_outlier_proportion(as.numeric(rowSums(raw_mat)))
  
  list(
    mean_expression = mean_val,
    variance = var_val,
    sd = sd_val,
    cv = cv_val,
    zero_fraction_feature = zero_fraction_feature,
    detection_freq_feature = detection_freq_feature,
    dispersion = dispersion_val,
    bcv = bcv_val,
    gene_cor = gene_cor,
    prop_outliers_feature = prop_outliers_feature
  )
}

# -----------------------------------------------------------------------------
# Splatter, zingeR & ZINB-WaVE Dropout and Zero-Inflation Modeling
# -----------------------------------------------------------------------------

#' Fit Zero-Probability Dropout Curve (Splatter & ZINB-WaVE)
#'
#' Fits an empirical logistic dropout curve: logit(P(Y=0)) = beta_0 + beta_1 * log(mu)
#' to model the dropout relationship with mean expression.
#'
#' @param counts Count matrix (genes x cells) or SingleCellExperiment.
#' @return A list containing intercept, slope, midpoint (inflection point), and R-squared.
#' @export
calc_zero_probability_curve <- function(counts) {
  if (inherits(counts, "SingleCellExperiment")) {
    counts <- SummarizedExperiment::assay(counts, "counts")
  }
  raw_mat <- as.matrix(counts)
  gene_means <- rowMeans(raw_mat)
  p0 <- rowMeans(raw_mat == 0)
  
  # Filter genes with non-trivial means and dropouts strictly between 0 and 1
  valid <- gene_means > 0 & p0 > 0 & p0 < 1
  if (sum(valid) < 10) {
    return(list(intercept = NA_real_, slope = NA_real_, midpoint = NA_real_, r_squared = NA_real_))
  }
  
  log_mu <- log(gene_means[valid])
  p0_sub <- p0[valid]
  # Clamp p0 for numerical stability in logit
  eps <- 1e-4
  p0_clamped <- pmin(pmax(p0_sub, eps), 1 - eps)
  logit_p0 <- log(p0_clamped / (1 - p0_clamped))
  
  fit <- tryCatch(stats::lm(logit_p0 ~ log_mu), error = function(e) NULL)
  if (is.null(fit)) {
    return(list(intercept = NA_real_, slope = NA_real_, midpoint = NA_real_, r_squared = NA_real_))
  }
  
  coefs <- stats::coef(fit)
  b0 <- coefs[1]
  b1 <- coefs[2]
  midpoint <- if (!is.na(b1) && abs(b1) > 1e-6) -b0 / b1 else NA_real_
  r2 <- summary(fit)$r.squared
  
  list(
    intercept = as.numeric(b0),
    slope = as.numeric(b1),
    midpoint = as.numeric(midpoint),
    r_squared = as.numeric(r2)
  )
}

#' Evaluate Zero-Probability Dropout Curve Fidelity (Splatter & ZINB-WaVE)
#'
#' Compares empirical dropout trends between reference and simulated datasets.
#'
#' @param ref_counts Reference count matrix.
#' @param sim_counts Simulated count matrix.
#' @return Named numeric vector of curve parameter differences.
#' @export
evaluate_zero_probability_curve <- function(ref_counts, sim_counts) {
  c_ref <- calc_zero_probability_curve(ref_counts)
  c_sim <- calc_zero_probability_curve(sim_counts)
  
  c(
    slope_error = abs(c_ref$slope - c_sim$slope),
    midpoint_error = abs(c_ref$midpoint - c_sim$midpoint),
    r2_discrepancy = abs(c_ref$r_squared - c_sim$r_squared)
  )
}

#' Calculate Posterior Excess Zero Weights (zingeR & ZINB-WaVE)
#'
#' Calculates the posterior probability of zero counts being technical dropouts
#' (excess zeros) versus biological sampling zeros under a Negative Binomial model.
#'
#' @param counts Count matrix (genes x cells).
#' @return A list with mean excess zero weight, gene-level weights, and estimated zero-inflation rate.
#' @export
calc_excess_zero_weights <- function(counts) {
  if (inherits(counts, "SingleCellExperiment")) {
    counts <- SummarizedExperiment::assay(counts, "counts")
  }
  raw_mat <- as.matrix(counts)
  gene_means <- rowMeans(raw_mat)
  gene_vars <- apply(raw_mat, 1, stats::var)
  p0_obs <- rowMeans(raw_mat == 0)
  
  # Method of moments estimate of NB size / dispersion: theta = mu^2 / (var - mu)
  theta <- ifelse(gene_vars > gene_means,
                  (gene_means^2) / pmax(1e-4, gene_vars - gene_means),
                  100) # Poisson-like when variance <= mean
  
  # Theoretical NB zero probability: P_NB(0) = (1 + mu / theta)^(-theta)
  p0_nb <- (1 + gene_means / theta)^(-theta)
  
  # Zero-inflation parameter pi_i: excess zeros beyond NB expectation
  pi_excess <- pmax(0, pmin(1, (p0_obs - p0_nb) / pmax(1e-4, 1 - p0_nb)))
  
  # Posterior probability of technical zero given observed zero:
  # w_i = P(excess zero | Y = 0) = pi / (pi + (1 - pi) * P_NB(0))
  denom <- pi_excess + (1 - pi_excess) * p0_nb
  w_excess <- ifelse(p0_obs > 0 & denom > 0, pi_excess / denom, 0)
  
  list(
    mean_excess_zero_weight = mean(w_excess, na.rm = TRUE),
    median_excess_zero_weight = stats::median(w_excess, na.rm = TRUE),
    mean_zero_inflation = mean(pi_excess, na.rm = TRUE),
    gene_excess_weights = w_excess
  )
}

# -----------------------------------------------------------------------------
# SymSim: Kinetic Transcription & Intrinsic vs. Extrinsic Noise Decomposition
# -----------------------------------------------------------------------------

#' Kinetic Noise Decomposition (SymSim & Elowitz et al.)
#'
#' Decomposes gene expression variance / squared coefficient of variation (CV^2)
#' into intrinsic transcriptional bursting noise and extrinsic cell-state noise.
#'
#' @param counts Matrix of expression counts (genes x cells) or SingleCellExperiment.
#' @param cell_states Optional factor or vector of cell states / subpopulation clusters.
#'   If provided, noise is partitioned into within-state (intrinsic) and between-state (extrinsic) components.
#'   If NULL, intrinsic noise is estimated via Poisson shot-noise expectation (1 / mean).
#'
#' @return A list containing mean intrinsic noise, mean extrinsic noise, noise ratio, and gene-level vectors.
#' @export
calc_kinetic_noise_decomposition <- function(counts, cell_states = NULL) {
  if (inherits(counts, "SingleCellExperiment")) {
    counts <- SummarizedExperiment::assay(counts, "counts")
  }
  raw_mat <- as.matrix(counts)
  gene_means <- rowMeans(raw_mat)
  gene_vars <- apply(raw_mat, 1, stats::var)
  
  valid <- which(gene_means > 0.05)
  if (length(valid) < 5) {
    return(list(mean_intrinsic_noise = NA_real_, mean_extrinsic_noise = NA_real_, noise_ratio = NA_real_))
  }
  
  sub_mat <- raw_mat[valid, , drop = FALSE]
  m_sub <- gene_means[valid]
  v_sub <- gene_vars[valid]
  cv2_total <- v_sub / (m_sub^2)
  
  if (!is.null(cell_states) && length(unique(cell_states)) >= 2) {
    states <- as.factor(cell_states)
    n_cells <- ncol(sub_mat)
    state_weights <- table(states) / n_cells
    
    # Within-state variance: intrinsic noise
    within_vars <- matrix(0, nrow = nrow(sub_mat), ncol = nlevels(states))
    state_means <- matrix(0, nrow = nrow(sub_mat), ncol = nlevels(states))
    
    lvls <- levels(states)
    for (idx in seq_along(lvls)) {
      c_idx <- which(states == lvls[idx])
      if (length(c_idx) >= 2) {
        state_means[, idx] <- rowMeans(sub_mat[, c_idx, drop = FALSE])
        within_vars[, idx] <- apply(sub_mat[, c_idx, drop = FALSE], 1, stats::var)
      } else {
        state_means[, idx] <- sub_mat[, c_idx]
        within_vars[, idx] <- 0
      }
    }
    
    # Expected within-state variance
    exp_within_var <- as.numeric(within_vars %*% as.numeric(state_weights))
    cv2_intrinsic <- exp_within_var / (m_sub^2)
    cv2_extrinsic <- pmax(0, cv2_total - cv2_intrinsic)
  } else {
    # Theoretical Poisson shot noise baseline: Var_poisson = mean => CV^2_poisson = 1 / mean
    cv2_intrinsic <- pmin(cv2_total, 1 / m_sub)
    cv2_extrinsic <- pmax(0, cv2_total - cv2_intrinsic)
  }
  
  ratio <- ifelse(cv2_total > 0, cv2_intrinsic / cv2_total, 0)
  
  list(
    mean_intrinsic_noise = mean(cv2_intrinsic, na.rm = TRUE),
    mean_extrinsic_noise = mean(cv2_extrinsic, na.rm = TRUE),
    mean_total_cv2 = mean(cv2_total, na.rm = TRUE),
    mean_intrinsic_fraction = mean(ratio, na.rm = TRUE),
    gene_intrinsic_noise = cv2_intrinsic,
    gene_extrinsic_noise = cv2_extrinsic
  )
}

#' Fit Chromatin Accessibility-Sparsity Polynomial Curve (simATAC)
#'
#' Fits a polynomial curve relating peak/bin mean accessibility to non-zero cell
#' proportion (NZP / detection frequency) as modeled in the simATAC framework
#' (Navidi et al., Genome Biology 2021): NZP = c0 + c1 * mean + c2 * mean^2.
#'
#' @param data Count matrix (features/peaks x cells) or data.frame with 'mean' and 'nzp'.
#' @param poly_degree Degree of polynomial (default: 2 for quadratic curve).
#'
#' @return A list with estimated coefficients (c0, c1, c2), R-squared, Spearman/Pearson
#'   correlations, and model fit summary.
#' @export
calc_accessibility_sparsity_curve <- function(data, poly_degree = 2) {
  if (is.matrix(data) || inherits(data, "Matrix")) {
    means <- rowMeans(data)
    nzp <- rowMeans(data > 0)
  } else if (is.data.frame(data) && all(c("mean", "nzp") %in% colnames(data))) {
    means <- data$mean
    nzp <- data$nzp
  } else {
    stop("Input 'data' must be a matrix/Matrix or a data.frame with 'mean' and 'nzp' columns.")
  }
  
  valid <- which(!is.na(means) & !is.na(nzp) & is.finite(means) & is.finite(nzp))
  means <- means[valid]
  nzp <- nzp[valid]
  
  if (length(means) < poly_degree + 2) {
    return(list(
      c0 = NA_real_, c1 = NA_real_, c2 = NA_real_,
      r_squared = NA_real_,
      spearman_cor = NA_real_,
      pearson_cor = NA_real_
    ))
  }
  
  df <- data.frame(y = nzp, x = means)
  fit <- stats::lm(y ~ stats::poly(x, degree = poly_degree, raw = TRUE), data = df)
  s_fit <- summary(fit)
  
  coefs <- stats::coef(fit)
  c0 <- unname(coefs[1])
  c1 <- ifelse(length(coefs) >= 2, unname(coefs[2]), 0)
  c2 <- ifelse(length(coefs) >= 3, unname(coefs[3]), 0)
  
  r2 <- s_fit$r.squared
  spearman_cor <- stats::cor(means, nzp, method = "spearman")
  pearson_cor <- stats::cor(means, nzp, method = "pearson")
  
  list(
    c0 = c0,
    c1 = c1,
    c2 = c2,
    r_squared = r2,
    spearman_cor = spearman_cor,
    pearson_cor = pearson_cor,
    poly_degree = poly_degree,
    model = fit
  )
}

#' Evaluate Accessibility-Sparsity Curve Concordance (simATAC)
#'
#' Compares the non-linear relationship between peak mean accessibility and non-zero
#' proportion (NZP) between reference and simulated scATAC-seq datasets.
#'
#' @param ref_data Reference count matrix (peaks x cells).
#' @param sim_data Simulated count matrix (peaks x cells).
#' @param poly_degree Degree of polynomial (default: 2).
#'
#' @return A list of reference and simulated curve parameters, absolute discrepancies,
#'   and curve prediction RMSE.
#' @export
evaluate_accessibility_sparsity_curve <- function(ref_data, sim_data, poly_degree = 2) {
  ref_res <- calc_accessibility_sparsity_curve(ref_data, poly_degree = poly_degree)
  sim_res <- calc_accessibility_sparsity_curve(sim_data, poly_degree = poly_degree)
  
  delta_c0 <- abs(sim_res$c0 - ref_res$c0)
  delta_c1 <- abs(sim_res$c1 - ref_res$c1)
  delta_c2 <- abs(sim_res$c2 - ref_res$c2)
  delta_r2 <- abs(sim_res$r_squared - ref_res$r_squared)
  delta_spearman <- abs(sim_res$spearman_cor - ref_res$spearman_cor)
  
  # Predict on common grid of means
  if (is.matrix(ref_data) || inherits(ref_data, "Matrix")) {
    ref_means <- rowMeans(ref_data)
  } else {
    ref_means <- ref_data$mean
  }
  
  max_m <- max(ref_means, na.rm = TRUE)
  if (is.finite(max_m) && max_m > 0) {
    grid_x <- seq(0, max_m, length.out = 100)
    ref_pred <- ref_res$c0 + ref_res$c1 * grid_x + ref_res$c2 * (grid_x^2)
    sim_pred <- sim_res$c0 + sim_res$c1 * grid_x + sim_res$c2 * (grid_x^2)
    curve_rmse <- sqrt(mean((ref_pred - sim_pred)^2, na.rm = TRUE))
  } else {
    curve_rmse <- NA_real_
  }
  
  list(
    ref_c0 = ref_res$c0,
    ref_c1 = ref_res$c1,
    ref_c2 = ref_res$c2,
    ref_r_squared = ref_res$r_squared,
    sim_c0 = sim_res$c0,
    sim_c1 = sim_res$c1,
    sim_c2 = sim_res$c2,
    sim_r_squared = sim_res$r_squared,
    delta_c0 = delta_c0,
    delta_c1 = delta_c1,
    delta_c2 = delta_c2,
    delta_r_squared = delta_r2,
    delta_spearman = delta_spearman,
    curve_rmse = curve_rmse
  )
}


