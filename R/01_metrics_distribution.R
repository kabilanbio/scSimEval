# -----------------------------------------------------------------------------
# 1D Univariate Distance & Divergence Metrics
# -----------------------------------------------------------------------------

#' Calculate Median Absolute Deviation (MAD)
#' @param ref Numeric vector of reference distribution.
#' @param sim Numeric vector of simulated distribution.
#' @param align Logical, whether to sort and quantile-align vectors. Default TRUE.
#' @return Numeric MAD value.
#' @export
calc_mad <- function(ref, sim, align = TRUE) {
  if (align) {
    aligned <- align_distributions(ref, sim)
    ref <- aligned$ref
    sim <- aligned$sim
  }
  stats::median(abs(ref - sim), na.rm = TRUE)
}

#' Calculate Kolmogorov-Smirnov Distance (1D KS)
#' @param ref Numeric vector of reference values.
#' @param sim Numeric vector of simulated values.
#' @return Maximum vertical distance between empirical CDFs (between 0 and 1).
#' @export
calc_ks <- function(ref, sim) {
  ref <- stats::na.omit(as.numeric(ref))
  sim <- stats::na.omit(as.numeric(sim))
  if (length(ref) == 0 || length(sim) == 0) return(NA_real_)
  
  if (requireNamespace("provenance", quietly = TRUE)) {
    as.numeric(provenance::KS.diss(ref, sim))
  } else {
    suppressWarnings(as.numeric(stats::ks.test(ref, sim)$statistic))
  }
}

#' Calculate Mean Absolute Error (MAE)
#' @param ref Numeric vector of reference values.
#' @param sim Numeric vector of simulated values.
#' @param align Logical, whether to sort and quantile-align vectors. Default TRUE.
#' @return Mean absolute difference.
#' @export
calc_mae <- function(ref, sim, align = TRUE) {
  if (align) {
    aligned <- align_distributions(ref, sim)
    ref <- aligned$ref
    sim <- aligned$sim
  }
  if (requireNamespace("MLmetrics", quietly = TRUE)) {
    MLmetrics::MAE(sim, ref)
  } else {
    mean(abs(sim - ref), na.rm = TRUE)
  }
}

#' Calculate Root Mean Squared Error (RMSE)
#' @param ref Numeric vector of reference values.
#' @param sim Numeric vector of simulated values.
#' @param align Logical, whether to sort and quantile-align vectors. Default TRUE.
#' @return Quadratic error penalizing large discrepancies.
#' @export
calc_rmse <- function(ref, sim, align = TRUE) {
  if (align) {
    aligned <- align_distributions(ref, sim)
    ref <- aligned$ref
    sim <- aligned$sim
  }
  if (requireNamespace("MLmetrics", quietly = TRUE)) {
    MLmetrics::RMSE(sim, ref)
  } else {
    sqrt(mean((sim - ref)^2, na.rm = TRUE))
  }
}

#' Calculate Distribution Overlapping Index (OV)
#' @param ref Numeric vector of reference values.
#' @param sim Numeric vector of simulated values.
#' @return Area of overlap under probability densities (0 to 1).
#' @export
calc_overlap <- function(ref, sim) {
  ref <- stats::na.omit(as.numeric(ref))
  sim <- stats::na.omit(as.numeric(sim))
  if (length(ref) < 2 || length(sim) < 2) return(NA_real_)
  
  if (requireNamespace("overlapping", quietly = TRUE)) {
    res <- tryCatch(
      as.numeric(overlapping::overlap(list(x = ref, y = sim))[["OV"]]),
      error = function(e) NA_real_
    )
    return(res)
  } else {
    from <- min(c(ref, sim))
    to <- max(c(ref, sim))
    d_ref <- stats::density(ref, from = from, to = to, n = 512)
    d_sim <- stats::density(sim, from = from, to = to, n = 512)
    dx <- d_ref$x[2] - d_ref$x[1]
    sum(pmin(d_ref$y, d_sim$y)) * dx
  }
}

#' Calculate Bhattacharyya Distance (BH)
#' @param ref Numeric vector of reference values.
#' @param sim Numeric vector of simulated values.
#' @param align Logical, whether to align lengths. Default TRUE.
#' @return Statistical divergence between discrete probability measures.
#' @export
calc_bhattacharyya <- function(ref, sim, align = TRUE) {
  if (align) {
    aligned <- align_distributions(ref, sim)
    ref <- aligned$ref
    sim <- aligned$sim
  }
  ref <- stats::na.omit(as.numeric(ref))
  sim <- stats::na.omit(as.numeric(sim))
  if (length(ref) == 0 || length(sim) == 0) return(NA_real_)
  
  min_val <- min(c(ref, sim), na.rm = TRUE)
  if (is.finite(min_val) && min_val <= 0) {
    shift <- abs(min_val) + 1e-6
    ref <- ref + shift
    sim <- sim + shift
  }
  s_ref <- sum(ref, na.rm = TRUE)
  s_sim <- sum(sim, na.rm = TRUE)
  if (s_ref <= 0 || s_sim <= 0) return(NA_real_)
  
  p_ref <- ref / s_ref
  p_sim <- sim / s_sim
  
  if (requireNamespace("philentropy", quietly = TRUE)) {
    tryCatch(
      as.numeric(philentropy::distance(rbind(p_ref, p_sim), method = "bhattacharyya")),
      error = function(e) {
        bc <- sum(sqrt(p_ref * p_sim), na.rm = TRUE)
        if (is.na(bc) || bc <= 0) Inf else -log(bc)
      }
    )
  } else {
    bc <- sum(sqrt(p_ref * p_sim), na.rm = TRUE)
    if (is.na(bc) || bc <= 0) return(Inf)
    -log(bc)
  }
}

#' Calculate 1D Wasserstein Metric / Earth Mover's Distance (WS)
#' Integrated from HelenaLC/simulation-comparison.
#'
#' @param ref Numeric vector of reference values.
#' @param sim Numeric vector of simulated values.
#' @param p Power of the Wasserstein metric (default 1 for standard earth mover's distance).
#' @return The 1D Wasserstein distance between the two empirical distributions.
#' @export
calc_wasserstein_1d <- function(ref, sim, p = 1) {
  ref <- stats::na.omit(as.numeric(ref))
  sim <- stats::na.omit(as.numeric(sim))
  if (length(ref) == 0 || length(sim) == 0) return(NA_real_)
  
  if (requireNamespace("waddR", quietly = TRUE)) {
    tryCatch(waddR::wasserstein_metric(ref, sim, p = p), error = function(e) NA_real_)
  } else {
    # Exact 1D Wasserstein-1 via integrated quantile difference
    aligned <- align_distributions(ref, sim, method = "quantile", n_points = 1000)
    mean(abs(aligned$ref - aligned$sim)^p)^(1/p)
  }
}

# -----------------------------------------------------------------------------
# 2D Bivariate Joint Relationship Metrics
# -----------------------------------------------------------------------------

#' Calculate Fasano-Franceschini 2D Kolmogorov-Smirnov Test Statistic
#' Integrated from simpipe.
#'
#' @param ref_mat 2-column numeric matrix for reference.
#' @param sim_mat 2-column numeric matrix for simulation.
#' @param threads CPU threads. Default 1.
#' @return Estimated 2D KS statistic.
#' @export
calc_fasano_franceschini <- function(ref_mat, sim_mat, threads = 1) {
  ref_mat <- as.matrix(stats::na.omit(ref_mat))
  sim_mat <- as.matrix(stats::na.omit(sim_mat))
  if (ncol(ref_mat) != 2 || ncol(sim_mat) != 2) stop("Requires 2-column matrices.")
  
  if (requireNamespace("fasano.franceschini.test", quietly = TRUE)) {
    res <- tryCatch(
      fasano.franceschini.test::fasano.franceschini.test(ref_mat, sim_mat, threads = threads),
      error = function(e) NULL
    )
    if (!is.null(res)) return(mean(res$estimate))
  }
  return(NA_real_)
}

#' Calculate Peacock 2D Kolmogorov-Smirnov Test Statistic
#' Integrated from HelenaLC/simulation-comparison.
#'
#' @param ref_mat 2-column numeric matrix for reference.
#' @param sim_mat 2-column numeric matrix for simulation.
#' @return Peacock test statistic.
#' @export
calc_peacock_2d <- function(ref_mat, sim_mat) {
  ref_mat <- as.matrix(stats::na.omit(ref_mat))
  sim_mat <- as.matrix(stats::na.omit(sim_mat))
  if (ncol(ref_mat) != 2 || ncol(sim_mat) != 2) stop("Requires 2-column matrices.")
  
  if (requireNamespace("Peacock.test", quietly = TRUE)) {
    res <- tryCatch(
      as.numeric(Peacock.test::peacock2(ref_mat, sim_mat)),
      error = function(e) NA_real_
    )
    return(res)
  }
  return(NA_real_)
}

#' Calculate 2D Bivariate Kernel Density Estimation (KDE) Test Statistic
#' Integrated from simpipe & SimBench.
#'
#' @param ref_mat 2-column numeric matrix for reference.
#' @param sim_mat 2-column numeric matrix for simulation.
#' @return z-statistic testing discrepancy between bivariate kernel density estimates.
#' @export
calc_kde_test <- function(ref_mat, sim_mat) {
  ref_mat <- as.matrix(stats::na.omit(ref_mat))
  sim_mat <- as.matrix(stats::na.omit(sim_mat))
  if (ncol(ref_mat) != 2 || ncol(sim_mat) != 2) stop("Requires 2-column matrices.")
  
  if (requireNamespace("ks", quietly = TRUE)) {
    res <- tryCatch(
      as.numeric(ks::kde.test(ref_mat, sim_mat)$zstat),
      error = function(e) NA_real_
    )
    return(res)
  }
  return(NA_real_)
}

#' Calculate 2D Earth Mover's Distance (2D EMD)
#' Integrated from HelenaLC/simulation-comparison.
#'
#' Computes 2D bivariate density via MASS::kde2d over a shared bounding box,
#' then measures optimal transport distance via emdist::emd2d.
#'
#' @param ref_mat 2-column numeric matrix for reference.
#' @param sim_mat 2-column numeric matrix for simulation.
#' @param n Grid resolution for 2D density estimation. Default is 25.
#' @return Normalized 2D Earth Mover's Distance.
#' @export
calc_emd_2d <- function(ref_mat, sim_mat, n = 25) {
  ref_mat <- as.matrix(stats::na.omit(ref_mat))
  sim_mat <- as.matrix(stats::na.omit(sim_mat))
  if (ncol(ref_mat) != 2 || ncol(sim_mat) != 2) stop("Requires 2-column matrices.")
  
  if (requireNamespace("MASS", quietly = TRUE) && requireNamespace("emdist", quietly = TRUE)) {
    tryCatch({
      lims <- c(
        range(c(ref_mat[, 1], sim_mat[, 1])),
        range(c(ref_mat[, 2], sim_mat[, 2]))
      )
      k_ref <- MASS::kde2d(ref_mat[, 1], ref_mat[, 2], n = n, lims = lims)
      k_sim <- MASS::kde2d(sim_mat[, 1], sim_mat[, 2], n = n, lims = lims)
      emdist::emd2d(k_ref$z, k_sim$z) / n
    }, error = function(e) NA_real_)
  } else {
    NA_real_
  }
}

# -----------------------------------------------------------------------------
# Metrics from countsimQC (Soneson & Robinson, 2018)
# -----------------------------------------------------------------------------

#' Calculate Area Between Empirical Cumulative Distribution Functions (eCDFs)
#'
#' Integrated from countsimQC (Soneson & Robinson). Measures the normalized area
#' between the eCDFs of reference and simulation across the shared support.
#'
#' @param ref Numeric vector of reference values.
#' @param sim Numeric vector of simulated values.
#' @return Normalized area between the two eCDFs.
#' @export
calc_ecdf_diffarea <- function(ref, sim) {
  ref <- stats::na.omit(as.numeric(ref))
  sim <- stats::na.omit(as.numeric(sim))
  if (length(ref) < 2 || length(sim) < 2) return(NA_real_)
  
  xmin <- min(c(ref, sim))
  xmax <- max(c(ref, sim))
  if (xmax == xmin) return(0)
  
  e1 <- stats::ecdf(ref)
  e2 <- stats::ecdf(sim)
  xv <- sort(unique(c(ref, sim)))
  ediff <- abs(e1(xv) - e2(xv))
  
  # Normalized trapezoidal integration matching caTools::trapz
  x_norm <- (xv - xmin) / (xmax - xmin)
  dx <- diff(x_norm)
  sum(dx * (ediff[-1] + ediff[-length(ediff)]) / 2)
}

#' Calculate Wald-Wolfowitz Runs Test on Pooled Distributions
#'
#' Integrated from countsimQC (Soneson & Robinson). Tests whether values from
#' the two datasets intermingle randomly when sorted. A significant left-sided
#' result indicates that identical values cluster together into fewer runs than
#' expected by chance, signaling distinct distributions.
#'
#' @param ref Numeric vector of reference values.
#' @param sim Numeric vector of simulated values.
#' @param alternative Alternative hypothesis: "left.sided" (default in countsimQC) or "two.sided".
#' @return A named list with the runs test statistic and p-value.
#' @export
calc_runs_test <- function(ref, sim, alternative = c("left.sided", "two.sided")) {
  alternative <- match.arg(alternative)
  ref <- stats::na.omit(as.numeric(ref))
  sim <- stats::na.omit(as.numeric(sim))
  if (length(ref) < 2 || length(sim) < 2) {
    return(list(runs_statistic = NA_real_, runs_pvalue = NA_real_))
  }
  
  df <- data.frame(
    val = c(ref, sim),
    ds = c(rep(1, length(ref)), rep(0, length(sim)))
  )
  df <- df[order(df$val), ]
  
  if (requireNamespace("randtests", quietly = TRUE)) {
    res <- tryCatch(
      randtests::runs.test(df$ds, threshold = 0.5, alternative = alternative, plot = FALSE),
      error = function(e) NULL
    )
    if (!is.null(res)) {
      return(list(runs_statistic = as.numeric(res$statistic), runs_pvalue = as.numeric(res$p.value)))
    }
  }
  
  # Base R exact Wald-Wolfowitz normal approximation
  r <- rle(df$ds)
  n_runs <- length(r$lengths)
  n1 <- sum(df$ds == 1)
  n0 <- sum(df$ds == 0)
  mu_r <- 1 + (2 * n1 * n0) / (n1 + n0)
  var_r <- (2 * n1 * n0 * (2 * n1 * n0 - n1 - n0)) / ((n1 + n0)^2 * (n1 + n0 - 1))
  
  if (is.na(var_r) || var_r <= 0) {
    return(list(runs_statistic = 0, runs_pvalue = 1))
  }
  
  z <- (n_runs - mu_r) / sqrt(var_r)
  p_val <- if (alternative == "left.sided") {
    stats::pnorm(z)
  } else {
    2 * stats::pnorm(-abs(z))
  }
  
  list(runs_statistic = as.numeric(z), runs_pvalue = as.numeric(p_val))
}

#' Calculate Nearest-Neighbor Label Mismatch Proportion
#'
#' Integrated from countsimQC (Soneson & Robinson). Evaluates whether the
#' dataset label composition in k-NN neighborhoods departs significantly from
#' the overall global dataset proportion using Chi-squared tests.
#' Works for both 1D numeric vectors and 2D matrices.
#'
#' @param ref Vector or matrix for reference dataset.
#' @param sim Vector or matrix for simulated dataset.
#' @param k Number of nearest neighbors (default max(5, 0.05 * N)).
#' @param subsample_size Number of subsampled points to test (default 300).
#' @return Fraction of points with significant neighbor composition mismatch (p <= 0.05).
#' @export
calc_nn_mismatch <- function(ref, sim, k = NULL, subsample_size = 300) {
  # Handle 1D or 2D
  if (is.matrix(ref) || is.data.frame(ref)) {
    ref_mat <- as.matrix(stats::na.omit(ref))
    sim_mat <- as.matrix(stats::na.omit(sim))
  } else {
    ref_mat <- matrix(stats::na.omit(as.numeric(ref)), ncol = 1)
    sim_mat <- matrix(stats::na.omit(as.numeric(sim)), ncol = 1)
  }
  
  n1 <- nrow(ref_mat)
  n2 <- nrow(sim_mat)
  if (n1 < 3 || n2 < 3) return(NA_real_)
  
  df_all <- rbind(ref_mat, sim_mat)
  labels <- c(rep("ref", n1), rep("sim", n2))
  n_total <- n1 + n2
  
  if (is.null(k)) {
    k <- max(5, round(0.05 * n_total))
  }
  k <- min(k, n_total - 1)
  
  # Subsample observations for speed
  idx <- sample(seq_len(n_total), min(n_total, subsample_size))
  overall_prop <- c(ref = n1 / n_total, sim = n2 / n_total)
  
  # Calculate distance and Chi-squared test for each sampled point
  mismatches <- vapply(idx, function(j) {
    diffs <- sweep(df_all, 2, df_all[j, ], "-")
    dists <- sqrt(rowSums(diffs^2))
    dists[j] <- Inf # exclude self
    nn_labels <- labels[order(dists)[seq_len(k)]]
    tab <- table(factor(nn_labels, levels = c("ref", "sim")))
    
    pval <- tryCatch(
      suppressWarnings(stats::chisq.test(tab, p = overall_prop)$p.value),
      error = function(e) 1
    )
    as.numeric(!is.na(pval) && pval <= 0.05)
  }, numeric(1))
  
  mean(mismatches, na.rm = TRUE)
}

#' Calculate Between-Dataset Silhouette Width
#'
#' Integrated from countsimQC (Soneson & Robinson). Treats dataset identity
#' (reference vs. simulation) as cluster labels to test if the two datasets
#' separate into distinct clusters or remain well-mixed.
#'
#' @param ref Vector or matrix for reference dataset.
#' @param sim Vector or matrix for simulated dataset.
#' @param subsample_size Number of subsampled points to evaluate (default 300).
#' @return A named list with global and local between-dataset silhouette widths.
#' @export
calc_between_dataset_silhouette <- function(ref, sim, subsample_size = 300) {
  if (is.matrix(ref) || is.data.frame(ref)) {
    ref_mat <- as.matrix(stats::na.omit(ref))
    sim_mat <- as.matrix(stats::na.omit(sim))
  } else {
    ref_mat <- matrix(stats::na.omit(as.numeric(ref)), ncol = 1)
    sim_mat <- matrix(stats::na.omit(as.numeric(sim)), ncol = 1)
  }
  
  n1 <- nrow(ref_mat)
  n2 <- nrow(sim_mat)
  if (n1 < 3 || n2 < 3) return(list(global_silh = NA_real_, local_silh = NA_real_))
  
  df_all <- rbind(ref_mat, sim_mat)
  ds_labels <- c(rep(1, n1), rep(2, n2))
  n_total <- n1 + n2
  
  idx <- sample(seq_len(n_total), min(n_total, subsample_size))
  k_local <- max(5, round(0.05 * n_total))
  
  silh_vals <- t(vapply(idx, function(j) {
    diffs <- sweep(df_all, 2, df_all[j, ], "-")
    dists <- sqrt(rowSums(diffs^2))
    
    dists_this <- dists[ds_labels == ds_labels[j] & seq_len(n_total) != j]
    dists_other <- dists[ds_labels != ds_labels[j]]
    
    a <- mean(dists_this, na.rm = TRUE)
    b <- mean(dists_other, na.rm = TRUE)
    s_global <- if (max(a, b) > 0) (b - a) / max(a, b) else 0
    
    # Local silhouette (k-NN)
    k_this <- min(k_local, length(dists_this))
    k_other <- min(k_local, length(dists_other))
    a_local <- mean(sort(dists_this)[seq_len(k_this)], na.rm = TRUE)
    b_local <- mean(sort(dists_other)[seq_len(k_other)], na.rm = TRUE)
    s_local <- if (max(a_local, b_local) > 0) (b_local - a_local) / max(a_local, b_local) else 0
    
    c(s_global = s_global, s_local = s_local)
  }, c(s_global = NA_real_, s_local = NA_real_)))
  
  list(
    between_dataset_silh_global = mean(silh_vals[, "s_global"], na.rm = TRUE),
    between_dataset_silh_local = mean(silh_vals[, "s_local"], na.rm = TRUE)
  )
}

# -----------------------------------------------------------------------------
# Combined Batch Calculators (Unifying simpipe, HelenaLC, SimBench, & countsimQC)
# -----------------------------------------------------------------------------

#' Compute All Univariate Distance & Accuracy Metrics for a Feature
#' @param ref Numeric vector of reference values.
#' @param sim Numeric vector of simulated values.
#' @param metric_prefix Optional prefix string for metric names.
#' @return Named list of univariate distance metrics.
#' @export
calc_all_univariate_metrics <- function(ref, sim, metric_prefix = "") {
  prefix <- if (nchar(metric_prefix) > 0) paste0(metric_prefix, "_") else ""
  
  runs_res <- calc_runs_test(ref, sim)
  silh_res <- calc_between_dataset_silhouette(ref, sim)
  
  res <- list(
    MAD = calc_mad(ref, sim),
    KS = calc_ks(ref, sim),
    MAE = calc_mae(ref, sim),
    RMSE = calc_rmse(ref, sim),
    OV = calc_overlap(ref, sim),
    Bhattacharyya = calc_bhattacharyya(ref, sim),
    Wasserstein = calc_wasserstein_1d(ref, sim),
    ECDF_DiffArea = calc_ecdf_diffarea(ref, sim),
    Runs_Statistic = runs_res$runs_statistic,
    Runs_PValue = runs_res$runs_pvalue,
    NN_Mismatch = calc_nn_mismatch(ref, sim),
    Between_Dataset_Silh_Global = silh_res$between_dataset_silh_global,
    Between_Dataset_Silh_Local = silh_res$between_dataset_silh_local
  )
  
  names(res) <- paste0(prefix, c(
    "MAD", "KS", "MAE", "RMSE", "OV", "Bhattacharyya", "Wasserstein",
    "ECDF_DiffArea", "Runs_Statistic", "Runs_PValue", "NN_Mismatch",
    "Between_Dataset_Silh_Global", "Between_Dataset_Silh_Local"
  ))
  return(res)
}

#' Compute All Bivariate (2D) Distance & Accuracy Metrics for Joint Distributions
#' @param ref_mat 2-column matrix for reference.
#' @param sim_mat 2-column matrix for simulation.
#' @param metric_prefix Optional prefix string for metric names.
#' @param threads CPU threads.
#' @return Named list of bivariate metrics.
#' @export
calc_all_bivariate_metrics <- function(ref_mat, sim_mat, metric_prefix = "", threads = 1) {
  prefix <- if (nchar(metric_prefix) > 0) paste0(metric_prefix, "_") else ""
  
  silh_2d <- calc_between_dataset_silhouette(ref_mat, sim_mat)
  
  res <- list(
    Fasano_Franceschini_2D_KS = calc_fasano_franceschini(ref_mat, sim_mat, threads = threads),
    Peacock_2D_KS = calc_peacock_2d(ref_mat, sim_mat),
    KDE_Bivariate_zstat = calc_kde_test(ref_mat, sim_mat),
    EMD_2D = calc_emd_2d(ref_mat, sim_mat),
    NN_Mismatch_2D = calc_nn_mismatch(ref_mat, sim_mat),
    Between_Dataset_Silh_2D_Global = silh_2d$between_dataset_silh_global,
    Between_Dataset_Silh_2D_Local = silh_2d$between_dataset_silh_local
  )
  
  names(res) <- paste0(prefix, c(
    "Fasano_Franceschini_2D_KS", "Peacock_2D_KS", "KDE_Bivariate_zstat",
    "EMD_2D", "NN_Mismatch_2D", "Between_Dataset_Silh_2D_Global",
    "Between_Dataset_Silh_2D_Local"
  ))
  return(res)
}

# -----------------------------------------------------------------------------
# Generative Model Selection & Goodness of Fit (scDesign3; Song et al., 2024)
# -----------------------------------------------------------------------------

#' Compute Model Information Criteria (AIC and BIC)
#'
#' Evaluates the statistical trade-off between model fit and parameter complexity
#' for single-cell generative models as benchmarked in scDesign3 (Song et al., Nat Biotechnol 2024).
#'
#' @param loglik Numeric vector or scalar of log-likelihood values.
#' @param n_params Numeric vector or scalar of number of estimated parameters (degrees of freedom).
#' @param n_obs Total number of independent observations (cells or cell-gene pairs).
#'
#' @return A named vector or data.frame containing loglik, n_params, n_obs, AIC, and BIC.
#' @export
calc_model_aic_bic <- function(loglik, n_params, n_obs) {
  if (length(loglik) != length(n_params)) {
    stop("loglik and n_params must have the same length.")
  }
  aic <- 2 * n_params - 2 * loglik
  bic <- log(n_obs) * n_params - 2 * loglik
  
  if (length(loglik) == 1) {
    return(c(loglik = as.numeric(loglik), n_params = as.numeric(n_params), n_obs = as.numeric(n_obs), AIC = aic, BIC = bic))
  } else {
    return(data.frame(
      loglik = as.numeric(loglik),
      n_params = as.numeric(n_params),
      n_obs = rep(n_obs, length.out = length(loglik)),
      AIC = aic,
      BIC = bic
    ))
  }
}

#' Marginal Model Goodness of Fit and Information Criteria for Single-Cell Simulators
#'
#' Evaluates gene-wise and aggregate marginal goodness-of-fit (Poisson, Negative Binomial, or Gaussian)
#' for simulated or fitted single-cell count matrices, computing total and mean AIC and BIC
#' across genes (scDesign3; Song et al., 2024).
#'
#' @param counts Matrix or data.frame of observed or simulated counts (genes x cells).
#' @param fitted_means Matrix or data.frame of model fitted means or expectations (genes x cells).
#' @param dispersions Optional vector of gene-level dispersion parameters for Negative Binomial model.
#'   If NULL, estimated via method-of-moments.
#' @param distribution Parametric distribution: "poisson", "nb" (Negative Binomial), or "gaussian".
#' @param n_params_per_gene Number of estimated parameters per gene. If NULL, defaults to 1 for Poisson,
#'   2 for NB, and 2 for Gaussian.
#'
#' @return A list containing:
#'   \item{total_loglik}{Sum of marginal log-likelihoods across all genes and cells.}
#'   \item{total_aic}{Aggregate AIC across all genes.}
#'   \item{total_bic}{Aggregate BIC across all genes.}
#'   \item{mean_gene_aic}{Mean AIC per gene.}
#'   \item{mean_gene_bic}{Mean BIC per gene.}
#'   \item{gene_loglik}{Vector of log-likelihoods per gene.}
#'   \item{gene_aic}{Vector of AIC per gene.}
#'   \item{gene_bic}{Vector of BIC per gene.}
#' @export
calc_marginal_aic_bic <- function(counts, fitted_means, dispersions = NULL, 
                                  distribution = c("poisson", "nb", "gaussian"),
                                  n_params_per_gene = NULL) {
  distribution <- match.arg(distribution)
  counts <- as.matrix(counts)
  fitted_means <- as.matrix(fitted_means)
  
  if (all(dim(counts) != dim(fitted_means))) {
    if (nrow(counts) == ncol(fitted_means) && ncol(counts) == nrow(fitted_means)) {
      fitted_means <- t(fitted_means)
    } else {
      stop("counts and fitted_means must have matching dimensions.")
    }
  }
  
  G <- nrow(counts)
  N <- ncol(counts)
  
  fitted_means[fitted_means <= 0] <- 1e-10
  
  if (is.null(n_params_per_gene)) {
    n_params_per_gene <- switch(distribution,
      poisson = 1,
      nb = 2,
      gaussian = 2
    )
  }
  
  if (distribution == "poisson") {
    ll_mat <- counts * log(fitted_means) - fitted_means - lgamma(counts + 1)
  } else if (distribution == "nb") {
    if (is.null(dispersions)) {
      dispersions <- vapply(seq_len(G), function(g) {
        y <- counts[g, ]
        mu <- fitted_means[g, ]
        var_y <- stats::var(y)
        mean_mu <- mean(mu)
        d <- (var_y - mean_mu) / max(mean_mu^2, 1e-4)
        max(1e-4, min(d, 100))
      }, numeric(1))
    }
    if (length(dispersions) == 1) dispersions <- rep(dispersions, G)
    theta_vec <- 1 / dispersions
    
    ll_mat <- matrix(0, nrow = G, ncol = N)
    for (g in seq_len(G)) {
      th <- theta_vec[g]
      y <- counts[g, ]
      mu <- fitted_means[g, ]
      ll_mat[g, ] <- lgamma(y + th) - lgamma(th) - lgamma(y + 1) + 
                     th * log(th / (th + mu)) + y * log(mu / (th + mu))
    }
  } else {
    var_res <- apply((counts - fitted_means)^2, 1, mean)
    var_res[var_res <= 0] <- 1e-6
    ll_mat <- matrix(0, nrow = G, ncol = N)
    for (g in seq_len(G)) {
      sig2 <- var_res[g]
      ll_mat[g, ] <- -0.5 * log(2 * pi * sig2) - ((counts[g, ] - fitted_means[g, ])^2) / (2 * sig2)
    }
  }
  
  gene_loglik <- rowSums(ll_mat)
  gene_aic <- 2 * n_params_per_gene - 2 * gene_loglik
  gene_bic <- log(N) * n_params_per_gene - 2 * gene_loglik
  
  total_loglik <- sum(gene_loglik)
  total_params <- sum(rep(n_params_per_gene, length.out = G))
  total_aic <- 2 * total_params - 2 * total_loglik
  total_bic <- log(N) * total_params - 2 * total_loglik
  
  list(
    total_loglik = as.numeric(total_loglik),
    total_aic = as.numeric(total_aic),
    total_bic = as.numeric(total_bic),
    mean_gene_aic = as.numeric(mean(gene_aic)),
    mean_gene_bic = as.numeric(mean(gene_bic)),
    gene_loglik = gene_loglik,
    gene_aic = gene_aic,
    gene_bic = gene_bic
  )
}

#' Perform Likelihood Ratio Test for Comparing Nested Single-Cell Simulation Models
#'
#' Direct port and generalization of \code{scDesign3::perform_lrt} (Song et al., Nat Biotechnol 2024).
#' Performs the likelihood ratio test to compare two nested simulation models (e.g., cell-type/covariate
#' model vs intercept-only null model, or spline trajectory vs linear model).
#'
#' @param alter_model Alternative model (more complex) or list of alternative models per gene, or numeric log-likelihoods.
#' @param null_model Null model (simpler, strictly nested) or list of null models per gene, or numeric log-likelihoods.
#' @param df_alter Degrees of freedom for alternative model (used if models are numeric log-likelihoods).
#' @param df_null Degrees of freedom for null model (used if models are numeric log-likelihoods).
#'
#' @return A data.frame containing:
#'   \item{LogLik_alter}{Log-likelihood under alternative model.}
#'   \item{LogLik_null}{Log-likelihood under null model.}
#'   \item{df_alter}{Degrees of freedom of alternative model.}
#'   \item{df_null}{Degrees of freedom of null model.}
#'   \item{LR_statistic}{Likelihood ratio statistic (-2 * (LL_null - LL_alter)).}
#'   \item{delta_df}{Difference in degrees of freedom (df_alter - df_null).}
#'   \item{p_value}{P-value from chi-squared test with delta_df degrees of freedom.}
#' @export
calc_likelihood_ratio_test <- function(alter_model, null_model, df_alter = NULL, df_null = NULL) {
  # Mode 1: Numeric log-likelihoods passed directly
  if (is.numeric(alter_model) && is.numeric(null_model)) {
    ll_alt <- alter_model
    ll_null <- null_model
    k_alt <- if (!is.null(df_alter)) df_alter else NA_real_
    k_null <- if (!is.null(df_null)) df_null else NA_real_
    
    lr <- 2 * (ll_alt - ll_null)
    delta_df <- k_alt - k_null
    
    p_val <- ifelse(!is.na(delta_df) & delta_df > 0 & lr >= 0,
                    stats::pchisq(lr, df = delta_df, lower.tail = FALSE),
                    NA_real_)
    
    return(data.frame(
      LogLik_alter = ll_alt,
      LogLik_null = ll_null,
      df_alter = k_alt,
      df_null = k_null,
      LR_statistic = lr,
      delta_df = delta_df,
      p_value = p_val
    ))
  }
  
  # Mode 2: Model objects (single model or list of models)
  if (inherits(alter_model, c("glm", "lm", "gam", "gamlss", "logLik")) || !is.list(alter_model)) {
    alter_model <- list(alter_model)
  }
  if (inherits(null_model, c("glm", "lm", "gam", "gamlss", "logLik")) || !is.list(null_model)) {
    null_model <- list(null_model)
  }
  
  n_genes <- min(length(alter_model), length(null_model))
  
  extract_df <- function(m) {
    if (inherits(m, "gamlss")) {
      m$df.fit
    } else if (inherits(m, "gam")) {
      sum(m$edf2 + m$edf1 - m$edf)
    } else if (inherits(m, "lm") || inherits(m, "glm")) {
      attr(stats::logLik(m), "df")
    } else if (!is.null(m$df)) {
      m$df
    } else {
      length(stats::coef(m))
    }
  }
  
  res_list <- lapply(seq_len(n_genes), function(i) {
    m1 <- alter_model[[i]]
    m2 <- null_model[[i]]
    
    k1 <- extract_df(m1)
    k2 <- extract_df(m2)
    
    ll1 <- as.numeric(stats::logLik(m1))
    ll2 <- as.numeric(stats::logLik(m2))
    
    lr <- 2 * (ll1 - ll2)
    delta_df <- k1 - k2
    
    p_val <- if (!is.na(delta_df) && delta_df > 0 && lr >= 0) {
      stats::pchisq(lr, df = delta_df, lower.tail = FALSE)
    } else {
      NA_real_
    }
    
    c(LogLik_alter = ll1, LogLik_null = ll2, df_alter = k1, df_null = k2,
      LR_statistic = lr, delta_df = delta_df, p_value = p_val)
  })
  
  df_res <- as.data.frame(do.call(rbind, res_list))
  if (!is.null(names(alter_model))) {
    rownames(df_res) <- names(alter_model)[seq_len(n_genes)]
  }
  return(df_res)
}

# -----------------------------------------------------------------------------
# Deep Generative Models: Maximum Mean Discrepancy & Fréchet Single-Cell Distance
# (scGAN, ACTIVA, scDiffusion, MichiGAN, scvi-tools, cfDiffusion)
# -----------------------------------------------------------------------------

#' Maximum Mean Discrepancy (MMD) with Gaussian RBF Kernel
#'
#' Evaluates the distributional discrepancy between empirical reference and
#' generated/simulated single-cell high-dimensional profiles or latent embeddings.
#'
#' @param ref_mat Matrix of reference cells (cells x features or features x cells).
#' @param sim_mat Matrix of simulated cells (cells x features or features x cells).
#' @param cells_as_cols Logical, whether cells are columns. If TRUE (default for scRNA-seq),
#'   the matrix is transposed so rows represent cells.
#' @param sigma Optional kernel bandwidth. If NULL, uses the median pairwise distance heuristic.
#' @param max_cells Maximum number of cells to subsample for speed (default 500).
#'
#' @return A list with MMD, squared MMD, and the kernel bandwidth sigma used.
#' @export
calc_mmd <- function(
  ref_mat,
  sim_mat,
  cells_as_cols = TRUE,
  sigma = NULL,
  max_cells = 500
) {
  if (inherits(ref_mat, "SingleCellExperiment")) {
    ref_mat <- SummarizedExperiment::assay(ref_mat, "counts")
  }
  if (inherits(sim_mat, "SingleCellExperiment")) {
    sim_mat <- SummarizedExperiment::assay(sim_mat, "counts")
  }
  
  x <- if (cells_as_cols) t(as.matrix(ref_mat)) else as.matrix(ref_mat)
  y <- if (cells_as_cols) t(as.matrix(sim_mat)) else as.matrix(sim_mat)
  
  # Subsample cells if needed
  if (nrow(x) > max_cells) x <- x[sample(nrow(x), max_cells), , drop = FALSE]
  if (nrow(y) > max_cells) y <- y[sample(nrow(y), max_cells), , drop = FALSE]
  
  m <- nrow(x)
  n <- nrow(y)
  
  # Pairwise squared Euclidean distance
  dist_xx <- as.matrix(stats::dist(x))^2
  dist_yy <- as.matrix(stats::dist(y))^2
  
  # Cross distance
  dist_xy <- matrix(0, nrow = m, ncol = n)
  for (i in seq_len(m)) {
    diffs <- sweep(y, 2, x[i, ], "-")
    dist_xy[i, ] <- rowSums(diffs^2)
  }
  
  # Median heuristic for bandwidth sigma
  if (is.null(sigma) || is.na(sigma) || sigma <= 0) {
    med_dist <- stats::median(dist_xx[upper.tri(dist_xx)])
    sigma <- if (!is.na(med_dist) && med_dist > 0) sqrt(0.5 * med_dist) else 1.0
  }
  
  gamma <- 1 / (2 * (sigma^2))
  
  k_xx <- exp(-gamma * dist_xx)
  k_yy <- exp(-gamma * dist_yy)
  k_xy <- exp(-gamma * dist_xy)
  
  # Unbiased MMD^2 estimator
  diag(k_xx) <- 0
  diag(k_yy) <- 0
  
  term_xx <- if (m > 1) sum(k_xx) / (m * (m - 1)) else 1
  term_yy <- if (n > 1) sum(k_yy) / (n * (n - 1)) else 1
  term_xy <- (2 * sum(k_xy)) / (m * n)
  
  mmd2 <- term_xx + term_yy - term_xy
  mmd_val <- sqrt(max(0, mmd2))
  
  list(
    mmd = as.numeric(mmd_val),
    mmd_squared = as.numeric(mmd2),
    bandwidth_sigma = as.numeric(sigma)
  )
}

#' Fréchet Single-Cell Distance (FSD)
#'
#' Computes the single-cell analogue of Fréchet Inception Distance (FID) on
#' low-dimensional PCA embeddings between reference and simulated cell populations.
#'
#' @param ref_mat Matrix of reference cells (cells x features or features x cells).
#' @param sim_mat Matrix of simulated cells (cells x features or features x cells).
#' @param cells_as_cols Logical, whether cells are columns (default TRUE).
#' @param n_pcs Number of principal components to evaluate (default 15).
#'
#' @return A list with Fréchet distance (FSD), mean discrepancy, and covariance trace discrepancy.
#' @export
calc_frechet_singlecell_distance <- function(
  ref_mat,
  sim_mat,
  cells_as_cols = TRUE,
  n_pcs = 15
) {
  x <- if (cells_as_cols) t(as.matrix(ref_mat)) else as.matrix(ref_mat)
  y <- if (cells_as_cols) t(as.matrix(sim_mat)) else as.matrix(sim_mat)
  
  # Joint PCA projection for consistent coordinate frame
  n_x <- nrow(x)
  n_y <- nrow(y)
  combined <- rbind(x, y)
  
  # Filter zero-variance features
  vars <- apply(combined, 2, stats::var)
  valid_feats <- which(!is.na(vars) & vars > 1e-8)
  if (length(valid_feats) < 2) {
    return(list(fsd = NA_real_, mean_discrepancy = NA_real_, cov_discrepancy = NA_real_))
  }
  combined <- combined[, valid_feats, drop = FALSE]
  
  rank_use <- min(n_pcs, ncol(combined) - 1, nrow(combined) - 1)
  pca <- tryCatch(stats::prcomp(combined, rank. = rank_use, scale. = TRUE), error = function(e) NULL)
  if (is.null(pca)) {
    return(list(fsd = NA_real_, mean_discrepancy = NA_real_, cov_discrepancy = NA_real_))
  }
  
  pts_ref <- pca$x[seq_len(n_x), , drop = FALSE]
  pts_sim <- pca$x[(n_x + 1):(n_x + n_y), , drop = FALSE]
  
  mu_ref <- colMeans(pts_ref)
  mu_sim <- colMeans(pts_sim)
  
  cov_ref <- stats::cov(pts_ref)
  cov_sim <- stats::cov(pts_sim)
  
  # Mean difference squared norm
  mean_diff_sq <- sum((mu_ref - mu_sim)^2)
  
  # Covariance square root product trace: Tr((cov_ref %*% cov_sim)^(1/2))
  # via eigenvalues of cov_ref %*% cov_sim
  cov_prod <- cov_ref %*% cov_sim
  eigen_vals <- eigen(cov_prod, only.values = TRUE)$values
  cov_tr_sqrt <- sum(sqrt(pmax(0, Re(eigen_vals))))
  
  cov_diff <- sum(diag(cov_ref)) + sum(diag(cov_sim)) - 2 * cov_tr_sqrt
  fsd_sq <- mean_diff_sq + max(0, cov_diff)
  
  list(
    fsd = sqrt(fsd_sq),
    fsd_squared = fsd_sq,
    mean_discrepancy = mean_diff_sq,
    cov_discrepancy = max(0, cov_diff)
  )
}


