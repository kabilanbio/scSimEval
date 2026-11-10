# -----------------------------------------------------------------------------
# Category II: Zero-Inflation, Dropout Kinetics (scRNA-seq), &
#              Chromatin Accessibility-Sparsity Curves (scATAC-seq)
# -----------------------------------------------------------------------------

#' Fit Zero-Probability Dropout Curve (Splatter & ZINB-WaVE)
#'
#' Fits an empirical logistic dropout curve: logit(P(Y=0)) = beta_0 + beta_1 * log(mu)
#' to model the scRNA-seq dropout relationship with mean expression.
#'
#' @param counts Count matrix (genes x cells) or SingleCellExperiment.
#' @return A list containing intercept, slope, midpoint (inflection point), and R-squared.
#' @examples
#' data(example_scrna, package = "scSimEval")
#' calc_zero_probability_curve(example_scrna$ref)
#' @export
calc_zero_probability_curve <- function(counts) {
  if (inherits(counts, "SingleCellExperiment")) {
    counts <- SummarizedExperiment::assay(counts, "counts")
  }
  gene_means <- as.numeric(Matrix::rowMeans(counts))
  p0 <- fast_sparse_zero_fraction_feature(counts)
  
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
#' Compares empirical logistic dropout trends between reference and simulated scRNA-seq datasets.
#'
#' @param ref_counts Reference count matrix.
#' @param sim_counts Simulated count matrix.
#' @return Named numeric vector of curve parameter differences (\code{slope_error}, \code{midpoint_error}, \code{r2_discrepancy}).
#' @examples
#' data(example_scrna, package = "scSimEval")
#' evaluate_zero_probability_curve(example_scrna$ref, example_scrna$sim)
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
#' @param counts Count matrix (features x cells).
#' @return A list with mean excess zero weight, gene-level weights, and estimated zero-inflation rate.
#' @examples
#' data(example_scrna, package = "scSimEval")
#' calc_excess_zero_weights(example_scrna$ref)
#' @export
calc_excess_zero_weights <- function(counts) {
  if (inherits(counts, "SingleCellExperiment")) {
    counts <- SummarizedExperiment::assay(counts, "counts")
  }
  gene_means <- as.numeric(Matrix::rowMeans(counts))
  gene_vars <- fast_row_vars(counts)
  p0_obs <- fast_sparse_zero_fraction_feature(counts)
  
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
#' @examples
#' data(example_scatac, package = "scSimEval")
#' calc_accessibility_sparsity_curve(example_scatac$ref)
#' @export
calc_accessibility_sparsity_curve <- function(data, poly_degree = 2) {
  if (inherits(data, "SingleCellExperiment")) {
    data <- SummarizedExperiment::assay(data, "counts")
  }
  if (is.matrix(data) || inherits(data, "Matrix")) {
    means <- as.numeric(Matrix::rowMeans(data))
    nzp <- 1 - fast_sparse_zero_fraction_feature(data)
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
  spearman_cor <- suppressWarnings(stats::cor(means, nzp, method = "spearman"))
  pearson_cor <- suppressWarnings(stats::cor(means, nzp, method = "pearson"))
  
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
#' @return A list containing \code{curve_rmse} and \code{r2_discrepancy} (along with fitted reference and simulated curve parameters).
#' @examples
#' data(example_scatac, package = "scSimEval")
#' evaluate_accessibility_sparsity_curve(example_scatac$ref, example_scatac$sim)
#' @export
evaluate_accessibility_sparsity_curve <- function(ref_data, sim_data, poly_degree = 2) {
  if (inherits(ref_data, "SingleCellExperiment")) {
    ref_data <- SummarizedExperiment::assay(ref_data, "counts")
  }
  if (inherits(sim_data, "SingleCellExperiment")) {
    sim_data <- SummarizedExperiment::assay(sim_data, "counts")
  }
  ref_res <- calc_accessibility_sparsity_curve(ref_data, poly_degree = poly_degree)
  sim_res <- calc_accessibility_sparsity_curve(sim_data, poly_degree = poly_degree)
  
  r2_discrepancy <- abs(sim_res$r_squared - ref_res$r_squared)
  
  # Predict on common grid of means
  if (is.matrix(ref_data) || inherits(ref_data, "Matrix")) {
    ref_means <- as.numeric(Matrix::rowMeans(ref_data))
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
    curve_rmse = curve_rmse,
    r2_discrepancy = r2_discrepancy,
    ref_r_squared = ref_res$r_squared,
    sim_r_squared = sim_res$r_squared
  )
}
