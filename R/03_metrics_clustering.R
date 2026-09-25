#' @title Cell Clustering and Population Separation Metrics
#' @description Comprehensive suite of unsupervised clustering quality metrics
#'   and supervised ground-truth concordance measures integrated from simpipe and HelenaLC.
#' @name metrics_clustering
NULL

#' Calculate Average Silhouette Width (ASW)
#' @param dist_mat Distance matrix among cells or numeric matrix (features x cells).
#' @param cluster_labels Vector of cluster or cell-type labels.
#' @return Mean silhouette width.
#' @export
calc_silhouette <- function(dist_mat, cluster_labels) {
  if (!inherits(dist_mat, "dist")) {
    dist_mat <- stats::dist(t(as.matrix(dist_mat)))
  }
  cluster_labels <- as.numeric(as.factor(cluster_labels))
  if (requireNamespace("cluster", quietly = TRUE)) {
    sil <- cluster::silhouette(cluster_labels, dist_mat)
    mean(sil[, 3], na.rm = TRUE)
  } else {
    NA_real_
  }
}

#' Calculate Dunn Index
#' @param dist_mat Distance matrix or numeric matrix.
#' @param cluster_labels Cluster assignments.
#' @return Dunn index (ratio of smallest inter-cluster to largest intra-cluster distance).
#' @export
calc_dunn <- function(dist_mat, cluster_labels) {
  if (!inherits(dist_mat, "dist")) {
    dist_mat <- stats::dist(t(as.matrix(dist_mat)))
  }
  cluster_labels <- as.numeric(as.factor(cluster_labels))
  if (requireNamespace("clValid", quietly = TRUE)) {
    clValid::dunn(dist_mat, cluster_labels)
  } else {
    NA_real_
  }
}

#' Calculate Cluster Connectivity
#' @param dist_mat Distance matrix or numeric matrix.
#' @param cluster_labels Cluster assignments.
#' @return Connectivity metric.
#' @export
calc_connectivity <- function(dist_mat, cluster_labels) {
  if (!inherits(dist_mat, "dist")) {
    dist_mat <- stats::dist(t(as.matrix(dist_mat)))
  }
  cluster_labels <- as.numeric(as.factor(cluster_labels))
  if (requireNamespace("clValid", quietly = TRUE)) {
    clValid::connectivity(dist_mat, cluster_labels)
  } else {
    NA_real_
  }
}

#' Calculate Davies-Bouldin Index (DB)
#' @param data Matrix with features in rows, cells in columns.
#' @param cluster_labels Cluster assignments.
#' @return Davies-Bouldin index (lower is better).
#' @export
calc_davies_bouldin <- function(data, cluster_labels) {
  cluster_labels <- as.numeric(as.factor(cluster_labels))
  if (requireNamespace("clusterSim", quietly = TRUE)) {
    clusterSim::index.DB(t(as.matrix(data)), cluster_labels)$DB
  } else {
    NA_real_
  }
}

#' Calculate Calinski-Harabasz Index (CH)
#' @param data Matrix with features in rows, cells in columns.
#' @param cluster_labels Cluster assignments.
#' @return Calinski-Harabasz index.
#' @export
calc_calinski_harabasz <- function(data, cluster_labels) {
  cluster_labels <- as.numeric(as.factor(cluster_labels))
  if (requireNamespace("fpc", quietly = TRUE)) {
    fpc::calinhara(t(as.matrix(data)), cluster_labels)
  } else {
    NA_real_
  }
}

#' Calculate Adjusted Rand Index (ARI)
#' @param pred Predicted cluster labels.
#' @param truth Ground truth cell type labels.
#' @return ARI value between -1 and 1.
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

#' Expected Mutual Information for Cluster Comparison
#'
#' Computes the exact expectation of mutual information under the generalized
#' hypergeometric model of randomness with fixed marginals (Vinh et al., 2010).
#'
#' @param a Integer vector of row sums (marginal cluster sizes of partition 1).
#' @param b Integer vector of column sums (marginal cluster sizes of partition 2).
#' @param N Total number of items / cells.
#' @return Expected mutual information in nats.
#' @keywords internal
calc_expected_mi <- function(a, b, N) {
  emi <- 0
  for (i in seq_along(a)) {
    a_i <- a[i]
    if (a_i == 0) next
    for (j in seq_along(b)) {
      b_j <- b[j]
      if (b_j == 0) next
      k_min <- max(1, a_i + b_j - N)
      k_max <- min(a_i, b_j)
      if (k_min > k_max) next
      
      k_vals <- k_min:k_max
      probs <- stats::dhyper(k_vals, a_i, N - a_i, b_j)
      term <- (k_vals / N) * log((N * k_vals) / (a_i * b_j))
      emi <- emi + sum(probs * term)
    }
  }
  emi
}

#' Calculate Adjusted Mutual Information (AMI)
#'
#' Computes the Adjusted Mutual Information (AMI) between two clusterings
#' or cell-type label assignments, correcting for chance expected mutual information
#' under the hypergeometric model with fixed marginals (Vinh et al., JMLR 2010),
#' as used in scMoMtF (Lan et al., PLOS Comput Biol 2024).
#'
#' @param pred Vector of predicted cluster labels.
#' @param truth Vector of ground truth cell type labels.
#' @param average_method Method to normalize mutual information: "arithmetic" (default, scikit-learn standard),
#'   "max", "min", or "geometric".
#' @return AMI score in [0, 1] (adjusted for chance).
#' @references Vinh, N. X., Epps, J., & Bailey, J. (2010). Information theoretic measures for clusterings comparison: Variants, properties, normalization and correction for chance. Journal of Machine Learning Research, 11, 2837-2854.
#' @references Lan, W., Ling, T., Chen, Q. et al. scMoMtF: An interpretable multitask learning framework for single-cell multi-omics data analysis. PLOS Comput Biol 20(12): e1012679 (2024).
#' @export
calc_ami <- function(pred, truth, average_method = c("arithmetic", "max", "min", "geometric")) {
  average_method <- match.arg(average_method)
  valid <- !is.na(pred) & !is.na(truth)
  pred <- as.character(pred[valid])
  truth <- as.character(truth[valid])
  
  if (length(pred) < 2) return(NA_real_)
  
  tab <- table(pred, truth)
  N <- sum(tab)
  if (N == 0) return(0)
  
  a <- rowSums(tab)
  b <- colSums(tab)
  
  if (length(a) <= 1 || length(b) <= 1) return(0)
  
  p_a <- a[a > 0] / N
  p_b <- b[b > 0] / N
  h_a <- -sum(p_a * log(p_a))
  h_b <- -sum(p_b * log(p_b))
  
  if (h_a == 0 || h_b == 0) return(0)
  
  p_tab <- tab / N
  mi <- 0
  for (i in seq_len(nrow(tab))) {
    for (j in seq_len(ncol(tab))) {
      if (tab[i, j] > 0) {
        mi <- mi + p_tab[i, j] * log((N * tab[i, j]) / (a[i] * b[j]))
      }
    }
  }
  
  emi <- calc_expected_mi(a, b, N)
  
  normalizer <- switch(
    average_method,
    "arithmetic" = (h_a + h_b) / 2,
    "max" = max(h_a, h_b),
    "min" = min(h_a, h_b),
    "geometric" = sqrt(h_a * h_b)
  )
  
  denom <- normalizer - emi
  if (denom <= 0) return(0)
  
  ami <- (mi - emi) / denom
  max(0, min(1, ami))
}

#' Calculate Clustering Accuracy (ACC)
#'
#' Integrated from scCluBench (Xu et al., AAAI 2026). Aligns predicted clusters
#' to ground truth classes using the Hungarian bipartite matching algorithm,
#' then computes the fraction of correctly mapped cells.
#'
#' @param pred Vector of predicted cluster labels.
#' @param truth Vector of ground truth cell type labels.
#' @return Clustering accuracy between 0 and 1.
#' @export
calc_clustering_accuracy <- function(pred, truth) {
  match_res <- hungarian_match(pred, truth)
  match_res$accuracy
}

#' Calculate Fowlkes-Mallows Index (FMI)
#'
#' Integrated from scCluBench (Xu et al., AAAI 2026). Computes the geometric mean
#' of pairwise cluster precision and recall.
#'
#' @param pred Vector of predicted cluster labels.
#' @param truth Vector of ground truth cell type labels.
#' @return FMI score between 0 and 1.
#' @export
calc_fmi <- function(pred, truth) {
  valid <- !is.na(pred) & !is.na(truth)
  pred <- as.character(pred[valid])
  truth <- as.character(truth[valid])
  if (length(pred) < 2) return(NA_real_)
  
  tab <- table(pred, truth)
  tp <- sum(choose(tab, 2))
  tp_fp <- sum(choose(rowSums(tab), 2))
  tp_fn <- sum(choose(colSums(tab), 2))
  if (tp_fp == 0 || tp_fn == 0) return(0)
  tp / sqrt(tp_fp * tp_fn)
}

#' Calculate Homogeneity, Completeness, and V-Measure
#'
#' Integrated from scCluBench (Xu et al., AAAI 2026).
#'
#' @param pred Vector of predicted cluster labels.
#' @param truth Vector of ground truth cell type labels.
#' @param beta Weight of completeness vs. homogeneity (default 1).
#' @return A named list containing homogeneity, completeness, and v_measure.
#' @export
calc_homogeneity_completeness_v_measure <- function(pred, truth, beta = 1) {
  valid <- !is.na(pred) & !is.na(truth)
  pred <- as.character(pred[valid])
  truth <- as.character(truth[valid])
  if (length(pred) == 0) {
    return(list(homogeneity = NA_real_, completeness = NA_real_, v_measure = NA_real_))
  }
  
  tab <- table(pred, truth)
  n <- sum(tab)
  if (n == 0) {
    return(list(homogeneity = 0, completeness = 0, v_measure = 0))
  }
  
  p_k <- rowSums(tab) / n
  p_c <- colSums(tab) / n
  h_c <- -sum(p_c[p_c > 0] * log(p_c[p_c > 0]))
  h_k <- -sum(p_k[p_k > 0] * log(p_k[p_k > 0]))
  
  # H(C|K)
  h_c_k <- 0
  for (k in seq_len(nrow(tab))) {
    n_k <- sum(tab[k, ])
    if (n_k > 0) {
      p_c_given_k <- tab[k, ] / n_k
      p_c_given_k <- p_c_given_k[p_c_given_k > 0]
      h_c_k <- h_c_k + (n_k / n) * (-sum(p_c_given_k * log(p_c_given_k)))
    }
  }
  
  # H(K|C)
  h_k_c <- 0
  for (c in seq_len(ncol(tab))) {
    n_c <- sum(tab[, c])
    if (n_c > 0) {
      p_k_given_c <- tab[, c] / n_c
      p_k_given_c <- p_k_given_c[p_k_given_c > 0]
      h_k_c <- h_k_c + (n_c / n) * (-sum(p_k_given_c * log(p_k_given_c)))
    }
  }
  
  hom <- if (h_c == 0) 1 else 1 - (h_c_k / h_c)
  com <- if (h_k == 0) 1 else 1 - (h_k_c / h_k)
  hom <- max(0, min(1, hom))
  com <- max(0, min(1, com))
  v <- if (hom + com == 0) 0 else (1 + beta) * (hom * com) / (beta * hom + com)
  v <- max(0, min(1, v))
  
  list(homogeneity = hom, completeness = com, v_measure = v)
}

#' Calculate Homogeneity Score
#' @param pred Vector of predicted cluster labels.
#' @param truth Vector of ground truth cell type labels.
#' @return Homogeneity score between 0 and 1.
#' @export
calc_homogeneity <- function(pred, truth) {
  calc_homogeneity_completeness_v_measure(pred, truth)$homogeneity
}

#' Calculate Completeness Score
#' @param pred Vector of predicted cluster labels.
#' @param truth Vector of ground truth cell type labels.
#' @return Completeness score between 0 and 1.
#' @export
calc_completeness <- function(pred, truth) {
  calc_homogeneity_completeness_v_measure(pred, truth)$completeness
}

#' Calculate V-Measure
#' @param pred Vector of predicted cluster labels.
#' @param truth Vector of ground truth cell type labels.
#' @param beta Weight parameter (default 1).
#' @return V-measure score between 0 and 1.
#' @export
calc_v_measure <- function(pred, truth, beta = 1) {
  calc_homogeneity_completeness_v_measure(pred, truth, beta = beta)$v_measure
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
#' Evaluates the deviation and goodness-of-fit of single-cell clustering labels directly
#' from raw or normalized count matrices without requiring ground-truth cell-type labels
#' (Fang et al., Genome Biology 2022; benchmarked in scDesign3, Song et al., Nat Biotechnol 2024).
#' CDI calculates the penalized log-likelihood (CDI-AIC and CDI-BIC) of count data
#' under a cluster-specific Poisson or Negative Binomial GLM with cell size factors.
#'
#' @param counts Matrix or data.frame of counts (genes x cells or cells x genes).
#' @param cluster_labels Vector of cluster assignments for each cell.
#' @param size_factors Optional numeric vector of cell-specific size factors. If NULL,
#'   computed as library size scaled by median library size.
#' @param model Distribution model to use: "poisson" (fast, default) or "nb" (Negative Binomial).
#' @param top_features Optional integer; if set, restricts CDI calculation to the top most
#'   variable features to accelerate computation on large matrices.
#'
#' @return A list containing:
#'   \item{cdi_aic}{Normalized CDI score based on AIC penalty (lower indicates better clustering).}
#'   \item{cdi_bic}{Normalized CDI score based on BIC penalty (favors parsimonious main clusters).}
#'   \item{loglik}{Total log-likelihood under the cluster-specific count model.}
#'   \item{deviance}{Model deviance relative to saturated count model.}
#'   \item{n_clusters}{Number of clusters evaluated.}
#'   \item{n_params}{Number of free parameters in the clustering model.}
#' @export
calc_cdi <- function(counts, cluster_labels, size_factors = NULL, 
                     model = c("poisson", "nb"), top_features = NULL) {
  model <- match.arg(model)
  counts <- as.matrix(counts)
  if (ncol(counts) != length(cluster_labels) && nrow(counts) == length(cluster_labels)) {
    counts <- t(counts) # Ensure features x cells
  }
  
  cluster_labels <- as.factor(cluster_labels)
  clusters <- levels(cluster_labels)
  K <- length(clusters)
  N <- ncol(counts)
  G <- nrow(counts)
  
  if (K < 2 || N < 2) {
    return(list(cdi_aic = NA_real_, cdi_bic = NA_real_, loglik = NA_real_, deviance = NA_real_, n_clusters = K, n_params = 0))
  }
  
  # Select top variable/expressed features if requested to accelerate
  if (!is.null(top_features) && top_features < G) {
    gene_vars <- apply(counts, 1, stats::var)
    keep_idx <- order(gene_vars, decreasing = TRUE)[seq_len(top_features)]
    counts <- counts[keep_idx, , drop = FALSE]
    G <- nrow(counts)
  }
  
  # Compute size factors if not provided
  if (is.null(size_factors)) {
    lib_sizes <- colSums(counts)
    med_lib <- stats::median(lib_sizes[lib_sizes > 0])
    if (is.na(med_lib) || med_lib == 0) med_lib <- 1
    size_factors <- lib_sizes / med_lib
    size_factors[size_factors <= 0] <- min(size_factors[size_factors > 0], 1e-4)
  }
  
  # Fit cluster model: mean for each gene in each cluster
  sum_s_by_k <- vapply(clusters, function(k) sum(size_factors[cluster_labels == k]), numeric(1))
  sum_s_by_k[sum_s_by_k == 0] <- 1e-6
  
  ind_mat <- stats::model.matrix(~ 0 + cluster_labels)
  colnames(ind_mat) <- clusters
  
  sum_y_by_k <- counts %*% ind_mat
  lambda_gk <- sweep(sum_y_by_k, 2, sum_s_by_k, "/")
  lambda_gk[lambda_gk < 1e-10] <- 1e-10
  
  mu_mat <- lambda_gk[, as.integer(cluster_labels), drop = FALSE]
  mu_mat <- sweep(mu_mat, 2, size_factors, "*")
  mu_mat[mu_mat < 1e-10] <- 1e-10
  
  if (model == "poisson") {
    ll_mat <- counts * log(mu_mat) - mu_mat - lgamma(counts + 1)
    ll_cluster <- sum(ll_mat)
    n_params <- K * G
  } else {
    disp_vec <- vapply(seq_len(G), function(g) {
      y <- counts[g, ]
      mu <- mu_mat[g, ]
      var_res <- mean((y - mu)^2)
      mean_mu <- mean(mu)
      disp <- (var_res - mean_mu) / max(mean_mu^2, 1e-4)
      max(1e-4, min(disp, 100))
    }, numeric(1))
    
    theta_vec <- 1 / disp_vec
    ll_cluster <- 0
    for (g in seq_len(G)) {
      th <- theta_vec[g]
      y <- counts[g, ]
      mu <- mu_mat[g, ]
      ll_g <- lgamma(y + th) - lgamma(th) - lgamma(y + 1) + 
              th * log(th / (th + mu)) + y * log(mu / (th + mu))
      ll_cluster <- ll_cluster + sum(ll_g)
    }
    n_params <- K * G + G
  }
  
  tot_obs <- N * G
  cdi_aic <- (2 * n_params - 2 * ll_cluster) / tot_obs
  cdi_bic <- (log(N) * n_params - 2 * ll_cluster) / tot_obs
  
  y_pos <- counts
  y_pos[y_pos == 0] <- 1
  ll_sat <- sum(counts * log(y_pos) - counts - lgamma(counts + 1))
  deviance <- 2 * (ll_sat - ll_cluster)
  
  list(
    cdi_aic = as.numeric(cdi_aic),
    cdi_bic = as.numeric(cdi_bic),
    loglik = as.numeric(ll_cluster),
    deviance = as.numeric(deviance),
    n_clusters = K,
    n_params = n_params
  )
}

#' Full Clustering Performance Evaluation
#'
#' Evaluates both unsupervised cluster separation (ASW, Dunn, Connectivity, DB, CH, Neighborhood Purity, CDI)
#' and supervised concordance against ground truth labels (Clustering Accuracy ACC,
#' Hungarian F1/Precision/Recall, ARI, NMI, AMI, FMI, Homogeneity, Completeness, V-measure).
#'
#' @param data Count or normalized matrix (features x cells).
#' @param cluster_info Cluster or cell type labels.
#' @param pred_clusters Optional predicted cluster labels (if different from ground truth). If NULL, k-means is automatically run.
#' @param dist_mat Optional precomputed distance matrix.
#' @param cell_types Optional alias for cluster_info.
#' @param ref_data Optional reference data to compute reference clustering quality.
#'
#' @return A named list of all clustering metrics.
#' @export
evaluate_clustering_metrics <- function(
  data,
  cluster_info = NULL,
  pred_clusters = NULL,
  dist_mat = NULL,
  cell_types = NULL,
  ref_data = NULL
) {
  if (is.null(cluster_info) && !is.null(cell_types)) {
    cluster_info <- cell_types
  }
  if (is.null(cluster_info)) {
    stop("cluster_info or cell_types must be provided.")
  }
  cluster_info <- as.factor(cluster_info)
  
  if (is.null(dist_mat)) {
    dist_mat <- stats::dist(t(as.matrix(data)))
  }
  
  res <- list(
    silhouette_sim = calc_silhouette(dist_mat, cluster_info),
    silhouette = calc_silhouette(dist_mat, cluster_info),
    dunn_sim = calc_dunn(dist_mat, cluster_info),
    dunn = calc_dunn(dist_mat, cluster_info),
    connectivity_sim = calc_connectivity(dist_mat, cluster_info),
    connectivity = calc_connectivity(dist_mat, cluster_info),
    davies_bouldin_sim = calc_davies_bouldin(data, cluster_info),
    davies_bouldin = calc_davies_bouldin(data, cluster_info),
    calinski_harabasz_sim = calc_calinski_harabasz(data, cluster_info),
    calinski_harabasz = calc_calinski_harabasz(data, cluster_info),
    neighborhood_purity = mean(calc_neighborhood_purity(dist_mat, cluster_info, is_distance = TRUE), na.rm = TRUE)
  )
  
  if (!is.null(ref_data)) {
    ref_dist <- stats::dist(t(as.matrix(ref_data)))
    res$silhouette_ref <- calc_silhouette(ref_dist, cluster_info)
    res$dunn_ref <- calc_dunn(ref_dist, cluster_info)
    res$davies_bouldin_ref <- calc_davies_bouldin(ref_data, cluster_info)
    res$calinski_harabasz_ref <- calc_calinski_harabasz(ref_data, cluster_info)
  }
  
  # Unsupervised Clustering Deviation Index (CDI; Fang et al., 2022 / scDesign3)
  cdi_res <- tryCatch(calc_cdi(data, cluster_info), error = function(e) NULL)
  if (!is.null(cdi_res) && !is.na(cdi_res$cdi_aic)) {
    res$cdi <- cdi_res$cdi_aic
    res$CDI_AIC <- cdi_res$cdi_aic
    res$CDI_BIC <- cdi_res$cdi_bic
  }
  
  # If pred_clusters is not provided, run k-means with k = nlevels(cluster_info)
  if (is.null(pred_clusters) && nlevels(cluster_info) >= 2) {
    k_centers <- nlevels(cluster_info)
    set.seed(42)
    km_fit <- tryCatch(stats::kmeans(t(as.matrix(data)), centers = k_centers, nstart = 5), error = function(e) NULL)
    if (!is.null(km_fit)) {
      pred_clusters <- factor(km_fit$cluster)
    }
  }
  
  # Supervised concordance if predicted clusters are supplied or inferred
  if (!is.null(pred_clusters)) {
    match_res <- hungarian_match(pred_clusters, cluster_info)
    hcv <- calc_homogeneity_completeness_v_measure(pred_clusters, cluster_info)
    
    res$clustering_accuracy <- match_res$accuracy
    res$accuracy <- match_res$accuracy
    res$hungarian_F1 <- match_res$mean_F1
    res$hungarian_precision <- match_res$mean_precision
    res$hungarian_recall <- match_res$mean_recall
    res$ari <- calc_ari(pred_clusters, cluster_info)
    res$ARI <- res$ari
    res$nmi <- calc_nmi(pred_clusters, cluster_info)
    res$NMI <- res$nmi
    res$ami <- calc_ami(pred_clusters, cluster_info)
    res$AMI <- res$ami
    res$fmi <- calc_fmi(pred_clusters, cluster_info)
    res$FMI <- res$fmi
    res$homogeneity <- hcv$homogeneity
    res$completeness <- hcv$completeness
    res$v_measure <- hcv$v_measure
    res$matched_pairs <- match_res$matched_pairs
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
#'
#' @return A list containing Generative Precision, Generative Recall, and Generative F1.
#' @export
calc_generative_precision_recall <- function(
  ref_mat,
  sim_mat,
  cells_as_cols = TRUE,
  k = 5,
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
  
  if (nrow(x) > max_cells) x <- x[sample(nrow(x), max_cells), , drop = FALSE]
  if (nrow(y) > max_cells) y <- y[sample(nrow(y), max_cells), , drop = FALSE]
  
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
  # Radius of k-NN sphere around each reference cell
  radius_x <- apply(dist_xx, 1, function(row) sort(row)[k_x])
  
  # Pairwise distance within simulation: dist_yy
  dist_yy <- as.matrix(stats::dist(y))
  diag(dist_yy) <- Inf
  # Radius of k-NN sphere around each simulated cell
  radius_y <- apply(dist_yy, 1, function(row) sort(row)[k_y])
  
  # Cross distance between reference and simulation: dist_xy (m x n)
  dist_xy <- matrix(0, nrow = m, ncol = n)
  for (i in seq_len(m)) {
    diffs <- sweep(y, 2, x[i, ], "-")
    dist_xy[i, ] <- sqrt(rowSums(diffs^2))
  }
  
  # Precision: For each simulated cell y_j, is it within at least one reference sphere?
  # dist_xy[i, j] <= radius_x[i] for some i
  precision_hits <- apply(dist_xy, 2, function(col_dists) {
    any(col_dists <= radius_x)
  })
  precision_val <- mean(precision_hits)
  
  # Recall: For each reference cell x_i, is it within at least one simulated sphere?
  # dist_xy[i, j] <= radius_y[j] for some j
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

# -----------------------------------------------------------------------------
# Epigenomic Cell-Type Annotation & Projection (EpiAnno, SCAN-ATAC-Sim)
# -----------------------------------------------------------------------------

#' Evaluate Epigenomic Supervised Cell-Type Annotation (EpiAnno / SCAN-ATAC-Sim)
#'
#' Trains a supervised classification model on reference single-cell chromatin accessibility
#' profiles (scATAC-seq / scCAS) and evaluates cell-type prediction / projection accuracy
#' on simulated cells against ground-truth cell-type labels (Chen et al., Nat Mach Intell 2022).
#'
#' @param ref_data Reference count/accessibility matrix (features x cells).
#' @param sim_data Simulated count/accessibility matrix (features x cells).
#' @param ref_celltypes Factor or character vector of cell types for reference cells.
#' @param sim_celltypes Factor or character vector of true cell types for simulated cells.
#' @param method Classification method: "knn" (k-nearest neighbors on PCA, default) or
#'   "centroid" (nearest centroid classifier).
#' @param n_pcs Number of principal components for dimensionality reduction (default 20).
#' @param k Number of nearest neighbors for k-NN (default 5).
#'
#' @return A list containing overall accuracy, balanced accuracy, macro F1, macro precision,
#'   macro recall, Cohen's kappa, and per-class performance metrics.
#' @export
evaluate_epigenomic_annotation <- function(
  ref_data,
  sim_data,
  ref_celltypes,
  sim_celltypes,
  method = c("knn", "centroid"),
  n_pcs = 20,
  k = 5
) {
  method <- match.arg(method)
  
  if (inherits(ref_data, "SingleCellExperiment")) {
    ref_data <- SummarizedExperiment::assay(ref_data, "counts")
  }
  if (inherits(sim_data, "SingleCellExperiment")) {
    sim_data <- SummarizedExperiment::assay(sim_data, "counts")
  }
  
  ref_mat <- as.matrix(ref_data)
  sim_mat <- as.matrix(sim_data)
  
  ref_celltypes <- as.character(ref_celltypes)
  sim_celltypes <- as.character(sim_celltypes)
  
  if (ncol(ref_mat) != length(ref_celltypes)) {
    stop("Number of columns in 'ref_data' must match length of 'ref_celltypes'.")
  }
  if (ncol(sim_mat) != length(sim_celltypes)) {
    stop("Number of columns in 'sim_data' must match length of 'sim_celltypes'.")
  }
  
  # Align common features
  common_feats <- intersect(rownames(ref_mat), rownames(sim_mat))
  if (length(common_feats) >= 10) {
    ref_sub <- ref_mat[common_feats, , drop = FALSE]
    sim_sub <- sim_mat[common_feats, , drop = FALSE]
  } else {
    n_feats <- min(nrow(ref_mat), nrow(sim_mat))
    ref_sub <- ref_mat[seq_len(n_feats), , drop = FALSE]
    sim_sub <- sim_mat[seq_len(n_feats), , drop = FALSE]
  }
  
  # Log-normalize (TF-IDF / log CPM equivalent)
  ref_cpm <- log2(sweep(ref_sub, 2, pmax(colSums(ref_sub), 1), "/") * 1e4 + 1)
  sim_cpm <- log2(sweep(sim_sub, 2, pmax(colSums(sim_sub), 1), "/") * 1e4 + 1)
  
  # Joint PCA projection
  comb_mat <- cbind(ref_cpm, sim_cpm)
  var_genes <- apply(comb_mat, 1, stats::var)
  keep_genes <- which(var_genes > 0)
  if (length(keep_genes) > 5) {
    top_g <- order(var_genes[keep_genes], decreasing = TRUE)[seq_len(min(1000, length(keep_genes)))]
    pca_input <- t(comb_mat[keep_genes[top_g], , drop = FALSE])
    pcs <- min(n_pcs, ncol(pca_input) - 1, nrow(pca_input) - 1)
    if (pcs >= 2) {
      pca_res <- stats::prcomp(pca_input, scale. = FALSE, rank. = pcs)
      ref_emb <- pca_res$x[seq_len(ncol(ref_mat)), , drop = FALSE]
      sim_emb <- pca_res$x[(ncol(ref_mat) + 1):nrow(pca_input), , drop = FALSE]
    } else {
      ref_emb <- t(ref_cpm)
      sim_emb <- t(sim_cpm)
    }
  } else {
    ref_emb <- t(ref_cpm)
    sim_emb <- t(sim_cpm)
  }
  
  # Supervised classification
  if (method == "knn") {
    k_use <- min(k, nrow(ref_emb) - 1)
    k_use <- max(1, k_use)
    pred_labels <- as.character(class::knn(
      train = ref_emb,
      test = sim_emb,
      cl = as.factor(ref_celltypes),
      k = k_use
    ))
  } else { # centroid
    classes <- unique(ref_celltypes)
    centroids <- t(vapply(classes, function(cl) {
      colMeans(ref_emb[ref_celltypes == cl, , drop = FALSE])
    }, numeric(ncol(ref_emb))))
    rownames(centroids) <- classes
    
    pred_labels <- vapply(seq_len(nrow(sim_emb)), function(i) {
      diffs <- sweep(centroids, 2, sim_emb[i, ], "-")
      dists <- rowSums(diffs^2)
      classes[which.min(dists)]
    }, character(1))
  }
  
  # Metrics calculation
  all_classes <- union(unique(sim_celltypes), unique(ref_celltypes))
  tbl <- table(factor(pred_labels, levels = all_classes), factor(sim_celltypes, levels = all_classes))
  
  total_cells <- length(sim_celltypes)
  acc <- sum(diag(tbl)) / total_cells
  
  # Macro metrics across all true classes present in test data
  true_classes <- unique(sim_celltypes)
  recalls <- numeric(length(true_classes))
  precisions <- numeric(length(true_classes))
  f1s <- numeric(length(true_classes))
  names(recalls) <- names(precisions) <- names(f1s) <- true_classes
  
  for (cl in true_classes) {
    tp <- tbl[cl, cl]
    fp <- sum(tbl[cl, ]) - tp
    fn <- sum(tbl[, cl]) - tp
    
    prec <- ifelse((tp + fp) > 0, tp / (tp + fp), 0)
    rec <- ifelse((tp + fn) > 0, tp / (tp + fn), 0)
    f1 <- ifelse((prec + rec) > 0, 2 * (prec * rec) / (prec + rec), 0)
    
    precisions[cl] <- prec
    recalls[cl] <- rec
    f1s[cl] <- f1
  }
  
  balanced_acc <- mean(recalls, na.rm = TRUE)
  macro_f1 <- mean(f1s, na.rm = TRUE)
  macro_prec <- mean(precisions, na.rm = TRUE)
  macro_rec <- mean(recalls, na.rm = TRUE)
  
  # Cohen's Kappa
  row_marginals <- rowSums(tbl) / total_cells
  col_marginals <- colSums(tbl) / total_cells
  p_e <- sum(row_marginals * col_marginals)
  kappa_val <- ifelse(p_e < 1, (acc - p_e) / (1 - p_e), 1)
  
  list(
    accuracy = as.numeric(acc),
    balanced_accuracy = as.numeric(balanced_acc),
    macro_f1 = as.numeric(macro_f1),
    macro_precision = as.numeric(macro_prec),
    macro_recall = as.numeric(macro_rec),
    cohen_kappa = as.numeric(kappa_val),
    per_cell_type_f1 = f1s,
    confusion_matrix = tbl,
    predicted_labels = pred_labels
  )
}


