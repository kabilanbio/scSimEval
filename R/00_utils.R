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
      n_points <- max(n_ref, n_sim)
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

