  if (is.null(df) || nrow(df) == 0) {
    if (is.data.frame(df)) df$Score_Norm <- numeric(0)
    return(df)
  }
  df$Score_Norm <- NA_real_
  metric_key <- if ("Property" %in% colnames(df)) {
    paste0(df$Property, "___", df$Metric)
  } else if ("Display_Category" %in% colnames(df)) {
    paste0(df$Display_Category, "___", df$Metric)
  } else {
    df$Metric
  }
  for (m_key in unique(metric_key)) {
    idx     <- metric_key == m_key
    met     <- df$Metric[idx][1]
    raw_val <- df$Score_Raw[idx]
    if (met %in% c("Runs_Statistic", "Between_Dataset_Silh_Global", "Between_Dataset_Silh_2D_Global")) {
      raw_val <- abs(raw_val)
    }
    non_na  <- raw_val[!is.na(raw_val) & is.finite(raw_val)]
    if (length(non_na) == 0) {
      df$Score_Norm[idx] <- NA_real_
      next
    }
    is_hib  <- met %in% .HIGHER_IS_BETTER_METRICS
    min_v   <- min(non_na)
    max_v   <- max(non_na)
    rng     <- max_v - min_v
    if (is.na(rng) || !is.finite(rng)) rng <- 0

    # 2D KS distances (Fasano/Peacock) are not naturally bounded to [0, 1].
    # They can be orders of magnitude larger than 1, so they must remain on the
    # regular min-max scale rather than being coerced into a 0-1 bounded metric.
    is_naturally_01 <- isTRUE(grepl("^(KS|OV|ari|nmi|ami|v_measure|ECDF_DiffArea|accuracy|precision|recall|f1_score|auc|pval_dist_score|F1|Top_DEG_Jaccard)", met, ignore.case = TRUE)) ||
                       isTRUE(min_v >= 0 && max_v <= 1 && rng <= 1)
    if (is_naturally_01) {
      if (is_hib) {
        normed <- pmax(0, pmin(1, raw_val))
      } else {
        normed <- pmax(0, pmin(1, 1 - raw_val))
      }
    } else if (rng < 1e-10) {
      if (is_hib) {
        normed <- if (max_v > 0) pmin(1, raw_val / max(1, max_v)) else rep(0.5, sum(idx))
      } else {
        normed <- 1 / (1 + pmax(0, raw_val))
      }
    } else if (is_hib) {
      normed <- (raw_val - min_v) / rng
    } else {
      normed <- 1 - (raw_val - min_v) / rng
    }
    df$Score_Norm[idx] <- pmax(0, pmin(1, normed))
  }
  still_na <- is.na(df$Score_Norm)
  if (any(still_na))
    df$Score_Norm[still_na] <- pmax(0, pmin(1, df$Score_Raw[still_na]))
  df
}
