test_that("plot_deg_bubble_matrix generates valid ggplot object from df and list", {
  # Mock data
  df <- data.frame(
    Method = rep(c("Method_A", "Method_B"), each = 3),
    Framework = rep(c("Simpipe (Duo et al., 2024)", "SimBench (Cao et al., 2021)", "Shaky Foundations (Crowell et al., 2023)"), 2),
    Metric = rep(c("DEG_Ratio", "Log2FC_Pearson_Corr", "Silhouette_Sim"), 2),
    Value = c(0.9, 0.85, 0.5, 1.1, 0.75, 0.45),
    stringsAsFactors = FALSE
  )
  
  p1 <- plot_deg_bubble_matrix(df)
  expect_s3_class(p1, "ggplot")
  
  # Test with evaluate_deg_fidelity mock output
  mock_res <- list(
    "Sim_1" = list(
      deg_summary_table = data.frame(
        Framework = c("Simpipe (Duo et al., 2024)", "SimBench (Cao et al., 2021)"),
        Metric = c("DEG_Ratio", "Log2FC_Pearson_Corr"),
        Value = c(0.95, 0.88),
        stringsAsFactors = FALSE
      )
    ),
    "Sim_2" = list(
      deg_summary_table = data.frame(
        Framework = c("Simpipe (Duo et al., 2024)", "SimBench (Cao et al., 2021)"),
        Metric = c("DEG_Ratio", "Log2FC_Pearson_Corr"),
        Value = c(0.80, 0.70),
        stringsAsFactors = FALSE
      )
    )
  )
  
  p2 <- plot_deg_bubble_matrix(mock_res)
  expect_s3_class(p2, "ggplot")
})
