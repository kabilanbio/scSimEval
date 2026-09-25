test_that("plot_distribution_qc returns valid ggplot with single and multiple simulators", {
  data("example_scrna", package = "scSimEval", envir = environment())
  
  p1 <- plot_distribution_qc(example_scrna$ref, example_scrna$sim)
  expect_s3_class(p1, "ggplot")

  sim_list <- list(
    "SimA" = example_scrna$sim,
    "SimB" = example_scrna$sim
  )
  p2 <- plot_distribution_qc(example_scrna$ref, sim_list)
  expect_s3_class(p2, "ggplot")
})

test_that("plot_metric_boxplots returns a valid ggplot object", {
  demo_df <- data.frame(
    Method = rep(c("Splatter", "scDesign3", "SymSim"), each = 4),
    Category = rep(c("(I) Distributional Properties", "(III) Cellular Structure & Concordance"), 6),
    Metric = rep(c("KS", "Wasserstein", "ARI", "NMI"), 3),
    Score = runif(12, 0.4, 0.95)
  )
  p <- plot_metric_boxplots(demo_df)
  expect_s3_class(p, "ggplot")
})

test_that("plot_scalability_benchmark returns valid ggplot and patchwork objects", {
  demo_sc <- data.frame(
    Method = rep(c("Splatter", "scDesign3", "SymSim"), each = 2),
    Category = "(VIII) Computational Scalability",
    Metric = rep(c("elapsed_time_seconds", "peak_memory_mb"), 3),
    Score = c(30, 700, 65, 1300, 45, 800)
  )
  p_run <- plot_scalability_benchmark(demo_sc, type = "runtime")
  expect_s3_class(p_run, "ggplot")

  p_mem <- plot_scalability_benchmark(demo_sc, type = "memory")
  expect_s3_class(p_mem, "ggplot")

  p_comp <- plot_scalability_benchmark(demo_sc, type = "composite", layout = "4panel")
  expect_true(inherits(p_comp, "ggplot") || inherits(p_comp, "patchwork"))
})

test_that("plot_metric_mds and plot_metric_pca return valid ggplot objects", {
  demo_df <- data.frame(
    Method = rep(c("Splatter", "scDesign3", "SymSim", "dyngen"), each = 3),
    Category = "(I) Distributional Properties",
    Metric = rep(c("KS", "Wasserstein", "MAD"), 4),
    Score = runif(12, 0.5, 0.95)
  )
  p_mds <- plot_metric_mds(demo_df)
  expect_s3_class(p_mds, "ggplot")

  p_pca <- plot_metric_pca(demo_df)
  expect_s3_class(p_pca, "ggplot")
})

test_that("plot_benchmark_summary returns a valid ggplot object", {
  data("example_multiomics", package = "scSimEval", envir = environment())
  
  res <- evaluate_multiomics_accuracy(
    ref_multi = example_multiomics$ref_multi,
    sim_multi = example_multiomics$sim_multi,
    verbose = FALSE
  )
  
  p <- plot_benchmark_summary(res$benchmark_summary_table, metric_type = "KS")
  expect_s3_class(p, "ggplot")
})

test_that("plot_trajectory_comparison returns a valid ggplot object", {
  data("example_scrna", package = "scSimEval", envir = environment())
  
  p <- plot_trajectory_comparison(example_scrna$ref, example_scrna$sim)
  expect_s3_class(p, "ggplot")
})

test_that("plot_cross_modal_coupling returns a valid ggplot object", {
  data("example_multiomics", package = "scSimEval", envir = environment())
  
  p <- plot_cross_modal_coupling(
    ref_rna = example_multiomics$ref_multi$rna,
    ref_atac = example_multiomics$ref_multi$atac,
    sim_rna = example_multiomics$sim_multi$rna,
    sim_atac = example_multiomics$sim_multi$atac,
    n_features = 10
  )
  expect_s3_class(p, "ggplot")
})

test_that("plot_metric_heatmap returns a valid ggplot object", {
  data("example_multiomics", package = "scSimEval", envir = environment())
  
  res <- evaluate_multiomics_accuracy(
    ref_multi = example_multiomics$ref_multi,
    sim_multi = example_multiomics$sim_multi,
    verbose = FALSE
  )
  
  sim_list <- list(
    Sim_Tool_A = res$benchmark_summary_table,
    Sim_Tool_B = res$benchmark_summary_table
  )
  
  p <- plot_metric_heatmap(sim_list, metric_name = "KS", top_n_properties = 8)
  expect_s3_class(p, "ggplot")
})

test_that("plot_benchmark_bubble_matrix returns a valid ggplot object with data.frame and list inputs", {
  demo_df <- data.frame(
    Method = rep(c("Splat", "scDesign3", "SymSim", "dyngen"), each = 4),
    Category = rep(c("Accuracy", "Accuracy", "Cellular Structure", "Scalability"), 4),
    Metric = rep(c("KS Distance", "Wasserstein", "Silhouette ASW", "CPU Time"), 4),
    Score = c(0.8, 0.7, 0.9, 0.6, 0.5, 0.8, 0.6, 0.5, 0.8, 0.7, 0.6, 0.5, 0.7, 0.6, 0.9, 0.5)
  )
  
  classes <- list(
    "Class 1" = c("Splat", "scDesign3"),
    "Class 2" = c("SymSim", "dyngen")
  )
  
  p1 <- plot_benchmark_bubble_matrix(demo_df, method_classes = classes)
  expect_s3_class(p1, "ggplot")
  
  p2 <- plot_bubble_matrix(demo_df)
  expect_s3_class(p2, "ggplot")
})

