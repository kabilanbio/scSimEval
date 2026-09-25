## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  warning = FALSE,
  message = FALSE
)


## ----workflow-diagram, echo=FALSE, fig.align='center', out.width='100%', fig.cap='**Figure 1: End-to-end benchmarking workflow of scSimEval.** The framework takes real and simulated single-cell multiomics data, evaluates 62 metrics across 8 foundational categories, and produces standardized scores, comprehensive figures, and ranking leaderboards.'----
knitr::include_graphics("figures/workflow_diagram.png")


## ----score-norm-workflow, echo=FALSE, fig.align='center', out.width='100%', fig.cap='**Figure 2: Score normalization and visual mapping workflow in scSimEval.** Box 1 shows the two-step normalization process: direction inversion for lower-is-better metrics followed by min-max scaling to a standard 0 to 1 range. Box 2 illustrates the visual mapping into bubble matrices, where bubble size reflects simulation quality and top performers are highlighted with bold squares.'----
knitr::include_graphics("figures/score_normalization_workflow.png")


## ----setup, include=FALSE-----------------------------------------------------
library(scSimEval)


## ----load-data----------------------------------------------------------------
# Load example datasets
data(example_scrna)
data(example_scatac)
data(example_multiomics)

# Check dimensions
cat("Real scRNA-seq dimensions:", dim(example_scrna$ref), "\n")
cat("Simulated scRNA-seq dimensions:", dim(example_scrna$sim), "\n")
cat("Cell types present:", levels(example_scrna$cell_types), "\n")


## ----unimodal-eval------------------------------------------------------------
# Evaluate unimodal accuracy on scRNA-seq
unimodal_res <- evaluate_simulation_accuracy(
  ref_data = example_scrna$ref,
  sim_data = example_scrna$sim,
  compute_bivariate = FALSE,
  verbose = FALSE
)

# Preview library size distribution metrics
lib_summary <- subset(unimodal_res$metrics_summary_table, Property == "library_size")
knitr::kable(head(lib_summary, 8), digits = 4, caption = "Library Size Statistical Distance Measures")


## ----plot-dist-qc-fig, echo=FALSE, fig.align='center', out.width='100%', fig.cap='**Figure 3: 14-panel comparative distribution QC layout.** Compares real biological data (warm brick red) and simulated data (steel blue) across single-cell properties, bivariate relationships, and biological signal retention.'----
knitr::include_graphics("figures/comparative_distribution_qc.png")


## ----clustering-eval----------------------------------------------------------
clust_res <- evaluate_clustering_metrics(
  data = example_scrna$sim,
  cell_types = example_scrna$cell_types,
  ref_data = example_scrna$ref
)

cat("Average Silhouette Width:", round(clust_res$silhouette, 4), "\n")
cat("Davies-Bouldin Index:", round(clust_res$davies_bouldin, 4), "\n")
cat("Adjusted Rand Index (ARI):", round(clust_res$ARI, 4), "\n")
cat("Normalized Mutual Information (NMI):", round(clust_res$NMI, 4), "\n")


## ----batch-eval---------------------------------------------------------------
batch_res <- evaluate_batch_metrics(
  data = example_scrna$sim,
  batch_info = example_scrna$batch_info,
  cell_types = example_scrna$cell_types,
  verbose = FALSE
)

cat("Batch Shannon Entropy:", round(batch_res$shannon_entropy, 4), "\n")
cat("PC Regression R2 (PCR):", round(batch_res$pcr_r2, 4), "\n")
cat("Cross-Batch Transfer Accuracy:", round(batch_res$cross_batch_accuracy, 4), "\n")


## ----trajectory-eval----------------------------------------------------------
# Automatically infer trajectories and evaluate fidelity
traj_res <- evaluate_trajectory_metrics(
  ref_data = example_scrna$ref,
  sim_data = example_scrna$sim,
  cell_types_ref = example_scrna$cell_types,
  cell_types_sim = example_scrna$cell_types
)

cat("Pseudotime Spearman Correlation:", round(traj_res$pseudotime_correlation, 4), "\n")
cat("Lineage Tree Branch Height RMSE:", round(as.numeric(traj_res$tree_height_rmse), 4), "\n")


## ----plot-scalability-fig, echo=FALSE, fig.align='center', out.width='100%', fig.cap='**Figure 4: Computational scalability benchmark suite.** 6-panel dashboard comparing execution runtime, peak RAM, runtime-memory tradeoff, CPU efficiency, cost footprint, and throughput across simulation methods.'----
knitr::include_graphics("figures/scalability_benchmark.png")


## ----master-pipeline----------------------------------------------------------
master_eval <- evaluate_multiomics_accuracy(
  ref_multi = example_multiomics$ref_multi,
  sim_multi = example_multiomics$sim_multi,
  cell_types = example_multiomics$cell_types,
  batch_info = example_multiomics$batch_info,
  cpu_time = example_multiomics$resource_stats$cpu_time,
  memory_mb = example_multiomics$resource_stats$memory_mb,
  system_time = example_multiomics$resource_stats$system_time,
  elapsed_time = example_multiomics$resource_stats$elapsed_time,
  verbose = FALSE
)

summary_df <- master_eval$benchmark_summary_table
cat("Total evaluated metric instances:", nrow(summary_df), "\n")
table(summary_df$Category)


## ----plot-bubble-matrix-fig, echo=FALSE, fig.align='center', out.width='100%', fig.cap='**Figure 5: Flagship multi-dimensional benchmarking bubble matrix.** Ranks simulation methods top-to-bottom across the 8 evaluation categories. Bubble size reflects standardized fidelity score ($0.00$ to $1.00$), and top performers are highlighted with bold squares.'----
knitr::include_graphics("figures/benchmark_bubble_matrix.png")


## ----eval-summary-bars, echo=FALSE, fig.align='center', out.width='100%', fig.cap='**Figure 6: Executive evaluation summary horizontal bar matrix.** Ranks simulators across the 8 categories and displays overall composite performance with exact score annotations.'----
knitr::include_graphics("figures/evaluation_summary_bars.png")


## ----plot-boxplots-fig, echo=FALSE, fig.align='center', out.width='100%', fig.cap='**Figure 7: Metric score distributions across the 8 categories.** Displays standardized scores ($[0, 1]$, higher is better) with individual data points and boxplots across candidate simulators.'----
knitr::include_graphics("figures/metric_boxplots.png")


## ----plot-heatmap-fig, echo=FALSE, fig.align='center', out.width='88%', fig.cap='**Figure 8: Complete benchmark metric heatmap.** Displays exact unnormalized raw scores in bold text inside every cell, grouped cleanly across the 8 evaluation categories.'----
knitr::include_graphics("figures/metric_heatmap.png")


## ----plot-pca-cat1, echo=FALSE, fig.align='center', out.width='85%', fig.cap='**Figure 9: Principal Component Analysis (PCA) ordination of simulation methods.** Projects simulators based on their metric profiles, showing global affinities and key discriminating metric vectors.'----
knitr::include_graphics("figures/individual_category_plots/pca_cat1_distribution.png")


## ----multi-dataset-eval-------------------------------------------------------
# Set up a collection of datasets
dataset_collection <- list(
  "Method 1-scRNA-seq" = list(ref = example_scrna$ref, sim = example_scrna$sim),
  "Method 1-scATAC-seq" = list(ref = example_scatac$ref, sim = example_scatac$sim),
  "Method 2-scRNA-seq" = list(ref = example_scrna$ref, sim = example_scrna$sim),
  "Method 2-scATAC-seq" = list(ref = example_scatac$ref, sim = example_scatac$sim)
)

# Run consolidated benchmarking in one step
consolidated_results <- evaluate_multiple_datasets(
  datasets = dataset_collection,
  pair_by_prefix = TRUE,
  compute_bivariate = FALSE,
  verbose = FALSE
)

# Preview consolidated results overview
knitr::kable(consolidated_results$dataset_overview, digits = 4, caption = "Consolidated Benchmark Overview")


## ----session-info-------------------------------------------------------------
sessionInfo()

