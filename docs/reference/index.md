# Package index

## Master Benchmarking Workflows

Unified end-to-end evaluation pipelines executing the 56 curated
evaluation metrics across unimodal and multiomics simulations.

- [`evaluate_simulation_accuracy()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_simulation_accuracy.md)
  : Evaluate Comprehensive Simulation Accuracy (Unimodal Omics Layer)
- [`evaluate_multiomics_accuracy()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_multiomics_accuracy.md)
  : Master Multiomics Simulation Benchmarking Pipeline
- [`evaluate_multiple_datasets()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_multiple_datasets.md)
  : Evaluate Multiple Single-Cell Datasets Simultaneously
- [`benchmark_simulation_suite()`](https://kabilanbio.github.io/scSimEval/reference/benchmark_simulation_suite.md)
  : Benchmark Simulation Methods Across Multiple Single-Cell Datasets
- [`combine_benchmark_results()`](https://kabilanbio.github.io/scSimEval/reference/combine_benchmark_results.md)
  [`append_benchmark_results()`](https://kabilanbio.github.io/scSimEval/reference/combine_benchmark_results.md)
  : Combine and Append Multiple Single-Cell Benchmark Results

## High-Resolution Visualization Suite (600 DPI)

Publication-grade visualization functions designed for multi-metric and
cross-simulator benchmarking.

- [`plot_benchmark_bubble_matrix()`](https://kabilanbio.github.io/scSimEval/reference/plot_benchmark_bubble_matrix.md)
  [`plot_bubble_matrix()`](https://kabilanbio.github.io/scSimEval/reference/plot_benchmark_bubble_matrix.md)
  : Plot Multi-Dimensional Benchmarking Bubble Matrix
- [`plot_evaluation_summary()`](https://kabilanbio.github.io/scSimEval/reference/plot_evaluation_summary.md)
  [`plot_benchmark_summary_bars()`](https://kabilanbio.github.io/scSimEval/reference/plot_evaluation_summary.md)
  [`plot_summary_bars()`](https://kabilanbio.github.io/scSimEval/reference/plot_evaluation_summary.md)
  : Plot Single-Cell Simulator Evaluation Summary
- [`plot_distribution_qc()`](https://kabilanbio.github.io/scSimEval/reference/plot_distribution_qc.md)
  : Plot Single-Cell Summary Distribution Quality (Comparative Overlays)
- [`plot_scalability_benchmark()`](https://kabilanbio.github.io/scSimEval/reference/plot_scalability_benchmark.md)
  : Plot Multi-Dimensional Computational Scalability Benchmark
- [`plot_metric_boxplots()`](https://kabilanbio.github.io/scSimEval/reference/plot_metric_boxplots.md)
  : Plot Comparison Boxplots Across Simulators Plot Benchmark Metric
  Score Distributions Across Simulators
- [`plot_category_metric_bars()`](https://kabilanbio.github.io/scSimEval/reference/plot_category_metric_bars.md)
  : Plot Benchmark Metric Performance Barplots by Category
- [`plot_metric_heatmap()`](https://kabilanbio.github.io/scSimEval/reference/plot_metric_heatmap.md)
  : Plot Multi-Simulator Comparative Metric Heatmap Across Canonical
  Categories
- [`plot_metric_pca()`](https://kabilanbio.github.io/scSimEval/reference/plot_metric_pca.md)
  : Principal Component Analysis (PCA) Dashboard of Benchmark Metrics
  and Methods
- [`plot_metric_mds()`](https://kabilanbio.github.io/scSimEval/reference/plot_metric_mds.md)
  : Multi-Dimensional Scaling (MDS) Ordination of Evaluation Metrics or
  Simulators
- [`plot_individual_metric_bar()`](https://kabilanbio.github.io/scSimEval/reference/plot_individual_metric_bar.md)
  : Plot Individual Metric Performance Barplot
- [`compute_dataset_embeddings()`](https://kabilanbio.github.io/scSimEval/reference/compute_dataset_embeddings.md)
  : Compute Low-Dimensional Embeddings for Reference and Simulated
  Single-Cell Datasets
- [`plot_dataset_embeddings()`](https://kabilanbio.github.io/scSimEval/reference/plot_dataset_embeddings.md)
  : Plot High-Dimensional Dataset Embeddings for Reference and Simulated
  Single-Cell Data
- [`compute_embedding_quality_metrics()`](https://kabilanbio.github.io/scSimEval/reference/compute_embedding_quality_metrics.md)
  : Compute Quantitative Quality Metrics for Low-Dimensional Cell
  Embeddings
- [`compute_method_leaderboard()`](https://kabilanbio.github.io/scSimEval/reference/compute_method_leaderboard.md)
  : Compute Simulator Method Performance Leaderboard
- [`zoom_plot()`](https://kabilanbio.github.io/scSimEval/reference/zoom_plot.md)
  : Interactive Step-by-Step Zoom In and Zoom Out for Plots

## Interactive Graphical Interface (Shiny App Studio)

Interactive Shiny web application providing point-and-click simulation
benchmarking, UMAP/t-SNE embedding explorations, and dynamic data
packaging.

- [`launch_scSimEval_app()`](https://kabilanbio.github.io/scSimEval/reference/launch_scSimEval_app.md)
  : Launch Interactive scSimEval Benchmarking Studio

## Category-Specific Evaluation Suites

Targeted evaluation functions executing specialized metric pipelines.

- [`evaluate_batch_metrics()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_batch_metrics.md)
  : Calculate Batch Integration Metrics
- [`evaluate_clustering_metrics()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_clustering_metrics.md)
  : Full Clustering Performance Evaluation
- [`evaluate_trajectory_metrics()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_trajectory_metrics.md)
  : Full Category VII Trajectory & Lineage Topology Evaluation (dynverse
  Workflow)
- [`evaluate_cross_batch_prediction()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_cross_batch_prediction.md)
  : Evaluate Cross-Batch Cell Classification Transfer F1
- [`evaluate_multiomics_integration()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_multiomics_integration.md)
  : Evaluate Multiomics Integration & Cross-Modal Alignment Fidelity
  (Category VI)
- [`evaluate_deg_fidelity()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_deg_fidelity.md)
  [`evaluate_rna_de_fidelity()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_deg_fidelity.md)
  : Evaluate scRNA-seq Differential Expression (DE) Fidelity
  (Sub-category V-A)
- [`evaluate_predictive_de_model()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_predictive_de_model.md)
  : Predictive Cell Identity Classification Using Simulated DEGs
  (Simpipe Approach 3)
- [`evaluate_zero_probability_curve()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_zero_probability_curve.md)
  : Evaluate Zero-Probability Dropout Curve Fidelity (Splatter &
  ZINB-WaVE)
- [`evaluate_accessibility_sparsity_curve()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_accessibility_sparsity_curve.md)
  : Evaluate Accessibility-Sparsity Curve Concordance (simATAC)
- [`benchmark_resource_usage()`](https://kabilanbio.github.io/scSimEval/reference/benchmark_resource_usage.md)
  : Benchmark Execution Runtime and Memory Usage

## Data Extraction & Empirical Properties

Extraction of feature-level and cell-level empirical properties.

- [`extract_cell_properties()`](https://kabilanbio.github.io/scSimEval/reference/extract_cell_properties.md)
  : Extract Comprehensive Cell-Level Properties
- [`extract_feature_properties()`](https://kabilanbio.github.io/scSimEval/reference/extract_feature_properties.md)
  : Extract Comprehensive Feature-Level (Gene / Peak) Properties
- [`extract_dataset_summary()`](https://kabilanbio.github.io/scSimEval/reference/extract_dataset_summary.md)
  : Extract Dataset Dimensional Properties and Summary Statistics

## Example Multiomics Datasets & Metadata

Pre-packaged synthetic and empirical single-cell multiomics benchmark
datasets and metadata annotation templates.

- [`example_scrna`](https://kabilanbio.github.io/scSimEval/reference/example_scrna.md)
  : Synthetic Example scRNA-seq Benchmarking Dataset
- [`example_scatac`](https://kabilanbio.github.io/scSimEval/reference/example_scatac.md)
  : Synthetic Example scATAC-seq Benchmarking Dataset
- [`example_multiomics`](https://kabilanbio.github.io/scSimEval/reference/example_multiomics.md)
  : Synthetic Example Multiomics Benchmarking Dataset
- [`example_rice_bud`](https://kabilanbio.github.io/scSimEval/reference/example_rice_bud.md)
  : Real Paired Single-Cell Multiomics Benchmark Dataset (Rice Bud,
  Seurat v5 Option B)
- [`example_soybean_root`](https://kabilanbio.github.io/scSimEval/reference/example_soybean_root.md)
  : Real Unpaired Single-Cell Multiomics Benchmark Dataset (Soybean
  Root, Seurat v5 Option B)
- [`example_cell_types`](https://kabilanbio.github.io/scSimEval/reference/example_cell_types.md)
  : Example Cell-Type Annotations Vector
- [`example_batch_info`](https://kabilanbio.github.io/scSimEval/reference/example_batch_info.md)
  : Example Batch Annotations Vector

## Computational Utilities & Fast Kernels

High-performance matrix calculation helpers, input validation, and
internal dimension-reduction utilities.

- [`standardize_benchmark_categories()`](https://kabilanbio.github.io/scSimEval/reference/standardize_benchmark_categories.md)
  : Standardize Category Names in Benchmark Tables
- [`filter_gold_standard_metrics()`](https://kabilanbio.github.io/scSimEval/reference/filter_gold_standard_metrics.md)
  : Filter Benchmark Data to the Curated Evaluation Measures Suite
- [`validate_sc_object()`](https://kabilanbio.github.io/scSimEval/reference/validate_sc_object.md)
  : Validate and Standardize Single-Cell Input Objects
- [`validate_multiomics_input()`](https://kabilanbio.github.io/scSimEval/reference/validate_multiomics_input.md)
  : Validate and Standardize Multiomics Input Objects
- [`calc_outlier_proportion()`](https://kabilanbio.github.io/scSimEval/reference/calc_outlier_proportion.md)
  : Calculate Proportion of Outliers in a Numeric Vector
- [`fast_row_vars()`](https://kabilanbio.github.io/scSimEval/reference/fast_row_vars.md)
  : Fast Row-Wise Variances for Dense or Sparse Single-Cell Matrices
- [`fast_col_vars()`](https://kabilanbio.github.io/scSimEval/reference/fast_col_vars.md)
  : Fast Column-Wise Variances for Dense or Sparse Single-Cell Matrices
- [`fast_sparse_zero_fraction_cell()`](https://kabilanbio.github.io/scSimEval/reference/fast_sparse_zero_fraction_cell.md)
  : Fast Zero Fraction per Cell (Column-wise)
- [`fast_sparse_zero_fraction_feature()`](https://kabilanbio.github.io/scSimEval/reference/fast_sparse_zero_fraction_feature.md)
  : Fast Zero Fraction per Feature (Row-wise)
- [`fast_cross_distance_sq()`](https://kabilanbio.github.io/scSimEval/reference/fast_cross_distance_sq.md)
  : Fast BLAS-Accelerated Cross-Distance Squared Matrix
- [`fast_pca()`](https://kabilanbio.github.io/scSimEval/reference/fast_pca.md)
  : Fast Truncated PCA on Top Highly Variable Features
- [`normalize_counts()`](https://kabilanbio.github.io/scSimEval/reference/normalize_counts.md)
  : Fast Count Normalization (Log2-CPM or TF-IDF)
- [`tfidf_transform()`](https://kabilanbio.github.io/scSimEval/reference/tfidf_transform.md)
  : Perform TF-IDF Normalization for scATAC-seq Peak Matrices
- [`subsample_dist_vector()`](https://kabilanbio.github.io/scSimEval/reference/subsample_dist_vector.md)
  : Subsample Large Vectors for Accelerated Distributional Distance
  Estimation

## Category 1: Global & Marginal Distributions

Evaluation functions assessing univariate and bivariate distribution
divergence across single-cell feature spaces.

- [`calc_ks()`](https://kabilanbio.github.io/scSimEval/reference/calc_ks.md)
  : Calculate Kolmogorov-Smirnov Distance (1D KS)
- [`calc_overlap()`](https://kabilanbio.github.io/scSimEval/reference/calc_overlap.md)
  : Calculate Kernel Density Overlap Value (OV)
- [`calc_wasserstein_1d()`](https://kabilanbio.github.io/scSimEval/reference/calc_wasserstein_1d.md)
  : Calculate 1D Wasserstein Metric / Earth Mover's Distance (WS)
- [`calc_mad()`](https://kabilanbio.github.io/scSimEval/reference/calc_mad.md)
  : Calculate Median Absolute Deviation (MAD)
- [`calc_mae()`](https://kabilanbio.github.io/scSimEval/reference/calc_mae.md)
  : Calculate Mean Absolute Error (MAE)
- [`calc_rmse()`](https://kabilanbio.github.io/scSimEval/reference/calc_rmse.md)
  : Calculate Root Mean Squared Error (RMSE)
- [`calc_bhattacharyya()`](https://kabilanbio.github.io/scSimEval/reference/calc_bhattacharyya.md)
  : Calculate Bhattacharyya Distance (BH)
- [`calc_runs_test()`](https://kabilanbio.github.io/scSimEval/reference/calc_runs_test.md)
  : Calculate Wald-Wolfowitz Runs Test on Pooled Distributions
- [`calc_nn_mismatch()`](https://kabilanbio.github.io/scSimEval/reference/calc_nn_mismatch.md)
  : Calculate Nearest-Neighbor Composition / Local Neighborhood
  Similarity (NN Mismatch)
- [`calc_between_dataset_silhouette()`](https://kabilanbio.github.io/scSimEval/reference/calc_between_dataset_silhouette.md)
  : Calculate Overall Clustering-Like Separation of Datasets via
  Silhouette-Based Measures
- [`calc_kde_ise()`](https://kabilanbio.github.io/scSimEval/reference/calc_kde_ise.md)
  : Calculate 1D Kernel Density Estimation Integrated Squared Error
  (KDE-ISE)
- [`calc_fasano_franceschini()`](https://kabilanbio.github.io/scSimEval/reference/calc_fasano_franceschini.md)
  : Calculate Fasano-Franceschini 2D Kolmogorov-Smirnov Test Statistic
- [`calc_kde_test()`](https://kabilanbio.github.io/scSimEval/reference/calc_kde_test.md)
  : Calculate 2D Bivariate Kernel Density Estimation (KDE) Test
  Statistic
- [`calc_emd_2d()`](https://kabilanbio.github.io/scSimEval/reference/calc_emd_2d.md)
  : Calculate 2D Earth Mover's Distance (2D EMD)
- [`calc_all_univariate_metrics()`](https://kabilanbio.github.io/scSimEval/reference/calc_all_univariate_metrics.md)
  : Compute All Univariate Distance & Accuracy Metrics for a Feature
- [`calc_all_bivariate_metrics()`](https://kabilanbio.github.io/scSimEval/reference/calc_all_bivariate_metrics.md)
  : Compute All Bivariate (2D) Distance & Accuracy Metrics for Joint
  Distributions

## Category 2: Sparsity & Zero-Inflation Dynamics

Assessment of zero-inflation dynamics, dropout probability curves, and
chromatin accessibility sparsity.

- [`calc_zero_probability_curve()`](https://kabilanbio.github.io/scSimEval/reference/calc_zero_probability_curve.md)
  : Fit Zero-Probability Dropout Curve (Splatter & ZINB-WaVE)
- [`calc_excess_zero_weights()`](https://kabilanbio.github.io/scSimEval/reference/calc_excess_zero_weights.md)
  : Calculate Posterior Excess Zero Weights (zingeR & ZINB-WaVE)
- [`calc_accessibility_sparsity_curve()`](https://kabilanbio.github.io/scSimEval/reference/calc_accessibility_sparsity_curve.md)
  : Fit Chromatin Accessibility-Sparsity Polynomial Curve (simATAC)

## Category 3: Cellular Manifold & Clustering Concordance

Concordance of unsupervised cluster definitions, cluster separation, and
latent manifold recovery.

- [`calc_ari()`](https://kabilanbio.github.io/scSimEval/reference/calc_ari.md)
  : Calculate Adjusted Rand Index (ARI)
- [`calc_nmi()`](https://kabilanbio.github.io/scSimEval/reference/calc_nmi.md)
  : Calculate Normalized Mutual Information (NMI)
- [`calc_silhouette()`](https://kabilanbio.github.io/scSimEval/reference/calc_silhouette.md)
  : Calculate Average Silhouette Width (ASW)
- [`calc_davies_bouldin()`](https://kabilanbio.github.io/scSimEval/reference/calc_davies_bouldin.md)
  : Calculate Davies-Bouldin Index (DB)
- [`calc_cdi()`](https://kabilanbio.github.io/scSimEval/reference/calc_cdi.md)
  : Compute Clustering Deviation Index (CDI) for Single-Cell Clustering
  Evaluation
- [`calc_neighborhood_purity()`](https://kabilanbio.github.io/scSimEval/reference/calc_neighborhood_purity.md)
  : Calculate Neighborhood Purity
- [`calc_mmd()`](https://kabilanbio.github.io/scSimEval/reference/calc_mmd.md)
  : Calculate Maximum Mean Discrepancy (MMD) on Single-Cell PCA
  Manifolds
- [`calc_generative_precision_recall()`](https://kabilanbio.github.io/scSimEval/reference/calc_generative_precision_recall.md)
  : Generative Precision and Recall for Single-Cell Manifolds

## Category 4: Technical Variation & Batch Mixing Dynamics

Quantification of technical variation, batch mixing, and local
neighborhood confounder integration.

- [`calc_cms()`](https://kabilanbio.github.io/scSimEval/reference/calc_cms.md)
  : Cell-Specific Mixing Score (CMS)
- [`calc_isi()`](https://kabilanbio.github.io/scSimEval/reference/calc_isi.md)
  : Inverse Simpson Index (ISI) for Batch Mixing

## Category 5: Biological Signal & Differential Fidelity

Evaluation of differentially expressed genes (Simpipe framework; Duo et
al., Genome Biology 2024) and differentially accessible regions via
pseudobulk edgeR-QLF (Libra framework; Squair et al., Nature
Communications 2021; Teo et al., Nature Communications 2024).

- [`run_edger_qlf_detrate()`](https://kabilanbio.github.io/scSimEval/reference/run_edger_qlf_detrate.md)
  : Run edgeR QLF Differential Expression with Cellular Detection Rate
  (edgeRQLFDetRate)
- [`calc_simpipe_deg_ratio()`](https://kabilanbio.github.io/scSimEval/reference/calc_simpipe_deg_ratio.md)
  : Calculate Simpipe Simulated DEG Ratio and True DEG Recovery
  (Approach 1)
- [`calc_simpipe_distribution_score()`](https://kabilanbio.github.io/scSimEval/reference/calc_simpipe_distribution_score.md)
  : Calculate Simpipe P-value Uniformity Distribution Score (Approach 2)
- [`evaluate_predictive_de_model()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_predictive_de_model.md)
  : Predictive Cell Identity Classification Using Simulated DEGs
  (Simpipe Approach 3)
- [`evaluate_deg_fidelity()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_deg_fidelity.md)
  [`evaluate_rna_de_fidelity()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_deg_fidelity.md)
  : Evaluate scRNA-seq Differential Expression (DE) Fidelity
  (Sub-category V-A)
- [`binarize_matrix()`](https://kabilanbio.github.io/scSimEval/reference/binarize_matrix.md)
  : Binarize Single-Cell Count Matrix (scATAC-seq Feature Binarization)
- [`to_pseudobulk()`](https://kabilanbio.github.io/scSimEval/reference/to_pseudobulk.md)
  : Construct Pseudobulk Count Matrices from Single-Cell Data (Libra
  Engine)
- [`run_edger_qlf_da()`](https://kabilanbio.github.io/scSimEval/reference/run_edger_qlf_da.md)
  : Run edgeR-QLF Differential Accessibility (DA) Analysis on Pseudobulk
  Data (Libra Engine)
- [`evaluate_atac_da_fidelity()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_atac_da_fidelity.md)
  : Evaluate scATAC-seq Differential Accessibility (DA) Fidelity
  (Sub-category V-B)

## Category 6: Multiomics Integration & Alignment

Gold-standard multiomics integration, batch mixing, and biological
conservation metrics (scib, scMultiSim, scCross, Matilda, scMoMtF).

- [`calc_batch_asw()`](https://kabilanbio.github.io/scSimEval/reference/calc_batch_asw.md)
  : Calculate Batch Average Silhouette Width (Batch ASW)
- [`calc_cell_type_asw()`](https://kabilanbio.github.io/scSimEval/reference/calc_cell_type_asw.md)
  : Calculate Cell-Type Average Silhouette Width (Cell-Type ASW)
- [`calc_graph_connectivity()`](https://kabilanbio.github.io/scSimEval/reference/calc_graph_connectivity.md)
  : Calculate KNN Graph Connectivity
- [`calc_integration_ari()`](https://kabilanbio.github.io/scSimEval/reference/calc_integration_ari.md)
  : Calculate Integration Adjusted Rand Index (ARI)
- [`calc_integration_nmi()`](https://kabilanbio.github.io/scSimEval/reference/calc_integration_nmi.md)
  : Calculate Integration Normalized Mutual Information (NMI)
- [`evaluate_multiomics_integration()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_multiomics_integration.md)
  : Evaluate Multiomics Integration & Cross-Modal Alignment Fidelity
  (Category VI)

## Category 7: Trajectory & Lineage Topology

HVG-filtered dynverse trajectory inference benchmarking (Saelens et al.,
Nature Biotechnology 2019) across milestone networks, geodesic
distances, branch assignments, and Random Forest feature importances.

- [`preprocess_trajectory_hvgs()`](https://kabilanbio.github.io/scSimEval/reference/preprocess_trajectory_hvgs.md)
  : Preprocess and Select Highly Variable Genes for Trajectory
  Evaluation (Saelens et al., 2019)
- [`infer_dynverse_trajectory()`](https://kabilanbio.github.io/scSimEval/reference/infer_dynverse_trajectory.md)
  : Infer Docker-Free dynwrap Milestone Trajectory from Normalized
  Expression
- [`calc_dynverse_him()`](https://kabilanbio.github.io/scSimEval/reference/calc_dynverse_him.md)
  : Compute dynverse Hamming-Ipsen-Mikhailov (HIM) Trajectory Topology
  Similarity
- [`calc_dynverse_edge_flip()`](https://kabilanbio.github.io/scSimEval/reference/calc_dynverse_edge_flip.md)
  : Compute dynverse Edge Flip Trajectory Topology Similarity
- [`calc_dynverse_geodesic_correlation()`](https://kabilanbio.github.io/scSimEval/reference/calc_dynverse_geodesic_correlation.md)
  : Compute dynverse Geodesic Distance Ordering Correlation
- [`calc_dynverse_f1_branches()`](https://kabilanbio.github.io/scSimEval/reference/calc_dynverse_f1_branches.md)
  : Compute dynverse Branch Assignment Overlap (F1_branches)
- [`calc_wcor_features()`](https://kabilanbio.github.io/scSimEval/reference/calc_wcor_features.md)
  : Compute Weighted Feature Importance Correlation (wcor_features,
  Saelens et al., 2019)
- [`evaluate_trajectory_metrics()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_trajectory_metrics.md)
  : Full Category VII Trajectory & Lineage Topology Evaluation (dynverse
  Workflow)

## Category 8: Computational Scalability & Efficiency

Core-normalized benchmarking of CPU execution work (Core-Seconds) and
per-core peak RAM footprint (MiB/core and GiB/core).

- [`benchmark_resource_usage()`](https://kabilanbio.github.io/scSimEval/reference/benchmark_resource_usage.md)
  : Benchmark Execution Runtime and Memory Usage
- [`compute_core_seconds()`](https://kabilanbio.github.io/scSimEval/reference/compute_core_seconds.md)
  : Compute Core-Seconds Runtime Metric
- [`compute_memory_per_core()`](https://kabilanbio.github.io/scSimEval/reference/compute_memory_per_core.md)
  : Compute Peak RAM per Allocated CPU Core
