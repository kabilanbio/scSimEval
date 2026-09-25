# Package index

## Master Benchmarking Workflows

Unified end-to-end evaluation pipelines executing the 62 canonical
metrics across unimodal and multiomics simulations.

- [`evaluate_simulation_accuracy()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_simulation_accuracy.md)
  : Evaluate Comprehensive Simulation Accuracy (Unimodal Omics Layer)
- [`evaluate_multiomics_accuracy()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_multiomics_accuracy.md)
  : Master Multiomics Simulation Benchmarking Pipeline
- [`evaluate_multiple_datasets()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_multiple_datasets.md)
  : Evaluate Multiple Single-Cell Datasets Simultaneously

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
- [`plot_metric_heatmap()`](https://kabilanbio.github.io/scSimEval/reference/plot_metric_heatmap.md)
  : Plot Multi-Simulator Comparative Metric Heatmap Across Canonical
  Categories
- [`plot_metric_pca()`](https://kabilanbio.github.io/scSimEval/reference/plot_metric_pca.md)
  : Principal Component Analysis (PCA) Dashboard of Benchmark Metrics
  and Methods
- [`plot_metric_mds()`](https://kabilanbio.github.io/scSimEval/reference/plot_metric_mds.md)
  : Multi-Dimensional Scaling (MDS) Ordination of Evaluation Metrics or
  Simulators

## Interactive Graphical Interface (Shiny App GUI)

Interactive Shiny web application providing point-and-click simulation
benchmarking and dynamic exploration.

- [`launch_scSimEval_app()`](https://kabilanbio.github.io/scSimEval/reference/launch_scSimEval_app.md)
  : Launch Interactive scSimEval Benchmarking Studio

## Category-Specific Evaluation Suites

Targeted evaluation functions executing specialized metric pipelines.

- [`evaluate_batch_metrics()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_batch_metrics.md)
  : Calculate Batch Integration Metrics
- [`evaluate_clustering_metrics()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_clustering_metrics.md)
  : Full Clustering Performance Evaluation
- [`evaluate_trajectory_metrics()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_trajectory_metrics.md)
  : Full Trajectory Accuracy Evaluation
- [`evaluate_cross_batch_prediction()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_cross_batch_prediction.md)
  : Evaluate Cross-Batch Cell Classification Transfer Accuracy
- [`evaluate_cross_modal_prediction()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_cross_modal_prediction.md)
  : Cross-Modal Cell-Type Label Transfer Accuracy
- [`evaluate_deg_fidelity()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_deg_fidelity.md)
  : Evaluate Differentially Expressed Gene (DEG) Fidelity
- [`evaluate_predictive_de_model()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_predictive_de_model.md)
  : Predictive Cell Identity Classification Using DE / Top Variable
  Features
- [`evaluate_simbench_signals()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_simbench_signals.md)
  : Evaluate the 5 SimBench Biological Signal Proportions
- [`evaluate_epigenomic_annotation()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_epigenomic_annotation.md)
  : Evaluate Epigenomic Supervised Cell-Type Annotation (EpiAnno /
  SCAN-ATAC-Sim)
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

## Example Multiomics Datasets

Pre-packaged synthetic and empirical single-cell multiomics benchmark
datasets.

- [`example_scrna`](https://kabilanbio.github.io/scSimEval/reference/example_scrna.md)
  : Synthetic Example scRNA-seq Benchmarking Dataset
- [`example_scatac`](https://kabilanbio.github.io/scSimEval/reference/example_scatac.md)
  : Synthetic Example scATAC-seq Benchmarking Dataset
- [`example_multiomics`](https://kabilanbio.github.io/scSimEval/reference/example_multiomics.md)
  : Synthetic Example Multiomics Benchmarking Dataset

## Category 1: Global & Marginal Distributions

Evaluation functions assessing univariate and bivariate distribution
divergence.

- [`calc_ks()`](https://kabilanbio.github.io/scSimEval/reference/calc_ks.md)
  : Calculate Kolmogorov-Smirnov Distance (1D KS)
- [`calc_wasserstein_1d()`](https://kabilanbio.github.io/scSimEval/reference/calc_wasserstein_1d.md)
  : Calculate 1D Wasserstein Metric / Earth Mover's Distance (WS)
  Integrated from HelenaLC/simulation-comparison.
- [`calc_mad()`](https://kabilanbio.github.io/scSimEval/reference/calc_mad.md)
  : Calculate Median Absolute Deviation (MAD)
- [`calc_mae()`](https://kabilanbio.github.io/scSimEval/reference/calc_mae.md)
  : Calculate Mean Absolute Error (MAE)
- [`calc_rmse()`](https://kabilanbio.github.io/scSimEval/reference/calc_rmse.md)
  : Calculate Root Mean Squared Error (RMSE)
- [`calc_overlap()`](https://kabilanbio.github.io/scSimEval/reference/calc_overlap.md)
  : Calculate Distribution Overlapping Index (OV)
- [`calc_bhattacharyya()`](https://kabilanbio.github.io/scSimEval/reference/calc_bhattacharyya.md)
  : Calculate Bhattacharyya Distance (BH)
- [`calc_ecdf_diffarea()`](https://kabilanbio.github.io/scSimEval/reference/calc_ecdf_diffarea.md)
  : Calculate Area Between Empirical Cumulative Distribution Functions
  (eCDFs)
- [`calc_runs_test()`](https://kabilanbio.github.io/scSimEval/reference/calc_runs_test.md)
  : Calculate Wald-Wolfowitz Runs Test on Pooled Distributions
- [`calc_nn_mismatch()`](https://kabilanbio.github.io/scSimEval/reference/calc_nn_mismatch.md)
  : Calculate Nearest-Neighbor Label Mismatch Proportion
- [`calc_between_dataset_silhouette()`](https://kabilanbio.github.io/scSimEval/reference/calc_between_dataset_silhouette.md)
  : Calculate Between-Dataset Silhouette Width
- [`calc_fasano_franceschini()`](https://kabilanbio.github.io/scSimEval/reference/calc_fasano_franceschini.md)
  : Calculate Fasano-Franceschini 2D Kolmogorov-Smirnov Test Statistic
  Integrated from simpipe.
- [`calc_peacock_2d()`](https://kabilanbio.github.io/scSimEval/reference/calc_peacock_2d.md)
  : Calculate Peacock 2D Kolmogorov-Smirnov Test Statistic Integrated
  from HelenaLC/simulation-comparison.
- [`calc_kde_test()`](https://kabilanbio.github.io/scSimEval/reference/calc_kde_test.md)
  : Calculate 2D Bivariate Kernel Density Estimation (KDE) Test
  Statistic Integrated from simpipe & SimBench.
- [`calc_emd_2d()`](https://kabilanbio.github.io/scSimEval/reference/calc_emd_2d.md)
  : Calculate 2D Earth Mover's Distance (2D EMD) Integrated from
  HelenaLC/simulation-comparison.
- [`calc_mmd()`](https://kabilanbio.github.io/scSimEval/reference/calc_mmd.md)
  : Maximum Mean Discrepancy (MMD) with Gaussian RBF Kernel
- [`calc_frechet_singlecell_distance()`](https://kabilanbio.github.io/scSimEval/reference/calc_frechet_singlecell_distance.md)
  : Fréchet Single-Cell Distance (FSD)
- [`calc_model_aic_bic()`](https://kabilanbio.github.io/scSimEval/reference/calc_model_aic_bic.md)
  : Compute Model Information Criteria (AIC and BIC)
- [`calc_marginal_aic_bic()`](https://kabilanbio.github.io/scSimEval/reference/calc_marginal_aic_bic.md)
  : Marginal Model Goodness of Fit and Information Criteria for
  Single-Cell Simulators
- [`calc_likelihood_ratio_test()`](https://kabilanbio.github.io/scSimEval/reference/calc_likelihood_ratio_test.md)
  : Perform Likelihood Ratio Test for Comparing Nested Single-Cell
  Simulation Models
- [`calc_outlier_proportion()`](https://kabilanbio.github.io/scSimEval/reference/calc_outlier_proportion.md)
  : Calculate Proportion of Outliers in a Numeric Vector
- [`calc_all_univariate_metrics()`](https://kabilanbio.github.io/scSimEval/reference/calc_all_univariate_metrics.md)
  : Compute All Univariate Distance & Accuracy Metrics for a Feature
- [`calc_all_bivariate_metrics()`](https://kabilanbio.github.io/scSimEval/reference/calc_all_bivariate_metrics.md)
  : Compute All Bivariate (2D) Distance & Accuracy Metrics for Joint
  Distributions

## Category 2: Cellular & Feature-Level Properties

Assessment of zero-inflation dynamics, library size distributions, and
kinetic noise decomposition.

- [`calc_zero_probability_curve()`](https://kabilanbio.github.io/scSimEval/reference/calc_zero_probability_curve.md)
  : Fit Zero-Probability Dropout Curve (Splatter & ZINB-WaVE)
- [`calc_excess_zero_weights()`](https://kabilanbio.github.io/scSimEval/reference/calc_excess_zero_weights.md)
  : Calculate Posterior Excess Zero Weights (zingeR & ZINB-WaVE)
- [`calc_kinetic_noise_decomposition()`](https://kabilanbio.github.io/scSimEval/reference/calc_kinetic_noise_decomposition.md)
  : Kinetic Noise Decomposition (SymSim & Elowitz et al.)
- [`calc_accessibility_sparsity_curve()`](https://kabilanbio.github.io/scSimEval/reference/calc_accessibility_sparsity_curve.md)
  : Fit Chromatin Accessibility-Sparsity Polynomial Curve (simATAC)

## Category 3: Cell Identity & Clustering Concordance

Concordance of unsupervised cluster definitions, cluster separation, and
latent manifold recovery.

- [`calc_ari()`](https://kabilanbio.github.io/scSimEval/reference/calc_ari.md)
  : Calculate Adjusted Rand Index (ARI)
- [`calc_nmi()`](https://kabilanbio.github.io/scSimEval/reference/calc_nmi.md)
  : Calculate Normalized Mutual Information (NMI)
- [`calc_ami()`](https://kabilanbio.github.io/scSimEval/reference/calc_ami.md)
  : Calculate Adjusted Mutual Information (AMI)
- [`calc_clustering_accuracy()`](https://kabilanbio.github.io/scSimEval/reference/calc_clustering_accuracy.md)
  : Calculate Clustering Accuracy (ACC)
- [`calc_fmi()`](https://kabilanbio.github.io/scSimEval/reference/calc_fmi.md)
  : Calculate Fowlkes-Mallows Index (FMI)
- [`calc_silhouette()`](https://kabilanbio.github.io/scSimEval/reference/calc_silhouette.md)
  : Calculate Average Silhouette Width (ASW)
- [`calc_dunn()`](https://kabilanbio.github.io/scSimEval/reference/calc_dunn.md)
  : Calculate Dunn Index
- [`calc_connectivity()`](https://kabilanbio.github.io/scSimEval/reference/calc_connectivity.md)
  : Calculate Cluster Connectivity
- [`calc_davies_bouldin()`](https://kabilanbio.github.io/scSimEval/reference/calc_davies_bouldin.md)
  : Calculate Davies-Bouldin Index (DB)
- [`calc_calinski_harabasz()`](https://kabilanbio.github.io/scSimEval/reference/calc_calinski_harabasz.md)
  : Calculate Calinski-Harabasz Index (CH)
- [`calc_homogeneity()`](https://kabilanbio.github.io/scSimEval/reference/calc_homogeneity.md)
  : Calculate Homogeneity Score
- [`calc_completeness()`](https://kabilanbio.github.io/scSimEval/reference/calc_completeness.md)
  : Calculate Completeness Score
- [`calc_v_measure()`](https://kabilanbio.github.io/scSimEval/reference/calc_v_measure.md)
  : Calculate V-Measure
- [`calc_homogeneity_completeness_v_measure()`](https://kabilanbio.github.io/scSimEval/reference/calc_homogeneity_completeness_v_measure.md)
  : Calculate Homogeneity, Completeness, and V-Measure
- [`calc_neighborhood_purity()`](https://kabilanbio.github.io/scSimEval/reference/calc_neighborhood_purity.md)
  : Calculate Neighborhood Purity
- [`calc_cdi()`](https://kabilanbio.github.io/scSimEval/reference/calc_cdi.md)
  : Compute Clustering Deviation Index (CDI) for Single-Cell Clustering
  Evaluation
- [`calc_generative_precision_recall()`](https://kabilanbio.github.io/scSimEval/reference/calc_generative_precision_recall.md)
  : Generative Precision and Recall for Single-Cell Manifolds

## Category 4: Batch Effect & Confounder Dynamics

Quantification of technical variation, batch mixing, and local density
factor conservation.

- [`calc_cms()`](https://kabilanbio.github.io/scSimEval/reference/calc_cms.md)
  : Cell-Specific Mixing Score (CMS)
- [`calc_isi()`](https://kabilanbio.github.io/scSimEval/reference/calc_isi.md)
  : Inverse Simpson Index (ISI) for Batch Mixing
- [`calc_seurat_mixing_metric()`](https://kabilanbio.github.io/scSimEval/reference/calc_seurat_mixing_metric.md)
  : Seurat Mixing Metric
- [`calc_ldf()`](https://kabilanbio.github.io/scSimEval/reference/calc_ldf.md)
  : Local Density Factor (LDF)
- [`calc_ldf_diff()`](https://kabilanbio.github.io/scSimEval/reference/calc_ldf_diff.md)
  : Local Density Differences (ldfDiff)
- [`calc_local_structure_metric()`](https://kabilanbio.github.io/scSimEval/reference/calc_local_structure_metric.md)
  : Local Structure Preservation Metric

## Category 5: Differential Expression & Signal Preservation

Recovery of biological contrast, differential gene detection accuracy,
and cell-cycle phase fidelity.

- [`calc_signal_de()`](https://kabilanbio.github.io/scSimEval/reference/calc_signal_de.md)
  : Detect Differential Expression (DE) Genes (Mean Shift via Limma or
  T-Test)
- [`calc_signal_dv()`](https://kabilanbio.github.io/scSimEval/reference/calc_signal_dv.md)
  : Detect Differential Variability (DV) Genes (Bartlett Test)
- [`calc_signal_dd()`](https://kabilanbio.github.io/scSimEval/reference/calc_signal_dd.md)
  : Detect Differential Distribution (DD) Genes (Kolmogorov-Smirnov
  Test)
- [`calc_signal_dp()`](https://kabilanbio.github.io/scSimEval/reference/calc_signal_dp.md)
  : Detect Differential Proportion / Zero-Inflation (DP) Genes (Chisq
  Test)
- [`calc_signal_bd()`](https://kabilanbio.github.io/scSimEval/reference/calc_signal_bd.md)
  : Detect Bimodal Distribution (BD) Genes (Bimodal Separation Index)
- [`calc_delta_variance()`](https://kabilanbio.github.io/scSimEval/reference/calc_delta_variance.md)
  : Calculate Delta Variance (Pseudoreplication Bias Metric)
- [`calc_cell_cycle_phase_fidelity()`](https://kabilanbio.github.io/scSimEval/reference/calc_cell_cycle_phase_fidelity.md)
  : Evaluate Cell Cycle Phase Distribution Fidelity
- [`calc_deconvolution_accuracy()`](https://kabilanbio.github.io/scSimEval/reference/calc_deconvolution_accuracy.md)
  : Evaluate Cell-Type Deconvolution & Mixture Proportion Accuracy
- [`calc_variance_decomposition()`](https://kabilanbio.github.io/scSimEval/reference/calc_variance_decomposition.md)
  : Variance Component Decomposition (BASiCS & muscat)
- [`calc_intraclass_correlation()`](https://kabilanbio.github.io/scSimEval/reference/calc_intraclass_correlation.md)
  : Calculate Intraclass Correlation Coefficient (ICC) (hierarchicell &
  rescueSim)
- [`calc_expected_mi()`](https://kabilanbio.github.io/scSimEval/reference/calc_expected_mi.md)
  : Expected Mutual Information for Cluster Comparison

## Category 6: Trajectory & Lineage Dynamics

Inference and benchmarking of differentiation pseudotime trajectories
and lineage branching topologies.

- [`infer_scrna_pseudotime()`](https://kabilanbio.github.io/scSimEval/reference/infer_scrna_pseudotime.md)
  : Automatically Infer Pseudotime Trajectory from scRNA-seq Counts
- [`infer_scrna_lineage_tree()`](https://kabilanbio.github.io/scSimEval/reference/infer_scrna_lineage_tree.md)
  : Automatically Infer Lineage Tree (hclust) from scRNA-seq Counts &
  Cell Types
- [`calc_pseudotime_correlation()`](https://kabilanbio.github.io/scSimEval/reference/calc_pseudotime_correlation.md)
  : Calculate Correlation of Geodesic Pseudotime Distances
- [`calc_tree_height_discrepancy()`](https://kabilanbio.github.io/scSimEval/reference/calc_tree_height_discrepancy.md)
  : Calculate Lineage Tree Branch Height Discrepancy (RMSE)

## Category 7: Multiomics Cross-Modality Coupling

Inter-modality regulatory coupling, peak-to-gene linkages, and joint
chromatin-expression modularity.

- [`calc_cross_modality_correlation()`](https://kabilanbio.github.io/scSimEval/reference/calc_cross_modality_correlation.md)
  : Calculate Cross-Modality Correlation Fidelity
- [`calc_foscttm()`](https://kabilanbio.github.io/scSimEval/reference/calc_foscttm.md)
  : Fraction of Samples Closer Than The True Match (FOSCTTM)
- [`calc_cross_modal_generation()`](https://kabilanbio.github.io/scSimEval/reference/calc_cross_modal_generation.md)
  : Cross-Modality In Silico Generation & Translation Fidelity
- [`calc_modality_alignment()`](https://kabilanbio.github.io/scSimEval/reference/calc_modality_alignment.md)
  : Modality Alignment & Omics Layer Mixing in Joint Latent Space
- [`calc_network_jaccard()`](https://kabilanbio.github.io/scSimEval/reference/calc_network_jaccard.md)
  : Jaccard Similarity of Predicted Regulatory Networks / Links
- [`calc_atac_rna_coupling()`](https://kabilanbio.github.io/scSimEval/reference/calc_atac_rna_coupling.md)
  : Chromatin Accessibility to RNA Regulatory Effect Coupling
- [`calc_coregulation_fidelity()`](https://kabilanbio.github.io/scSimEval/reference/calc_coregulation_fidelity.md)
  : Multi-Omics Co-Regulation & Modularity Fidelity
- [`calc_peak_coaccessibility_fidelity()`](https://kabilanbio.github.io/scSimEval/reference/calc_peak_coaccessibility_fidelity.md)
  : Evaluate Single-Cell ATAC Peak Co-Accessibility Fidelity (SCRIP)
- [`calc_coexpression_module_fidelity()`](https://kabilanbio.github.io/scSimEval/reference/calc_coexpression_module_fidelity.md)
  : Evaluate Gene Co-Expression Module Fidelity (ESCO)
- [`calc_accessibility_profile_concordance()`](https://kabilanbio.github.io/scSimEval/reference/calc_accessibility_profile_concordance.md)
  : Evaluate Chromatin Accessibility Profile Concordance (DiTSim)

## Category 8: Computational Scalability & Resource Footprint

Empirical benchmarking of CPU runtime, peak RAM consumption, and
execution efficiency.

- [`benchmark_resource_usage()`](https://kabilanbio.github.io/scSimEval/reference/benchmark_resource_usage.md)
  : Benchmark Execution Runtime and Memory Usage
