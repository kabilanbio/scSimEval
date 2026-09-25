# =============================================================================
# Demonstration: Ground-Truth-Free Multiomics Benchmarking (56 Selected Metrics)
# Operating strictly on:
#   1. Real & simulated single-cell multiomics (scRNA-seq + scATAC-seq)
#   2. Cell type annotations
#   3. Batch / donor annotations
#   4. Logged computational scalability statistics (CPU time, Peak RAM)
#   5. Trajectory & lineage features auto-inferred from scRNA-seq counts
# Package: scSimEval
# =============================================================================

# Dynamically source all R scripts
r_dir <- if (dir.exists("scSimEval/R")) {
  "scSimEval/R"
} else if (dir.exists("R")) {
  "R"
} else {
  file.path("..", "R")
}

scripts <- sort(list.files(r_dir, pattern = "\\.R$", full.names = TRUE))
for (s in scripts) source(s)

set.seed(42)

# Load packaged example datasets or generate on the fly
if (file.exists("data/example_scrna.rda") && file.exists("data/example_scatac.rda")) {
  message("Loading packaged example datasets (example_scrna, example_scatac)...")
  load("data/example_scrna.rda")
  load("data/example_scatac.rda")
  ref_rna <- example_scrna$ref
  sim_rna <- example_scrna$sim
  ref_atac <- example_scatac$ref
  sim_atac <- example_scatac$sim
  cell_types <- example_scrna$cell_types
  batches <- example_scrna$batch_info
  n_genes <- nrow(ref_rna)
  n_cells <- ncol(ref_rna)
} else {
  message("Generating mock multiomics reference and simulated data (scRNA-seq + scATAC-seq)...")
  n_genes <- 60
  n_cells <- 80
  cell_types <- factor(rep(c("TypeA", "TypeB"), each = n_cells / 2))
  batches <- factor(rep(c("Batch1", "Batch2"), length.out = n_cells))
  base_mu_A <- rnorm(n_genes, mean = 4, sd = 1)
  base_mu_B <- base_mu_A
  base_mu_B[1:15] <- base_mu_B[1:15] + 2.5
  mu_ref <- cbind(matrix(rep(base_mu_A, n_cells / 2), nrow = n_genes),
                  matrix(rep(base_mu_B, n_cells / 2), nrow = n_genes))
  ref_rna <- matrix(rnbinom(n_genes * n_cells, mu = pmax(0.1, mu_ref), size = 1.2),
                    nrow = n_genes, ncol = n_cells,
                    dimnames = list(paste0("Gene_", 1:n_genes), paste0("Cell_", 1:n_cells)))
  sim_rna <- matrix(rnbinom(n_genes * n_cells, mu = pmax(0.1, mu_ref * 0.95), size = 1.1),
                    nrow = n_genes, ncol = n_cells,
                    dimnames = list(paste0("Gene_", 1:n_genes), paste0("Cell_", 1:n_cells)))
  ref_atac <- matrix(rbinom(n_genes * n_cells, size = 2, prob = 0.18),
                     nrow = n_genes, ncol = n_cells,
                     dimnames = list(paste0("Peak_", 1:n_genes), paste0("Cell_", 1:n_cells)))
  sim_atac <- matrix(rbinom(n_genes * n_cells, size = 2, prob = 0.17),
                     nrow = n_genes, ncol = n_cells,
                     dimnames = list(paste0("Peak_", 1:n_genes), paste0("Cell_", 1:n_cells)))
}

# =============================================================================
# Category 1 & 2: Distributional & Summary Property Distances
# =============================================================================
message("\n--- [Category 1 & 2] Evaluating Statistical Distribution Accuracies ---")
unimodal_res <- evaluate_simulation_accuracy(
  ref_data = ref_rna,
  sim_data = sim_rna,
  compute_bivariate = FALSE,
  verbose = FALSE
)

message("Summary Table (Unimodal Distribution Metrics):")
lib_subset <- subset(unimodal_res$metrics_summary_table, Property == "library_size")
print(lib_subset, row.names = FALSE)

# Generative Model Goodness of Fit (AIC/BIC & LRT)
message("\nEvaluating Generative Model Goodness of Fit (AIC/BIC & LRT)...")
fitted_means <- matrix(rowMeans(sim_rna), nrow = n_genes, ncol = n_cells)
m_aic_bic <- calc_marginal_aic_bic(sim_rna, fitted_means, distribution = "nb")
message("  - Total Marginal AIC: ", round(m_aic_bic$total_aic, 2))
message("  - Total Marginal BIC: ", round(m_aic_bic$total_bic, 2))

lrt_res <- calc_likelihood_ratio_test(
  alter_model = m_aic_bic$total_loglik,
  null_model = m_aic_bic$total_loglik - 25,
  df_alter = 2 * n_genes,
  df_null = n_genes
)
message("  - Likelihood Ratio Statistic (LRT): ", round(lrt_res$LR_statistic, 2))
message("  - LRT Chi-squared p-value: ", format.pval(lrt_res$p_value, digits = 3))

# Zero-Probability Dropout Curve & Excess Zero Weights
message("\nEvaluating Zero-Probability Dropout Curve & Excess Zero Weights...")
zero_eval <- evaluate_zero_probability_curve(ref_rna, sim_rna)
message("  - Dropout Curve Slope Discrepancy: ", round(zero_eval["slope_error"], 4))
message("  - Dropout Curve Midpoint Discrepancy: ", round(zero_eval["midpoint_error"], 4))

excess_zeros <- calc_excess_zero_weights(ref_rna)
message("  - Mean Posterior Excess Zero Weight (zingeR): ", round(excess_zeros$mean_excess_zero_weight, 4))

# Deep Generative Manifold Distances (MMD & Fréchet Distance)
message("\nEvaluating Deep Generative Manifold Distances (MMD & Fréchet Distance)...")
mmd_res <- calc_mmd(ref_rna, sim_rna)
message("  - Maximum Mean Discrepancy (MMD): ", round(mmd_res$mmd, 4))

fsd_res <- calc_frechet_singlecell_distance(ref_rna, sim_rna, n_pcs = 10)
message("  - Fréchet Single-Cell Distance (FSD): ", round(fsd_res$fsd, 4))

# Kinetic Noise Decomposition
message("\nEvaluating Kinetic Noise Decomposition (SymSim)...")
sym_noise <- calc_kinetic_noise_decomposition(sim_rna, cell_states = cell_types)
message("  - Mean Intrinsic Noise: ", round(sym_noise$mean_intrinsic_noise, 4))
message("  - Mean Extrinsic Noise: ", round(sym_noise$mean_extrinsic_noise, 4))

# Accessibility-Sparsity Polynomial Curve
message("\nEvaluating Accessibility-Sparsity Polynomial Curve (simATAC)...")
atac_curve_res <- evaluate_accessibility_sparsity_curve(ref_atac, sim_atac)
message("  - Reference Accessibility Curve R2: ", round(atac_curve_res$ref_r_squared, 4))
message("  - Simulation Accessibility Curve R2: ", round(atac_curve_res$sim_r_squared, 4))
message("  - Delta Curve R2: ", round(atac_curve_res$delta_r_squared, 4))

# =============================================================================
# Category 3: Cell Clustering & Population Concordance
# =============================================================================
message("\n--- [Category 3] Evaluating Cell Clustering & Concordance ---")
pred_clusters <- factor(kmeans(t(ref_rna), centers = 2)$cluster)

clust_metrics <- evaluate_clustering_metrics(
  data = sim_rna,
  cell_types = cell_types,
  ref_data = ref_rna
)

message("Silhouette Width (ASW): ", round(clust_metrics$silhouette_sim, 4))
message("Dunn Index: ", round(clust_metrics$dunn_sim, 4))
message("Davies-Bouldin (DB) Index: ", round(clust_metrics$davies_bouldin_sim, 4))
message("Adjusted Rand Index (ARI): ", round(clust_metrics$ari, 4))
message("Normalized Mutual Information (NMI): ", round(clust_metrics$nmi, 4))
message("Adjusted Mutual Information (AMI): ", round(clust_metrics$ami, 4))
message("Clustering Deviation Index (CDI): ", round(clust_metrics$cdi, 4))

# =============================================================================
# Category 4: Batch Effect Mixing & Confounder Preservation
# =============================================================================
message("\n--- [Category 4] Evaluating Batch Effect Mixing ---")
batch_metrics <- evaluate_batch_metrics(
  data = sim_rna,
  batch_info = batches,
  pre_data = ref_rna,
  cell_types = cell_types,
  verbose = FALSE
)

message("Cell-Specific Mixing Score (CMS): ", round(batch_metrics$cms, 4))
message("Inverse Simpson Index (ISI / LISI): ", round(batch_metrics$isi, 4))
message("Seurat Mixing Metric: ", round(batch_metrics$seurat_mixing_metric, 4))
message("Local Density Factor Difference (ldfDiff): ", round(batch_metrics$ldf_diff, 4))
message("PC Regression Batch R2 (PCR): ", round(batch_metrics$pcr_r2, 4))
message("Batch Shannon Entropy: ", round(batch_metrics$shannon_entropy, 4))
message("Cross-Batch Classification Accuracy: ", round(batch_metrics$cross_batch_accuracy, 4))

# Hierarchical Intraclass Correlation (ICC)
icc_res <- calc_intraclass_correlation(sim_rna, donor_labels = batches, n_genes = 30)
message("Hierarchical Intraclass Correlation (ICC): ", round(icc_res$mean_icc, 4))

# =============================================================================
# Category 5: Biological Signal & Downstream Performance
# =============================================================================
message("\n--- [Category 5] Evaluating Biological Signal Fidelity ---")
signal_res <- evaluate_simbench_signals(
  ref_mat = ref_rna,
  sim_mat = sim_rna,
  ref_celltypes = cell_types,
  sim_celltypes = cell_types
)
print(signal_res, row.names = FALSE)

# Predictive cell identity classifier
pred_model <- evaluate_predictive_de_model(
  data = sim_rna,
  group = cell_types
)
message("Predictive Classification Accuracy: ", round(pred_model$accuracy, 4))
message("Predictive Classification Macro F1: ", round(pred_model$F1, 4))

# Cell-Type Deconvolution & Mixture Accuracy
deconv_res <- calc_deconvolution_accuracy(
  true_proportions = cell_types,
  estimated_proportions = cell_types,
  batch_ref = batches,
  batch_sim = batches
)
message("Cell-Type Mixture Deconvolution RMSE: ", round(deconv_res$deconvolution_rmse, 4))
message("Cell-Type Mixture Deconvolution JSD: ", round(deconv_res$deconvolution_jsd, 4))

# Delta Variance (Pseudoreplication bias)
dv_res <- calc_delta_variance(sim_rna, replicates = batches, conditions = cell_types)
message("Mean Delta Variance: ", round(dv_res$mean_delta_variance, 4))

# =============================================================================
# Category 6: Trajectory & Lineage Dynamics (Inferred from scRNA-seq)
# =============================================================================
message("\n--- [Category 6] Evaluating Trajectory Dynamics Inferred from scRNA-seq ---")
pt_ref <- infer_scrna_pseudotime(ref_rna)
pt_sim <- infer_scrna_pseudotime(sim_rna)
message("Inferred reference pseudotime range: [", round(min(pt_ref), 2), ", ", round(max(pt_ref), 2), "]")

traj_res <- evaluate_trajectory_metrics(
  ref_data = ref_rna,
  sim_data = sim_rna,
  cell_types_ref = cell_types,
  cell_types_sim = cell_types
)
message("Pseudotime Correlation: ", round(traj_res$pseudotime_correlation, 4))
message("Lineage Tree Height Discrepancy (RMSE): ", round(traj_res$tree_height_rmse, 4))

# =============================================================================
# Category 7: Cross-Modality Coupling & Modularity
# =============================================================================
message("\n--- [Category 7] Evaluating Cross-Modality Regulatory Linkage & Alignment ---")
ref_multi <- list(rna = ref_rna, atac = ref_atac)
sim_multi <- list(rna = sim_rna, atac = sim_atac)

cross_cor_res <- calc_cross_modality_correlation(
  ref_mod1 = ref_rna,
  ref_mod2 = ref_atac,
  sim_mod1 = sim_rna,
  sim_mod2 = sim_atac
)
message("Cross-Modality Correlation KS Distance: ", round(cross_cor_res$cross_modality_cor_KS, 4))

cpred <- evaluate_cross_modal_prediction(sim_atac, cell_types)
message("Cross-Modal Label Transfer Accuracy (scATAC -> scRNA cell types): ", round(cpred$cross_modal_accuracy, 4))

# FOSCTTM Cell Alignment
p1 <- stats::prcomp(t(sim_rna), rank. = 5)$x
p2 <- stats::prcomp(t(sim_atac), rank. = 5)$x
fos_res <- calc_foscttm(p1, p2)
message("FOSCTTM Alignment Score: ", round(fos_res$foscttm, 4))
message("Match@1 Rate: ", round(fos_res$match_at_1, 4))

# In Silico Translation & Generation Fidelity
gen_fid <- calc_cross_modal_generation(ref_atac, sim_atac)
message("In Silico Generation Mean Cell PCC: ", round(gen_fid$mean_cell_pcc, 4))
message("In Silico Generation RMSE: ", round(gen_fid$rmse, 4))

# ATAC-RNA Regulatory Coupling
coupling_res <- calc_atac_rna_coupling(sim_atac, sim_rna)
message("ATAC-RNA Regulatory Coupling Mean Correlation: ", round(coupling_res$mean_coupling_cor, 4))

# Peak Co-Accessibility Fidelity
peak_coacc <- calc_peak_coaccessibility_fidelity(ref_atac, sim_atac, top_n_peaks = 50)
message("Peak Co-Accessibility RV Coefficient: ", round(peak_coacc$rv_coefficient, 4))

# Co-Expression Module Fidelity
coexpr_fid <- calc_coexpression_module_fidelity(ref_rna, sim_rna, n_modules = 3, top_genes = 60)
message("Co-Expression Module Correlation R: ", round(coexpr_fid$module_correlation_r, 4))
message("Co-Expression Modularity Fidelity: ", round(coexpr_fid$modularity_fidelity, 4))

# Chromatin Accessibility Profile Concordance
acc_prof <- calc_accessibility_profile_concordance(ref_atac, sim_atac, cell_types = cell_types)
message("Accessibility Profile Global PCC: ", round(acc_prof$global_pcc, 4))
message("Accessibility Profile KL Divergence: ", round(acc_prof$kl_divergence, 4))

# =============================================================================
# Category 8: Computational Scalability & Master Pipeline Execution
# =============================================================================
message("\n--- [Category 8 & Master Pipeline] Executing evaluate_multiomics_accuracy ---")

res_master <- evaluate_multiomics_accuracy(
  ref_multi = ref_multi,
  sim_multi = sim_multi,
  cell_types = cell_types,
  batch_info = batches,
  cpu_time = 12.4,
  memory_mb = 285.6,
  system_time = 1.8,
  elapsed_time = 14.5,
  compute_bivariate = FALSE,
  verbose = FALSE
)

message("Total Metrics in Master Summary Table: ", nrow(res_master$benchmark_summary_table))
message("\nMaster Summary Table Head:")
print(head(res_master$benchmark_summary_table, 15), row.names = FALSE)

message("\n====================================================================")
message("ALL 56 CURATED GROUND-TRUTH-FREE METRICS EXECUTED CLEANLY!")
message("====================================================================")
