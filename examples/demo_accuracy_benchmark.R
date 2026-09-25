# =============================================================================
# Demonstration: Benchmarking Single-Cell Multiomics Simulation Accuracy
# Project: scSimEval
# =============================================================================

# Dynamically locate R/ directory whether run from repo root or project root
r_dir <- if (dir.exists("scSimEval/R")) {
  "scSimEval/R"
} else if (dir.exists("R")) {
  "R"
} else {
  file.path("..", "R")
}

scripts <- sort(list.files(r_dir, pattern = "\\.R$", full.names = TRUE))
for (s in scripts) {
  source(s)
}

set.seed(42)

message("Generating mock reference and simulated scRNA-seq and scATAC-seq data...")

# Mock scRNA-seq (100 genes x 200 cells)
n_genes <- 100
n_cells <- 200

# Reference RNA: Negative Binomial counts
ref_rna <- matrix(rnbinom(n_genes * n_cells, mu = 4, size = 1.2),
                  nrow = n_genes, ncol = n_cells,
                  dimnames = list(paste0("Gene_", 1:n_genes), paste0("Cell_", 1:n_cells)))

# Simulated RNA: slightly perturbed Negative Binomial
sim_rna <- matrix(rnbinom(n_genes * n_cells, mu = 3.8, size = 1.1),
                  nrow = n_genes, ncol = n_cells,
                  dimnames = list(paste0("Gene_", 1:n_genes), paste0("Cell_", 1:n_cells)))

# Mock scATAC-seq (100 peaks x 200 cells, sparse binary/low counts)
ref_atac <- matrix(rbinom(n_genes * n_cells, size = 2, prob = 0.15),
                   nrow = n_genes, ncol = n_cells,
                   dimnames = list(paste0("Peak_", 1:n_genes), paste0("Cell_", 1:n_cells)))

sim_atac <- matrix(rbinom(n_genes * n_cells, size = 2, prob = 0.14),
                   nrow = n_genes, ncol = n_cells,
                   dimnames = list(paste0("Peak_", 1:n_genes), paste0("Cell_", 1:n_cells)))

# -----------------------------------------------------------------------------
# Part 1: Evaluate a Single Modality (scRNA-seq) with the 8 simpipe Metrics
# -----------------------------------------------------------------------------
message("\n--- Running Unimodal Accuracy Evaluation (scRNA-seq) ---")
rna_eval <- evaluate_simulation_accuracy(
  ref_data = ref_rna,
  sim_data = sim_rna,
  compute_bivariate = FALSE, # Set TRUE when fasano.franceschini.test & ks installed
  verbose = TRUE
)

message("\nTop 15 Univariate Accuracy Metrics (scRNA-seq):")
print(head(rna_eval$metrics_summary_table, 15))

# -----------------------------------------------------------------------------
# Part 2: Evaluate Multiomics Simulation (scRNA-seq + scATAC-seq)
# -----------------------------------------------------------------------------
message("\n--- Running Multiomics Accuracy Evaluation (scRNA + scATAC) ---")
ref_multi <- list(mod1 = ref_rna, mod2 = ref_atac)
sim_multi <- list(mod1 = sim_rna, mod2 = sim_atac)

multi_eval <- evaluate_multiomics_accuracy(
  ref_multi = ref_multi,
  sim_multi = sim_multi,
  mod1_name = "scRNA",
  mod2_name = "scATAC",
  verbose = FALSE
)

message("\nCross-Modality Correlation Accuracy Metrics:")
print(multi_eval$cross_modality_correlation)

message("\nMultiomics Summary Table (Sample Rows):")
print(head(multi_eval$multimodal_summary_table, 15))

message("\nDone! All 8 simpipe accuracy metrics executed successfully.")
