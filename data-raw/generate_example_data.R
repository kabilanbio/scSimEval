# =============================================================================
# Script to Generate Synthetic Example Datasets for scSimEval
# Generates:
#   - example_scrna: list(ref, sim, cell_types, batch_info)
#   - example_scatac: list(ref, sim, cell_types, batch_info)
#   - example_multiomics: list(ref_multi, sim_multi, cell_types, batch_info, resource_stats)
# =============================================================================

set.seed(42)

n_genes <- 60
n_peaks <- 60
n_cells <- 80

# 1. Annotations
cell_names <- paste0("Cell_", sprintf("%02d", 1:n_cells))
gene_names <- paste0("Gene_", sprintf("%02d", 1:n_genes))
peak_names <- paste0("Peak_", sprintf("%02d", 1:n_peaks))

cell_types <- factor(rep(c("TypeA", "TypeB"), each = n_cells / 2),
                     levels = c("TypeA", "TypeB"))
names(cell_types) <- cell_names

batch_info <- factor(rep(c("Batch1", "Batch2"), length.out = n_cells),
                     levels = c("Batch1", "Batch2"))
names(batch_info) <- cell_names

# 2. scRNA-seq Count Matrices (Negative Binomial with cell-type shifts)
base_mu_A <- rnorm(n_genes, mean = 4.0, sd = 0.8)
base_mu_B <- base_mu_A
base_mu_B[1:15] <- base_mu_B[1:15] + 2.5 # Differential expression signal

mu_ref_rna <- cbind(
  matrix(rep(base_mu_A, n_cells / 2), nrow = n_genes),
  matrix(rep(base_mu_B, n_cells / 2), nrow = n_genes)
)

ref_rna <- matrix(
  rnbinom(n_genes * n_cells, mu = pmax(0.1, mu_ref_rna), size = 1.2),
  nrow = n_genes, ncol = n_cells,
  dimnames = list(gene_names, cell_names)
)

sim_rna <- matrix(
  rnbinom(n_genes * n_cells, mu = pmax(0.1, mu_ref_rna * 0.95), size = 1.1),
  nrow = n_genes, ncol = n_cells,
  dimnames = list(gene_names, cell_names)
)

example_scrna <- list(
  ref = ref_rna,
  sim = sim_rna,
  cell_types = cell_types,
  batch_info = batch_info
)

# 3. scATAC-seq Peak Accessibility Matrices (Binomial/Poisson peak counts)
base_prob_A <- runif(n_peaks, min = 0.1, max = 0.3)
base_prob_B <- base_prob_A
base_prob_B[1:15] <- pmin(0.85, base_prob_B[1:15] + 0.3) # Differential accessibility signal

prob_ref_atac <- cbind(
  matrix(rep(base_prob_A, n_cells / 2), nrow = n_peaks),
  matrix(rep(base_prob_B, n_cells / 2), nrow = n_peaks)
)

ref_atac <- matrix(
  rbinom(n_peaks * n_cells, size = 2, prob = prob_ref_atac),
  nrow = n_peaks, ncol = n_cells,
  dimnames = list(peak_names, cell_names)
)

sim_atac <- matrix(
  rbinom(n_peaks * n_cells, size = 2, prob = prob_ref_atac * 0.92),
  nrow = n_peaks, ncol = n_cells,
  dimnames = list(peak_names, cell_names)
)

example_scatac <- list(
  ref = ref_atac,
  sim = sim_atac,
  cell_types = cell_types,
  batch_info = batch_info
)

# 4. Master Unified Multiomics Bundle
example_multiomics <- list(
  ref_multi = list(rna = ref_rna, atac = ref_atac),
  sim_multi = list(rna = sim_rna, atac = sim_atac),
  cell_types = cell_types,
  batch_info = batch_info,
  resource_stats = list(
    cpu_time = 14.2,
    memory_mb = 385.4,
    system_time = 1.1,
    elapsed_time = 15.3
  )
)

# 5. Save to data/ with xz compression
if (!dir.exists("data")) dir.create("data", recursive = TRUE)

save(example_scrna, file = "data/example_scrna.rda", compress = "xz")
save(example_scatac, file = "data/example_scatac.rda", compress = "xz")
save(example_multiomics, file = "data/example_multiomics.rda", compress = "xz")

message("Successfully generated and saved example_scrna, example_scatac, and example_multiomics to data/.")
