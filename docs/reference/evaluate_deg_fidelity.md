# Evaluate Differentially Expressed Gene (DEG) Fidelity

Evaluates the preservation of differentially expressed genes (DEGs) and
biological signal between empirical reference and simulated datasets
across three landmark single-cell benchmarking frameworks:

- **Simpipe (Duo et al., 2024):** True DEG ratio, Distribution score
  (p-value uniformity via Pearson Chi-Square goodness-of-fit test on
  remaining genes after DEG removal), and supervised machine learning
  classification (Accuracy, Precision, Recall, F1) using simulated DEGs
  to predict cell identity.

- **SimBench (Cao et al., 2021):** Symmetric Mean Absolute Percentage
  Error (SMAPE) on DEG proportions, SimBench DE fidelity score, log2
  fold-change effect size Pearson and Spearman correlation, and top-N
  DEG Jaccard overlap.

- **Shaky Foundations (Crowell et al., 2023):** Group silhouette
  separation width, silhouette discrepancy, and Percent Variance
  Explained (PVE) by group assignment.

## Usage

``` r
evaluate_deg_fidelity(
  ref_data,
  sim_data,
  ref_celltypes,
  sim_celltypes,
  group1 = NULL,
  group2 = NULL,
  fdr_cutoff = 0.05,
  logfc_cutoff = 0.5,
  top_n_de = 100,
  classifier = c("knn", "svm", "rf"),
  run_pvalue_uniformity = TRUE
)
```

## Arguments

- ref_data:

  Reference single-cell expression or count matrix (features x cells).

- sim_data:

  Simulated single-cell expression or count matrix (features x cells).

- ref_celltypes:

  Factor or character vector of cell types for reference cells.

- sim_celltypes:

  Factor or character vector of cell types for simulated cells.

- group1:

  Optional character string specifying group 1 name. Default NULL
  (auto-selects top abundant type).

- group2:

  Optional character string specifying group 2 name. Default NULL
  (auto-selects second most abundant type).

- fdr_cutoff:

  Adjusted p-value significance threshold. Default is 0.05.

- logfc_cutoff:

  Absolute log2 fold-change cutoff. Default is 0.5.

- top_n_de:

  Number of top ranked DEGs to evaluate for Jaccard overlap and
  classification. Default is 100.

- classifier:

  Supervised classifier: "knn" (default), "svm", or "rf".

- run_pvalue_uniformity:

  Logical, whether to run Simpipe's Chi-square p-value uniformity test.
  Default is TRUE.

## Value

A list containing:

- `deg_summary_table`: Comprehensive data.frame uniting all 15 metrics
  across Simpipe, SimBench, and Shaky Foundations.

- `deg_genes_ref`: Character vector of significant DEGs in reference.

- `deg_genes_sim`: Character vector of significant DEGs in simulation.

- `contrast`: String describing the evaluated group contrast.

- `logfc_ref`: Named numeric vector of reference log2 fold-changes.

- `logfc_sim`: Named numeric vector of simulated log2 fold-changes.

- `ml_classification`: Detailed classifier performance metrics.

- `pvalue_uniformity`: Chi-square goodness-of-fit test results.
