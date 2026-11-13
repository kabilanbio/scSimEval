# Calculate Simpipe Simulated DEG Ratio and True DEG Recovery (Approach 1)

Evaluates the ability of simulation methods to simulate Differentially
Expressed Genes (DEGs) according to the Simpipe benchmarking framework
(Duo et al., Genome Biology 2024). Uses `edgeRQLFDetRate` on reference
data as ground truth, and compares candidate DEGs returned by the
simulation method.

## Usage

``` r
calc_simpipe_deg_ratio(
  ref_data,
  sim_data,
  ref_celltypes,
  sim_celltypes,
  sim_degs = NULL,
  fdr_cutoff = 0.05,
  logfc_cutoff = 0,
  min_cells = 5
)
```

## Arguments

- ref_data:

  Reference count matrix (genes x cells).

- sim_data:

  Simulated count matrix (genes x cells).

- ref_celltypes:

  Factor or character vector of reference cell groups.

- sim_celltypes:

  Factor or character vector of simulated cell groups.

- sim_degs:

  Optional list of character vectors containing candidate DEGs already
  returned by the simulation method for each pairwise group contrast. If
  NULL, edgeRQLFDetRate is executed on sim_data to detect candidate
  DEGs.

- fdr_cutoff:

  FDR threshold for true and candidate DEGs (default 0.05).

- logfc_cutoff:

  Absolute log2 fold change threshold (default 0).

- min_cells:

  Minimum number of cells per group to perform DE (default 5).

## Value

A list containing:

- P:

  Overall weighted Simpipe true DEG recovery proportion (P)

- deg_ratio:

  Overall ratio of candidate/simulated DEGs to true DEGs (N_given /
  N_true)

- pairwise_results:

  Data frame with per-pair details: Group1, Group2, N_true, N_given,
  N_overlap, P_true, DEG_ratio, Weight

- true_degs:

  Named list of true DEG gene vectors per pair

- candidate_degs:

  Named list of candidate/simulated DEG gene vectors per pair

## Details

For 2 groups (\\n = 2\\): \$\$P = P\_{\text{true}}\$\$ For 3 or more
groups (\\n \> 2\\): \$\$P = \sum\_{i=1}^M P\_{\text{true}, i} \cdot
\frac{N\_{\text{given}, i}}{N\_{\text{given}}}\$\$ where
\\P\_{\text{true}, i}\\ is the proportion of true DEGs identified from
two cell groups through `edgeRQLFDetRate`, \\N\_{\text{given}, i}\\
denotes the DEG number returned by the simulation method from the
\\i\\-th group pair, and \\N\_{\text{given}}\\ denotes the total DEG
number across all group pairs.

Also calculates the candidate-to-true DEG ratio: \$\$\text{DEG Ratio} =
\frac{N\_{\text{given}}}{N\_{\text{true}}}\$\$

## Examples

``` r
data(example_scrna, package = "scSimEval")
res <- calc_simpipe_deg_ratio(
  ref_data = example_scrna$ref,
  sim_data = example_scrna$sim,
  ref_celltypes = example_scrna$cell_types,
  sim_celltypes = example_scrna$cell_types
)
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
res$P
#> [1] 0
res$deg_ratio
#> [1] 0
```
