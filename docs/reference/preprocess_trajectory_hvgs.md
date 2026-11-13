# Preprocess and Select Highly Variable Genes for Trajectory Evaluation (Saelens et al., 2019)

Implements the standardized scRNA-seq quality control and HVG selection
workflow from the *dynverse* benchmarking study (*Saelens et al., Nature
Biotechnology, 2019*):

1.  Harsh cell outlier filtering (3 MADs from median on library size,
    detected genes, and mitochondrial percentage).

2.  Low-expression gene filtering (expressed in \\\ge 5\\\\ of cells and
    average expression \\\> 0.02\\).

3.  HVG selection via mean-variance decomposition
    ([`scran::modelGeneVar`](https://rdrr.io/pkg/scran/man/modelGeneVar.html)
    with `FDR < 0.05` and biological variance `> 0.5`, with adaptive
    top-HVG fallback for sparse/small benchmark sets).

4.  Iterative removal of all-zero cells and uninformative zero features
    until convergence.

## Usage

``` r
preprocess_trajectory_hvgs(
  ref_data,
  sim_data = NULL,
  cell_types_ref = NULL,
  cell_types_sim = NULL,
  nmads = 3,
  min_cell_frac = 0.05,
  min_mean_expr = 0.02,
  fdr_threshold = 0.05,
  bio_var_threshold = 0.5,
  min_hvgs = 30L,
  max_hvgs = 1500L
)
```

## Arguments

- ref_data:

  Reference count matrix, `SingleCellExperiment`, or `Seurat` object.

- sim_data:

  Optional simulated count matrix, `SingleCellExperiment`, or `Seurat`
  object.

- cell_types_ref:

  Optional vector of reference cell-type labels.

- cell_types_sim:

  Optional vector of simulated cell-type labels.

- nmads:

  MAD threshold for outlier cell removal (default `3`).

- min_cell_frac:

  Minimum fraction of cells expressing a retained gene (default `0.05`).

- min_mean_expr:

  Minimum average expression count for a retained gene (default `0.02`).

- fdr_threshold:

  FDR cutoff for HVG biological variance (default `0.05`).

- bio_var_threshold:

  Minimum biological variance component for strict HVG selection
  (default `0.5`).

- min_hvgs:

  Minimum HVG floor with adaptive fallback if fewer pass strict
  thresholds (default `30L`).

- max_hvgs:

  Maximum HVG cap for tractable Random Forest importance fitting
  (default `1500L`).

## Value

A list containing HVG-filtered log-normalized expression matrices
(`ref_expr`, `sim_expr`), selected HVG names (`hvgs`), synced cell-type
vectors (`cell_types_ref`, `cell_types_sim`), HVG statistics table
(`hvg_stats`), and filtering counts (`summary`).

## Examples

``` r
data(example_scrna, package = "scSimEval")
prep <- preprocess_trajectory_hvgs(example_scrna$ref, example_scrna$sim,
  cell_types_ref = example_scrna$cell_types, cell_types_sim = example_scrna$cell_types)
length(prep$hvgs)
#> [1] 60
```
