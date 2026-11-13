# Evaluate Multiomics Integration & Cross-Modal Alignment Fidelity (Category VI)

Evaluates the multiomics integration and cross-modal alignment fidelity
of simulated single-cell multiomics data across 5 gold-standard
integration metrics extracted from key single-cell benchmarking suites
(`scib`, `scMultiSim`, `scCross`, `Matilda`, and `scMoMtF`):

1.  **Batch ASW (`batch_asw`):** Average Silhouette Width evaluating
    batch mixing within cell types (scib, scMultiSim).

2.  **Cell-Type ASW (`cell_type_asw`):** Average Silhouette Width
    evaluating cell-type bio-conservation in the multiomics space
    (scCross, scib).

3.  **Graph Connectivity (`graph_connectivity`):** Relative size of the
    largest connected component in the k-NN graph for each cell type
    (scCross, scMultiSim, scib).

4.  **ARI (`ari`):** Adjusted Rand Index comparing joint multiomics
    clustering against ground-truth cell types (scib, scCross,
    scMultiSim, Matilda, scMoMtF).

5.  **NMI (`nmi`):** Normalized Mutual Information quantifying
    information concordance between joint clustering and cell types
    (scib, scCross, scMultiSim, Matilda, scMoMtF).

## Usage

``` r
evaluate_multiomics_integration(
  ref_multi = NULL,
  sim_multi,
  cell_types = NULL,
  batch_info = NULL,
  embed = NULL,
  reduction_name = NULL,
  n_pcs = 15,
  k = 30,
  threads = 1,
  verbose = TRUE
)
```

## Arguments

- ref_multi:

  Optional reference multiomics list of matrices (e.g.
  `list(rna = ..., atac = ...)`).

- sim_multi:

  Simulated multiomics list of matrices or a matrix representing joint
  embedding coordinates.

- cell_types:

  Optional cell-type vector or list/data.frame of annotations for each
  cell.

- batch_info:

  Optional batch/donor/replicate vector or list/data.frame of
  annotations for each cell.

- embed:

  Optional precomputed joint embedding matrix (cells x dimensions). If
  NULL, computed automatically from `sim_multi`.

- reduction_name:

  Optional name of the precomputed dimensional reduction to extract from
  Seurat/SCE objects.

- n_pcs:

  Number of principal components per modality to retain for joint
  embedding (default 15).

- k:

  Number of nearest neighbors for graph-based metrics (default 30).

- threads:

  Number of CPU threads for parallel processing (default 1).

- verbose:

  Logical; whether to print progress messages (default TRUE).

## Value

A list of class `scSimEval_integration` containing:

- metrics:

  Named list of the 5 calculated integration metric values (batch_asw,
  cell_type_asw, graph_connectivity, ari, nmi).

- integration_summary_table:

  A tidy data.frame containing the 5 metrics with columns `Category`,
  `Property`, `Metric`, `Value`, and `Modality`.

- embedding:

  The joint multiomics embedding matrix used for evaluation.

## Examples

``` r
data(example_multiomics, package = "scSimEval")
ct <- rep(c("A", "B"), each = 50)
res <- evaluate_multiomics_integration(sim_multi = example_multiomics$sim_multi, cell_types = ct)
#> Using precomputed dimensional reduction from sim_multi for Category VI.
```
