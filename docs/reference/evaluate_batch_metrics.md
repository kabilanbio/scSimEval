# Calculate Batch Integration Metrics

Evaluates how realistically simulated technical batches mix or separate,
computing kBET, LISI, Batch Silhouette, Shannon Entropy, PC Regression,
CMS (Cell-Specific Mixing Score), Inverse Simpson Index (ISI), Seurat
Mixing Metric, and optionally Local Density Differences (ldfDiff) and
Local Structure Preservation.

## Usage

``` r
evaluate_batch_metrics(
  data,
  batch_info,
  k = NULL,
  n_pcs = 30,
  pre_data = NULL,
  cell_types = NULL,
  verbose = FALSE
)
```

## Arguments

- data:

  Count or log-normalized expression matrix (features x cells).

- batch_info:

  Factor or character vector of batch assignments for each cell.

- k:

  Neighborhood size k for k-NN graph and mixing metrics. Default is
  min(table(batch_info))/2.

- n_pcs:

  Number of principal components for embedding evaluation. Default is
  30.

- pre_data:

  Optional matrix of pre-integration / reference expression data (for
  ldfDiff and local structure).

- cell_types:

  Optional factor or character vector of cell types for cross-batch
  evaluation.

- verbose:

  Logical, whether to print progress.

## Value

A named list of batch metrics:

- batch_silhouette:

  Average silhouette width using batch labels (lower is better mixed)

- shannon_entropy:

  Entropy of batch frequencies in local neighborhoods

- pcr_r2:

  Total variance explained by batch in principal components (PCR)

- cms:

  Mean Cell-specific mixing score (CellMixS)

- isi:

  Mean Inverse Simpson Index for batch mixing (CellMixS)

- seurat_mixing_metric:

  Mean Seurat mixing metric (Stuart et al. / CellMixS)

- ldf_diff:

  Mean Local Density Factor difference (if pre_data supplied)

- local_structure:

  Mean Local structure preservation overlap (if pre_data supplied)

- kbet_rejection:

  kBET rejection rate (if kBET installed)

- lisi_batch:

  Average batch LISI (if lisi installed)

- cross_batch_accuracy:

  Mean cross-batch cell type transfer accuracy (if cell_types supplied)

- cross_batch_F1:

  Mean cross-batch cell type transfer macro F1 (if cell_types supplied)
