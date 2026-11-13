# Calculate Batch Integration Metrics

Evaluates how realistically simulated technical batches mix or separate,
computing the five canonical Category IV (Technical Variation & Batch
Mixing Dynamics) metrics: Batch Silhouette, PC Regression R2, CMS
(Cell-Specific Mixing Score), Inverse Simpson Index (ISI), and
cross-batch cell type classification macro F1.

## Usage

``` r
evaluate_batch_metrics(
  data,
  batch_info,
  k = NULL,
  n_pcs = 30,
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

- cell_types:

  Optional factor or character vector of cell types for cross-batch
  evaluation.

- verbose:

  Logical, whether to print progress.

## Value

A named list of batch metrics:

- batch_silhouette:

  Average silhouette width using batch labels (lower is better mixed)

- pcr_r2:

  Total variance explained by batch in principal components (PCR)

- cms:

  Mean Cell-specific mixing score (CellMixS)

- isi:

  Mean Inverse Simpson Index for batch mixing (CellMixS)

- cross_batch_F1:

  Mean cross-batch cell type transfer macro F1 (if cell_types supplied)

## Examples

``` r
data(example_scrna, package = "scSimEval")
evaluate_batch_metrics(example_scrna$ref, example_scrna$batch_info,
                       cell_types = example_scrna$cell_types)
#> $batch_silhouette
#> [1] -0.001626347
#> 
#> $pcr_r2
#> [1] 0.01058094
#> 
#> $cms
#> [1] 0
#> 
#> $isi
#> [1] 1.923783
#> 
#> $cross_batch_F1
#> [1] 0.4460483
#> 
```
