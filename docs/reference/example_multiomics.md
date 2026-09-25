# Synthetic Example Multiomics Benchmarking Dataset

A unified single-cell multiomics dataset pairing scRNA-seq and
scATAC-seq profiles across the same 80 cells, accompanied by cell-type
labels, batch annotations, and logged computational scalability metrics.

## Usage

``` r
data(example_multiomics)
```

## Format

A list with 5 elements:

- ref_multi:

  Named list of reference matrices: `rna` (60 x 80) and `atac` (60 x
  80).

- sim_multi:

  Named list of simulated matrices: `rna` (60 x 80) and `atac` (60 x
  80).

- cell_types:

  A factor of cell-type identities ("TypeA", "TypeB") for the 80 cells.

- batch_info:

  A factor of batch annotations ("Batch1", "Batch2") for the 80 cells.

- resource_stats:

  A list of computational scalability metrics: `cpu_time`, `memory_mb`,
  `system_time`, `elapsed_time`.

## Source

Synthetic multiomics benchmark data generated for testing `scSimEval`.

## Examples

``` r
data(example_multiomics)
names(example_multiomics$ref_multi)
#> [1] "rna"  "atac"
names(example_multiomics$sim_multi)
#> [1] "rna"  "atac"
```
