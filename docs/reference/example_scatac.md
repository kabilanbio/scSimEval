# Synthetic Example scATAC-seq Benchmarking Dataset

A synthetic single-cell ATAC sequencing peak accessibility dataset
designed for testing, benchmarking, and demonstrating evaluation metrics
in `scSimEval`. Contains paired empirical reference accessibility counts
and simulated accessibility counts across 60 peaks and 80 cells matching
the cells in `example_scrna`.

## Usage

``` r
data(example_scatac)
```

## Format

A list with 4 elements:

- ref:

  A numeric matrix of reference peak accessibility counts (60 peaks x 80
  cells).

- sim:

  A numeric matrix of simulated peak accessibility counts (60 peaks x 80
  cells).

- cell_types:

  A factor of cell-type identities ("TypeA", "TypeB") for the 80 cells.

- batch_info:

  A factor of batch annotations ("Batch1", "Batch2") for the 80 cells.

## Source

Synthetic benchmark data generated via binomial sampling with cell-type
specific accessibility shifts.

## Examples

``` r
data(example_scatac)
dim(example_scatac$ref)
#> [1] 60 80
dim(example_scatac$sim)
#> [1] 60 80
table(example_scatac$batch_info)
#> 
#> Batch1 Batch2 
#>     40     40 
```
