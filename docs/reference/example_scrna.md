# Synthetic Example scRNA-seq Benchmarking Dataset

A synthetic single-cell RNA sequencing dataset designed for testing,
benchmarking, and demonstrating evaluation metrics in `scSimEval`.
Contains paired empirical reference counts and simulated counts across
60 genes and 80 cells from two cell types and two batches.

## Usage

``` r
data(example_scrna)
```

## Format

A list with 4 elements:

- ref:

  A numeric matrix of reference count data (60 genes x 80 cells).

- sim:

  A numeric matrix of simulated count data (60 genes x 80 cells).

- cell_types:

  A factor of cell-type identities ("TypeA", "TypeB") for the 80 cells.

- batch_info:

  A factor of batch annotations ("Batch1", "Batch2") for the 80 cells.

## Source

Synthetic benchmark data generated via negative binomial sampling with
cell-type specific expression shifts.

## Examples

``` r
data(example_scrna)
dim(example_scrna$ref)
#> [1] 60 80
dim(example_scrna$sim)
#> [1] 60 80
table(example_scrna$cell_types)
#> 
#> TypeA TypeB 
#>    40    40 
```
