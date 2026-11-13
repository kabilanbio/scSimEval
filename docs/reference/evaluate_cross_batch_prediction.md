# Evaluate Cross-Batch Cell Classification Transfer F1

Evaluates whether simulated or integrated multi-batch single-cell data
preserve cell-type identity across distinct experimental batches or
donors. Trains a k-NN classifier on each individual batch and tests the
macro F1 score of predictions on all other remaining batches.

## Usage

``` r
evaluate_cross_batch_prediction(coords, cell_types, batch_info, k = 5)
```

## Arguments

- coords:

  Embedding coordinates (cells x dimensions).

- cell_types:

  Factor or character vector of ground-truth cell type labels.

- batch_info:

  Factor or character vector of batch assignments.

- k:

  Integer number of nearest neighbors for the k-NN classifier (default
  5).

## Value

A list containing:

- mean_cross_batch_F1:

  Overall mean macro-averaged F1 score across all directed batch pairs

## Examples

``` r
coords <- matrix(stats::rnorm(100), 50, 2)
batch <- factor(rep(c("B1", "B2"), length.out = 50))
ct <- factor(rep(c("TypeA", "TypeB"), length.out = 50))
evaluate_cross_batch_prediction(coords, ct, batch)
#> $mean_cross_batch_F1
#> [1] NaN
#> 
```
