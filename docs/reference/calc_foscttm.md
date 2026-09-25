# Fraction of Samples Closer Than The True Match (FOSCTTM)

Computes the Fraction of Samples Closer Than The True Match (FOSCTTM) to
quantify cell alignment error across single-cell multi-omics modalities
in a shared latent or integrated space (Zhai et al., Genome Biology
2024; Liu et al., Nat Biotechnol 2023). A score of 0 represents perfect
alignment (where the true paired cell is the nearest neighbor), while
0.5 corresponds to random chance.

## Usage

``` r
calc_foscttm(x, y, metric = c("euclidean", "cosine"))
```

## Arguments

- x:

  Coordinates matrix for Modality 1 (cells x dimensions).

- y:

  Coordinates matrix for Modality 2 (cells x dimensions), where row i of
  y corresponds to the true match of row i of x.

- metric:

  Distance metric to evaluate: "euclidean" (default) or "cosine".

## Value

A list containing:

- foscttm:

  Bidirectional mean FOSCTTM score across all cells (lower is better, 0
  to 0.5).

- foscttm_xy:

  Directional FOSCTTM from Modality 1 to Modality 2.

- foscttm_yx:

  Directional FOSCTTM from Modality 2 to Modality 1.

- match_at_1:

  Top-1 match rate (proportion of cells where the true match is rank 1).

- match_at_5:

  Top-5 match rate (proportion of cells where the true match is in top
  5).

- cell_foscttm:

  Vector of bidirectional FOSCTTM scores for individual cells.
