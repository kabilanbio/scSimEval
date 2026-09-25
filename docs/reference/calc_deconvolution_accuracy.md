# Evaluate Cell-Type Deconvolution & Mixture Proportion Accuracy

Evaluates the accuracy of cell-type mixture proportion estimation /
deconvolution methods between empirical reference and simulated cell
populations (Chen et al., Bioinformatics 2021). Can accept precomputed
proportion matrices OR cell-type vectors (and optional batch labels) to
calculate proportion preservation across batches.

## Usage

``` r
calc_deconvolution_accuracy(
  true_proportions,
  estimated_proportions,
  batch_ref = NULL,
  batch_sim = NULL
)
```

## Arguments

- true_proportions:

  Vector, matrix, or data.frame of reference cell-type proportions, or
  factor/vector of reference cell types.

- estimated_proportions:

  Vector, matrix, or data.frame of estimated / simulated cell-type
  proportions, or factor/vector of simulated cell types.

- batch_ref:

  Optional batch labels for reference cells (if true_proportions is a
  cell-type vector).

- batch_sim:

  Optional batch labels for simulated cells (if estimated_proportions is
  a cell-type vector).

## Value

A list containing RMSE, MAE, Pearson correlation (r), Spearman
correlation (rho), Jensen-Shannon Divergence (JSD), and Total Variation
Distance (TVD).
