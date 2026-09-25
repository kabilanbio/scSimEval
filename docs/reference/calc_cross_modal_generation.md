# Cross-Modality In Silico Generation & Translation Fidelity

Evaluates the biological accuracy and reconstruction fidelity of
cross-modality generation (e.g. predicting scATAC from scRNA or vice
versa; Zhai et al., Genome Biology 2024). Computes cell-wise and
feature-wise Pearson and Spearman correlations, cosine similarity, RMSE,
and MAE between predicted and measured multi-omics profiles.

## Usage

``` r
calc_cross_modal_generation(true_data, pred_data)
```

## Arguments

- true_data:

  Matrix or data.frame of measured features across cells (features x
  cells).

- pred_data:

  Matrix or data.frame of in-silico generated / predicted features
  (features x cells).

## Value

A list containing cell-wise and feature-wise correlation, cosine
similarity, RMSE, and MAE.
