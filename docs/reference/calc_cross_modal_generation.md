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

## Examples

``` r
data(example_multiomics, package = "scSimEval")
r_rna <- example_multiomics$ref_multi$rna
r_atac <- example_multiomics$ref_multi$atac
calc_cross_modal_generation(r_rna, r_atac)
#> $mean_cell_pcc
#> [1] 0.05812026
#> 
#> $median_cell_pcc
#> [1] 0.06043748
#> 
#> $mean_cell_scc
#> [1] 0.07055168
#> 
#> $median_cell_scc
#> [1] 0.07361193
#> 
#> $mean_feat_pcc
#> [1] 0.04695087
#> 
#> $median_feat_pcc
#> [1] 0.05451325
#> 
#> $mean_feat_scc
#> [1] 0.05410635
#> 
#> $median_feat_scc
#> [1] 0.06122589
#> 
#> $mean_cell_cosine
#> [1] 0.5955659
#> 
#> $rmse
#> [1] 10.3143
#> 
#> $mae
#> [1] 8.019169
#> 
```
