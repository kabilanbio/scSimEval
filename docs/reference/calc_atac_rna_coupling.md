# Chromatin Accessibility to RNA Regulatory Effect Coupling

Evaluates whether simulated paired scATAC-seq and scRNA-seq profiles
capture the true regulatory coupling where regional chromatin opening
gates target gene transcription (scMultiSim; Li et al., Nat Methods
2023).

## Usage

``` r
calc_atac_rna_coupling(atac_data, rna_data, linked_pairs = NULL)
```

## Arguments

- atac_data:

  Matrix of chromatin peak accessibilities (peaks x cells).

- rna_data:

  Matrix of gene expression counts (genes x cells).

- linked_pairs:

  Optional 2-column data.frame of known linked peak-gene pairs. If NULL,
  assumes 1-to-1 matching by row index.

## Value

A list containing mean coupling correlation, positive coupling ratio,
and mean R-squared.

## Examples

``` r
data(example_multiomics, package = "scSimEval")
r_rna <- example_multiomics$ref_multi$rna
r_atac <- example_multiomics$ref_multi$atac
calc_atac_rna_coupling(r_rna, r_atac)
#> $mean_coupling_cor
#> [1] 0.04695087
#> 
#> $median_coupling_cor
#> [1] 0.05451325
#> 
#> $positive_coupling_ratio
#> [1] 0.7333333
#> 
#> $mean_r_squared
#> [1] 0.01642647
#> 
#> $coupling_correlations
#>  [1]  0.143080984  0.061246595  0.024165441  0.063673642 -0.130450929
#>  [6]  0.054906216  0.040003354  0.284619033 -0.161589078  0.092782956
#> [11]  0.056827722 -0.150152159  0.065508545  0.072124391  0.091422074
#> [16]  0.186972993  0.084977247  0.220297687 -0.051546027  0.007535943
#> [21] -0.151189660  0.263398545  0.138251376 -0.179604629 -0.150785742
#> [26]  0.009952578  0.091184234  0.005294338  0.121019590 -0.030406654
#> [31] -0.145036850  0.223321070  0.116178412  0.017882860  0.250277479
#> [36] -0.084308396 -0.199028063  0.105880280  0.114525236  0.083386721
#> [41]  0.247704631  0.172547802 -0.107081839 -0.020605772  0.015767972
#> [46]  0.054120277 -0.053465517  0.037649296  0.032179824  0.011030644
#> [51]  0.043811622  0.024359015  0.093292879  0.013863964 -0.039400714
#> [56]  0.094970084  0.320797165  0.107757941 -0.022786788  0.133940401
#> 
```
