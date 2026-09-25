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
