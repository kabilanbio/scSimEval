# Modality Alignment & Omics Layer Mixing in Joint Latent Space

Quantifies how well different single-cell omics layers mix in a unified
latent embedding space without sacrificing biological clustering (Zhai
et al., Genome Biology 2024; Luecken et al., Nat Methods 2022).

## Usage

``` r
calc_modality_alignment(embedding, modalities, cell_types = NULL, k = 15)
```

## Arguments

- embedding:

  Joint coordinates matrix (cells x latent_dimensions).

- modalities:

  Factor or character vector indicating the modality of each cell (e.g.
  "RNA" vs "ATAC").

- cell_types:

  Optional factor or character vector of cell type labels.

- k:

  Number of nearest neighbors for neighborhood connectivity calculation
  (default 15).

## Value

A list containing modality_asw, modality_mixing_score, and
mean_cross_modality_neighbor_frac.

## Examples

``` r
data(example_multiomics, package = "scSimEval")
r_rna <- example_multiomics$ref_multi$rna
r_atac <- example_multiomics$ref_multi$atac
joint_emb <- rbind(t(r_rna[seq_len(10), ]), t(r_atac[seq_len(10), ]))
calc_modality_alignment(joint_emb, rep(c("RNA", "ATAC"), each = ncol(r_rna)))
#> $modality_asw
#> [1] 0.3818172
#> 
#> $modality_mixing_score
#> [1] 0.6181828
#> 
#> $mean_cross_modality_neighbor_frac
#> [1] 0.0825
#> 
```
