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
