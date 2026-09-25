# Detect Differential Expression (DE) Genes (Mean Shift via Limma or T-Test)

Detect Differential Expression (DE) Genes (Mean Shift via Limma or
T-Test)

## Usage

``` r
calc_signal_de(exprs_mat, cell_types, p_sig = 0.05)
```

## Arguments

- exprs_mat:

  Log-normalized matrix (genes x cells).

- cell_types:

  Factor or binary vector of 2 cell types.

- p_sig:

  Significance threshold (default 0.05).

## Value

Vector of adjusted p-values.
