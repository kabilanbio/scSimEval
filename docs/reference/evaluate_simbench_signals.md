# Evaluate the 5 SimBench Biological Signal Proportions

Compares proportions of genes exhibiting DE, DV, DD, DP, and BD between
empirical reference and simulated datasets.

## Usage

``` r
evaluate_simbench_signals(
  ref_mat,
  sim_mat,
  ref_celltypes,
  sim_celltypes,
  p_sig = 0.05,
  bi_cutoff = 0.3
)
```

## Arguments

- ref_mat:

  Reference matrix.

- sim_mat:

  Simulated matrix.

- ref_celltypes:

  Reference cell type labels (subsets to top 2 abundant types).

- sim_celltypes:

  Simulated cell type labels.

- p_sig:

  Significance cutoff (default 0.05).

- bi_cutoff:

  Bimodal index cutoff (default 0.3).

## Value

A tidy data.frame comparing biological signal proportions.
