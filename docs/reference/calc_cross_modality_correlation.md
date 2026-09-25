# Calculate Cross-Modality Correlation Fidelity

Evaluates whether cross-modality feature relationships (e.g.
peak-to-gene links, promoter accessibility vs. gene expression, or mRNA
vs. surface protein) are accurately captured by the simulation method
compared to the empirical reference.

## Usage

``` r
calc_cross_modality_correlation(
  ref_mod1,
  ref_mod2,
  sim_mod1,
  sim_mod2,
  feature_pairs = NULL,
  method = c("spearman", "pearson")
)
```

## Arguments

- ref_mod1:

  Reference matrix for Modality 1 (features x cells).

- ref_mod2:

  Reference matrix for Modality 2 (features x cells).

- sim_mod1:

  Simulated matrix for Modality 1.

- sim_mod2:

  Simulated matrix for Modality 2.

- feature_pairs:

  Optional 2-column data.frame/matrix of paired feature names/indices.

- method:

  Correlation method: "spearman" (default) or "pearson".

## Value

A named list of the 7 univariate accuracy metrics comparing the
reference vs. simulated cross-modality correlation distributions.
