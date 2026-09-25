# Compute Clustering Deviation Index (CDI) for Single-Cell Clustering Evaluation

Evaluates the deviation and goodness-of-fit of single-cell clustering
labels directly from raw or normalized count matrices without requiring
ground-truth cell-type labels (Fang et al., Genome Biology 2022;
benchmarked in scDesign3, Song et al., Nat Biotechnol 2024). CDI
calculates the penalized log-likelihood (CDI-AIC and CDI-BIC) of count
data under a cluster-specific Poisson or Negative Binomial GLM with cell
size factors.

## Usage

``` r
calc_cdi(
  counts,
  cluster_labels,
  size_factors = NULL,
  model = c("poisson", "nb"),
  top_features = NULL
)
```

## Arguments

- counts:

  Matrix or data.frame of counts (genes x cells or cells x genes).

- cluster_labels:

  Vector of cluster assignments for each cell.

- size_factors:

  Optional numeric vector of cell-specific size factors. If NULL,
  computed as library size scaled by median library size.

- model:

  Distribution model to use: "poisson" (fast, default) or "nb" (Negative
  Binomial).

- top_features:

  Optional integer; if set, restricts CDI calculation to the top most
  variable features to accelerate computation on large matrices.

## Value

A list containing:

- cdi_aic:

  Normalized CDI score based on AIC penalty (lower indicates better
  clustering).

- cdi_bic:

  Normalized CDI score based on BIC penalty (favors parsimonious main
  clusters).

- loglik:

  Total log-likelihood under the cluster-specific count model.

- deviance:

  Model deviance relative to saturated count model.

- n_clusters:

  Number of clusters evaluated.

- n_params:

  Number of free parameters in the clustering model.
