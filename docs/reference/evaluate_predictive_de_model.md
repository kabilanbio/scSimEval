# Predictive Cell Identity Classification Using Simulated DEGs (Simpipe Approach 3)

Constructs classification models (Support Vector Machine with radial
kernel, Random Forest, or k-NN) using simulated gene expression matrices
and candidate DEGs (Duo et al., Genome Biology 2024 / Simpipe). Genes
with zero variance are filtered out, and DEGs are used to train on 80%
of cells and evaluate on 20% test cells (or cross-validation). Computes
Accuracy, Macro-Precision, Macro-Recall, and Macro-F1 score.

## Usage

``` r
evaluate_predictive_de_model(
  data,
  group,
  de_features = NULL,
  n_top = 50,
  method = c("svm", "knn", "rf")
)
```

## Arguments

- data:

  Matrix of counts or normalized expressions (genes x cells).

- group:

  Cell type / cluster labels.

- de_features:

  Optional character vector of candidate DEG feature names. If NULL, top
  variable features are used.

- n_top:

  Number of top features to select if de_features is NULL (default 50).

- method:

  Classifier algorithm: "svm" (default in Simpipe), "knn", or "rf".

## Value

A list containing accuracy, precision, recall, and F1 score.

## Examples

``` r
data(example_scrna, package = "scSimEval")
evaluate_predictive_de_model(example_scrna$ref, group = example_scrna$cell_types)
#> $accuracy
#> [1] 0.5625
#> 
#> $precision
#> [1] 0.5833333
#> 
#> $recall
#> [1] 0.5793651
#> 
#> $F1
#> [1] 0.5607843
#> 
```
