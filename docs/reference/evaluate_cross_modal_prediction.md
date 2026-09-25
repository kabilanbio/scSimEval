# Cross-Modal Cell-Type Label Transfer Accuracy

Trains a model using ONLY Modality 1 (e.g., scATAC-seq) and predicts
cell-type identities defined in Modality 2 (e.g., scRNA-seq) to test
cross-modal biological coherence.

## Usage

``` r
evaluate_cross_modal_prediction(mod1_data, cell_types)
```

## Arguments

- mod1_data:

  Matrix for Modality 1 (features x cells).

- cell_types:

  Ground truth cell-type vector from Modality 2.

## Value

A list with cross-modal classification accuracy and macro F1 score.
