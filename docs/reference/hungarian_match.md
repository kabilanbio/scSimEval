# Hungarian Maximum Weight Matching for Cluster Assignment

Matches predicted cluster labels to ground truth labels using the
Hungarian algorithm (Kuhn-Munkres) on the contingency/F1 matrix. Adapted
from HelenaLC/simulation-comparison and simpipe.

## Usage

``` r
hungarian_match(pred, truth)
```

## Arguments

- pred:

  Vector of predicted cluster labels.

- truth:

  Vector of ground truth cell-type labels.

## Value

A list with matched cluster pairs, precision, recall, and macro F1
score.
