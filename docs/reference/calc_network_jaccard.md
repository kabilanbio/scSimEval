# Jaccard Similarity of Predicted Regulatory Networks / Links

Computes the Jaccard similarity index of top-k edge predictions between
two networks or top correlated feature pairs between reference and
simulated datasets.

## Usage

``` r
calc_network_jaccard(edges1, edges2, k = NULL, directed = TRUE)
```

## Arguments

- edges1:

  Data frame or character vector representing the first network.

- edges2:

  Data frame or character vector representing the second network.

- k:

  Integer cutoff for top edge selection (optional).

- directed:

  Logical, whether edges are directed (default TRUE) or undirected
  (FALSE).

## Value

Jaccard similarity index in \[0, 1\] (\|E1 cap E2\| / \|E1 cup E2\|).
