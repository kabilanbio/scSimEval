# Filter Benchmark Data to the Curated Evaluation Measures Suite

Filters any benchmark results table to the curated 56 evaluation
measures across all 8 canonical categories (while preserving legacy
benchmark table compatibility), eliminating unrequested or collinear
duplicate metrics to prevent distortion of simulator fidelity
leaderboards.

## Usage

``` r
filter_gold_standard_metrics(df)
```

## Arguments

- df:

  A data.frame of benchmark metrics or a scSimEval benchmark session
  object.

## Value

Cleaned data.frame containing the curated evaluation measures.

## Examples

``` r
data(example_scrna, package = "scSimEval")
```
