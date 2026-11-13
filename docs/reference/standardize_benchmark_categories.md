# Standardize Category Names in Benchmark Tables

Standardizes legacy, informal, or variant evaluation category labels to
the canonical 8 scSimEval categories across benchmark data frames.

## Usage

``` r
standardize_benchmark_categories(df)
```

## Arguments

- df:

  Data frame containing a `Category` and/or `Metric` column.

## Value

Data frame with standardized `Category` assignments.

## Examples

``` r
df <- data.frame(Category = "Distribution", Metric = "KS", Value = 0.1)
standardize_benchmark_categories(df)
#>                        Category Metric Value
#> 1 (I) Distributional Properties     KS   0.1
```
