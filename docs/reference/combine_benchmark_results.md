# Combine and Append Multiple Single-Cell Benchmark Results

Harmonizes, standardizes, and consolidates multiple single-cell
simulator evaluation datasets or saved benchmark sessions into a unified
multi-dataset benchmark object.

## Usage

``` r
combine_benchmark_results(
  ...,
  list_of_benchmarks = NULL,
  dataset_names = NULL,
  filter_gold_standard = TRUE,
  recompute_normalized = TRUE
)

append_benchmark_results(
  ...,
  list_of_benchmarks = NULL,
  dataset_names = NULL,
  filter_gold_standard = TRUE,
  recompute_normalized = TRUE
)
```

## Arguments

- ...:

  One or more benchmark evaluation objects, tidy data frames, or file
  paths (ending in `.rds` or `.csv`).

- list_of_benchmarks:

  Optional list of benchmark objects or data frames.

- dataset_names:

  Optional character vector specifying custom names for each dataset
  being combined. If `NULL`, dataset names are automatically extracted
  from the input objects or generated cleanly (e.g., `"Dataset 1"`,
  `"Dataset 2"`).

- filter_gold_standard:

  Logical. If `TRUE` (default), filters all combined datasets to the 48
  curated evaluation metrics.

- recompute_normalized:

  Logical. If `TRUE` (default), re-normalizes raw metric scores across
  all evaluated simulators and datasets using direction-aware scaling.

## Value

An S3 object of class
`c("scSimEval_multidataset_benchmark", "scSimEval_consolidated", "list")`
containing:

- `benchmark_summary_table`: Tidy data frame of cross-dataset average
  fidelity scores per simulator and metric (ready for
  [`plot_benchmark_bubble_matrix`](https://kabilanbio.github.io/scSimEval/reference/plot_benchmark_bubble_matrix.md)).

- `granular_benchmark_table`: Complete combined data frame preserving
  individual dataset rows.

- `dataset_summary_table`: Summary table of all combined datasets and
  simulator counts.

- `dataset_summary_df`: Alias of `dataset_summary_table`.

- `cross_dataset_stats`: Data frame of mean, SD, SE, and sample count
  per simulator-metric pair.

- `methods`: Unique simulator names evaluated across all datasets.

- `datasets`: Vector of distinct dataset names.

- `is_multidataset`: Logical indicating whether more than one dataset is
  included.

- `n_datasets`: Integer number of datasets.

- `n_simulators`: Integer number of unique simulators.

## Details

When the same simulator (e.g., `"scDesign3"` or `"scMultiSim"`) is
evaluated across multiple empirical datasets,
`combine_benchmark_results` maintains each distinct dataset evaluation
and provides cross-dataset aggregate fidelity scores across the 48
curated evaluation metrics. This enables rigorous comparative evaluation
and bubble matrix visualization of average simulator performance across
diverse biological contexts.

## Examples

``` r
# \donttest{
data(example_scrna, package = "scSimEval")
res1 <- evaluate_simulation_accuracy(example_scrna$ref, example_scrna$sim)
#> [1/5] Extracting cell-level properties...
#> [2/5] Extracting feature-level properties...
#> [3/5] Computing univariate accuracy metrics...
#> [4/5] Computing bivariate metrics (Fasano-Franceschini, KDE zstat, 2D EMD)...
#> [5/5] Computing zero-probability and manifold distances...
#> Evaluating batch integration & confounder mixing...
#> Evaluating scRNA-seq biological signal fidelity (Libra Pseudobulk edgeR-QLF)...
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> [6/7] Evaluating trajectory and lineage topology (HVG + dynverse, scRNA-seq only)...
#> Unimodal accuracy evaluation complete.
res2 <- evaluate_simulation_accuracy(example_scrna$ref, example_scrna$sim)
#> [1/5] Extracting cell-level properties...
#> [2/5] Extracting feature-level properties...
#> [3/5] Computing univariate accuracy metrics...
#> [4/5] Computing bivariate metrics (Fasano-Franceschini, KDE zstat, 2D EMD)...
#> [5/5] Computing zero-probability and manifold distances...
#> Evaluating batch integration & confounder mixing...
#> Evaluating scRNA-seq biological signal fidelity (Libra Pseudobulk edgeR-QLF)...
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> [6/7] Evaluating trajectory and lineage topology (HVG + dynverse, scRNA-seq only)...
#> Unimodal accuracy evaluation complete.
comb <- combine_benchmark_results(res1, res2, dataset_names = c("Cohort A", "Cohort B"))
#> Warning: Item 1 could not be parsed as a valid benchmark table; skipping.
#> Warning: Item 2 could not be parsed as a valid benchmark table; skipping.
#> Error in combine_benchmark_results(res1, res2, dataset_names = c("Cohort A",     "Cohort B")): None of the supplied inputs contained valid benchmark data.
print(comb$datasets)
#> Error: object 'comb' not found
# }
```
