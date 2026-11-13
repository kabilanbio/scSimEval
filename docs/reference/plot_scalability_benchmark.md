# Plot Multi-Dimensional Computational Scalability Benchmark

Evaluates computational resource efficiency across single-cell
simulators using core-normalized Category VIII benchmarks:
**Core-Seconds** (\\\text{Elapsed Time (s)} \times N\_{\text{cores}}\\;
Lannelongue et al., 2021; Wratten et al., 2021; Grealey et al., 2022)
and **Peak RAM per Core** (\\\text{Peak RAM} / N\_{\text{cores}}\\;
Holtgrewe et al., 2015; Boisvert et al., 2012; Chen et al., 2022),
alongside resource tradeoff frontiers, cost footprints, and simulation
throughput.

## Usage

``` r
plot_scalability_benchmark(
  benchmark_data,
  type = c("composite", "runtime", "memory", "tradeoff", "cost", "throughput",
    "efficiency"),
  layout = c("4panel", "6panel"),
  modality = c("total", "rna", "atac"),
  cell_count = 1000,
  palette = NULL,
  base_size = 11,
  time_unit = c("core_seconds", "core_hours", "seconds", "minutes"),
  memory_unit = c("MiB/core", "GiB/core", "MB/core", "GB/core", "MiB", "MB")
)
```

## Arguments

- benchmark_data:

  Benchmark summary table or named list containing Category VIII
  metrics.

- type:

  Character. Visualization type: `"composite"` (multi-panel dashboard),
  `"runtime"` (core-seconds CPU work barplot), `"memory"` (peak RAM per
  core barplot), `"tradeoff"` (core-seconds vs memory-per-core scatter
  curve), `"cost"` (combined resource footprint index), `"throughput"`
  (cells/sec), or `"efficiency"` (CPU/elapsed concurrency ratio if CPU
  time available). Default `"composite"`.

- layout:

  Character. For composite dashboards: `"4panel"` (runtime, memory,
  tradeoff, cost footprint) or `"6panel"` (adds throughput). Default
  `"4panel"`.

- modality:

  Character string specifying the modality layer to benchmark: `"total"`
  (default, combined runtime and memory), `"rna"` (scRNA-seq layer
  only), or `"atac"` (scATAC-seq layer only).

- cell_count:

  Integer. Number of cells processed in the benchmark, used to calculate
  simulation throughput (`cells / second`). Default `1000`.

- palette:

  Optional named character vector of simulator colors.

- base_size:

  Numeric. Base font size. Default `11`.

- time_unit:

  Character. Display unit for runtime: `"core_seconds"` (default),
  `"core_hours"`, `"seconds"`, or `"minutes"`.

- memory_unit:

  Character. Display unit for peak RAM per core: `"MiB/core"` (default),
  `"GiB/core"`, `"MB/core"`, or `"GB/core"`.

## Value

A `ggplot` or `patchwork` object.

## Details

![Computational Scalability and Resource
Footprint](figures/scalability_benchmark.png)

## References

Lannelongue, L., et al. (2021). Green Algorithms: Quantifying the carbon
footprint of computation. *Advanced Science*, 8(12), 2100707.

Holtgrewe, M., et al. (2015). Methods for the detection and assembly of
novel sequence. *Bioinformatics*, 31(12), 1904–1912.

Boisvert, S., et al. (2012). Ray Meta: scalable de novo metagenome
assembly. *Genome Biology*, 13(12), R122.

Chen, W., et al. (2022). STEPS 4.0: Fast and memory-efficient
simulations. *Frontiers in Neuroinformatics*, 16, 883742.

## Examples

``` r
data(example_scrna, package = "scSimEval")
res <- evaluate_simulation_accuracy(example_scrna$ref, example_scrna$sim)
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
```
