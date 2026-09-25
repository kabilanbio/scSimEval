# Plot Multi-Dimensional Computational Scalability Benchmark

Evaluates computational resource efficiency across single-cell
simulators using raw benchmarking measurements (elapsed real wall-clock
time, total CPU time, peak resident Plot Multi-Dimensional Computational
Scalability Benchmark

## Usage

``` r
plot_scalability_benchmark(
  benchmark_data,
  type = c("composite", "runtime", "memory", "tradeoff", "cost", "throughput",
    "efficiency"),
  layout = c("4panel", "6panel"),
  cell_count = 1000,
  palette = NULL,
  base_size = 11
)
```

## Arguments

- benchmark_data:

  Benchmark summary table or named list containing Category VIII
  metrics.

- type:

  Character. Visualization type: `"composite"` (multi-panel dashboard),
  `"runtime"` (elapsed wall-clock time barplot), `"memory"` (peak RAM
  barplot), `"tradeoff"` (runtime vs memory scatter/line curve),
  `"cost"` (combined resource footprint index), `"throughput"`
  (cells/sec), or `"efficiency"` (CPU/elapsed concurrency ratio if CPU
  time available). Default `"composite"`.

- layout:

  Character. For composite dashboards: `"4panel"` (runtime, memory,
  tradeoff, cost footprint) or `"6panel"` (adds throughput). Default
  `"4panel"`.

- cell_count:

  Integer. Number of cells processed in the benchmark, used to calculate
  simulation throughput (`cells / second`). Default `1000`.

- palette:

  Optional named character vector of simulator colors.

- base_size:

  Numeric. Base font size. Default `11`.

## Value

A `ggplot` or `patchwork` object.

## Details

Evaluates computational resource efficiency across single-cell
simulators using raw benchmarking measurements (elapsed real wall-clock
time in seconds and peak resident memory consumption in MiB).

![Computational Scalability and Resource
Footprint](figures/scalability_benchmark.png)
