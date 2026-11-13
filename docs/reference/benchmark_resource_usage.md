# Benchmark Execution Runtime and Memory Usage

Times an expression and measures change in memory allocation, computing
both raw runtime/memory metrics and core-normalized Category VIII
scalability indicators (`core_seconds` and `memory_per_core_mib`).

## Usage

``` r
benchmark_resource_usage(expr, sim_cores = NULL)
```

## Arguments

- expr:

  Expression or function to benchmark.

- sim_cores:

  Optional positive integer specifying CPU cores used during execution
  (default `1L`).

## Value

A list with elapsed seconds/minutes, peak memory allocated (MiB and MB),
`core_seconds`, `core_hours`, `memory_per_core_mib`,
`memory_per_core_gib`, and `sim_cores`.

## References

Lannelongue, L., et al. (2021). Green Algorithms. *Advanced Science*,
8(12), 2100707. Holtgrewe, M., et al. (2015). *Bioinformatics*, 31(12),
1904–1912.

## Examples

``` r
benchmark_resource_usage({ Sys.sleep(0.01); 1 + 1 })
#> $result
#> [1] 2
#> 
#> $core_seconds
#> [1] 0.03
#> 
#> $core_hours
#> [1] 8.333333e-06
#> 
#> $memory_per_core_mib
#> [1] 0.00440979
#> 
#> $memory_per_core_gib
#> [1] 4.306436e-06
#> 
#> $memory_per_core_mb
#> [1] 0.004624
#> 
#> $memory_per_core_gb
#> [1] 4.624e-06
#> 
#> $elapsed_seconds
#> [1] 0.03
#> 
#> $elapsed_minutes
#> [1] 5e-04
#> 
#> $memory_mb
#> [1] 0.00440979
#> 
#> $memory_mib
#> [1] 0.00440979
#> 
#> $memory_mb_si
#> [1] 0.004624
#> 
```
