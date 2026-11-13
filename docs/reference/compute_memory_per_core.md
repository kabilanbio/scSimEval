# Compute Peak RAM per Allocated CPU Core

Normalizes peak memory consumption by the number of allocated processor
cores (\\\text{Peak RAM} / N\_{\text{cores}}\\; Holtgrewe et al., 2015;
Boisvert et al., 2012; Chen et al., 2022).

## Usage

``` r
compute_memory_per_core(
  memory_val,
  memory_unit = "MiB",
  sim_cores = 1L,
  target_unit = c("MiB", "GiB", "MB", "GB")
)
```

## Arguments

- memory_val:

  Numeric peak memory consumption.

- memory_unit:

  Memory unit of `memory_val`: `"MiB"` (default) or `"MB"`.

- sim_cores:

  Positive integer specifying the number of CPU cores allocated (default
  `1L`).

- target_unit:

  Target memory unit: `"MiB"` (MiB/core, default), `"GiB"` (GiB/core),
  `"MB"` (MB/core), or `"GB"` (GB/core).

## Value

Numeric Peak RAM per allocated core in requested target unit.

## References

Holtgrewe, M., Kuchenbecker, L., & Reinert, K. (2015). Methods for the
detection and assembly of novel sequence in high-throughput sequencing
data. *Bioinformatics*, 31(12), 1904–1912.

Boisvert, S., Raymond, F., Godzaridis, E., Laviolette, F., & Corbeil, J.
(2012). Ray Meta: scalable de novo metagenome assembly and profiling.
*Genome Biology*, 13(12), R122.

Chen, W., et al. (2022). STEPS 4.0: Fast and memory-efficient molecular
simulations of neurons at the nanoscale. *Frontiers in
Neuroinformatics*, 16, 883742.
