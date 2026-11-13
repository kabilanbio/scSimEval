# Compute Core-Seconds Runtime Metric

Quantifies cumulative CPU work across all allocated processing cores
(\\\text{Elapsed Time (s)} \times N\_{\text{cores}}\\; Lannelongue et
al., 2021; Wratten et al., 2021; Grealey et al., 2022).

## Usage

``` r
compute_core_seconds(elapsed_time, time_unit = "seconds", sim_cores = 1L)
```

## Arguments

- elapsed_time:

  Numeric elapsed wall-clock runtime.

- time_unit:

  Time unit of `elapsed_time`: `"seconds"` (default) or `"minutes"`.

- sim_cores:

  Positive integer specifying the number of CPU cores allocated (default
  `1L`).

## Value

Numeric Core-Seconds (\\\text{core}\cdot\text{s}\\).

## References

Lannelongue, L., Grealey, J., & Inouye, M. (2021). Green Algorithms:
Quantifying the carbon footprint of computation. *Advanced Science*,
8(12), 2100707.

Wratten, L., Wilm, A., & Göke, J. (2021). Reproducible, scalable, and
shareable analysis pipelines with bioinformatics workflow managers.
*Nature Methods*, 18, 1161–1168.

Grealey, J., et al. (2022). The Carbon Footprint of Bioinformatics.
*Molecular Biology and Evolution*, 39(3), msac034.
