#' Synthetic Example scRNA-seq Benchmarking Dataset
#'
#' A synthetic single-cell RNA sequencing dataset designed for testing, benchmarking,
#' and demonstrating evaluation metrics in \code{scSimEval}. Contains paired
#' empirical reference counts and simulated counts across 60 genes and 80 cells
#' from two cell types and two batches.
#'
#' @docType data
#' @name example_scrna
#' @usage data(example_scrna)
#' @format A list with 4 elements:
#' \describe{
#'   \item{ref}{A numeric matrix of reference count data (60 genes x 80 cells).}
#'   \item{sim}{A numeric matrix of simulated count data (60 genes x 80 cells).}
#'   \item{cell_types}{A factor of cell-type identities ("TypeA", "TypeB") for the 80 cells.}
#'   \item{batch_info}{A factor of batch annotations ("Batch1", "Batch2") for the 80 cells.}
#' }
#' @source Synthetic benchmark data generated via negative binomial sampling with
#'   cell-type specific expression shifts.
#' @examples
#' data(example_scrna)
#' dim(example_scrna$ref)
#' dim(example_scrna$sim)
#' table(example_scrna$cell_types)
"example_scrna"

#' Synthetic Example scATAC-seq Benchmarking Dataset
#'
#' A synthetic single-cell ATAC sequencing peak accessibility dataset designed for testing,
#' benchmarking, and demonstrating evaluation metrics in \code{scSimEval}. Contains paired
#' empirical reference accessibility counts and simulated accessibility counts across 60 peaks
#' and 80 cells matching the cells in \code{example_scrna}.
#'
#' @docType data
#' @name example_scatac
#' @usage data(example_scatac)
#' @format A list with 4 elements:
#' \describe{
#'   \item{ref}{A numeric matrix of reference peak accessibility counts (60 peaks x 80 cells).}
#'   \item{sim}{A numeric matrix of simulated peak accessibility counts (60 peaks x 80 cells).}
#'   \item{cell_types}{A factor of cell-type identities ("TypeA", "TypeB") for the 80 cells.}
#'   \item{batch_info}{A factor of batch annotations ("Batch1", "Batch2") for the 80 cells.}
#' }
#' @source Synthetic benchmark data generated via binomial sampling with
#'   cell-type specific accessibility shifts.
#' @examples
#' data(example_scatac)
#' dim(example_scatac$ref)
#' dim(example_scatac$sim)
#' table(example_scatac$batch_info)
"example_scatac"

#' Synthetic Example Multiomics Benchmarking Dataset
#'
#' A unified single-cell multiomics dataset pairing scRNA-seq and scATAC-seq
#' profiles across the same 80 cells, accompanied by cell-type labels, batch annotations,
#' and logged computational scalability metrics.
#'
#' @docType data
#' @name example_multiomics
#' @usage data(example_multiomics)
#' @format A list with 5 elements:
#' \describe{
#'   \item{ref_multi}{Named list of reference matrices: \code{rna} (60 x 80) and \code{atac} (60 x 80).}
#'   \item{sim_multi}{Named list of simulated matrices: \code{rna} (60 x 80) and \code{atac} (60 x 80).}
#'   \item{cell_types}{A factor of cell-type identities ("TypeA", "TypeB") for the 80 cells.}
#'   \item{batch_info}{A factor of batch annotations ("Batch1", "Batch2") for the 80 cells.}
#'   \item{resource_stats}{A list of computational scalability metrics: \code{cpu_time}, \code{memory_mb}, \code{system_time}, \code{elapsed_time}.}
#' }
#' @source Synthetic multiomics benchmark data generated for testing \code{scSimEval}.
#' @examples
#' data(example_multiomics)
#' names(example_multiomics$ref_multi)
#' names(example_multiomics$sim_multi)
"example_multiomics"
