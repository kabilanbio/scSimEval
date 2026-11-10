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
#'   \item{ref}{A \code{SingleCellExperiment} object containing reference count data (60 genes x 80 cells), cell types and batch annotations in \code{colData}, and precomputed PCA reductions in \code{reducedDims}.}
#'   \item{sim}{A \code{SingleCellExperiment} object containing simulated count data (60 genes x 80 cells), cell types and batch annotations in \code{colData}, and precomputed PCA reductions in \code{reducedDims}.}
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
#'   \item{ref}{A \code{SingleCellExperiment} object containing reference peak accessibility counts (60 peaks x 80 cells), cell types and batch annotations in \code{colData}, and precomputed PCA reductions in \code{reducedDims}.}
#'   \item{sim}{A \code{SingleCellExperiment} object containing simulated peak accessibility counts (60 peaks x 80 cells), cell types and batch annotations in \code{colData}, and precomputed PCA reductions in \code{reducedDims}.}
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
#'   \item{ref_multi}{Named list of \code{SingleCellExperiment} objects: \code{rna} (60 x 80) and \code{atac} (60 x 80), complete with precomputed PCA reductions in \code{reducedDims}.}
#'   \item{sim_multi}{Named list of \code{SingleCellExperiment} objects: \code{rna} (60 x 80) and \code{atac} (60 x 80), complete with precomputed PCA reductions in \code{reducedDims}.}
#'   \item{cell_types}{A factor of cell-type identities ("TypeA", "TypeB") for the 80 cells.}
#'   \item{batch_info}{A factor of batch annotations ("Batch1", "Batch2") for the 80 cells.}
#'   \item{resource_stats}{A list of computational scalability metrics: \code{cpu_time}, \code{memory_mb}, \code{system_time}, \code{elapsed_time}.}
#' }
#' @source Synthetic multiomics benchmark data generated for testing \code{scSimEval}.
#' @examples
#' data(example_multiomics)
#' names(example_multiomics$ref_multi)
#' names(example_multiomics$sim_multi)
#' @export
"example_multiomics"

#' Example Cell-Type Annotations Vector
#'
#' A named factor vector containing discrete cell-type labels ("TypeA", "TypeB")
#' for 80 single cells, matching the column names of \code{\link{example_scrna}} and
#' \code{\link{example_scatac}}. Useful as a reference template for formatting cell-type
#' metadata inputs in \code{scSimEval}.
#'
#' @docType data
#' @name example_cell_types
#' @usage data(example_cell_types)
#' @format A named factor with 80 elements across 2 levels: "TypeA" (40 cells) and "TypeB" (40 cells).
#'   Names correspond to individual cell barcodes (\code{"Cell_01"} to \code{"Cell_80"}).
#' @source Derived from \code{example_scrna$cell_types}.
#' @examples
#' data(example_cell_types)
#' head(example_cell_types)
#' table(example_cell_types)
#' @export
"example_cell_types"

#' Example Batch Annotations Vector
#'
#' A named factor vector containing technical batch annotations ("Batch1", "Batch2")
#' for 80 single cells, matching the column names of \code{\link{example_scrna}} and
#' \code{\link{example_scatac}}. Useful as a reference template for formatting batch
#' metadata inputs in \code{scSimEval}.
#'
#' @docType data
#' @name example_batch_info
#' @usage data(example_batch_info)
#' @format A named factor with 80 elements across 2 levels: "Batch1" (40 cells) and "Batch2" (40 cells).
#'   Names correspond to individual cell barcodes (\code{"Cell_01"} to \code{"Cell_80"}).
#' @source Derived from \code{example_scrna$batch_info}.
#' @examples
#' data(example_batch_info)
#' head(example_batch_info)
#' table(example_batch_info)
#' @export
"example_batch_info"

#' Real Paired Single-Cell Multiomics Benchmark Dataset (Rice Bud)
#'
#' A real, pre-processed paired single-cell RNA-seq and ATAC-seq benchmark dataset of
#' rice bud tissue (\emph{Oryza sativa}, GEO accession GSE232863). Contains paired
#' empirical reference profiles (200 genes x 100 cells for RNA, 200 peaks x 100 cells for ATAC)
#' across 8 biological cell types and 2 batches, alongside simulated outputs from 6 major
#' single-cell multiomics simulators (\code{scMultiSim}, \code{scDesign3}, \code{MOSim},
#' \code{scCross}, \code{scMoMtF}, and \code{Matilda}) and observed scalability benchmarks
#' (runtime in seconds and peak RAM in MiB).
#'
#' @docType data
#' @name example_rice_bud
#' @usage data(example_rice_bud)
#' @format A list with 4 elements:
#' \describe{
#'   \item{metadata}{List containing organism, tissue, GEO accession (GSE232863), data type (Paired), and description.}
#'   \item{reference}{List containing empirical reference sparse matrices \code{rna} (200 x 100) and \code{atac} (200 x 100), \code{cell_types} (8 cell types), and \code{batches} (2 batches).}
#'   \item{simulated}{Named list of simulated outputs for 6 simulators: \code{scMultiSim}, \code{scDesign3}, \code{MOSim}, \code{scCross}, \code{scMoMtF}, and \code{Matilda}, each with \code{rna} and \code{atac} count matrices.}
#'   \item{scalability_benchmarks}{A data.frame with benchmarked \code{Method}, \code{Runtime_Seconds}, and \code{PeakRAM_MiB}.}
#' }
#' @source Gene Expression Omnibus (GEO) GSE232863. Preprocessed via standard rice tissue pipeline.
#' @examples
#' data(example_rice_bud)
#' dim(example_rice_bud$reference$rna)
#' dim(example_rice_bud$reference$atac)
#' names(example_rice_bud$simulated)
"example_rice_bud"

#' Real Unpaired Single-Cell Multiomics Benchmark Dataset (Soybean Root)
#'
#' A real, pre-processed unpaired single-cell RNA-seq and ATAC-seq benchmark dataset of
#' soybean root tissue (\emph{Glycine max}). Contains independent empirical reference profiles
#' (200 genes x 100 cells for RNA, 200 peaks x 100 cells for ATAC) with modality-specific
#' cell-type and batch annotations, alongside simulated outputs from 3 multiomics simulators
#' (\code{scMultiSim}, \code{scDesign3}, and \code{scCross}) and observed
#' scalability benchmarks.
#'
#' @docType data
#' @name example_soybean_root
#' @usage data(example_soybean_root)
#' @format A list with 4 elements:
#' \describe{
#'   \item{metadata}{List containing organism, tissue, data type (Unpaired), and description.}
#'   \item{reference}{List containing empirical reference matrices \code{rna} (200 x 100), \code{atac} (200 x 100), \code{rna_cell_types}, \code{rna_batches}, \code{atac_cell_types}, and \code{atac_batches}.}
#'   \item{simulated}{Named list of simulated outputs for 3 simulators: \code{scMultiSim}, \code{scDesign3}, and \code{scCross}, each with \code{rna} and \code{atac} count matrices.}
#'   \item{scalability_benchmarks}{A data.frame with benchmarked runtime and peak RAM for RNA, ATAC, and combined execution.}
#' }
#' @source Preprocessed from soybean root single-cell atlas datasets.
#' @examples
#' data(example_soybean_root)
#' dim(example_soybean_root$reference$rna)
#' dim(example_soybean_root$reference$atac)
#' names(example_soybean_root$simulated)
"example_soybean_root"

