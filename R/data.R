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

#' Real Paired Single-Cell Multiomics Benchmark Dataset (Rice Bud, Seurat v5 Option B)
#'
#' A real, pre-processed paired single-cell RNA-seq and ATAC-seq benchmark dataset of
#' rice bud tissue (\emph{Oryza sativa}, GEO accession GSE232863), packaged in
#' \strong{Option B (Strict S4 Standardization)} using \code{Seurat} v5 objects.
#' Contains paired empirical reference profiles (200 shared top-expressed genes x 100 cells for RNA,
#' 200 shared top-accessible peaks x 100 cells for ATAC) across 4 major biological cell types
#' (\code{Bud_primordium}, \code{Epidermis}, \code{Parenchyma_cell}, \code{Root_meristem}, 25 cells each)
#' and 2 biological batches (\code{Bud_21}, \code{Bud_22}), alongside simulated outputs from 6 major
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
#'   \item{reference}{List containing Option B standardized \code{Seurat} v5 objects:
#'     \code{combined} (multi-modal \code{Seurat} with \code{RNA} 200 x 100, \code{ATAC} 200 x 100, \code{cell_type}, \code{batch}, and \code{integrated}/\code{pca} reductions),
#'     \code{rna} (unimodal \code{Seurat} object for scRNA-seq / Shiny upload), and
#'     \code{atac} (unimodal \code{Seurat} object for scATAC-seq / Shiny upload).}
#'   \item{simulated}{Named list of simulated outputs for 6 simulators: \code{scMultiSim}, \code{scDesign3}, \code{MOSim}, \code{scCross}, \code{scMoMtF}, and \code{Matilda}. Each simulator entry contains \code{combined}, \code{rna}, and \code{atac} \code{Seurat} objects matching the reference structure.}
#'   \item{scalability_benchmarks}{A data.frame with benchmarked \code{Method}, \code{Runtime_Seconds}, and \code{PeakRAM_MiB}.}
#' }
#' @source Gene Expression Omnibus (GEO) GSE232863. Preprocessed via standard rice tissue pipeline.
#' @examples
#' data(example_rice_bud)
#' example_rice_bud$reference$combined
#' dim(example_rice_bud$reference$rna)
#' dim(example_rice_bud$reference$atac)
#' names(example_rice_bud$simulated)
#' \donttest{
#' evaluate_multiomics_accuracy(
#'   ref_multi = example_rice_bud$reference$combined,
#'   sim_multi = example_rice_bud$simulated$scMultiSim$combined,
#'   pairing   = "paired"
#' )
#' }
"example_rice_bud"

#' Real Unpaired Single-Cell Multiomics Benchmark Dataset (Soybean Root, Seurat v5 Option B)
#'
#' A real, pre-processed unpaired single-cell RNA-seq and ATAC-seq benchmark dataset of
#' soybean root tissue (\emph{Glycine max}), packaged in \strong{Option B (Strict S4 Standardization)}
#' using \code{Seurat} v5 objects. Contains independent empirical reference profiles
#' (200 shared top-expressed genes x 100 RNA cells, 200 shared top-accessible peaks x 100 ATAC cells)
#' across 4 major biological cell types (\code{Cortex}, \code{Endodermis}, \code{Epidermis}, \code{Procambium},
#' 25 RNA + 25 ATAC cells each) and 2 biological replicates (\code{rep1}, \code{rep2}), alongside
#' simulated outputs from 3 unpaired multiomics simulators (\code{scMultiSim}, \code{scDesign3},
#' and \code{scCross}) and observed scalability benchmarks.
#'
#' @docType data
#' @name example_soybean_root
#' @usage data(example_soybean_root)
#' @format A list with 4 elements:
#' \describe{
#'   \item{metadata}{List containing organism, tissue, data type (Unpaired), and description.}
#'   \item{reference}{List containing Option B standardized \code{Seurat} v5 objects:
#'     \code{combined} (multi-modal \code{Seurat} spanning 200 cells with \code{RNA} 200 x 100, \code{ATAC} 200 x 100, \code{cell_type}, \code{batch}, \code{modality}, and \code{integrated}/\code{joint.cca.umap} reductions),
#'     \code{rna} (unimodal \code{Seurat} object for scRNA-seq / Shiny upload), and
#'     \code{atac} (unimodal \code{Seurat} object for scATAC-seq / Shiny upload).}
#'   \item{simulated}{Named list of simulated outputs for 3 simulators: \code{scMultiSim}, \code{scDesign3}, and \code{scCross}. Each simulator entry contains \code{combined}, \code{rna}, and \code{atac} \code{Seurat} objects matching the reference structure.}
#'   \item{scalability_benchmarks}{A data.frame with benchmarked runtime and peak RAM for RNA, ATAC, and combined execution.}
#' }
#' @source Preprocessed from soybean root single-cell atlas datasets.
#' @examples
#' data(example_soybean_root)
#' example_soybean_root$reference$combined
#' dim(example_soybean_root$reference$rna)
#' dim(example_soybean_root$reference$atac)
#' names(example_soybean_root$simulated)
#' \donttest{
#' evaluate_multiomics_accuracy(
#'   ref_multi = list(
#'     "scRNA-seq"  = example_soybean_root$reference$rna,
#'     "scATAC-seq" = example_soybean_root$reference$atac
#'   ),
#'   sim_multi = list(
#'     "scRNA-seq"  = example_soybean_root$simulated$scMultiSim$rna,
#'     "scATAC-seq" = example_soybean_root$simulated$scMultiSim$atac
#'   ),
#'   pairing = "unpaired"
#' )
#' }
"example_soybean_root"

