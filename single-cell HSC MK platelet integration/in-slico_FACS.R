# In-silico FACS for HSCs, megakaryocytes and platelets.
# Source this file from HSC_Megakaryocyte_Platelet_altas_integration.R.

hsc_markers <- c("SPINK2", "CYTL1", "EGFL7", "GATA1", "GATA2", "CD34")
mk_platelet_markers <- c("PF4", "PPBP", "MAST1", "ITGA2B", "GP9", "GP1BA")

# Pure marker gate; usable on a small matrix without Seurat.
classify_marker_counts <- function(counts) {
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("Count matrix needs gene row names and cell column names.")
  }
  missing <- setdiff(c(hsc_markers, mk_platelet_markers), rownames(counts))
  if (length(missing)) {
    stop("Marker genes are missing after gene-symbol mapping: ",
         paste(missing, collapse = ", "))
  }
  detected <- function(gene, threshold) as.numeric(counts[gene, ]) > threshold
  hsc_positive <- Reduce(`+`, lapply(hsc_markers, detected, threshold = 0))
  hsc <- hsc_positive > 3L
  mk_platelet <- Reduce(`|`, lapply(mk_platelet_markers, detected, threshold = 1))
  label <- rep("Neither", ncol(counts))
  label[hsc] <- "HSC"
  label[mk_platelet] <- "MK_or_Platelet"
  label[hsc & mk_platelet] <- "HSC_and_MK_or_Platelet"
  names(label) <- colnames(counts)
  label
}

# Keep cells with at least 250 raw UMIs, then apply the marker gate.
# A sample is retained only when at least 10 target cells remain.
label_and_filter_sc_object <- function(sc_obj, filter_cells = TRUE,
                                       min_umi = 250L, min_target_cells = 10L,
                                       sample_col = "sample_id") {
  if (inherits(sc_obj, "Seurat")) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("SeuratObject is required.")
    counts <- SeuratObject::GetAssayData(sc_obj, assay = "RNA", layer = "counts")
    metadata <- sc_obj[[]]
  } else if (inherits(sc_obj, "SingleCellExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("SummarizedExperiment is required for SingleCellExperiment inputs.")
    }
    counts <- SummarizedExperiment::assay(sc_obj, "counts")
    metadata <- as.data.frame(SummarizedExperiment::colData(sc_obj))
  } else {
    stop("Input must be a Seurat or SingleCellExperiment object.")
  }
  if (!identical(colnames(counts), colnames(sc_obj))) {
    stop("RNA count columns must match object cell names and order.")
  }
  if (!sample_col %in% names(metadata)) stop("Missing metadata column: ", sample_col)
  sample_id <- as.character(metadata[[sample_col]])
  if (anyNA(sample_id) || any(!nzchar(sample_id))) {
    stop("Each cell needs a nonempty sample_id.")
  }

  umi <- colSums(counts)
  passes_umi <- umi >= min_umi
  label <- rep("Below_UMI_threshold", ncol(counts))
  if (any(passes_umi)) {
    label[passes_umi] <- classify_marker_counts(counts[, passes_umi, drop = FALSE])
  }
  names(label) <- colnames(counts)
  if (inherits(sc_obj, "Seurat")) {
    sc_obj$marker_condition <- label
  } else {
    SummarizedExperiment::colData(sc_obj)$marker_condition <- label
  }
  if (!filter_cells) return(sc_obj)

  target <- passes_umi & label %in% c("HSC", "MK_or_Platelet", "HSC_and_MK_or_Platelet")
  target_by_sample <- table(sample_id[target])
  retained_samples <- names(target_by_sample)[target_by_sample >= min_target_cells]
  keep <- target & sample_id %in% retained_samples
  message("Target cells after UMI and marker gates: ", sum(target),
          "; retained after sample minimum: ", sum(keep),
          " across ", length(retained_samples), " samples.")
  if (!any(keep)) return(NULL)
  if (inherits(sc_obj, "Seurat")) {
    base::subset(sc_obj, cells = colnames(sc_obj)[keep])
  } else {
    sc_obj[, keep]
  }
}
