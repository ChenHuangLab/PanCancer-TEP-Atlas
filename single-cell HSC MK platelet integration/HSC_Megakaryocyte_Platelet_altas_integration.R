#!/usr/bin/env Rscript
# HSC / megakaryocyte / platelet atlas integration from prepared Seurat RDS files.
# The filename follows the exploratory source script for traceability.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2L != 0L ||
    (length(args) > 0L && any(!args[seq(1L, length(args), by = 2L)] %in%
                              c("--input-dir", "--output-dir")))) {
  stop("Usage: Rscript HSC_Megakaryocyte_Platelet_altas_integration.R ",
       "[--input-dir PATH] [--output-dir PATH]")
}
options <- list("--input-dir" = file.path(script_dir, "input"),
                "--output-dir" = file.path(script_dir, "results"))
if (length(args)) {
  for (i in seq(1L, length(args), by = 2L)) options[[args[i]]] <- args[i + 1L]
}
input_dir <- options[["--input-dir"]]
output_dir <- options[["--output-dir"]]

required_packages <- c("Seurat", "SeuratObject", "harmony")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install the required R packages: ", paste(missing_packages, collapse = ", "))
}
library(Seurat)
source(file.path(script_dir, "in-slico_FACS.R"), local = TRUE)

disease_folders <- c("BloodCancer", "Malignant Lymphoma", "Multiple Myeloma",
                     "Osteosarcoma", "other Cancer")
rds_filename <- "first_qc_Megakaryocyte.rds"
rds_paths <- file.path(input_dir, disease_folders, rds_filename)
missing_inputs <- rds_paths[!file.exists(rds_paths)]
if (length(missing_inputs)) {
  stop("Prepared Seurat RDS inputs are missing:\n", paste(missing_inputs, collapse = "\n"))
}

seurat_list <- list()
for (i in seq_along(rds_paths)) {
  dataset_name <- disease_folders[i]
  message("Reading ", dataset_name, ": ", rds_paths[i])
  obj <- readRDS(rds_paths[i])
  if (!inherits(obj, "Seurat")) stop(dataset_name, ": input RDS must contain a Seurat object.")
  if (!"RNA" %in% SeuratObject::Assays(obj)) stop(dataset_name, ": RNA assay is missing.")
  SeuratObject::DefaultAssay(obj) <- "RNA"
  input_layers <- SeuratObject::Layers(obj[["RNA"]])
  if (sum(startsWith(input_layers, "counts")) > 1L) {
    obj[["RNA"]] <- SeuratObject::JoinLayers(obj[["RNA"]])
  }
  if (!"sample_id" %in% colnames(obj[[]])) {
    stop(dataset_name, ": sample_id metadata is required for per-sample QC and Harmony.")
  }
  obj$cancer_type <- dataset_name
  obj <- label_and_filter_sc_object(obj, filter_cells = TRUE, min_umi = 250L,
                                    min_target_cells = 10L, sample_col = "sample_id")
  if (is.null(obj)) {
    warning(dataset_name, ": no samples met the minimum of 10 target cells; skipped.")
    next
  }
  # Pairing source group with sample ID prevents collisions between studies.
  obj$integration_batch <- paste(dataset_name, obj$sample_id, sep = "::")
  seurat_list[[dataset_name]] <- obj
}
if (!length(seurat_list)) stop("No samples passed QC and marker gating.")

if (length(seurat_list) == 1L) {
  combined <- seurat_list[[1L]]
} else {
  cell_prefixes <- gsub("[^A-Za-z0-9]+", "_", names(seurat_list))
  combined <- merge(x = seurat_list[[1L]], y = seurat_list[-1L],
                    add.cell.ids = cell_prefixes)
}
# Seurat v5 merge can retain a separate counts layer per input. Harmony here
# uses one RNA assay, so rejoin layers before sample-level feature selection.
rna_layers <- SeuratObject::Layers(combined[["RNA"]])
if (sum(startsWith(rna_layers, "counts")) > 1L) {
  combined[["RNA"]] <- SeuratObject::JoinLayers(combined[["RNA"]])
}
if (ncol(combined) <= 30L) stop("At least 31 target cells are needed for 30 PCs.")
sample_counts <- table(combined$integration_batch)
if (length(sample_counts) < 2L) stop("Harmony needs at least two retained samples.")
message("Retained ", ncol(combined), " cells from ", length(sample_counts), " samples.")

# Select HVGs within each sample, then select the top 2,000 across samples.
cell_list <- SplitObject(combined, split.by = "integration_batch")
cell_list <- lapply(cell_list, function(obj) {
  obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                       scale.factor = 1e4, verbose = FALSE)
  FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000L,
                       verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = cell_list, nfeatures = 2000L)
features <- intersect(features, rownames(combined))
if (length(features) <= 30L) stop("At least 31 shared HVGs are needed for 30 PCs.")

combined <- NormalizeData(combined, normalization.method = "LogNormalize",
                          scale.factor = 1e4, verbose = FALSE)
VariableFeatures(combined) <- features
combined <- ScaleData(combined, features = features, verbose = FALSE)
combined <- RunPCA(combined, features = features, npcs = 30L, verbose = FALSE)
set.seed(1314L)
combined <- harmony::RunHarmony(combined, group.by.vars = "integration_batch",
                                reduction = "pca", dims.use = 1:30,
                                theta = 6, lambda = 1, verbose = FALSE)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30,
                          verbose = FALSE)
combined <- FindClusters(combined, resolution = 0.1, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30,
                    seed.use = 1314L, verbose = FALSE)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(combined, file.path(output_dir, "HSC_MK_Platelet_integrated.rds"))
sample_meta <- unique(combined[[]][, c("integration_batch", "sample_id", "cancer_type")])
sample_meta$retained_target_cells <- as.integer(
  sample_counts[as.character(sample_meta$integration_batch)])
utils::write.csv(sample_meta,
                 file.path(output_dir, "retained_target_cells_per_sample.csv"),
                 row.names = FALSE)

# Plot only metadata columns actually present in these prepared inputs.
plot_by <- intersect(c("sample_id", "sample_source", "marker_condition", "tissue",
                       "disease", "disease_state"), colnames(combined[[]]))
for (field in plot_by) {
  plot <- DimPlot(combined, reduction = "umap", group.by = field)
  plot_file <- paste0("UMAP_", field, ".pdf")
  ggplot2::ggsave(file.path(output_dir, plot_file), plot = plot,
                  width = 9, height = 7, units = "in")
}
message("Wrote integrated object and plots to: ", normalizePath(output_dir))
