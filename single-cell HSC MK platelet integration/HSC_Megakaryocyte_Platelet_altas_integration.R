#!/usr/bin/env Rscript
# Publication version of 2026_07_05_single_cell_atlas_construction.R.
# The exploratory source file is preserved outside this repository.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
args <- commandArgs(trailingOnly = TRUE)
allowed <- c("--input-dir", "--output-dir", "--gene-map", "--cluster-map")
if (length(args) %% 2L || any(!args[seq(1L, length(args), by = 2L)] %in% allowed)) {
  stop("Usage: Rscript HSC_Megakaryocyte_Platelet_altas_integration.R ",
       "[--input-dir PATH] [--output-dir PATH] [--gene-map PATH] ",
       "[--cluster-map PATH]")
}
opts <- list("--input-dir" = file.path(script_dir, "input"),
             "--output-dir" = file.path(script_dir, "results"))
if (length(args)) {
  for (i in seq(1L, length(args), by = 2L)) opts[[args[i]]] <- args[i + 1L]
}
input_dir <- opts[["--input-dir"]]
output_dir <- opts[["--output-dir"]]
gene_map_path <- if ("--gene-map" %in% args) opts[["--gene-map"]] else
  file.path(input_dir, "mapped_genes.csv")
cluster_map_path <- if ("--cluster-map" %in% args) opts[["--cluster-map"]] else
  file.path(input_dir, "cluster_celltype_map.csv")

packages <- c("Seurat", "SeuratObject", "Matrix", "harmony")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install the required R packages: ", paste(missing_packages, collapse = ", "))
}
library(Seurat)

groups <- c("BloodCancer", "Malignant Lymphoma", "Multiple Myeloma",
            "Osteosarcoma", "other Cancer", "other disease")
inputs <- file.path(input_dir, groups, "first_qc_cleaned.rds")
missing_inputs <- inputs[!file.exists(inputs)]
if (length(missing_inputs)) {
  stop("Missing prepared Seurat inputs:\n", paste(missing_inputs, collapse = "\n"))
}
objects <- setNames(vector("list", length(groups)), groups)
for (i in seq_along(groups)) {
  message("Reading ", groups[i], ": ", inputs[i])
  obj <- readRDS(inputs[i])
  if (!inherits(obj, "Seurat")) stop(inputs[i], " does not contain a Seurat object.")
  if (!"RNA" %in% SeuratObject::Assays(obj)) stop(inputs[i], " has no RNA assay.")
  if (!"sample_id" %in% colnames(obj[[]])) stop(inputs[i], " has no sample_id column.")
  SeuratObject::DefaultAssay(obj) <- "RNA"
  obj$cancer_type <- groups[i]
  objects[[i]] <- obj
}
prefixes <- gsub("[^A-Za-z0-9]+", "_", names(objects))
combined <- merge(objects[[1L]], y = objects[-1L], add.cell.ids = prefixes)
if (sum(startsWith(SeuratObject::Layers(combined[["RNA"]]), "counts")) > 1L) {
  combined[["RNA"]] <- SeuratObject::JoinLayers(combined[["RNA"]])
}
counts <- SeuratObject::LayerData(combined, assay = "RNA", layer = "counts")
if (!inherits(counts, "sparseMatrix")) counts <- methods::as(counts, "dgCMatrix")
# Apply the manuscript's UMI threshold to raw counts before gene mapping.
cells_250 <- colnames(counts)[Matrix::colSums(counts) >= 250]
if (!length(cells_250)) stop("No cells satisfy the 250 raw UMI threshold.")
combined <- subset(combined, cells = cells_250)
counts <- counts[, cells_250, drop = FALSE]

# The mapping table is produced by the upstream MyGene.info annotation step.
# Like the source script, it determines which original gene IDs are retained.
if (!file.exists(gene_map_path)) stop("Missing upstream gene map: ", gene_map_path)
gene_map <- utils::read.csv(gene_map_path, header = FALSE, stringsAsFactors = FALSE)
if (ncol(gene_map) < 3L) stop("Gene map needs original_id, entrezgene, symbol columns.")
gene_map <- gene_map[, 1:3]
names(gene_map) <- c("original_id", "entrezgene", "symbol")
if (nrow(gene_map) && gene_map$original_id[1] == "original_id") gene_map <- gene_map[-1L, ]
gene_map$original_id <- trimws(gene_map$original_id)
gene_map$symbol <- trimws(gene_map$symbol)
gene_map <- gene_map[!is.na(gene_map$original_id) & nzchar(gene_map$original_id) &
                       !is.na(gene_map$symbol) & nzchar(gene_map$symbol), ]
gene_map <- gene_map[!duplicated(gene_map$original_id), ]
gene_ids <- rownames(counts)
lookup <- match(gene_ids, gene_map$original_id)
keep <- !is.na(lookup)
if (!any(keep)) stop("No count-matrix genes match mapped_genes.csv.")
symbols <- gene_ids[keep]
is_ensg <- startsWith(symbols, "ENSG")
symbols[is_ensg] <- gene_map$symbol[lookup[keep][is_ensg]]
counts <- counts[keep, , drop = FALSE]

# Sparse summation handles duplicate mappings and existing symbol collisions.
symbol_levels <- unique(symbols)
if (anyDuplicated(symbols)) {
  incidence <- Matrix::sparseMatrix(
    i = seq_along(symbols), j = match(symbols, symbol_levels), x = 1,
    dims = c(length(symbols), length(symbol_levels)))
  counts <- Matrix::t(incidence) %*% counts
}
rownames(counts) <- symbol_levels
if (anyDuplicated(rownames(counts))) stop("Duplicate gene symbols remain.")
message("Retained ", nrow(counts), " mapped genes and ", ncol(counts), " cells.")

unmapped_path <- file.path(input_dir, "unmapped_genes.csv")
if (file.exists(unmapped_path)) {
  unmapped <- utils::read.csv(unmapped_path, header = FALSE, stringsAsFactors = FALSE)[[1L]]
  original_counts <- SeuratObject::LayerData(combined, assay = "RNA", layer = "counts")
  unmatched <- intersect(unmapped, rownames(original_counts))
  denominator <- sum(original_counts)
  fraction <- if (denominator > 0) sum(original_counts[unmatched, , drop = FALSE]) /
    denominator else NA_real_
  message("Unmapped gene UMI fraction before filtering: ", signif(fraction, 4L))
}

metadata <- combined[[]]
metadata$nCount_RNA <- NULL
metadata$nFeature_RNA <- NULL
combined <- CreateSeuratObject(counts = counts, meta.data = metadata)
meta <- combined[[]]
if (anyNA(meta$sample_id) || any(!nzchar(as.character(meta$sample_id)))) {
  stop("sample_id must be defined for every cell.")
}
combined$integration_batch <- paste(meta$cancer_type, meta$sample_id, sep = "::")
sample_sizes <- table(combined$integration_batch)
retained_samples <- names(sample_sizes[sample_sizes >= 10L])
if (!length(retained_samples)) stop("No sample has at least 10 cells after UMI filtering.")
combined <- subset(combined, cells = colnames(combined)[
  combined$integration_batch %in% retained_samples])
if (ncol(combined) <= 30L) stop("At least 31 cells are needed for 30 PCs.")
if (length(retained_samples) < 2L) stop("Harmony needs at least two retained samples.")
message("Retained ", ncol(combined), " cells from ", length(retained_samples), " samples.")

sample_objects <- SplitObject(combined, split.by = "integration_batch")
sample_objects <- lapply(sample_objects, function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize",
                     scale.factor = 1e4, verbose = FALSE)
  FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000L,
                       verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = sample_objects, nfeatures = 2000L)
features <- intersect(features, rownames(combined))
if (length(features) <= 30L) stop("At least 31 HVGs are needed for 30 PCs.")
combined <- NormalizeData(combined, normalization.method = "LogNormalize",
                          scale.factor = 1e4, verbose = FALSE)
VariableFeatures(combined) <- features
combined <- ScaleData(combined, features = features, verbose = FALSE)
combined <- RunPCA(combined, features = features, npcs = 30L, verbose = FALSE)
set.seed(1314L)
combined <- harmony::RunHarmony(combined, group.by.vars = "integration_batch",
                                theta = 6, lambda = 1, verbose = FALSE)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30,
                          verbose = FALSE)
combined <- FindClusters(combined, resolution = 0.1, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30,
                    min.dist = 0.2, seed.use = 1314L, verbose = FALSE)

# The exploratory mapping assigns conflicting types to clusters 18 and 19,
# and refers to clusters beyond its initialized range. Supply a reviewed CSV.
if (file.exists(cluster_map_path)) {
  annotation <- utils::read.csv(cluster_map_path, stringsAsFactors = FALSE)
  if (!all(c("ClusterID", "celltype") %in% names(annotation))) {
    stop("Cluster mapping must contain ClusterID and celltype columns.")
  }
  ids <- as.character(annotation$ClusterID)
  labels <- trimws(as.character(annotation$celltype))
  if (anyNA(ids) || anyNA(labels) || any(!nzchar(ids)) || any(!nzchar(labels)) ||
      anyDuplicated(ids)) {
    stop("Cluster mapping must have unique, nonempty ClusterID and celltype values.")
  }
  observed <- as.character(sort(unique(as.integer(as.character(combined$seurat_clusters)))))
  if (!setequal(ids, observed)) {
    stop("Cluster mapping must cover exactly the observed clusters: ",
         paste(observed, collapse = ", "))
  }
  combined$celltype <- setNames(labels, ids)[as.character(combined$seurat_clusters)]
  message("Applied cell-type labels from ", cluster_map_path)
} else {
  message("No cluster mapping supplied; output retains unannotated seurat_clusters.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(combined, file.path(output_dir, "first_qc_final_harmony.rds"))
utils::write.csv(data.frame(gene = features),
                 file.path(output_dir, "integration_features.csv"), row.names = FALSE)
meta <- combined[[]]
samples <- unique(meta[, c("integration_batch", "sample_id", "cancer_type")])
samples$retained_cells <- as.integer(table(meta$integration_batch)[samples$integration_batch])
utils::write.csv(samples, file.path(output_dir, "retained_cells_per_sample.csv"),
                 row.names = FALSE)
for (field in intersect(c("sample_id", "sample_source", "marker_condition", "tissue",
                          "disease", "disease_state", "seurat_clusters", "celltype"),
                        names(meta))) {
  p <- DimPlot(combined, reduction = "umap", group.by = field)
  ggplot2::ggsave(file.path(output_dir, paste0("UMAP_", field, ".pdf")), plot = p,
                  width = 9, height = 7, units = "in")
}
message("Wrote integrated object and summaries to: ", normalizePath(output_dir))
