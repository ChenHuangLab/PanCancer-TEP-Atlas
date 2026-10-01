#### verified HSC-megakaryocyte-platelet single-cell atlas constructions ####
library(data.table)
library(Seurat)
library(Matrix)
library(R.matlab)
rm(list = ls())
library(tidyverse)
message("--> 步骤 1: 正在配置路径...")

# Project root directory.
base_dir <- "~/huangjianxiang/Megakaryocyte/MK_dataset_supplementary/"

# Disease-group folders.
disease_folders <- c(
  "BloodCancer",
  "Malignant Lymphoma",
  "Multiple Myeloma",
  "Osteosarcoma",
  "other Cancer",
  "other disease"
)

# Input RDS filename.
rds_filename <- "first_qc_cleaned.rds"

# Build input file paths.
rds_file_paths <- file.path(base_dir, disease_folders,rds_filename)

# Inspect input paths.
print(rds_file_paths)


# ---------------------------------------------------------------- #
# Step 2: Read Seurat objects.
# ---------------------------------------------------------------- #
message("\n--> 步骤 2: 正在循环读取 ", length(rds_file_paths), " 个 Seurat RDS 文件...")

# Initialize the object list.
seurat_list <- list()

for (i in 1:length(rds_file_paths)) {
  
  file_path <- rds_file_paths[i]
  dataset_name <- disease_folders[i] # Use the folder name as the dataset name.
  
  message("  -> 正在读取: ", dataset_name)
  
  # Check that the input exists.
  if (file.exists(file_path)) {
    
    # Read the RDS object.
    seurat_obj <- readRDS(file_path)
    
    # Record the source dataset.
    seurat_obj$cancer_type <- dataset_name
    
    # Store the object under its dataset name.
    seurat_list[[dataset_name]] <- seurat_obj
    
  } else {
    # Warn and skip missing inputs.
    warning("文件未找到，已跳过: ", file_path)
  }
}


# ---------------------------------------------------------------- #
# Step 3: Merge Seurat objects.
# ---------------------------------------------------------------- #
# Require at least one valid object.
if (length(seurat_list) > 0) {
  
  message("\n--> 步骤 3: 正在合并 ", length(seurat_list), " 个 Seurat 对象...")
  
  # Use a single object directly.
  if (length(seurat_list) == 1) {
    combined_seurat <- seurat_list[[1]]
  } else {
    # Merge all objects.
    # Prefix cell names with dataset names to keep them unique.
    combined_seurat <- merge(x = seurat_list[[1]], 
                             y = seurat_list[-1], 
                             add.cell.ids = names(seurat_list))
  }
  
  message("--> 所有样本成功合并!")
  
} else {
  stop("错误：未能读取任何有效的 Seurat 对象，无法进行合并。")
}


# ---------------------------------------------------------------- #
# Step 4: Inspect the merged object.
# ---------------------------------------------------------------- #
message("\n--- 最终验证 ---")

# Print the merged object.
print(combined_seurat)
columns_to_check <- c(
  "sample_source", 
  "disease", 
  "disease_state", 
  "tissue",
  "sample_id"
)
if (!is.null(combined_seurat@meta.data)) {
  
  message("\n--- 使用 dplyr 生成的综合汇总表 ---")
  
  # Extract cell metadata.
  metadata <- combined_seurat@meta.data
  
  # Summarize cells and samples by group.
  summary_table <- metadata %>%
    group_by(disease, disease_state, tissue, sample_source) %>%
    summarise(
      Cell_Count = n(), # Cells per group.
      Unique_Sample_IDs = n_distinct(sample_id), # Unique samples per group.
      .groups = 'drop' # Drop grouping.
    ) %>%
    arrange(disease, disease_state) # Sort by disease and disease state.
  
  # Print the summary table.
  print(as.data.frame(summary_table))
  
} else {
  
  stop("错误：无法访问 Seurat 对象的元数据。")
  
}
sum(is.na(combined_seurat$sample_id))
sum(is.na(combined_seurat$sample_source))
sum(is.na(combined_seurat$disease))
sum(is.na(combined_seurat$disease_state))
sum(is.na(combined_seurat$tissue))
#### Harmonize gene identifiers ####
library(biomaRt)
library(dplyr)
library(stringr)
all_genes_in_seurat <- rownames(combined_seurat)
# Count Ensembl IDs starting with ENSG.
ensg_genes <- sum(str_starts(all_genes_in_seurat, "ENSG"))
message(paste("Ensembl ID (以ENSG开头) 的数量:", ensg_genes))
ac_genes <- sum(str_starts(all_genes_in_seurat,"AC"))
#write.csv(all_genes_in_seurat,file = "seurat_genes.csv")
write.table(all_genes_in_seurat, 
            file = "seurat_genes.csv", 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)
#### Quantify UMIs from unmapped genes ####
unmapped_file_path <- "unmapped_genes.csv"
if (!file.exists(unmapped_file_path)) {
  stop("错误: 未找到 'unmapped_gene_ids.csv' 文件。请先运行Python脚本导出该文件。")
}

unmapped_genes_df <- read.csv(unmapped_file_path, header = FALSE)
unmapped_genes_vector <- unmapped_genes_df$V1

message(paste("成功加载了", length(unmapped_genes_vector), "个未能映射的基因ID。"))


# Calculate the unmapped UMI fraction.
message("\n--- 正在计算表达量比例 ---")

# Extract raw counts.
# Inspect the available counts layers.
# List all RNA assay layers.
all_assay_layers <- Layers(combined_seurat, assay = "RNA")

# Select layers starting with counts.
count_layers <- grep("^counts", all_assay_layers, value = TRUE)

# Inspect the selected layers.
message(paste("在 'RNA' assay 中共发现", length(all_assay_layers), "个层。"))
message(paste("经过精确筛选后，实际用于计算的 'counts' 层有", length(count_layers), "个:"))
print(count_layers) # Confirm the layer names.
# Initialize totals.
# Initialize count accumulators.
total_counts_all_layers <- 0
unmapped_counts_all_layers <- 0

# Iterate over counts layers.
for (layer in count_layers) {
  message(paste("  -> 真正处理的层:", layer)) # Report the current layer.
  
  counts_matrix_layer <- GetAssayData(combined_seurat, assay = "RNA", layer = layer)
  print(dim(counts_matrix_layer))
  total_counts_all_layers <- total_counts_all_layers + sum(Matrix::colSums(counts_matrix_layer))
  
  genes_to_check_in_layer <- intersect(unmapped_genes_vector, rownames(counts_matrix_layer))
  
  if (length(genes_to_check_in_layer) > 0) {
    unmapped_counts_layer <- sum(Matrix::colSums(counts_matrix_layer[genes_to_check_in_layer, , drop = FALSE]))
    unmapped_counts_all_layers <- unmapped_counts_all_layers + unmapped_counts_layer
  }
}

# Calculate the fraction and percentage.
if (total_counts_all_layers > 0) {
  proportion <- unmapped_counts_all_layers / total_counts_all_layers
  percentage <- proportion * 100
} else {
  proportion <- 0
  percentage <- 0
}

# Report the results.
message("\n--- 计算完成 ---")
message(paste("所有 'counts' 层的总 UMI 数量:", format(total_counts_all_layers, big.mark = ",", scientific = FALSE)))
message(paste("所有 'counts' 层中未能映射基因的总 UMI 数量:", format(unmapped_counts_all_layers, big.mark = ",", scientific = FALSE)))
message(paste("未能映射基因的表达量占比: ", round(percentage, 4), "%"))

# Optional filtering threshold.
threshold <- 2.0 
#### Retain mapped genes and convert Ensembl IDs to symbols ####
mapped_genes_df <- read.csv("mapped_genes.csv", header = FALSE)
mapped_genes_vector <- mapped_genes_df$V1
# Keep mapped features.
# Subset by gene identifier.
seurat_integrated_filtered <- subset(combined_seurat, features = mapped_genes_vector)

# Equivalent indexing alternative.
# seurat_integrated_filtered <- seurat_integrated[mapped_genes_to_keep, ]


# Inspect the filtered object.
message(paste("过滤后，新对象包含", nrow(seurat_integrated_filtered), "个基因。"))

# Check object dimensions.
#print(dim(seurat_integrated))
print(dim(seurat_integrated_filtered))
ensg_genes <- sum(str_starts(rownames(seurat_integrated_filtered), "ENSG"))
#### Remap Ensembl IDs ####
names(mapped_genes_df) <- c("original_id", "entrezgene", "symbol")
# Separate Ensembl IDs and gene symbols.
message("正在从Seurat对象中分离ENSG基因和Symbol基因...")

# Extract feature names and counts using Seurat v5 methods.
all_genes <- rownames(seurat_integrated_filtered)
seurat_integrated_filtered <- JoinLayers(seurat_integrated_filtered)
counts_matrix <- LayerData(seurat_integrated_filtered, layer = "counts")
all_genes_from_matrix <- rownames(counts_matrix)
# Identify Ensembl IDs.
is_ensg <- str_starts(all_genes, "ENSG")

# Split features and counts into two groups.
ensg_genes <- all_genes_from_matrix[is_ensg]
symbol_genes <- all_genes_from_matrix[!is_ensg]


counts_ensg_part <- counts_matrix[ensg_genes, ]
counts_symbol_part <- counts_matrix[symbol_genes, ]

message(paste("分离完成:", length(ensg_genes), "个ENSG基因,", length(symbol_genes), "个Symbol基因。"))


# Map and clean Ensembl IDs.
message("正在对ENSG基因进行映射...")

# Select mappings for the current genes.
relevant_mapping <- mapped_genes_df[mapped_genes_df$original_id %in% ensg_genes, ]

# Keep the first mapping for each Ensembl ID.
# Assign one target symbol per Ensembl ID.
unique_mapping <- relevant_mapping[!duplicated(relevant_mapping$original_id), ]

# Create a named Ensembl-to-symbol lookup.
ensg_to_symbol_map <- setNames(unique_mapping$symbol, unique_mapping$original_id)

# Look up symbols for Ensembl rows.
new_symbol_names <- ensg_to_symbol_map[rownames(counts_ensg_part)]

# Identify successfully mapped genes.
mapped_indices <- !is.na(new_symbol_names)

# Retain mapped rows and assign symbol names.
counts_ensg_mapped <- counts_ensg_part[mapped_indices, ]
rownames(counts_ensg_mapped) <- new_symbol_names[mapped_indices]

message(paste("成功将", sum(mapped_indices), "个ENSG ID映射到Gene Symbol。"))


# Aggregate duplicate gene symbols.
message("正在检查并聚合重复的Gene Symbol...")
if (any(duplicated(rownames(counts_ensg_mapped)))) {
  message("检测到重复行名，开始使用 Matrix.utils 进行高效聚合...")
  library(Matrix.utils)
  # Sparse aggregation: start.
  # Sum duplicate rows without densifying the matrix.
  final_counts_ensg_part <- aggregate.Matrix(
    x = counts_ensg_mapped,
    groupings = rownames(counts_ensg_mapped),
    fun = 'sum'
  )
  # Sparse aggregation: end.
  
  message("ENSG部分聚合完成。")
} else {
  message("未检测到重复行名，无需聚合。")
  final_counts_ensg_part <- counts_ensg_mapped
}


# Rebuild the counts matrix and Seurat object.
message("正在合并矩阵...")

# Combine original symbols and mapped Ensembl rows.
# Bind the sparse matrices.
final_counts_matrix_pre <- rbind(counts_symbol_part, final_counts_ensg_part)

message("正在进行最终的全局聚合检查...")
# Check for duplicate symbols after merging.
if (any(duplicated(rownames(final_counts_matrix_pre)))) {
  message("检测到全局重复行名，开始最终聚合...")
  
  # Global aggregation: start.
  # Sum duplicate symbols across both groups.
  final_counts_matrix <- aggregate.Matrix(
    x = final_counts_matrix_pre,
    groupings = rownames(final_counts_matrix_pre),
    fun = 'sum'
  )
  # Global aggregation: end.
  
  message("全局聚合完成。")
} else {
  message("未产生新的重复，无需全局聚合。")
  final_counts_matrix <- final_counts_matrix_pre
}



# Create a Seurat object with harmonized counts and original metadata.
### Use the aggregated gene-symbol matrix.
seurat_final_harmonized <- CreateSeuratObject(
  counts = final_counts_matrix,
  meta.data = seurat_integrated_filtered@meta.data
)

message("🎉 任务完成！")
# Verify the rebuilt object.
message("最终对象的维度信息:")
print(dim(seurat_final_harmonized))
message("最终对象中是否还存在ENSG ID?")
print(paste("剩余ENSG数量:", sum(str_starts(rownames(seurat_final_harmonized), "ENSG"))))
### Inspect the Seurat object.
seurat_final_harmonized
# Inspect assay information.
# Confirm the RNA assay.
seurat_final_harmonized@assays
# Confirm the counts layer.
Layers(seurat_final_harmonized, assay = "RNA")

# Inspect cell metadata.
# Check that metadata rows match the cells.
head(seurat_final_harmonized@meta.data)
tail(seurat_final_harmonized@meta.data)
nrow(seurat_final_harmonized@meta.data) # Should equal the number of cells.
#### Inspect dotted gene names ####
all_genes <- rownames(seurat_final_harmonized)
genes_with_dots_mask <- grepl("\\.", all_genes)
sum(genes_with_dots_mask)
genes_to_remove <- all_genes[genes_with_dots_mask]
message("待移除基因名示例:")
print(head(genes_to_remove, 10))
# Retain genes without dots in their names.
# Invert the dotted-name mask.
genes_to_keep_mask <- !grepl("\\.", rownames(seurat_final_harmonized))

# Count retained genes.
message(paste("原始基因数:", nrow(seurat_final_harmonized)))
message(paste("将要保留的基因数:", sum(genes_to_keep_mask)))

# Subset the Seurat object.
# Filter the feature dimension.
# Retain all cells.
seurat_cleaned <- seurat_final_harmonized[genes_to_keep_mask, ]

message("已创建新的、不含带小数点基因名的Seurat对象: 'seurat_cleaned'")
seurat_cleaned@assays
# Confirm the counts layer.
Layers(seurat_cleaned, assay = "RNA")
combined_seurat <- seurat_cleaned
save(combined_seurat,file = "~/huangjianxiang/Megakaryocyte/MK_dataset_supplementary/first_qc_merged_cleaned_genes.rds")
#### Exclude samples with fewer than 10 cells ####
sample_counts <- table(combined_seurat$sample_id)
# Identify source datasets.
message("--- 步骤1: 正在从细胞名称中提取高级别 dataset_id ---")


# Count cells by dataset_id before filtering.
message("\n--- 步骤2: 正在按 dataset_id 统计细胞总数 (筛选前) ---")
dataset_counts <- table(combined_seurat$dataset_id)
print(dataset_counts)


# Select datasets to retain.
message("\n--- 步骤3: 正在找出需要保留的数据集 ---")
datasets_to_keep <- names(dataset_counts[dataset_counts >= 20])
message("将保留 ", length(datasets_to_keep), " 个细胞数不少于10的数据集。")


# Filter by dataset_id.
message("\n--- 步骤4: 正在根据 dataset_id 进行筛选 ---")
# Select dataset IDs in datasets_to_keep.
# Retain cells from those datasets.
sample_counts <- table(combined_seurat$sample_id)
message("筛选前的样本总数: ", length(sample_counts))

samples_to_keep <- names(sample_counts[sample_counts >= 10])
message("将保留 ", length(samples_to_keep), " 个细胞数不少于10的样本。")

# Filter samples.
# Subset the Seurat object.
# Select sample IDs in samples_to_keep.
# Retain cells from those samples.
filtered_seurat <- subset(combined_seurat,subset = dataset_id %in% datasets_to_keep)
filtered_seurat <- subset(combined_seurat, subset = sample_id %in% samples_to_keep)


# Verify filtering results.
message("\n--- 验证筛选结果 ---")
message("筛选前的总细胞数: ", ncol(combined_seurat))
message("筛选后的总细胞数: ", ncol(filtered_seurat))

message("\n筛选后剩余的样本及其细胞数:")
# Check that retained samples have at least 10 cells.
print(table(filtered_seurat$sample_id))
sce.all <- filtered_seurat
#sce.all <- JoinLayers(sce.all)
#### Find highly variable genes from counts ####
cell.list <- SplitObject(sce.all, split.by = "sample_id")
library(future)
library(future.apply)
# Configure parallel processing.
# Use multicore on Linux/macOS or multisession on Windows.
plan(multisession, workers = 12)  # Adjust workers to available CPU cores.
cell.list <- future_lapply(cell.list, function(seurat_obj) {
  seurat_obj <- FindVariableFeatures(seurat_obj,
                                     selection.method = "vst",
                                     nfeatures = 2000)
  return(seurat_obj)
})
plan()
features <-  SelectIntegrationFeatures(cell.list)
save(features,file = "first_qc_high_variable_2000.Rdata")
### Normalize, reduce dimensions and cluster.
library(future)

# Switch to sequential processing.
plan(sequential)
sce.all <- NormalizeData(sce.all)
sce.all <- ScaleData(sce.all,features = features)
sce.all <- RunPCA(sce.all, features = features)
sce.all <- RunUMAP(sce.all, reduction = "pca", dims = 1:30)
DimPlot(sce.all,reduction = "umap",group.by = "sample_id",raster=FALSE)
DimPlot(sce.all,reduction = "umap",group.by = "sample_source",raster=FALSE)
DimPlot(sce.all,reduction = "umap",group.by = "marker_condition",raster=FALSE)
#### Correct batch effects with Harmony ####
library(harmony)
sce.all <- RunHarmony(sce.all, group.by.vars = c("sample_id"), theta = 6,lambda = 1)
### Memory settings.

sce.all <- FindNeighbors(sce.all, reduction = 'harmony', dims = 1:30)
sce.all <- FindClusters(sce.all, resolution = 0.1)
sce.all <- RunUMAP(sce.all, reduction = "harmony", dims = 1:30,min.dist = 0.2)
#ElbowPlot(sce.all)
DimPlot(sce.all, reduction = "umap", group.by = "sample_id")
DimPlot(sce.all,reduction = "umap",group.by = "sample_source",raster=FALSE)
DimPlot(sce.all,reduction = "umap",group.by = "marker_condition",raster=FALSE)
DimPlot(sce.all,reduction = "umap",group.by = "tissue",raster=FALSE)
DimPlot(sce.all,reduction = "umap",split.by = "tissue",group.by = "marker_condition",raster=FALSE)

DimPlot(sce.all,reduction = "umap",group.by = "seurat_clusters",raster=FALSE,label = T)
SaveSeuratRds(sce.all,file = "/mnt/3/huang/huangjianxiang/Megakaryocyte/MK_dataset_supplementary/first_qc_final_harmony.rds")

#### Annotate immune-like megakaryocytes, HSCs and platelets ####
celltype=data.frame(ClusterID=c(0:24),
                    celltype= 'Unknown') # Initialize the annotation table.
celltype[celltype$ClusterID %in% c( 28,7,15,27,3,29,30,31),2]='Hematopoietic stem cell' 
celltype[celltype$ClusterID %in% c( 19,1,2,33,24, ),2]='T cell like Megakaryocyte' 
celltype[celltype$ClusterID %in% c( 17,18,19, ),2]='B cell like Megakaryocyte' 
celltype[celltype$ClusterID %in% c( 6,21 ),2]='Platelet' 
celltype[celltype$ClusterID %in% c( 0,5,11,12,14,18 ),2]='Megakaryocyte' 

head(celltype)
celltype
table(celltype$celltype)
sce.all@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce.all@meta.data[which(sce.all@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce.all@meta.data$celltype)
DimPlot(sce.all,reduction = "umap",group.by = "celltype",raster = FALSE,label = F)
SaveSeuratRds(sce.all,file = "/mnt/3/huang/huangjianxiang/Megakaryocyte/MK_dataset_supplementary/first_qc_final_harmony.rds")


