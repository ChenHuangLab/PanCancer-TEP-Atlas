#### verified HSC-megakaryocyte-platelet single-cell atlas constructions ####
library(data.table)
library(Seurat)
library(Matrix)
library(R.matlab)
rm(list = ls())
library(tidyverse)
message("--> 步骤 1: 正在配置路径...")

# 您的项目根目录
base_dir <- "~/huangjianxiang/Megakaryocyte/MK_dataset_supplementary/"

# 您提到的五个主要疾病类别文件夹
disease_folders <- c(
  "BloodCancer",
  "Malignant Lymphoma",
  "Multiple Myeloma",
  "Osteosarcoma",
  "other Cancer",
  "other disease"
)

# 目标 RDS 文件的固定名称
rds_filename <- "first_qc_cleaned.rds"

# 使用 file.path 智能地构建每个 rds 文件的完整路径
rds_file_paths <- file.path(base_dir, disease_folders,rds_filename)

# 打印出来检查一下路径是否正确
print(rds_file_paths)


# ---------------------------------------------------------------- #
# 步骤 2: 循环读取所有 Seurat 对象
# ---------------------------------------------------------------- #
message("\n--> 步骤 2: 正在循环读取 ", length(rds_file_paths), " 个 Seurat RDS 文件...")

# 初始化一个空列表，用于存储读入的 Seurat 对象
seurat_list <- list()

for (i in 1:length(rds_file_paths)) {
  
  file_path <- rds_file_paths[i]
  dataset_name <- disease_folders[i] # 使用文件夹名作为数据集的名字
  
  message("  -> 正在读取: ", dataset_name)
  
  # 检查文件是否存在，防止因某个文件缺失而报错
  if (file.exists(file_path)) {
    
    # 读取 RDS 文件
    seurat_obj <- readRDS(file_path)
    
    # 【最佳实践】为每个对象添加一个新的元数据列，标记其来源数据集
    seurat_obj$cancer_type <- dataset_name
    
    # 将读入并标记好的对象存入列表，并用数据集名称为其命名
    seurat_list[[dataset_name]] <- seurat_obj
    
  } else {
    # 如果文件不存在，则发出警告并跳过
    warning("文件未找到，已跳过: ", file_path)
  }
}


# ---------------------------------------------------------------- #
# 步骤 3: 合并列表中的所有 Seurat 对象
# ---------------------------------------------------------------- #
# 确保列表中至少有一个对象可以合并
if (length(seurat_list) > 0) {
  
  message("\n--> 步骤 3: 正在合并 ", length(seurat_list), " 个 Seurat 对象...")
  
  # 如果只有一个对象，直接使用它
  if (length(seurat_list) == 1) {
    combined_seurat <- seurat_list[[1]]
  } else {
    # 使用 merge 函数合并列表中的所有对象
    # add.cell.ids 会使用列表的名字 (我们设置的数据集名) 作为细胞名的前缀，确保唯一性
    combined_seurat <- merge(x = seurat_list[[1]], 
                             y = seurat_list[-1], 
                             add.cell.ids = names(seurat_list))
  }
  
  message("--> 所有样本成功合并!")
  
} else {
  stop("错误：未能读取任何有效的 Seurat 对象，无法进行合并。")
}


# ---------------------------------------------------------------- #
# 步骤 4: 验证合并结果
# ---------------------------------------------------------------- #
message("\n--- 最终验证 ---")

# 打印合并后的 Seurat 对象信息
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
  
  # 从Seurat对象中提取元数据
  metadata <- combined_seurat@meta.data
  
  # 使用dplyr进行分组和计数
  summary_table <- metadata %>%
    group_by(disease, disease_state, tissue, sample_source) %>%
    summarise(
      Cell_Count = n(), # 计算每个分组的细胞数
      Unique_Sample_IDs = n_distinct(sample_id), # 计算每个分组有多少个独特的样本ID
      .groups = 'drop' # 取消分组
    ) %>%
    arrange(disease, disease_state) # 按疾病排序
  
  # 打印漂亮的汇总表
  print(as.data.frame(summary_table))
  
} else {
  
  stop("错误：无法访问 Seurat 对象的元数据。")
  
}
sum(is.na(combined_seurat$sample_id))
sum(is.na(combined_seurat$sample_source))
sum(is.na(combined_seurat$disease))
sum(is.na(combined_seurat$disease_state))
sum(is.na(combined_seurat$tissue))
#### 对怪异的基因名进行一个校正的过程 ####
library(biomaRt)
library(dplyr)
library(stringr)
all_genes_in_seurat <- rownames(combined_seurat)
# 1. 检查以 "ENSG" 开头的基因 (Ensembl IDs)
ensg_genes <- sum(str_starts(all_genes_in_seurat, "ENSG"))
message(paste("Ensembl ID (以ENSG开头) 的数量:", ensg_genes))
ac_genes <- sum(str_starts(all_genes_in_seurat,"AC"))
#write.csv(all_genes_in_seurat,file = "seurat_genes.csv")
write.table(all_genes_in_seurat, 
            file = "seurat_genes.csv", 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)
#### 未能成功映射的基因的一个表达量的占比
unmapped_file_path <- "unmapped_genes.csv"
if (!file.exists(unmapped_file_path)) {
  stop("错误: 未找到 'unmapped_gene_ids.csv' 文件。请先运行Python脚本导出该文件。")
}

unmapped_genes_df <- read.csv(unmapped_file_path, header = FALSE)
unmapped_genes_vector <- unmapped_genes_df$V1

message(paste("成功加载了", length(unmapped_genes_vector), "个未能映射的基因ID。"))


# --- 2. 计算表达量比例 ---
message("\n--- 正在计算表达量比例 ---")

# 从Seurat对象中提取原始counts矩阵 (这是最准确的)
# 您的对象有多个layers，我们使用默认的 'counts' layer
# a. 首先获取 'RNA' assay 中的【所有】层
all_assay_layers <- Layers(combined_seurat, assay = "RNA")

# b. 然后，使用 grep 手动、精确地筛选出以 "counts" 开头的层
count_layers <- grep("^counts", all_assay_layers, value = TRUE)

# c. 验证我们筛选出的层是否正确
message(paste("在 'RNA' assay 中共发现", length(all_assay_layers), "个层。"))
message(paste("经过精确筛选后，实际用于计算的 'counts' 层有", length(count_layers), "个:"))
print(count_layers) # 打印出来确认一下
# b. 初始化累加器
# d. 初始化累加器
total_counts_all_layers <- 0
unmapped_counts_all_layers <- 0

# e. 循环遍历【已正确筛选的】counts 层列表
for (layer in count_layers) {
  message(paste("  -> 真正处理的层:", layer)) # 修改了日志信息以示区别
  
  counts_matrix_layer <- GetAssayData(combined_seurat, assay = "RNA", layer = layer)
  print(dim(counts_matrix_layer))
  total_counts_all_layers <- total_counts_all_layers + sum(Matrix::colSums(counts_matrix_layer))
  
  genes_to_check_in_layer <- intersect(unmapped_genes_vector, rownames(counts_matrix_layer))
  
  if (length(genes_to_check_in_layer) > 0) {
    unmapped_counts_layer <- sum(Matrix::colSums(counts_matrix_layer[genes_to_check_in_layer, , drop = FALSE]))
    unmapped_counts_all_layers <- unmapped_counts_all_layers + unmapped_counts_layer
  }
}

# f. 计算最终的比例和百分比 (后续代码与之前相同)
if (total_counts_all_layers > 0) {
  proportion <- unmapped_counts_all_layers / total_counts_all_layers
  percentage <- proportion * 100
} else {
  proportion <- 0
  percentage <- 0
}

# 打印最终结果
message("\n--- 计算完成 ---")
message(paste("所有 'counts' 层的总 UMI 数量:", format(total_counts_all_layers, big.mark = ",", scientific = FALSE)))
message(paste("所有 'counts' 层中未能映射基因的总 UMI 数量:", format(unmapped_counts_all_layers, big.mark = ",", scientific = FALSE)))
message(paste("未能映射基因的表达量占比: ", round(percentage, 4), "%"))

# --- 3. (可选步骤) 过滤操作 (与之前相同，subset函数可直接用于v5对象) ---
threshold <- 2.0 
#### 我的策略是对这些不能映射的基因进行去除，然后可以映射的映射为gene_symbol，然后在每一个数据集中执行这样的操作 ####
mapped_genes_df <- read.csv("mapped_genes.csv", header = FALSE)
mapped_genes_vector <- mapped_genes_df$V1
# Seurat 提供了非常方便的 subset 函数来筛选 features (基因)
# 这是最推荐、最清晰的方式
seurat_integrated_filtered <- subset(combined_seurat, features = mapped_genes_vector)

# 您也可以使用R基础的方括号语法，效果是相同的
# seurat_integrated_filtered <- seurat_integrated[mapped_genes_to_keep, ]


# --- 步骤 3: 验证结果 ---
message(paste("过滤后，新对象包含", nrow(seurat_integrated_filtered), "个基因。"))

# 检查一下维度变化
#print(dim(seurat_integrated))
print(dim(seurat_integrated_filtered))
ensg_genes <- sum(str_starts(rownames(seurat_integrated_filtered), "ENSG"))
#### 将这部分的ensemble_id重新映射 ####
names(mapped_genes_df) <- c("original_id", "entrezgene", "symbol")
# --- 步骤 2: 分离Seurat对象中的数据 ---
message("正在从Seurat对象中分离ENSG基因和Symbol基因...")

# 获取当前的基因名和counts矩阵 (使用兼容v5的函数)
all_genes <- rownames(seurat_integrated_filtered)
seurat_integrated_filtered <- JoinLayers(seurat_integrated_filtered)
counts_matrix <- LayerData(seurat_integrated_filtered, layer = "counts")
all_genes_from_matrix <- rownames(counts_matrix)
# 识别哪些是ENSG ID
is_ensg <- str_starts(all_genes, "ENSG")

# 将基因名和矩阵分成两部分
ensg_genes <- all_genes_from_matrix[is_ensg]
symbol_genes <- all_genes_from_matrix[!is_ensg]


counts_ensg_part <- counts_matrix[ensg_genes, ]
counts_symbol_part <- counts_matrix[symbol_genes, ]

message(paste("分离完成:", length(ensg_genes), "个ENSG基因,", length(symbol_genes), "个Symbol基因。"))


# --- 步骤 3: 执行映射与清理 ---
message("正在对ENSG基因进行映射...")

# 筛选映射文件中与我们数据相关的行
relevant_mapping <- mapped_genes_df[mapped_genes_df$original_id %in% ensg_genes, ]

# 处理一对多问题：如果一个ENSG ID对应多个Symbol，我们只保留第一个
# 这确保了每个ENSG ID只有一个唯一的映射目标
unique_mapping <- relevant_mapping[!duplicated(relevant_mapping$original_id), ]

# 创建一个从ENSG到Symbol的“字典”（命名向量）
ensg_to_symbol_map <- setNames(unique_mapping$symbol, unique_mapping$original_id)

# 将ENSG行名替换为Symbol行名
new_symbol_names <- ensg_to_symbol_map[rownames(counts_ensg_part)]

# 找出哪些ENSG成功找到了对应的Symbol
mapped_indices <- !is.na(new_symbol_names)

# 更新counts矩阵，只保留成功映射的基因，并赋予新的Symbol行名
counts_ensg_mapped <- counts_ensg_part[mapped_indices, ]
rownames(counts_ensg_mapped) <- new_symbol_names[mapped_indices]

message(paste("成功将", sum(mapped_indices), "个ENSG ID映射到Gene Symbol。"))


# --- 步骤 4: 聚合重复的Gene Symbol (核心步骤) ---
message("正在检查并聚合重复的Gene Symbol...")
if (any(duplicated(rownames(counts_ensg_mapped)))) {
  message("检测到重复行名，开始使用 Matrix.utils 进行高效聚合...")
  library(Matrix.utils)
  # --- 适配部分 1: START ---
  # 使用 aggregate.Matrix 直接在稀疏矩阵上操作，高效且内存友好
  final_counts_ensg_part <- aggregate.Matrix(
    x = counts_ensg_mapped,
    groupings = rownames(counts_ensg_mapped),
    fun = 'sum'
  )
  # --- 适配部分 1: END ---
  
  message("ENSG部分聚合完成。")
} else {
  message("未检测到重复行名，无需聚合。")
  final_counts_ensg_part <- counts_ensg_mapped
}


# --- 步骤 5 (已优化): 重建最终的counts矩阵和Seurat对象 ---
message("正在合并矩阵...")

# 将原始的Symbol部分和处理好的ENSG部分合并
# 两个部分都应该是稀疏矩阵，rbind会高效地合并它们
final_counts_matrix_pre <- rbind(counts_symbol_part, final_counts_ensg_part)

message("正在进行最终的全局聚合检查...")
# 检查合并后是否产生了新的重复行名
if (any(duplicated(rownames(final_counts_matrix_pre)))) {
  message("检测到全局重复行名，开始最终聚合...")
  
  # --- 适配部分 2: START ---
  # 再次使用 aggregate.Matrix 进行全局聚合
  final_counts_matrix <- aggregate.Matrix(
    x = final_counts_matrix_pre,
    groupings = rownames(final_counts_matrix_pre),
    fun = 'sum'
  )
  # --- 适配部分 2: END ---
  
  message("全局聚合完成。")
} else {
  message("未产生新的重复，无需全局聚合。")
  final_counts_matrix <- final_counts_matrix_pre
}



# 使用完全统一的counts矩阵，以及原始的meta.data，创建一个新的、干净的Seurat对象
### 使用聚合后的一个矩阵，现在就是所有的基因名都是一个可以映射的状态
seurat_final_harmonized <- CreateSeuratObject(
  counts = final_counts_matrix,
  meta.data = seurat_integrated_filtered@meta.data
)

message("🎉 任务完成！")
# --- 验证 ---
message("最终对象的维度信息:")
print(dim(seurat_final_harmonized))
message("最终对象中是否还存在ENSG ID?")
print(paste("剩余ENSG数量:", sum(str_starts(rownames(seurat_final_harmonized), "ENSG"))))
### 检查seurat对象是否正常
seurat_final_harmonized
# 3. 检查 Assay 信息
# 确认 Assay 名称是 "RNA"
seurat_final_harmonized@assays
# 确认 counts layer 存在
Layers(seurat_final_harmonized, assay = "RNA")

# 4. 快速浏览一下元数据 (meta.data)
# 确保它看起来是完整的，并且行数与细胞数一致
head(seurat_final_harmonized@meta.data)
tail(seurat_final_harmonized@meta.data)
nrow(seurat_final_harmonized@meta.data) # 这个数字应该等于 dim() 输出的第二个数字
#### 还有一些重复基因名的基因
all_genes <- rownames(seurat_final_harmonized)
genes_with_dots_mask <- grepl("\\.", all_genes)
sum(genes_with_dots_mask)
genes_to_remove <- all_genes[genes_with_dots_mask]
message("待移除基因名示例:")
print(head(genes_to_remove, 10))
# 我们想要保留的是不包含小数点的基因
# 所以我们使用 '!' 来反转上面的逻辑向量
genes_to_keep_mask <- !grepl("\\.", rownames(seurat_final_harmonized))

# 计算一下将要保留的基因数量
message(paste("原始基因数:", nrow(seurat_final_harmonized)))
message(paste("将要保留的基因数:", sum(genes_to_keep_mask)))

# 对Seurat对象进行子集操作
# 第一个维度是基因（features），我们只保留'genes_to_keep_mask'为TRUE的
# 第二个维度是细胞（cells），我们保留所有细胞，所以留空
seurat_cleaned <- seurat_final_harmonized[genes_to_keep_mask, ]

message("已创建新的、不含带小数点基因名的Seurat对象: 'seurat_cleaned'")
seurat_cleaned@assays
# 确认 counts layer 存在
Layers(seurat_cleaned, assay = "RNA")
combined_seurat <- seurat_cleaned
save(combined_seurat,file = "~/huangjianxiang/Megakaryocyte/MK_dataset_supplementary/first_qc_merged_cleaned_genes.rds")
#### 开始进一步质控，去除细胞数小于10的样本 ####
sample_counts <- table(combined_seurat$sample_id)
# 这是最稳健的方法，我们再次从细胞名称中提取数据集ID
message("--- 步骤1: 正在从细胞名称中提取高级别 dataset_id ---")


# --- 步骤 2: 按新的 'dataset_id' 统计每个数据集的总细胞数 ---
message("\n--- 步骤2: 正在按 dataset_id 统计细胞总数 (筛选前) ---")
dataset_counts <- table(combined_seurat$dataset_id)
print(dataset_counts)


# --- 步骤 3: 找出细胞数大于或等于10的数据集ID ---
message("\n--- 步骤3: 正在找出需要保留的数据集 ---")
datasets_to_keep <- names(dataset_counts[dataset_counts >= 20])
message("将保留 ", length(datasets_to_keep), " 个细胞数不少于10的数据集。")


# --- 步骤 4: 【核心操作】使用 subset() 函数进行最终筛选 ---
message("\n--- 步骤4: 正在根据 dataset_id 进行筛选 ---")
# subset = dataset_id %in% datasets_to_keep 的意思是：
# 只保留元数据中 'dataset_id' 的值存在于 'datasets_to_keep' 向量中的那些细胞
sample_counts <- table(combined_seurat$sample_id)
message("筛选前的样本总数: ", length(sample_counts))

samples_to_keep <- names(sample_counts[sample_counts >= 10])
message("将保留 ", length(samples_to_keep), " 个细胞数不少于10的样本。")

# --- 步骤3: 【核心操作】使用 subset() 函数进行筛选 ---
# subset() 是Seurat中非常方便的筛选函数
# subset = sample_id %in% samples_to_keep 的意思是：
# 只保留元数据中 'sample_id' 这一列的值存在于 'samples_to_keep' 向量中的那些细胞
filtered_seurat <- subset(combined_seurat,subset = dataset_id %in% datasets_to_keep)
filtered_seurat <- subset(combined_seurat, subset = sample_id %in% samples_to_keep)


# --- 步骤4: 验证结果 ---
message("\n--- 验证筛选结果 ---")
message("筛选前的总细胞数: ", ncol(combined_seurat))
message("筛选后的总细胞数: ", ncol(filtered_seurat))

message("\n筛选后剩余的样本及其细胞数:")
# 再次运行table()，您会看到所有样本的细胞数都 >= 10
print(table(filtered_seurat$sample_id))
sce.all <- filtered_seurat
#sce.all <- JoinLayers(sce.all)
#### 直接在counts矩阵上找高变基因，加快速度 ####
cell.list <- SplitObject(sce.all, split.by = "sample_id")
library(future)
library(future.apply)
# 设置并行策略
# multicore 在 Linux/macOS 下可用，windows 下用 multisession
plan(multisession, workers = 12)  # 根据 CPU 核心数调整 workers
cell.list <- future_lapply(cell.list, function(seurat_obj) {
  seurat_obj <- FindVariableFeatures(seurat_obj,
                                     selection.method = "vst",
                                     nfeatures = 2000)
  return(seurat_obj)
})
plan()
features <-  SelectIntegrationFeatures(cell.list)
save(features,file = "first_qc_high_variable_2000.Rdata")
### 找完高变之后进行一个标准化，降维聚类分析 ####
library(future)

# 1️⃣ 设置单线程（顺序执行）
plan(sequential)
sce.all <- NormalizeData(sce.all)
sce.all <- ScaleData(sce.all,features = features)
sce.all <- RunPCA(sce.all, features = features)
sce.all <- RunUMAP(sce.all, reduction = "pca", dims = 1:30)
DimPlot(sce.all,reduction = "umap",group.by = "sample_id",raster=FALSE)
DimPlot(sce.all,reduction = "umap",group.by = "sample_source",raster=FALSE)
DimPlot(sce.all,reduction = "umap",group.by = "marker_condition",raster=FALSE)
#### harmony去除批次效应 ####
library(harmony)
sce.all <- RunHarmony(sce.all, group.by.vars = c("sample_id"), theta = 6,lambda = 1)
### 加大内存设置

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
# 3. 现在，在这个经过处理的、毫无歧义的 sce.all 对象上运行 IntegrateLayers
#    这个对象结构完整，并且 scVI 只能选择 counts 数据
# sce.all.integrated.scvi <- IntegrateLayers(
#   object = sce.all,           # <--- 使用我们修改后的原始对象
#   features = features,
#   method = scVIIntegration,
#   batch = "sample_id",
#   new.reduction = "integrated.scvi",
#   
#   # 其他训练参数保持不变
#   max_epochs = 400,
#   gene_likelihood = "nb",
#   early_stopping = TRUE,
#   num_workers=111,
#   conda_env = "/home/hc/anaconda3/envs/r_443/bin/python",
#   verbose = TRUE
# )
#### 细分免疫like-MK，HSC，以及血小板 ####
celltype=data.frame(ClusterID=c(0:24),
                    celltype= 'Unknown') #随便构建一个空的数据框
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


