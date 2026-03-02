######################################################################
library(scmap)
library(Seurat)
library(SingleR)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(matrixStats)
library(mclust)
library(cowplot)
library(ggplot2)
library(dplyr)
library(scrapper)
library(batchelor)

#--------------------------bulk data
load("./CUPM/data/RNAseq/GSE183202/GSE183202_Exp_Cli.RData")
exp1=log2(fpkm1+1)
#------------------filter out 低表达的gene
index1=rowMeans(exp1<1)
exp1=exp1[which(index1<0.75),]
cli1[cli1==""]=NA
ref_exp=exp1[,!is.na(cli1$groups)]
ref_label=cli1[!is.na(cli1$groups),"groups"]

######################################################################
#------------------scRNA-seq data
load("./CUPM/data/scRNAseq/scaled_GSE183916.RData")
sct_exp <- GetAssayData(data1, assay = "SCT", layer = "data")
index2=rowMeans(sct_exp==0)
sct_exp=sct_exp[which(index2<0.75),]
common_genes <- intersect(rownames(sct_exp), rownames(ref_exp))
length(common_genes)
test_mat <- sct_exp[common_genes, , drop=FALSE]         # genes x cells
ref_mat  <- ref_exp[common_genes, , drop=FALSE]    # genes x samples
bulk_labels <- ref_label
seurat_obj <- data1

# 1. 构建 reference
ref_sce <- SingleCellExperiment(list(logcounts = as.matrix(ref_mat)))
ref_sce$group <- bulk_labels
rowData(ref_sce)$feature_symbol <- rownames(ref_mat)
# 2. 构建 query
test_sce <- SingleCellExperiment(list(logcounts = as.matrix(test_mat)))
rowData(test_sce)$feature_symbol <- rownames(test_mat)
#特征选择与索引构建
ref_sce <- selectFeatures(ref_sce, n_features = 1000, suppress_plot = TRUE)
ref_sce <- indexCluster(ref_sce, cluster_col = "group")
test_sce <- selectFeatures(test_sce, n_features = 1000, suppress_plot = TRUE)
# 3. 映射
scmap_res <- scmapCluster(
  projection = test_sce,
  index_list = list(ref = metadata(ref_sce)$scmap_cluster_index),
  threshold = 0.05
)
str(scmap_res)
table(scmap_res$combined_labs)
summary(scmap_res$scmap_cluster_siml)
# 4. 结果
scmap_label <- scmap_res$combined_labs
###############################################################
#--------B_cell
# 指定目标细胞类型
cell_type_name <- "B_cell"
# 找到满足“该类型 + unassigned分组”的细胞索引
idx_unassigned <- which(seurat_obj$celltype == cell_type_name & scmap_label == "unassigned")
# 检查数量
n <- length(idx_unassigned)
if (n == 0) {
  stop("没有需要处理的 unassigned 细胞")
}

# 计算分配比例（7:2:1）
n1 <- floor(n * 0.7)  # k1
n2 <- floor(n * 0.2)  # k2
n3 <- n - n1 - n2     # k3

# 随机划分
set.seed(123)  # 可选，保证可重复
idx_shuffled <- sample(idx_unassigned)

# 将结果写回 group_vec
scmap_label[idx_shuffled[1:n1]]="k1"
scmap_label[idx_shuffled[(n1 + 1):(n1 + n2)]]="k2"
scmap_label[idx_shuffled[(n1 + n2 + 1):n]]="k3"





seurat_obj$scmap_label <- scmap_res$combined_labs
DimPlot(seurat_obj, reduction="umap", group.by = "scmap_label") + ggtitle("scmap")



