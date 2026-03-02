######################################################################
library(Seurat)
library(SingleR)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(matrixStats)
library(mclust)
library(cowplot)
library(ggplot2)
library(dplyr)
# ---- 参数（可调整阈值） ----
singleR_margin_thresh <- 0.15   # SingleR 最佳分数与第二最佳分数的差阈值（可调）
corr_score_thresh <- 0.2        # Spearman 相关性阈值（按情况调）

#--------------------------bulk data
load("./CUPM/data/RNAseq/GSE183202/GSE183202_Exp_Cli.RData")
exp1=log2(fpkm1+1)
#------------------filter out 低表达的gene
index1=rowMeans(exp1<1)
exp1=exp1[which(index1<0.5),]
cli1[cli1==""]=NA
ref_exp=exp1[,!is.na(cli1$groups)]
ref_label=cli1[!is.na(cli1$groups),"groups"]

######################################################################
#------------------scRNA-seq data
load("./CUPM/data/scRNAseq/scaled_GSE183916.RData")
sct_exp <- GetAssayData(data1, assay = "SCT", layer = "data")
common_genes <- intersect(rownames(sct_exp), rownames(ref_exp))
test_mat <- sct_exp[common_genes, , drop=FALSE]         # genes x cells
ref_mat  <- ref_exp[common_genes, , drop=FALSE]    # genes x samples
bulk_labels <- ref_label
seurat_obj <- data1
# ---- 2) SingleR（每细胞 & pseudo-bulk） ----
# 2A. 直接对每个 cell 运行（genes x cells）
singleR_res_cells <- SingleR(test = test_mat, ref = ref_mat, labels = bulk_labels, prune = TRUE)
# 得到每个细胞的 scores (matrix), labels, 注：pruned.labels 可作为更保守输出

# 2B. 伪-bulk（按 Seurat cluster 或 Idents 聚合）
if(is.null(Idents(seurat_obj))) {
  warning("Seurat 对象尚未设置 Idents (cluster)。将按细胞全部合并作为单一 pseudo-bulk。建议先运行 FindClusters.")
}
clusters <- as.character(Idents(seurat_obj))
unique_clusters <- unique(clusters)
pseudo_bulk <- sapply(unique_clusters, function(cl) {
  rowMeans(as.matrix(test_mat[, clusters == cl, drop=FALSE]))
})
# pseudo_bulk: genes x cluster
singleR_res_pseudo <- SingleR(test = pseudo_bulk, ref = ref_mat, labels = bulk_labels, prune = TRUE)

# ---- 4) 相关性比对（Spearman）: 每个细胞与每个 group 质心做相关性 ----
# 确保 centroids 为 genes x groups，与 test_mat genes x cells
groups <- unique(bulk_labels)
centroids <- sapply(groups, function(g) {
  rowMeans(ref_mat[, bulk_labels == g, drop=FALSE])
})
centroids <- as.matrix(centroids) # genes x groups (注意 orientation)
# 计算 groups x cells 的相关性矩阵
test_dense <- as.matrix(test_mat)
cors <- cor(centroids, test_dense, method = "spearman") # groups x cells
corr_label <- apply(cors, 2, function(x) {
  idx <- which.max(x)
  names(idx) <- NULL
  return(rownames(cors)[idx])
})
corr_score <- apply(cors, 2, max)

# ---- 5) 置信度度量与过滤 ----
# SingleR: 计算每细胞 best score 与 second best score 的差（margin）
single_scores <- singleR_res_cells$scores  # matrix: test x ref.labels (可能 genes? 实际为 cells x ref)
# 注意 SingleR 返回结构：scores 的行与 test 对应
best_idx <- apply(single_scores, 1, which.max)
best_scores <- apply(single_scores, 1, max)
second_scores <- sapply(1:nrow(single_scores), function(i) {
  ord <- order(single_scores[i, ], decreasing = TRUE)
  if(length(ord) >= 2) single_scores[i, ord[2]] else NA
})
single_margin <- best_scores - second_scores
single_label <- singleR_res_cells$labels
# 如果你想使用 pruned.labels 的更保守结果：
single_pruned <- singleR_res_cells$pruned.labels

# 把结果写入 meta.data
meta <- seurat_obj@meta.data
meta$SingleR_label <- single_label
meta$SingleR_pruned <- single_pruned
meta$SingleR_bestScore <- best_scores
meta$SingleR_secondScore <- second_scores
meta$SingleR_margin <- single_margin
meta$Corr_label <- corr_label
meta$Corr_score <- corr_score
seurat_obj@meta.data <- meta

# ---- 6) 共识 (majority vote) 以及置信度过滤示例 ----
all_methods <- data.frame(cell = colnames(test_mat),
                          SingleR = meta$SingleR_label,
                          SingleR_pruned = meta$SingleR_pruned,
                          Corr = meta$Corr_label,
                          Corr_score = meta$Corr_score,
                          SingleR_margin = meta$SingleR_margin,
                          SingleR_score = meta$SingleR_bestScore,
                          stringsAsFactors = FALSE)
# 简单多数投票（若三法一致则标 consensus，否则 "ambiguous"）
get_consensus <- function(s1, s2) {
  labs <- c(s1, s2)
  labs <- labs[!is.na(labs)]
  if(length(labs) == 0) return(NA)
  tb <- table(labs)
  if(max(tb) >= 1) {
    return(names(tb)[which.max(tb)])
  } else {
    return("ambiguous")
  }
}
all_methods$Consensus_label <- mapply(get_consensus, all_methods$SingleR, all_methods$Corr)

# 置信度过滤示例：保守策略
all_methods$SingleR_confident <- all_methods$SingleR_margin >= singleR_margin_thresh & !is.na(all_methods$SingleR)
all_methods$Corr_confident <- all_methods$Corr_score >= corr_score_thresh

# 最终保守注释：至少两法“高置信”且一致；否则标为 "low_confidence"
all_methods$Final_label <- mapply(function(cons, sconf, cconf, s, c) {
  # if majority consensus and at least two confident methods agree -> accept, else low_confidence
  if(cons != "ambiguous") {
    # count how many methods (that are confident) equal consensus
    nconf_agree <- sum(c(sconf & (s==cons),
                         cconf & (c==cons)), na.rm=TRUE)
    if(nconf_agree >= 1) return(cons)
  }
  return("low_confidence")
}, all_methods$Consensus_label, all_methods$SingleR_confident, all_methods$Corr_confident,
   all_methods$SingleR, all_methods$Corr)

# 把最终结果写回 Seurat meta
seurat_obj@meta.data$Consensus_label <- all_methods$Consensus_label
seurat_obj@meta.data$Final_label <- all_methods$Final_label
seurat_obj@meta.data$SingleR_confident <- all_methods$SingleR_confident
seurat_obj@meta.data$Corr_confident <- all_methods$Corr_confident

# ---- 7) 评估方法一致性（confusion matrix，ARI） ----
# 用 ARI 衡量方法间的一致性（只计算非 NA 的细胞）
labels_df <- all_methods %>% select(cell, SingleR, Corr, Final_label)
# ARI pairwise
ari_SC <- adjustedRandIndex(labels_df$SingleR, labels_df$Corr)

ari_summary <- data.frame(pair = c("SingleR_vs_Corr"),
                          ARI = c(ari_SC))
print(ari_summary)

# 显示混淆矩阵（SingleR vs Transfer）
conf_SC <- table(labels_df$SingleR, labels_df$Corr)

# ---- 8) 可视化：UMAP 图（并列三法 + Final_label + QC 指标） ----
# 如果 Seurat 对象没有 UMAP embedding，先跑 PCA + UMAP
if(is.null(seurat_obj@reductions$umap)) {
  if(is.null(seurat_obj@reductions$pca)) {
    seurat_obj <- ScaleData(seurat_obj, verbose=FALSE)
    seurat_obj <- RunPCA(seurat_obj, verbose=FALSE)
  }
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose=FALSE)
}
p1 <- DimPlot(seurat_obj, reduction="umap", group.by = "SingleR_label") + ggtitle("SingleR")
p2 <- DimPlot(seurat_obj, reduction="umap", group.by = "Corr_label") + ggtitle("Correlation")
p3 <- DimPlot(seurat_obj, reduction="umap", group.by = "Final_label") + ggtitle("Final_label (conservative)")

# QC 小提琴图（示例：nFeature_RNA, nCount_RNA, percent.mt）
if(!"percent.mt" %in% colnames(seurat_obj@meta.data)) {
  # 计算 percent.mt 如果未计算
  if(grepl("MT-", rownames(seurat_obj)[1], ignore.case=TRUE) || TRUE) {
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  } else {
    seurat_obj[["percent.mt"]] <- 0
  }
}
p_qc1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA"), group.by = "Final_label", pt.size = 0.1) + ggtitle("nFeature_RNA by Final_label")
p_qc2 <- VlnPlot(seurat_obj, features = c("nCount_RNA"), group.by = "Final_label", pt.size = 0.1) + ggtitle("nCount_RNA by Final_label")
p_qc3 <- VlnPlot(seurat_obj, features = c("percent.mt"), group.by = "Final_label", pt.size = 0.1) + ggtitle("percent.mt by Final_label")

# 合并图
umap_grid <- plot_grid(p1, p2, p3, ncol=2)
qc_grid <- plot_grid(p_qc1, p_qc2, p_qc3, ncol=3)
final_plot <- plot_grid(umap_grid, qc_grid, ncol=1, rel_heights = c(2,1))
final_plot
# 保存图像
ggsave("./CUPM/results/V2_202507/CMS_PM_Subgroup/UMAP_comparison_methods.pdf", final_plot, width = 14, height = 12, dpi = 300)

# ---- 9) 输出表格/结果 ----
write.csv(all_methods, file = "./CUPM/results/V2_202507/CMS_PM_Subgroup/single_cell_annotation_comparison.csv", row.names = FALSE)
write.csv(ari_summary, file = "./CUPM/results/V2_202507/CMS_PM_Subgroup/annotation_ARI_summary.csv", row.names = FALSE)
saveRDS(seurat_obj, file = "./CUPM/results/V2_202507/CMS_PM_Subgroup/seurat_obj_with_annotations.rds")





