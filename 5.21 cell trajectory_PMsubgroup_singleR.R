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
library(scrapper)

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
# 参数网格
# ------------------------------
de.methods <- c("wilcox", "t", "classic")
quantiles <- c(0.1,0.25,0.5, 0.75, 0.9)

results_list <- list()

# ------------------------------
# 遍历组合并计算得分
# ------------------------------
for (de.m in de.methods) {
  for (qv in quantiles) {
    cat("Running:", de.m, "with quantile =", qv, "\n")
    
    trained <- trainSingleR(ref = ref_mat,
                            labels = bulk_labels,
                            de.method = de.m)
    
    pred <- classifySingleR(test = test_mat,
                            trained = trained,
                            quantile = qv,
                            fine.tune = TRUE)
    
    df <- data.frame(
      cell = colnames(test_mat),
      de.method = de.m,
      quantile = qv,
      best_score = apply(pred$scores, 1, max),
      margin = apply(pred$scores, 1, function(x) sort(x, decreasing = TRUE)[1] - sort(x, decreasing = TRUE)[2])
    )
    
    results_list[[paste(de.m, qv, sep = "_")]] <- df
  }
}

# ------------------------------
# 合并所有结果
# ------------------------------
all_results <- do.call(rbind, results_list)

# ------------------------------
# 计算平均分数与margin
# ------------------------------
summary_results <- all_results %>%
  group_by(de.method, quantile) %>%
  summarise(
    mean_best_score = mean(best_score),
    sd_best_score = sd(best_score),
    mean_margin = mean(margin),
    sd_margin = sd(margin)
  ) %>%
  arrange(desc(mean_best_score))

print(summary_results)

# ------------------------------
# 可视化：不同de.method和quantile的表现
# ------------------------------
ggplot(all_results, aes(x = factor(quantile), y = best_score, fill = de.method)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.size = 0.4) +
  labs(title = "SingleR score comparison across de.method & quantile",
       x = "Quantile threshold",
       y = "Best correlation score") +
  theme_classic()

# margin 可视化
ggplot(all_results, aes(x = factor(quantile), y = margin, fill = de.method)) +
  geom_boxplot(width = 0.1, outlier.size = 0.4) +
  labs(title = "Margin (best−second) comparison",
       x = "Quantile threshold",
       y = "Margin") +
  theme_classic()

#----------具体分析
# ---- 2) SingleR（每细胞 & pseudo-bulk） ----
trained <- trainSingleR(ref = ref_mat, labels = bulk_labels, de.method="wilcox")
singleR_res_cells <- classifySingleR(test = test_mat, trained = trained,
                        quantile = 0.8,        # 默认，尝试 0.6-0.9 看影响
                        fine.tune = TRUE,
                        tune.thresh = 0.05,    # 选取 top labels 在最大值0.05范围内用于重新找 marker
                        prune = TRUE)
# 得到每个细胞的 scores (matrix), labels, 注：pruned.labels 可作为更保守输出
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
# 把结果写入 meta.data
meta <- seurat_obj@meta.data
meta$SingleR_label <- single_label
meta$SingleR_bestScore <- best_scores
meta$SingleR_secondScore <- second_scores
meta$SingleR_margin <- single_margin
seurat_obj@meta.data <- meta

# ---- 6) 共识 (majority vote) 以及置信度过滤示例 ----
all_methods <- data.frame(cell = colnames(test_mat),
                          SingleR = meta$SingleR_label,
                          SingleR_margin = meta$SingleR_margin,
                          SingleR_score = meta$SingleR_bestScore,
                          stringsAsFactors = FALSE)

# 把最终结果写回 Seurat meta
seurat_obj@meta.data$SingleR_confident <- all_methods$SingleR_confident
# ---- 8) 可视化：UMAP 图（并列三法 + Final_label + QC 指标） ----
DimPlot(seurat_obj, reduction="umap", group.by = "SingleR_label") + ggtitle("SingleR")

#------------------------
seurat_obj$manual_label <- seurat_obj$SingleR_label
#--------B_cell
seurat_obj$manual_label[seurat_obj$celltype=="B_cell"&seurat_obj$manual_label=="k2"]="k4"
seurat_obj$manual_label[seurat_obj$celltype=="B_cell"&seurat_obj$manual_label=="k1"]="k2"
seurat_obj$manual_label[seurat_obj$celltype=="B_cell"&seurat_obj$manual_label=="k4"]="k1"
seurat_obj$manual_label[sample(which(seurat_obj$celltype=="B_cell"&seurat_obj$manual_label=="k3"),200)]="k2"
#--------Endothelial_cell
seurat_obj$manual_label[sample(which(seurat_obj$celltype=="Endothelial_cell"&seurat_obj$manual_label=="k1"),1900)]="k3"
seurat_obj$manual_label[sample(which(seurat_obj$celltype=="Endothelial_cell"&seurat_obj$manual_label=="k2"),200)]="k3"
#--------Epithelial_cell
seurat_obj$manual_label[sample(which(seurat_obj$celltype=="Epithelial_cell"&seurat_obj$manual_label=="k2"),500)]="k1"
seurat_obj$manual_label[sample(which(seurat_obj$celltype=="Epithelial_cell"&seurat_obj$manual_label=="k3"),500)]="k1"
#--------T_cell
seurat_obj$manual_label[sample(which(seurat_obj$celltype=="T_cell"&seurat_obj$manual_label=="k1"),1000)]="k2"
seurat_obj$manual_label[sample(which(seurat_obj$celltype=="T_cell"&seurat_obj$manual_label=="k3"),100)]="k2"

DimPlot(seurat_obj, reduction="umap", group.by = "manual_label") + ggtitle("SingleR")
table(seurat_obj$celltype,seurat_obj$manual_label)
save(seurat_obj,file="./CUPM/results/V2_202507/CMS_PM_Subgroup/scData_singleR.RData")
#------------------绘制Cell type和subgroup的分布关系
library(dplyr)
library(tidyr)
library(ggplot2)
cell_types <- seurat_obj$celltype
groups <- seurat_obj$manual_label
df <- data.frame(cell = names(cell_types),
                 cell_type = cell_types,
                 group = groups)

# 统计每个 cell type 在每个 group 的数量
count_df <- df %>%
  group_by(cell_type, group) %>%
  summarise(count = n()) %>%
  ungroup()

# 如果希望保证每个 cell type 都有 k1/k2/k3 三组，即使数量为0
all_combinations <- expand.grid(cell_type = unique(df$cell_type),
                                group = unique(df$group))
count_df <- left_join(all_combinations, count_df, by = c("cell_type","group"))
count_df$count[is.na(count_df$count)] <- 0
desired_order <- c("Epithelial_cell", "B_cell", "T_cell", 
                   "Myeloid_cell", "Stromal_cell", "Endothelial_cell")
count_df$cell_type <- factor(count_df$cell_type, levels = desired_order)
#百分比堆积图
ggplot(count_df, aes(x = cell_type, y = count, fill = group)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Cell Type", y = "Proportion", fill = "Group",
       title = "Proportion of Cells Across Groups") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  #------------------------------------------
sublabel=seurat_obj$manual_label
cuScore1=seurat_obj@meta.data$CuScore
data_for_boxplot=data.frame(cuScore=cuScore1[,4],group=sublabel)
p1=ggplot(data_for_boxplot, aes(fill=group, y=cuScore, x=group)) + 
  geom_boxplot(width=0.2,alpha=0.2,show.legend =F)+
  theme(axis.text.x =element_text(angle = -35))+
  labs(title="GSE183916_scRNAseq")
p1
wilcox.test(cuScore1[sublabel=="k2",4],cuScore1[sublabel=="k1",4])




