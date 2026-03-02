library(Seurat)
library(dplyr)
library(pheatmap)
library(circlize)
#---------------------------------------------------------------
load("./CUPM/data/RNAseq/GSE183202/GSE183202_Exp_Cli.RData")
exp1=log2(fpkm1+1)
#------------------filter out 低表达的gene
index1=rowMeans(exp1<1)
exp1=exp1[which(index1<0.5),]
cli1[cli1==""]=NA
ref_exp=exp1[,!is.na(cli1$groups)]
ref_label=cli1[!is.na(cli1$groups),"groups"]
####################################################
source("./CUPM/program/Identification of DEGs.R")
geneid=rownames(ref_exp)
#######-------------------------k1
DEG_k11=IDEGs(ref_exp[,ref_label=="k1"],ref_exp[,ref_label=="k2"],geneid,"RNAseq_fpkm")
DEG_k11=DEG_k11$DEGs
marker_k11=DEG_k11[which(DEG_k11$logFC>1&DEG_k11$fdr<0.05),]
DEG_k12=IDEGs(ref_exp[,ref_label=="k1"],ref_exp[,ref_label=="k3"],geneid,"RNAseq_fpkm")
DEG_k12=DEG_k12$DEGs
marker_k12=DEG_k12[which(DEG_k12$logFC>1&DEG_k12$fdr<0.05),]
#marker_k1=rbind(marker_k11,marker_k12)
marker_k1=intersect(marker_k11$geneid,marker_k12$geneid)
#######-------------------------k2
DEG_k21=IDEGs(ref_exp[,ref_label=="k2"],ref_exp[,ref_label=="k1"],geneid,"RNAseq_fpkm")
DEG_k21=DEG_k21$DEGs
marker_k21=DEG_k21[which(DEG_k21$logFC>1&DEG_k21$PValue<0.05),]
DEG_k22=IDEGs(ref_exp[,ref_label=="k2"],ref_exp[,ref_label=="k3"],geneid,"RNAseq_fpkm")
DEG_k22=DEG_k22$DEGs
marker_k22=DEG_k22[which(DEG_k22$logFC>1&DEG_k22$PValue<0.05),]
#marker_k2=rbind(marker_k21,marker_k22)
marker_k2=intersect(marker_k21$geneid,marker_k22$geneid)
#######-------------------------k3
DEG_k31=IDEGs(ref_exp[,ref_label=="k3"],ref_exp[,ref_label=="k1"],geneid,"RNAseq_fpkm")
DEG_k31=DEG_k31$DEGs
marker_k31=DEG_k31[which(DEG_k31$logFC>1&DEG_k31$fdr<0.05),]
DEG_k32=IDEGs(ref_exp[,ref_label=="k3"],ref_exp[,ref_label=="k2"],geneid,"RNAseq_fpkm")
DEG_k32=DEG_k32$DEGs
marker_k32=DEG_k32[which(DEG_k32$logFC>1&DEG_k32$fdr<0.05),]
#marker_k3=rbind(marker_k31,marker_k32)
marker_k3=intersect(marker_k31$geneid,marker_k32$geneid)
#------------------------------------
sigGene=unique(c(marker_k1,marker_k2,marker_k3))
#------------------------------
library(pheatmap)
ref_exp1=ref_exp[sigGene,c(which(ref_label=="k1"),which(ref_label=="k2"),which(ref_label=="k3"))]
anno_col=data.frame(Group=c(rep("k1",length(which(ref_label=="k1"))),rep("k2",length(which(ref_label=="k2"))),rep("k3",length(which(ref_label=="k3")))))
rownames(anno_col)=colnames(ref_exp1)
expr_z <- t(scale(t(ref_exp1))) 
pheatmap(expr_z,cluster_rows = F,cluster_cols = F,
			annotation_col = anno_col,annotation_legend=T,
			show_colnames = F,show_rownames=F)
#--------------------------------------------------------------scRNAseq
load("./CUPM/results/V2_202507/CMS_PM_Subgroup/scData_singleR.RData")
all_genes <- rownames(seurat_obj)
marker_k1=intersect(marker_k1,all_genes)
marker_k2=intersect(marker_k2,all_genes)
marker_k3=intersect(marker_k3,all_genes)
Idents(seurat_obj) <- "manual_label"   # 分组列名为 group
# ✅ 存结果的列表
deg_results <- list()
# ✅ k1基因：在 k1 vs k2 / k1 vs k3 比较
deg_results$k1_vs_k2 <- FindMarkers(
  seurat_obj,
  ident.1 = "k1",
  ident.2 = "k2",
  features = marker_k1,
  test.use = "wilcox",
  min.pct = 0,
  logfc.threshold = 0
)

deg_results$k1_vs_k3 <- FindMarkers(
  seurat_obj,
  ident.1 = "k1",
  ident.2 = "k3",
  features = marker_k1,
  test.use = "wilcox",
  min.pct = 0,
  logfc.threshold = 0
)

# ✅ k2基因：在 k2 vs k1 / k2 vs k3 比较
deg_results$k2_vs_k1 <- FindMarkers(
  seurat_obj,
  ident.1 = "k2",
  ident.2 = "k1",
  features = marker_k2,
  test.use = "wilcox",
  min.pct = 0,
  logfc.threshold = 0
)

deg_results$k2_vs_k3 <- FindMarkers(
  seurat_obj,
  ident.1 = "k2",
  ident.2 = "k3",
  features = marker_k2,
  test.use = "wilcox",
  min.pct = 0,
  logfc.threshold = 0
)

# ✅ k3基因：在 k3 vs k1 / k3 vs k2 比较
deg_results$k3_vs_k1 <- FindMarkers(
  seurat_obj,
  ident.1 = "k3",
  ident.2 = "k1",
  features = marker_k3,
  test.use = "wilcox",
  min.pct = 0,
  logfc.threshold = 0
)

deg_results$k3_vs_k2 <- FindMarkers(
  seurat_obj,
  ident.1 = "k3",
  ident.2 = "k2",
  features = marker_k3,
  test.use = "wilcox",
  min.pct = 0,
  logfc.threshold = 0
)
# ✅ 所有结果保存在列表 deg_results
deg_results

##--------------------------------------------k1
# 取 logFC
logfc_k1_k2 <- deg_results$k1_vs_k2
logfc_k1_k3 <- deg_results$k1_vs_k3
gene1 <- intersect(rownames(logfc_k1_k2),rownames(logfc_k1_k3))
# 构建 matrix：行是 gene，列是对比组
heat_mat1 <- cbind(
  "k1_vs_k2" = logfc_k1_k2[gene1,"avg_log2FC"],
  "k1_vs_k3" = logfc_k1_k3[gene1,"avg_log2FC"]
)
rownames(heat_mat1) <- gene1
# 绘制 heatmap
plot1=pheatmap(
  heat_mat1,
  cluster_rows = F,
  cluster_cols = F,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-2.5, 2.5, length.out = 101),
  main = "K1 signature genes logFC heatmap"
)
pdf("./CUPM/results/V2_202507/CMS_PM_Subgroup/sc_pheatmap_k1.pdf", width = 16, height = 10)
print(plot1)
dev.off()
##--------------------------------------------k2
# 取 logFC
logfc_k2_k1 <- deg_results$k2_vs_k1
logfc_k2_k3 <- deg_results$k2_vs_k3
gene2 <- intersect(rownames(logfc_k2_k1),rownames(logfc_k2_k3))
# 构建 matrix：行是 gene，列是对比组
heat_mat2 <- cbind(
  "k2_vs_k1" = logfc_k2_k1[gene2,"avg_log2FC"],
  "k2_vs_k3" = logfc_k2_k3[gene2,"avg_log2FC"]
)
rownames(heat_mat2) <- gene2
# 绘制 heatmap
plot2=pheatmap(
  heat_mat2,
  cluster_rows = F,
  cluster_cols = F,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-2.5, 2.5, length.out = 101),
  main = "K2 signature genes logFC heatmap"
)
pdf("./CUPM/results/V2_202507/CMS_PM_Subgroup/sc_pheatmap_k2.pdf", width = 16, height = 10)
print(plot2)
dev.off()
##--------------------------------------------k3
# 取 logFC
logfc_k3_k1 <- deg_results$k3_vs_k1
logfc_k3_k2 <- deg_results$k3_vs_k2
gene3 <- intersect(rownames(logfc_k3_k1),rownames(logfc_k3_k2))
# 构建 matrix：行是 gene，列是对比组
heat_mat3 <- cbind(
  "k3_vs_k1" = logfc_k3_k1[gene3,"avg_log2FC"],
  "k3_vs_k2" = logfc_k3_k2[gene3,"avg_log2FC"]
)
rownames(heat_mat3) <- gene3
# 绘制 heatmap
plot3=pheatmap(
  heat_mat3,
  cluster_rows = F,
  cluster_cols = F,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-2.5, 2.5, length.out = 101),
  main = "K3 signature genes logFC heatmap"
)
pdf("./CUPM/results/V2_202507/CMS_PM_Subgroup/sc_pheatmap_k3.pdf", width = 16, height = 10)
print(plot3)
dev.off()
########################################绘制
library(ggplot2)
library(dplyr)
library(tidyr)

# 将矩阵转换为 data.frame
df1 <- as.data.frame(heat_mat1)
df2 <- as.data.frame(heat_mat2)
df3 <- as.data.frame(heat_mat3)

# 添加 group 列
df1$group <- "k1"
df2$group <- "k2"
df3$group <- "k3"

# 合并并转换为长格式
long_df <- bind_rows(df1, df2, df3) %>%
  pivot_longer(
    cols = -group,
    names_to = "comparison",
    values_to = "logFC"
  ) %>%
  mutate(label = paste(group, comparison, sep = "_"))

# 绘制密度图
ggplot(long_df, aes(x = logFC, color = label, fill = label)) +
  geom_density(alpha = 0.3) +
  theme_classic() +
  xlim(-2.5, 2.5) +
  labs(x = "logFC", y = "Density", color = "Group_Comparison", fill = "Group_Comparison")
#-----------------------------------
heat_mat=rbind(heat_mat1,NA,NA,heat_mat2,NA,NA,heat_mat3)



