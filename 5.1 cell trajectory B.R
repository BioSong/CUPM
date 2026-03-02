library(Seurat)
library(monocle)
#---------scRNAseq annotated with cell type
#load("./CUPM/results/1_UmapAnalysis/scaled_GSE183916.RData")#---------data1
load("./CUPM/data/scRNAseq/scaled_GSE183916.RData")
cells_with_B <- rownames(data1@meta.data[data1@meta.data$pat_id == "B", ])
data1=subset(data1,cells=cells_with_B)
#--------------------Cuproptosis Related Genes
geneset=read.table("./CUPM/data/CuproptosisRelatedGenes/CRGenes.txt",sep="\t",header=T,as.is=T,quote="\"") 
cuGeneSet <- geneset$SYMBOL[c(which(geneset$Group=="group1"),which(geneset$Group=="group3"))]
#Set these genes in the scRNAseq data object, because the next several functions will depend on them.
data2 <- as.CellDataSet(data1, assay = "SCT", reduction = "umap")
data2 <- estimateSizeFactors(data2)
data2 <- estimateDispersions(data2)
data21 <- setOrderingFilter(data2, cuGeneSet)
plot_ordering_genes(data21)
data21 <- reduceDimension(data21, max_components = 2,
    method = 'DDRTree')
data21 <- orderCells(data21)
plot_cell_trajectory(data21, color_by = "loc")
plot_cell_trajectory(data21, color_by = "celltype")
plot_cell_trajectory(data21, color_by = "CMS")
plot_cell_trajectory(data21, color_by = "CuScore$Final_Score")+
  scale_color_viridis_c(option = "inferno")


#----------------------------
pData(data21)$celltype_highlight <- ifelse(
  pData(data21)$celltype %in% c("Epithelial_cell", "Stromal_cell"),
  pData(data21)$celltype,
  "Other"
)
plot_cell_trajectory(data21, color_by = "celltype_highlight") +
  scale_color_manual(
    values = c(
      "Epithelial_cell" = "#1f77b4",     # 蓝色
      "Stromal_cell" = "#ff7f0e",        # 橙色
      "Other" = "lightgray"              # 灰色
    )
  ) +
  ggtitle("Highlight of Epithelial and Stromal Cells in Full Trajectory") +
  theme_minimal()
#--------------------------------------------
exp_MSN <- FetchData(data1, vars = "MSN", slot = "data", assay = "SCT")
plot_cell_trajectory(data21, color_by = exp_MSN[,1])+
  scale_color_viridis_c(option = "inferno")
#----------------------------------------------------
library(ggplot2)
library(dplyr)

# 提取 MSN 表达、伪时间、状态
df_msn <- data.frame(
  Expression = exprs(data21)["MSN", ],
  Pseudotime = pData(data21)$Pseudotime,
  Branch = as.factor(pData(data21)$State)
)
# log10 转换表达（避免大量 0）
df_msn$log_expr <- log10(df_msn$Expression + 1e-2)

df_msn$sct_expr <- exp_MSN[,1]

# 绘图
ggplot(df_msn, aes(x = Pseudotime, y = sct_expr, color = Branch)) +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE, color = "black") +
  theme_bw(base_size = 13) +
  labs(title = "MSN", x = "Pseudo-time", y = "Relative Expression (log10)") +
  theme(plot.title = element_text(face = "italic", hjust = 0.5),
        legend.position = "right")
