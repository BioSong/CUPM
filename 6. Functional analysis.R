library(Seurat)
library(cluster)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(factoextra)
#------------------------------------------------------
load("./CUPM/results/V2_202507/CMS_PM_Subgroup/scData_singleR.RData")
#-------绘制primary和PM CUscore的密度分布曲线
cuscore1 <- seurat_obj@meta.data$CuScore[,4]
loc_labels <- seurat_obj@meta.data$loc  # primary 和 metastasis 标签
# 创建数据框，方便绘图
plot_data <- data.frame(Cu_score = cuscore1, loc = loc_labels)
# 绘制密度分布图
ggplot(plot_data, aes(x = Cu_score, fill = loc)) +
  geom_density(alpha = 0.5) +  # alpha 控制透明度
  scale_fill_manual(values = c("primary" = "blue", "metastasis" = "red")) +  # 设置颜色
  theme_minimal() +  # 主题
  labs(title = "Density Plot of Cu Score by Loc", x = "Cu Score", y = "Density") +
  theme(legend.title = element_blank())  # 去除图例标题
#-----------------利用kmeans对cells进行分类
geneset=read.table("./CUPM/data/CuproptosisRelatedGenes/CRGenes.txt",sep="\t",header=T,as.is=T,quote="\"")
Cu_genes=unique(geneset$SYMBOL[c(which(geneset$Group=="group1"),which(geneset$Group=="group3"))])
Cu_expression_data <- GetAssayData(seurat_obj, assay = "SCT", slot = "data")
Cu_genes=intersect(Cu_genes,rownames(Cu_expression_data))
Cu_expression_data=as.matrix(Cu_expression_data[Cu_genes,])
summary(apply(Cu_expression_data,1,sd))
#--------------------------对Cell聚类需把Cell作为行
data_for_clustering <- t(Cu_expression_data)
data_scaled <- scale(data_for_clustering)
dim(data_scaled)      # 应为：Cell × Gene
head(data_scaled[,1:5])
summary(apply(data_scaled, 2, sd))  # 每个基因方差应 ~1
data_scaled=data_scaled[,which(colMeans(is.na(data_scaled))==0)]
#_---------------------------------------利用Seurat进行聚类
Cu_genes=colnames(data_scaled)
DefaultAssay(seurat_obj) <- "RNA"
#---------------------挑选metastasis的样本
met_cells <- WhichCells(seurat_obj, expression = loc == "metastasis")
seurat_meta <- subset(seurat_obj, cells = met_cells)
# 3. 数据标准化 & 缩放
seurat_meta <- ScaleData(seurat_meta, features = Cu_genes)
# 4. PCA降维
seurat_meta <- RunPCA(seurat_meta, features = Cu_genes)
# 5. 构建近邻图
seurat_meta <- FindNeighbors(seurat_meta, dims = 1:10)
# 6. 聚类
seurat_meta <- FindClusters(seurat_meta, resolution = 0.1)
# 7. 可视化检查
DimPlot(seurat_meta, reduction="umap", group.by = "seurat_clusters") + ggtitle("seurat_clusters")
#-----------------------------------------比较CU score的差别
seurat_meta$Cu_score <- seurat_meta@meta.data$CuScore[, 4] 
FeaturePlot(seurat_meta, reduction="umap", features = "Cu_score") + 
  ggtitle("Cu score in Metastasis Cells")
ggplot(seurat_meta@meta.data, aes(x = Cu_score, fill = seurat_clusters)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  labs(title = "Cu score Density in Metastasis Cells",
       x = "Cu score (singscore)", y = "Density") +
  scale_fill_brewer(palette = "Set2", name = "Cluster")
