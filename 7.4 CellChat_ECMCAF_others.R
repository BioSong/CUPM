library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(CellChat)

#---------scRNAseq annotated with cell type
load("./CUPM/data/scRNAseq/scaled_GSE183916.RData")
celltype=data1@meta.data$celltype
celltype[celltype=="Endothelial_cell"]="Stromal_cell"
data1@meta.data$celltype=celltype
#-----------------------------------stromal
load('./CUPM/results/V2_202507/5_CuScore_Stromal/Seurat_Stromal.RData')
#-------------------------------CAF
load('./CUPM/results/V2_202507/5_CuScore_Stromal/Seurat_CAF.RData')
#------------------------------
# 1) 把 CAF 子群注回到 stromal
caf_subtypes <- caf@meta.data$subtype
names(caf_subtypes) <- rownames(caf@meta.data)  # cell barcodes
# 初始化列（避免覆盖已有重要信息）
stromal@meta.data$CAF_subtype <- NA
stromal@meta.data[names(caf_subtypes), "CAF_subtype"] <- caf_subtypes
# 2) 构建 interaction_group：CAF 用 subtype 名，其它细胞保留 subtype 名
stromal@meta.data$interaction_group <- ifelse(
  stromal@meta.data$subtype == "CAF",
  as.character(stromal@meta.data$CAF_subtype),   # 含 ECM-CAF, iCAF, apCAF...
  as.character(stromal@meta.data$subtype)       # Epithelial, T cell, Myeloid...
)
stromal@meta.data[is.na(stromal$interaction_group),"interaction_group"]="CAF"
#------------------------------------------------
# 1) 把 stromal 子群注回到 data1
stromal_subtype <- stromal@meta.data$interaction_group
names(stromal_subtype) <- rownames(stromal@meta.data)  # cell barcodes
# 初始化列（避免覆盖已有重要信息）
data1@meta.data$stromal_subtype <- NA
data1@meta.data[names(stromal_subtype), "stromal_subtype"] <- stromal_subtype

# 2) 构建 interaction_group：CAF 用 subtype 名，其它细胞保留 celltype 名
data1@meta.data$interaction_group <- ifelse(
  data1@meta.data$celltype == "Stromal_cell",
  as.character(data1@meta.data$stromal_subtype),   # 含 ECM-CAF, iCAF, apCAF...
  as.character(data1@meta.data$celltype)       # Epithelial, T cell, Myeloid...
)
#------------------------------------------------------
library(CellChat)
# Load the ligand-receptor interaction database for human
data("CellChatDB.human")
# Check the contents of the database
str(CellChatDB.human)
# Load the ligand-receptor interaction database for human
CellChatDB <- CellChatDB.human
#------------------------------------------------------
load("./CUPM/results/V2_202507/5_CuScore_Stromal/cellchat_caf/scdata_cellchat.RData")
data1$CuScore_Final=data1@meta.data$CuScore$Final_Score
data1=subset(data1, subset = loc == "metastasis")
count1=data1@assays$RNA@counts
data1@meta.data$samples=factor(colnames(count1))
data1@meta.data$orig.ident=colnames(count1)
# Convert the Seurat object to a CellChat object
cellchat <- createCellChat(object = data1, assay = "SCT",group.by ="interaction_group")
# Set the CellChatDB as the ligand-receptor interaction database for the analysis
cellchat@DB <- CellChatDB

# Preprocess data to compute the cell communication network
cellchat <- subsetData(cellchat) # Subset the data to focus on the relevant cell types or pathways
cellchat <- identifyOverExpressedGenes(cellchat) # Identify overexpressed genes
cellchat <- identifyOverExpressedInteractions(cellchat) # Identify overexpressed interactions

# Compute the cell-cell communication network
cellchat <- computeCommunProb(cellchat) # Compute the communication probabilities
cellchat <- computeCommunProbPathway(cellchat) # Compute the pathways
cellchat <- aggregateNet(cellchat) # Aggregate the network for visualization

cellchat <- computeNetSimilarity(cellchat)
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")

# Visualize the results
# 假设 cellchat 已完成 computeCommunProb() 和 aggregateNet() 步骤
groupSize <- as.numeric(table(cellchat@idents))

# 使用 netVisual_circle 绘制聚合后的通信网络图（圆形布局）
netVisual_circle(cellchat@net$weight, 
                 vertex.weight = groupSize,
                 weight.scale = TRUE,
                 vertex.label.cex = 1.0,        # 节点标签字体大小
                 edge.width.max = 10,          # 最大边宽
                 title.name = "Aggregated cell-cell communication")
netVisual_circle(cellchat@net$weight, 
                 sources.use="ECMCAF",
                 vertex.weight = groupSize,
                 weight.scale = TRUE,
                 vertex.label.cex = 1.0,        # 节点标签字体大小
                 edge.width.max = 10,          # 最大边宽
                 title.name = "Aggregated cell-cell communication")
netVisual_circle(cellchat@net$weight, 
                 targets.use="Epithelial_cell",
                 vertex.weight = groupSize,
                 weight.scale = TRUE,
                 vertex.label.cex = 1.0,        # 节点标签字体大小
                 edge.width.max = 10,          # 最大边宽
                 title.name = "Aggregated cell-cell communication")
#-----------------------------------------------细胞类型对之间的通信强度热图
netVisual_heatmap(cellchat, 
                  measure = "weight", 
                  color.heatmap = c("lightblue", "navy"), 
                  width = 5, height = 5, 
                  font.size = 8, font.size.title = 10, 
                  title.name = "Communication strength (weight)")
#可视化各细胞群体在信号传导网络中的角色
# 如果未运行，先计算中心性
cellchat <- netAnalysis_computeCentrality(cellchat)

# 然后再画 outgoing role 热图
netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern = "outgoing",           # or "incoming"
  color.heatmap = "Oranges",
  font.size = 10,
  font.size.title = 12,
  width = 6, height = 4
)





