library(CellChat)
library(Seurat)
library(dplyr)
# ---------------------------
load("./CUPM/results/V2_202507/FunctionalAnalysis/Seurat_meta_cluster.RData")
seurat_meta$Cu_group=paste("k",seurat_meta$RNA_snn_res.0.1,sep="")

# prepare data: use normalized data slot (SCT data preferred if present)
data.input <- GetAssayData(seurat_meta, assay = "SCT", slot = "data")
meta_input <- seurat_meta@meta.data
# build CellChat object
meta_input$labels <- meta_input$Cu_group  # 或 Idents(seurat_meta)
# 5. 构建 CellChat 对象（不要再 as.matrix/as.data.frame）
cellchat <- createCellChat(object = data.input, meta = meta_input, group.by = "labels")
# set database (human)
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
# preprocess
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
# plot circle of interactions (top)
#可视化细胞间通信数量 & 强度
netVisual_circle(cellchat@net$count, vertex.weight = as.numeric(table(cellchat@idents)),
                 weight.scale = TRUE, label.edge = FALSE)
netVisual_circle(cellchat@net$weight, vertex.weight = as.numeric(table(cellchat@idents)),
                 weight.scale = TRUE, label.edge = FALSE)
#查看不同cluster之间的主导信号通路
pathways.show <- cellchat@netP$pathways
head(pathways.show)
#可视化某条信号通路层次结构
netVisual_aggregate(cellchat, signaling = "PARs")
#########################################################
#对比 Cu_high (k0) vs Cu_low (k1)
# 1. 准备数据：只保留 k0 和 k1
cells_use <- rownames(seurat_meta@meta.data)[seurat_meta$Cu_group %in% c("k0","k1")]
data_input_sub <- GetAssayData(seurat_meta, assay = "SCT", slot = "data")[, cells_use]
meta_input_sub <- seurat_meta@meta.data[cells_use, ]

# 2. 设置分组信息
meta_input_sub$labels <- meta_input_sub$Cu_group  # CellChat 分组依据

# 3. 构建 CellChat 对象
cellchat_sub <- createCellChat(object = data_input_sub, meta = meta_input_sub, group.by = "labels")

# 4. 设置数据库（人类）
cellchat_sub@DB <- CellChatDB.human

# 5. 数据预处理
cellchat_sub <- subsetData(cellchat_sub)
cellchat_sub <- identifyOverExpressedGenes(cellchat_sub)
cellchat_sub <- identifyOverExpressedInteractions(cellchat_sub)

# 6. 计算通信概率
cellchat_sub <- computeCommunProb(cellchat_sub)
cellchat_sub <- filterCommunication(cellchat_sub, min.cells = 10)
cellchat_sub <- computeCommunProbPathway(cellchat_sub)
cellchat_sub <- aggregateNet(cellchat_sub)

# 7. 可视化整体通信网络
groupSize <- as.numeric(table(cellchat_sub@idents))
netVisual_circle(cellchat_sub@net$count, vertex.weight = groupSize,
                 weight.scale = TRUE, label.edge = FALSE,
                 title.name = "Number of interactions")

# 8. 挑选特定信号通路（可选）
# pathways.show <- c("VEGF", "FGF")  # 示例
# netVisual_aggregate(cellchat_sub, signaling = pathways.show)
library(CellChat)

# 1. 拆分对象为两个组
cellchat_k0 <- subsetCellChat(cellchat_sub, idents.use = "k0")
cellchat_k1 <- subsetCellChat(cellchat_sub, idents.use = "k1")

# 2. 设置数据库（可重复确保一致）
cellchat_k0@DB <- CellChatDB.human
cellchat_k1@DB <- CellChatDB.human

# 3. 数据预处理
cellchat_k0 <- subsetData(cellchat_k0)
cellchat_k1 <- subsetData(cellchat_k1)

cellchat_k0 <- identifyOverExpressedGenes(cellchat_k0)
cellchat_k1 <- identifyOverExpressedGenes(cellchat_k1)

cellchat_k0 <- identifyOverExpressedInteractions(cellchat_k0)
cellchat_k1 <- identifyOverExpressedInteractions(cellchat_k1)

# 4. 计算通信概率
cellchat_k0 <- computeCommunProb(cellchat_k0)
cellchat_k1 <- computeCommunProb(cellchat_k1)

cellchat_k0 <- filterCommunication(cellchat_k0, min.cells = 10)
cellchat_k1 <- filterCommunication(cellchat_k1, min.cells = 10)

cellchat_k0 <- computeCommunProbPathway(cellchat_k0)
cellchat_k1 <- computeCommunProbPathway(cellchat_k1)

cellchat_k0 <- aggregateNet(cellchat_k0)
cellchat_k1 <- aggregateNet(cellchat_k1)

# 5. 差异通信分析（k0 vs k1）
cellchat_list <- list("Cu_high" = cellchat_k0, "Cu_low" = cellchat_k1)
cellchat_diff <- compareInteractions(cellchat_list, show.legend = F)

# 6. 可视化差异通信
# 整体网络强度对比
netVisual_diffInteraction(cellchat_diff, weight.scale = T)

# 差异信号通路显示
pathways.show <- cellchat_diff@netP$pathways
netVisual_diffInteraction(cellchat_diff, signaling = pathways.show)
