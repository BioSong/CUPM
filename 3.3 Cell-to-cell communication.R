#-------------------------------------------------终端
# 创建一个独立的环境，使用 Python 3.10
conda create -n cellchat_env python=3.10 -y

# 激活环境
conda activate cellchat_env

# 安装需要的 Python 包（推荐使用 conda 源）
conda install conda-forge::numpy
conda install conda-forge::pip
install conda-forge::scikit-learn
conda install conda-forge::umap-learn

#------------------------------------R
library(reticulate)
Sys.unsetenv("RETICULATE_PYTHON")
# 设置使用刚刚创建的 conda 环境
use_condaenv("cellchat_env", required = TRUE)
# 检查当前 Python 配置
py_config()
# 验证模块是否可用
py_module_available("numpy")       # TRUE
py_module_available("umap")        # TRUE
py_module_available("sklearn")     # TRUE

library(CellChat)
# Load the ligand-receptor interaction database for human
data("CellChatDB.human")
# Check the contents of the database
str(CellChatDB.human)
# Load the ligand-receptor interaction database for human
CellChatDB <- CellChatDB.human

load("./CUPM/data/scRNAseq/scaled_GSE183916.RData")
count1=data1@assays$RNA@counts
data1@meta.data$samples=factor(colnames(count1))
data1@meta.data$orig.ident=colnames(count1)
cells_with_B <- rownames(data1@meta.data[data1@meta.data$pat_id == "B", ])
data1_B=subset(data1,cells=cells_with_B)
# Convert the Seurat object to a CellChat object
cellchat <- createCellChat(object = data1, assay = "SCT",group.by ="celltype")
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
                 sources.use="T_cell",
                 vertex.weight = groupSize,
                 weight.scale = TRUE,
                 vertex.label.cex = 1.0,        # 节点标签字体大小
                 edge.width.max = 10,          # 最大边宽
                 title.name = "Aggregated cell-cell communication")
netVisual_circle(cellchat@net$weight, 
                 targets.use="T_cell",
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

# 或者画 incoming role 热图
netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern = "incoming",
  color.heatmap = "Blues"
)
# 画出所有通路中，细胞作为 signal sender 的活跃程度
# 获取前10个最强通信通路
infoflow <- cellchat@netP$pathways.info
top_pathways <- head(infoflow$pathway[order(infoflow$info.flow, decreasing = TRUE)], 10)

netVisual_bubble(
  cellchat, 
  sources.use = NULL,  # 使用所有源细胞类型
  targets.use = NULL,   # 使用所有目标细胞类型
  signaling = NULL,     # 显示所有信号通路
  remove.isolate = TRUE # 移除无交互的条目
) + 
coord_flip()  # 翻转坐标轴使通路在y轴

