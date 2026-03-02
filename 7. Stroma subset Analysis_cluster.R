library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

#---------scRNAseq annotated with cell type
load("./CUPM/data/scRNAseq/scaled_GSE183916.RData")
celltype=data1@meta.data$celltype
celltype[celltype=="Endothelial_cell"]="Stromal_cell"
data1@meta.data$celltype=celltype
#exp1=data1@assays$SCT@scale.data
#loc=data1@meta.data$loc#-----------primary, metastasis
#cuScore1=data1@meta.data$CuScore
#pat_id=data1@meta.data$pat_id
#-------------------------------------------------------------
# 1. subset stromal cells
stromal <- subset(data1, subset = celltype == "Stromal_cell")
# 2. re-SCT
stromal <- SCTransform(stromal, verbose = FALSE)
# 3. PCA
stromal <- RunPCA(stromal, verbose = FALSE)
# 4. graph-based clustering
stromal <- FindNeighbors(stromal, dims = 1:30)
stromal <- FindClusters(stromal, resolution = 0.8)
# 5. UMAP (visualization only)
stromal <- RunUMAP(stromal, dims = 1:30)
DimPlot(stromal, group.by = "seurat_clusters", label = TRUE, repel = TRUE)

# 自定义marker基因集（按需添加或修改）
markers_stromal <- FindAllMarkers(
  stromal,
  only.pos = TRUE,          # 仅保留上调基因
  min.pct = 0.1,           # 在至少1%的细胞中表达
  logfc.threshold = 0.1,    # 对数倍变化阈值
  test.use = "wilcox"       # 使用Wilcoxon秩和检验
)
head(markers_stromal)
#Tumor-associated endothelial cells （TECs)
#“tumor-associated adipocytes (CAAs)”
stromal_markers_list <- list(
  CAF=c("COL1A2", "COL3A1", "COL1A1","VIM"),
  MSCs=c("NT5E", "THY1", "ENG","SOX4","PROM1","ZFP36L2","TFCP2L1"),
  CAAs=c("ADIPOQ","FABP4","PLIN1"),
  TECs=c("FLT1","THY1","VWF","IGFBP4","CALCRL", "GRB10", "HSPG2","VEGFA"),
  pericyte = c("CSPG4", "PDGFRB","MCAM"),
  others=c("PLVAP", "PECAM1", "CDH5","EPB41L1","PLEKHG6","RAB11FIP1","GLG1")
)
#检查marker是否包含在表达谱中
marker_check <- lapply(stromal_markers_list, function(markers){
  present <- markers[markers %in% rownames(stromal)]
  missing <- markers[!markers %in% rownames(stromal)]
  list(
    present = present,
    missing = missing
  )
})

marker_check
#---------------------------------------------------
for (i in names(stromal_markers_list)) {
  stromal <- AddModuleScore(
    stromal,
    features = list(stromal_markers_list[[i]]),
    name = paste0(i, "_score")
  )
}

#判定每个 cluster 与subtype 关联最强
fc_col <- "avg_log2FC"
clusters <- sort(unique(markers_stromal$cluster))
# overlap count
count_mat <- sapply(stromal_markers_list, function(gs){
  sapply(clusters, function(cl){
    sum(markers_stromal$gene[markers_stromal$cluster == cl] %in% gs)
  })
})
rownames(count_mat) <- clusters

# logFC sum
fc_mat <- sapply(stromal_markers_list, function(gs){
  sapply(clusters, function(cl){
    sel <- markers_stromal$cluster == cl & markers_stromal$gene %in% gs
    if(!any(sel)) return(0)
    sum(markers_stromal[[fc_col]][sel], na.rm = TRUE)
  })
})
rownames(fc_mat) <- clusters
choose_subtype <- sapply(seq_along(clusters), function(i){
  cl <- clusters[i]
  cnts <- count_mat[i, ]
  
  if(all(cnts == 0)) return("Other")
  
  best <- names(cnts)[cnts == max(cnts)]
  
  if(length(best) == 1){
    return(best)
  } else {
    # 破局：使用 logFC 总和
    fc_vals <- fc_mat[i, best]
    return(best[which.max(fc_vals)])
  }
})
names(choose_subtype) <- clusters
choose_subtype
annotation_df <- data.frame(
  cluster = clusters,
  subtype = unname(choose_subtype),
  count_mat[as.character(clusters), ],
  row.names = NULL,
  check.names = FALSE
)
print(annotation_df)
#-----------------------------endothelial and caf markers
stromal$seurat_clusters <- as.character(stromal$seurat_clusters)
stromal$subtype <- choose_subtype[stromal$seurat_clusters]
stromal$subtype[is.na(stromal$subtype)] <- "Other"
my_colors <- c(
  "#1b589eff",  # subtype1
  "#d95f02",  # subtype2
  "#d96683ff",  # subtype3
  "#454344ff",  # subtype4
  "#66a61e",  # subtype5
  "#e6ab02"   # subtype6
)
DimPlot(stromal, group.by = "subtype", label = TRUE, repel = TRUE, cols = my_colors)

#score 验证，任选一个
FeaturePlot(stromal, features = c("CAAs_score1","CAF_score1", "MSCs_score1","others_score1","pericyte_score1","TECs_score1"))
save(stromal,file="./CUPM/results/V2_202507/5_CuScore_Stromal/Seurat_Stromal.RData")

FeaturePlot(stromal, features = c("others_score1","TECs_score1"))


