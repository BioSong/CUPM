library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

#------------------------------------
load('./CUPM/results/V2_202507/5_CuScore_Stromal/Seurat_Stromal.RData')
caf <- subset(stromal, subset = subtype == "CAF")
caf@meta.data$subtype=NULL
caf@meta.data <- caf@meta.data[, !grepl("_score1$", colnames(caf@meta.data))]
caf$CuScore_Final=caf@meta.data$CuScore$Final_Score
#----------------------------------每种cluster比较primary和PM间的CuScore差异
cluster=stromal@meta.data$SCT_snn_res.0.8
loc[which(loc=="primary")]="1_primary";loc[which(loc=="metastasis")]="2_metastasis"
df_caf <- caf@meta.data %>%
  select(cluster = SCT_snn_res.0.8, 
         CuScore = CuScore_Final,
         loc)
# 创建分组变量
df_caf <- df_caf %>%
  mutate(
    group = ifelse(loc == "primary",
                   "Primary",
                   paste0("PM_CAF_", cluster))
  )
df_caf$group <- factor(
  df_caf$group,
  levels = c("Primary",
             "PM_CAF_1", "PM_CAF_5", "PM_CAF_10", "PM_CAF_11", "PM_CAF_14", "PM_CAF_17")
)
# 箱线图绘制
p <- ggplot(df_caf, aes(x = group, y = CuScore, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.6, alpha = 0.5) +
  theme_classic() +
  labs(x = "", y = "CuScore") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

print(p)
#-----------------------------------------------
t.test(df_caf[df_caf$group=="PM_CAF_1","CuScore"],df_caf[df_caf$group=="Primary","CuScore"])
t.test(df_caf[df_caf$group=="PM_CAF_5","CuScore"],df_caf[df_caf$group=="Primary","CuScore"])
t.test(df_caf[df_caf$group=="PM_CAF_10","CuScore"],df_caf[df_caf$group=="Primary","CuScore"])
t.test(df_caf[df_caf$group=="PM_CAF_11","CuScore"],df_caf[df_caf$group=="Primary","CuScore"])
t.test(df_caf[df_caf$group=="PM_CAF_14","CuScore"],df_caf[df_caf$group=="Primary","CuScore"])
t.test(df_caf[df_caf$group=="PM_CAF_17","CuScore"],df_caf[df_caf$group=="Primary","CuScore"])
#----------------------------------------------
caf <- subset(caf, subset = loc == "metastasis")
markers_caf <- FindAllMarkers(
  caf,
  only.pos = TRUE,          # 仅保留上调基因
  min.pct = 0.1,           # 在至少1%的细胞中表达
  logfc.threshold = 0.1,    # 对数倍变化阈值
  test.use = "wilcox"       # 使用Wilcoxon秩和检验
)
head(markers_caf)
# CAF subtype signatures
caf_markers_list <- list(
  myCAF = c("ACTA2", "TAGLN", "COL1A1", "COL1A2","CTGF","TNC","PDGFRB"),
  iCAF = c("IL6","TIMP1","AEBP1","LIF", "CXCL12","COL14A1","C1R"),
  apCAF = c("SLPI", "CD74","NKAIN4","IRF4"),
  ECMCAF=c( "COL11A1", "FN1", "SPARC", "POSTN", "MMP11", "PXDN","THBS2")
)

for (i in names(caf_markers_list)) {
  caf <- AddModuleScore(
    caf,
    features = list(caf_markers_list[[i]]),
    name = paste0(i, "_score")
  )
}

#判定每个 cluster 与subtype 关联最强
fc_col <- "avg_log2FC"
clusters <- sort(unique(markers_caf$cluster))
# overlap count
count_mat <- sapply(caf_markers_list, function(gs){
  sapply(clusters, function(cl){
    sum(markers_caf$gene[markers_caf$cluster == cl] %in% gs)
  })
})
rownames(count_mat) <- clusters

# logFC sum
fc_mat <- sapply(caf_markers_list, function(gs){
  sapply(clusters, function(cl){
    sel <- markers_caf$cluster == cl & markers_caf$gene %in% gs
    if(!any(sel)) return(0)
    sum(markers_caf[[fc_col]][sel], na.rm = TRUE)
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
#-------------------------
index=as.character(caf$seurat_clusters)
subtype <- choose_subtype[index]
names(subtype)=names(caf$seurat_clusters)
caf$subtype=subtype
my_colors <- c(
  "#1b589eff",  # subtype1
  "#d95f02",  # subtype2
  "#d96683ff",  # subtype3
  "#454344ff",  # subtype4
  "#66a61e",  # subtype5
  "#e6ab02"   # subtype6
)
DimPlot(caf, group.by = "subtype", label = TRUE, repel = TRUE, cols = my_colors)
#---------------------------------
cluster_col <- "seurat_clusters"
#DotPlot：展示 ECM vs myCAF marker 表达差异
DotPlot(caf, features = c(caf_markers_list$myCAF, caf_markers_list$ECMCAF),
        group.by = cluster_col) +
  scale_color_gradient(low = "#e7eaf6", high = "#273c75") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
VlnPlot(caf, features = c("ECMCAF_score1", "myCAF_score1"),
        group.by = cluster_col) +
  theme_classic() +
  ylab("Module score")

FeaturePlot(caf, features = "ECMCAF_score1", reduction = "umap",
            pt.size = 0.6, min.cutoff="q10", max.cutoff="q90") &
  scale_color_gradient(low = "#d6e4ff", high = "#0a2351")

FeaturePlot(caf, features = "myCAF_score1", reduction = "umap",
            pt.size = 0.6, min.cutoff="q10", max.cutoff="q90") &
  scale_color_gradient(low = "#ffe6e6", high = "#7a0a0a")
save(caf,file="./CUPM/results/V2_202507/5_CuScore_Stromal/Seurat_CAF.RData")
