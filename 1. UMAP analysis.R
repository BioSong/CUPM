load("./CUPM/data/scRNAseq/scaled_GSE183916.RData")
# 5. 可视化UMAP结果
# 5.1 根据patient_id
DimPlot(data1, reduction = "umap", group.by = "pat_id2") + ggtitle("UMAP Visualization")
DimPlot(data1, label = FALSE,  group.by="pat_id2", 
          cols = c('B-primary' = "#ff75a0", "B-metastasis" = "#93329e",'C-metastasis' = "#f0a500", 'D-metastasis' = "#1687a7", 'E-metastasis' = "#666666", "F-metastasis" = "#95e1d3"))
#---------------只展示patient B
selected_samples <- c("B-metastasis", "B-primary")
selected_cells <- rownames(data1@meta.data[data1@meta.data$pat_id2 %in% selected_samples, ])
DimPlot(data1, reduction = "umap",cells = selected_cells, group.by = "pat_id2")
# 5.2 根据sample location
DimPlot(data1, reduction = "umap", group.by = "loc") + ggtitle("UMAP Visualization")

# 6. 如果想自定义群体分组，也可以指定group.by参数
# 例如：DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type")
#--------------------进行细胞类型注释
#-----------------------------
library(future)
plan("multicore", workers = 8)
cluster_markers <- FindAllMarkers(
  data1,
  only.pos = TRUE,          # 仅保留上调基因
  min.pct = 0.1,           # 在至少1%的细胞中表达
  logfc.threshold = 0.1,    # 对数倍变化阈值
  test.use = "wilcox"       # 使用Wilcoxon秩和检验
)
plan(sequential)
# 查看每个cluster的前5个标记基因
top5_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

print(top5_markers)

# 5.2 根据sample location
DimPlot(data1, reduction = "umap", group.by = "celltype") + ggtitle("UMAP_Visualization")
save(data1,file="./CUPM/results/1_UmapAnalysis/GSE183916_scaled_MAC.RData")
# 绘制标记基因点图
DotPlot(data1, 
        features = unique(unlist(canonical_markers)),
        cols = c("blue", "red"),
        dot.scale = 6,
        group.by = "pat_id2") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

# 展示关键标记的UMAP
FeaturePlot(data1, 
            features = c("IGHA1", "CD40", "CCL7", "MKI67", "SPP1", "CD4"),
            reduction = "umap",
            order = TRUE,
            combine = FALSE)

#---------------------
