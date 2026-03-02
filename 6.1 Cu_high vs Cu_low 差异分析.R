library(Seurat)
library(dplyr)
library(stringr)  # 用于换行处理

load("./CUPM/results/V2_202507/FunctionalAnalysis/Seurat_meta_cluster.RData")
seurat_meta$Cu_group=paste("K",seurat_meta$RNA_snn_res.0.1,sep="")
meta <- seurat_meta@meta.data
cu_col="Cu_score"
cluster_means <- meta %>%
  as.data.frame() %>%
  group_by(Cu_group) %>%
  summarise(mean_Cu = mean(.data[[cu_col]], na.rm = TRUE),
            median_Cu = median(.data[[cu_col]], na.rm = TRUE),
            n_cells = n()) %>%
  arrange(desc(mean_Cu))

print(cluster_means)
# 3) 选择平均 Cu score 最高与最低的 cluster（各 1 个）
high_cluster <- as.character(cluster_means$Cu_group[1])
low_cluster  <- as.character(cluster_means$Cu_group[nrow(cluster_means)])

message("Selected high cluster: ", high_cluster)
message("Selected low cluster : ", low_cluster)

# 4) 标注每个 cell 为 Cu_high / Cu_low / other
seurat_meta$Cu_group_highlow <- ifelse(seurat_meta$Cu_group == high_cluster, "Cu_high",
                                ifelse(seurat_meta$Cu_group == low_cluster, "Cu_low", "other"))

# 简单统计
table(seurat_meta$Cu_group_highlow)

# 5) 提取子集对象（方便后续DE）
seurat_high <- subset(seurat_meta, subset = Cu_group_highlow == "Cu_high")
seurat_low  <- subset(seurat_meta, subset = Cu_group_highlow == "Cu_low")

# 6) 快速可视化：boxplot 显示三组 cluster mean/median（可选）
library(ggplot2)
p <- ggplot(seurat_meta@meta.data, aes_string(x = "Cu_group", y = cu_col)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_jitter(width = 0.2, size = 0.3, alpha = 0.4) +
  theme_minimal() +
  labs(title = "Cu score per Cu_group (metastasis)", y = cu_col, x = "Cu_group")
ggsave("./CUPM/results/V2_202507/FunctionalAnalysis/clustering_results/CuScore_per_cluster_boxplot.pdf", p, width = 6, height = 4)
#---------------------------K0和K1，差异分析
Idents(seurat_meta) <- "Cu_group"
deg_results <- FindMarkers(
  object = seurat_meta,
  ident.1 = "K0",
  ident.2 = "K1",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
deg_filtered <- deg_results %>%
  dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)

# 分别提取上调和下调的基因
up_genes <- rownames(deg_filtered[deg_filtered$avg_log2FC > 0.25, ])
down_genes <- rownames(deg_filtered[deg_filtered$avg_log2FC < -0.25, ])
#----------------------------------火山图
geneset=read.table("./CUPM/data/CuproptosisRelatedGenes/CRGenes.txt",sep="\t",header=T,as.is=T,quote="\"")
Cu_genes=unique(geneset$SYMBOL[c(which(geneset$Group=="group1"),which(geneset$Group=="group3"))])
volcano_df <- deg_filtered %>%
  mutate(
    gene = rownames(deg_filtered),
    logFC = avg_log2FC,
    negLogFDR = -log10(p_val_adj),
    regulation = case_when(
      avg_log2FC > 0.25 & p_val_adj < 0.05 ~ "Up",
      avg_log2FC < -0.25 & p_val_adj < 0.05 ~ "Down",
      TRUE ~ "NotSig"
    ),
    Cu_related = ifelse(rownames(deg_filtered) %in% Cu_genes, "Cu_gene", "Other")
  )
volcano_df$label <- ifelse(
  volcano_df$regulation != "NotSig" & volcano_df$Cu_related == "Cu_gene",
  volcano_df$gene,
  ""
)
# 绘制火山图
p <- ggplot(volcano_df, aes(x = logFC, y = negLogFDR)) +
  geom_point(aes(color = regulation), alpha = 0.6, size = 1.5) +
  geom_point(
    data = subset(volcano_df, Cu_related == "Cu_gene" & regulation != "NotSig"),
    aes(x = logFC, y = negLogFDR),
    color = "black", size = 2.8, shape = 21, stroke = 0.8, fill = "yellow"
  ) +
  ggrepel::geom_text_repel(
    data = subset(volcano_df, Cu_related == "Cu_gene" & regulation != "NotSig"),
    aes(label = label),
    size = 3,
    max.overlaps = 50
  ) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NotSig" = "grey")) +
  coord_cartesian(xlim = c(-10, 10)) +  # ✅ 控制x轴范围
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_classic() +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P value",
    color = "Regulation",
    title = "Volcano Plot Highlighting Significant Cuproptosis-related Genes"
  )

p
#---------------------------功能分析
library(org.Hs.eg.db)
library(clusterProfiler)
upgene_ids <- bitr(
  up_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
upgene_ids <- upgene_ids$ENTREZID
#------------------------
kk_upgene <- enrichKEGG(
  gene         = upgene_ids,
  organism     = "hsa",
  pvalueCutoff = 0.05
)
kk_upgene@result <- kk_upgene@result[
  !grepl("^hsa05", kk_upgene@result$ID),
]
# Dotplot
dotplot(kk_upgene)
# Barplot
barplot(kk_upgene)
#---------------------------------------
library(dplyr)
library(ggalluvial)
library(ggplot2)

# 转成 data.frame
df <- as.data.frame(kk_upgene)
df <- df %>% arrange(p.adjust)
dotplot_order <- df %>% arrange(p.adjust) %>% pull(Description)

df$Description <- factor(df$Description, levels = dotplot_order)

df_sankey <- df %>%
    mutate(Width = -log10(p.adjust))

ggplot(df_sankey,
              aes(axis1 = category, axis2 = subcategory, axis3 = Description,
                  y = Width)) +
    geom_alluvium(aes(fill = category), width = 1/12, alpha = 0.7) +
    geom_stratum(width = 1/12, fill = "grey80", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("Category", "SubCategory", "Description"),
                     expand = c(0.15, 0.05)) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 12, face = "bold")) +
    ggtitle("KEGG Pathway Classification Sankey Plot")


