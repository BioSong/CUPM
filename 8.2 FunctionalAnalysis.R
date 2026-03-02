# install.packages("Seurat")
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("matrixStats")
# install.packages("tibble")
# install.packages("tidyr")

# BiocManager::install(c("decoupleR", "progeny", "dorothea", "SingleCellExperiment"))

library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(matrixStats)
library(SingleCellExperiment)
library(decoupleR)
library(progeny)
library(dorothea)
#-------------------------------------------
load("./CUPM/results/V2_202507/6_CuScore_upstream/scdata_copykat.RData")
load("./CUPM/results/V2_202507/6_CuScore_upstream/cds_pseudotime.RData")

#2) monocle3：只用 malignant tumor（上皮+aneuploid）做轨迹
library(monocle3)
library(SeuratWrappers)
library(Seurat)

tumor <- subset(data1, subset = interaction_group == "Epithelial_cell" & malignant_status == "malignant" & pat_id=="B")

# 在 tumor 子集内重做标准流程（避免之前整合/全局结构影响轨迹）
DefaultAssay(tumor) <- "RNA"
tumor <- NormalizeData(tumor)
tumor <- FindVariableFeatures(tumor, nfeatures = 3000)
tumor <- ScaleData(tumor)
tumor <- RunPCA(tumor)
tumor <- RunUMAP(tumor, dims = 1:30)
tumor <- FindNeighbors(tumor, dims = 1:30)
tumor <- FindClusters(tumor, resolution = 0.8)
tumor@meta.data$pseudotime=colData(cds)$pseudotime_manual
tumor@meta.data$CAF_module=colData(cds)$ECMCAF_Response

# ========= 你需要按自己对象修改的列名 =========
PT_COL      <- "pseudotime"
COPYKAT_COL <- "copykat_pred"      # aneuploid/diploid
CUSER_COL   <- "Final_CuScore"
CAF_COL     <- "CAF_module"   # 你的CAF模块分数（可选）
# ============================================
seu=tumor
stopifnot(PT_COL %in% colnames(seu@meta.data))
stopifnot(COPYKAT_COL %in% colnames(seu@meta.data))

seu_sub <- subset(seu, subset =  copykat_pred == "aneuploid")

# 可选：只保留表达较好的基因（后续 pseudobulk 更稳）
DefaultAssay(seu_sub) <- "RNA"
seu_sub <- NormalizeData(seu_sub, verbose = FALSE)
seu_sub <- FindVariableFeatures(seu_sub, selection.method = "vst", nfeatures = 4000, verbose = FALSE)
#pseudotime 分箱（30 bins）并做 pseudobulk（平均表达）
set.seed(1)
NBINS <- 30

md <- seu_sub@meta.data %>%
  mutate(
    pt = .data[[PT_COL]]
  ) %>%
  filter(!is.na(pt)) %>%
  arrange(pt) %>%
  mutate(
    pt_bin = ntile(pt, NBINS),
    pt_bin = factor(pt_bin, levels = as.character(1:NBINS))
  )

seu_sub@meta.data$pt_bin <- md$pt_bin[match(rownames(seu_sub@meta.data), rownames(md))]

table(seu_sub$pt_bin, useNA = "ifany")
expr <- GetAssayData(seu_sub, slot = "data")  # log-normalized

# 分组均值：genes x bins
bins <- levels(seu_sub$pt_bin)

pb_mat <- sapply(bins, function(b) {
  cells <- WhichCells(seu_sub, expression = pt_bin == b)
  if (length(cells) == 0) return(rep(NA_real_, nrow(expr)))
  Matrix::rowMeans(expr[, cells, drop = FALSE])
})

pb_mat <- as.matrix(pb_mat)
rownames(pb_mat) <- rownames(expr)
colnames(pb_mat) <- paste0("bin_", bins)

dim(pb_mat)
#准备每个 bin 的元信息（pt 均值、CuScore 均值、CAF 均值）
bin_meta <- seu_sub@meta.data %>%
  as_tibble(rownames = "cell") %>%
  filter(!is.na(pt_bin)) %>%
  group_by(pt_bin) %>%
  summarise(
    n_cells = dplyr::n(),
    pt_mean = mean(.data[[PT_COL]], na.rm = TRUE),
    CuScore_mean = if (CUSER_COL %in% colnames(seu_sub@meta.data)) mean(.data[[CUSER_COL]], na.rm = TRUE) else NA_real_,
    CAF_mean = if (CAF_COL %in% colnames(seu_sub@meta.data)) mean(.data[[CAF_COL]], na.rm = TRUE) else NA_real_
  ) %>%
  ungroup() %>%
  mutate(bin = paste0("bin_", as.character(pt_bin)))

bin_meta
#PROGENy 通路活性（基于 pseudobulk）
ORG <- "Human"  # 或 "Mouse"
TOP <- 500      # progeny 推荐常用 500
# progeny 要求输入：基因(行) × 样本(列)
pro_act <- progeny(
  pb_mat,
  scale = TRUE,
  organism = ORG,
  top = TOP
)
pro_act1 <- progeny(
  exp1,
  scale = TRUE,
  organism = ORG,
  top = TOP
)

# pro_act: pathways x bins
pro_act[1:5, 1:5]
#---------------------------相关分析
# pro_act: bins x pathways （你打印出来就是这样）
stopifnot(all(rownames(pro_act) %in% bin_meta$bin))

pro_df <- as.data.frame(pro_act) %>%
  tibble::rownames_to_column("bin") %>%
  dplyr::left_join(
    bin_meta %>% dplyr::select(bin, pt_mean, CuScore_mean, CAF_mean, n_cells),
    by = "bin"
  ) %>%
  tidyr::pivot_longer(
    cols = -c(bin, pt_mean, CuScore_mean, CAF_mean, n_cells),
    names_to = "pathway",
    values_to = "activity"
  )

str(pro_df)

pro_stats <- pro_df %>%
  group_by(pathway) %>%
  summarise(
    rho_pt = suppressWarnings(cor(activity, pt_mean, method = "spearman", use = "complete.obs")),
    p_pt   = suppressWarnings(cor.test(activity, pt_mean, method = "spearman")$p.value),
    rho_Cu = suppressWarnings(cor(activity, CuScore_mean, method = "spearman", use = "complete.obs")),
    p_Cu   = if (all(is.na(CuScore_mean))) NA_real_ else suppressWarnings(cor.test(activity, CuScore_mean, method = "spearman")$p.value),
    rho_CAF= suppressWarnings(cor(activity, CAF_mean, method = "spearman", use = "complete.obs")),
    p_CAF  = if (all(is.na(CAF_mean))) NA_real_ else suppressWarnings(cor.test(activity, CAF_mean, method = "spearman")$p.value)
  ) %>%
  mutate(
    fdr_pt = p.adjust(p_pt, method = "BH"),
    fdr_Cu = p.adjust(p_Cu, method = "BH"),
    fdr_CAF= p.adjust(p_CAF, method = "BH")
  ) %>%
  arrange(fdr_pt)

pro_stats %>% head(10)

#DoRothEA TF 活性
# 选择 TF regulons（按置信度过滤）
if (ORG == "Human") {
  regulon <- dorothea_hs
} else {
  regulon <- dorothea_mm
}
# 常用只取 A/B/C 级（更可靠）
regulon <- regulon %>%
  filter(confidence %in% c("A", "B", "C")) %>%
  dplyr::select(tf, target, mor, confidence) %>%
  distinct()

regulon %>% head()
# decoupleR 通常期望：features(行) x samples(列)
# run_ulm 输出长表：sample, statistic, source(tf)
tf_ulm <- run_ulm(
  mat = pb_mat,
  network = regulon,
  .source = "tf",
  .target = "target",
  .mor = "mor",
)

tf_ulm %>% head()
#整理为 TF × bin 的活性矩阵，并计算与 pseudotime / CuScore / CAF 的相关
tf_df <- tf_ulm %>%
  select(bin = condition, tf = source, activity = score) %>%
  left_join(bin_meta %>% select(bin, pt_mean, CuScore_mean, CAF_mean, n_cells), by = "bin")

tf_stats <- tf_df %>%
  group_by(tf) %>%
  summarise(
    rho_pt = suppressWarnings(cor(activity, pt_mean, method = "spearman", use = "complete.obs")),
    p_pt   = suppressWarnings(cor.test(activity, pt_mean, method = "spearman")$p.value),
    rho_Cu = suppressWarnings(cor(activity, CuScore_mean, method = "spearman", use = "complete.obs")),
    p_Cu   = if (all(is.na(CuScore_mean))) NA_real_ else suppressWarnings(cor.test(activity, CuScore_mean, method = "spearman")$p.value),
    rho_CAF= suppressWarnings(cor(activity, CAF_mean, method = "spearman", use = "complete.obs")),
    p_CAF  = if (all(is.na(CAF_mean))) NA_real_ else suppressWarnings(cor.test(activity, CAF_mean, method = "spearman")$p.value)
  ) %>%
  mutate(
    fdr_pt = p.adjust(p_pt, method = "BH"),
    fdr_Cu = p.adjust(p_Cu, method = "BH"),
    fdr_CAF= p.adjust(p_CAF, method = "BH")
  ) %>%
  arrange(fdr_pt)

tf_stats %>% head(10)

#------------------------------------------plot
#作图：沿 pseudotime 的通路/TF 趋势
#PROGENy：Top pathways
topP <- pro_stats %>% slice_head(n = 6) %>% pull(pathway)
# 过滤数据以包含所需的通路
df_plot <- pro_df %>%
  filter(pathway %in% topP) %>%
  mutate(pathway = factor(pathway, levels = topP))

# 创建图形
p <- ggplot(df_plot, aes(x = pt_mean, y = activity)) +
  geom_point(color = "black", alpha = 0.8, size = 1.5) +  # 固定点大小
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.8, color = "blue") +
  facet_wrap(~ pathway, scales = "free_y") +  # 按通路分面
  theme_bw() +
  labs(x = "Mean pseudotime (bin)", y = "PROGENy activity (scaled)")

# 显示图形
print(p)


p <- ggplot(df_plot, aes(x = pt_mean, y = activity, color = pathway)) +
  geom_point(alpha = 0.6, size = 1.6) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.9) +  # 关键：se = FALSE
  theme_bw() +
  labs(
    x = "Mean pseudotime (bin)",
    y = "PROGENy activity (scaled)",
    color = "Pathway"
  ) +
  theme(
    legend.title = element_text(face = "bold"),
    legend.key.height = grid::unit(0.5, "lines")
  )

p

#DoRothEA：Top TFs
topTF <- tf_stats %>% slice_head(n = 6) %>% pull(tf)
ggplot(
  tf_df %>% filter(tf %in% topTF),
  aes(x = pt_mean, y = activity)
) +
  geom_point(aes(size = n_cells), alpha = 0.8) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.8) +
  facet_wrap(~ tf, scales = "free_y") +
  theme_bw() +
  labs(x = "Mean pseudotime (bin)", y = "DoRothEA TF activity (ULM statistic)", size = "Cells/bin")
#-----------------------合并在一张图
df_top <- tf_df %>%
  filter(tf %in% topTF) %>%
  group_by(tf) %>%
  mutate(activity_z = as.numeric(scale(activity))) %>%
  ungroup() %>%
  mutate(tf = factor(tf, levels = topTF))

ggplot(df_top, aes(pt_mean, activity_z, color = tf, group = tf)) +
  geom_point(alpha = 0.35, size = 1.2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.9) +
  theme_bw() +
  labs(x = "Mean pseudotime (bin)", y = "TF activity (z-score)", color = "TF")

# CAF module / CuScore 也用分箱趋势
bin_meta_long <- bin_meta %>%
  select(bin, pt_mean, CuScore_mean, CAF_mean, n_cells) %>%
  pivot_longer(cols = c(CuScore_mean, CAF_mean), names_to = "metric", values_to = "value")

ggplot(bin_meta_long, aes(x = pt_mean, y = value)) +
  geom_point(aes(size = n_cells), alpha = 0.8) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.8) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_bw() +
  labs(x = "Mean pseudotime (bin)", y = "Mean score (bin)")

#------------------
#boxplot展示primary和PM的pseudotime差异
# 确保 df2 的 loc 列为因子
df2$loc <- factor(df2$loc, levels = c("primary", "metastasis"))

# 创建 boxplot
p <- ggplot(df2, aes(x = loc, y = pt, fill = loc)) +
  geom_boxplot(outlier.shape = NA, width = 0.3) +  # 去掉箱形图的离群点
  labs(x = "Location", y = "Pseudotime") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(face = "bold", size = 16),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14)
  )

# 显示图形
print(p)

# 保存图形
ggsave("pseudotime_boxplot.pdf", p, width = 8, height = 6)
ggsave("pseudotime_boxplot.png", p, width = 8, height = 6, dpi = 300)






write.csv(pro_stats, "./CUPM/results/V2_202507/6_CuScore_upstream/PROGENy_pseudotime_stats.csv", row.names = FALSE)
write.csv(tf_stats,  "./CUPM/results/V2_202507/6_CuScore_upstream/DoRothEA_TF_pseudotime_stats.csv", row.names = FALSE)

# 也保存长表，方便画更复杂的图
write.csv(pro_df, "./CUPM/results/V2_202507/6_CuScore_upstream/PROGENy_activity_long.csv", row.names = FALSE)
write.csv(tf_df,  "./CUPM/results/V2_202507/6_CuScore_upstream/DoRothEA_activity_long.csv", row.names = FALSE)



