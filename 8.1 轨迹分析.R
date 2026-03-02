load("./CUPM/results/V2_202507/6_CuScore_upstream/scdata_copykat.RData")

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

# 转 monocle3 cds
cds <- as.cell_data_set(tumor)
reducedDims(cds)$UMAP <- Embeddings(tumor, "umap")
# 关键 meta
colData(cds)$pat_id   <- tumor$pat_id
colData(cds)$loc      <- tumor$loc
colData(cds)$CuScore  <- tumor$Final_CuScore
colData(cds)$cluster  <- factor(tumor$seurat_clusters)

cds <- preprocess_cds(cds, num_dim = 50)
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE)

# Root：优先用 B-primary 的 malignant 作为起点
root_cells <- colnames(cds)[colData(cds)$pat_id == "B" & colData(cds)$loc == "primary"]
if (length(root_cells) > 50) {
  cds <- order_cells(cds, root_cells = root_cells)
} else {
  warning("Too few B-primary malignant cells for robust rooting; ordering without explicit root.")
  cds <- order_cells(cds)
}
#--------------------------------
library(dplyr)
reduc <- "UMAP"
df <- data.frame(
  cell = colnames(cds),
  pt   = monocle3::pseudotime(cds),
  loc  = colData(cds)$loc,
  CuScore = colData(cds)$Final_CuScore  # 你的 CuScore 列名若不同改这里
) %>%
  filter(!is.na(pt), !is.na(loc), !is.na(CuScore))

colData(cds)$pseudotime_manual <- NA_real_
colData(cds)$pseudotime_manual[match(res1$cell, colnames(cds))] <- res1$pt_manual
save(cds,file="./CUPM/results/V2_202507/6_CuScore_upstream/cds_pseudotime.RData")
#5A：轨迹图（按 loc / pseudotime 着色）
p1 <- plot_cells(cds, color_cells_by = "loc", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
p2 <- plot_cells(cds, color_cells_by = "pseudotime_manual", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
p1; p2
#-------------------------
df2 <- data.frame(
  cell = colnames(cds),
  pt = colData(cds)$pseudotime_manual,
  pt2=pseudotime(cds),
  loc=colData(cds)$loc,
  CuScore=colData(cds)$Final_CuScore,
  S = colData(cds)$S.Score,
  G2M = colData(cds)$G2M.Score,
  pct_mt = colData(cds)$percent.m
)
cor(df2$pt, df2$CuScore, use="complete.obs", method="spearman")
cor(df2$pt, df2$G2M, use="complete.obs", method="spearman")
cor(df2$pt, df2$pct_mt, use="complete.obs", method="spearman")
cor(df2$pt, df2$S, use="complete.obs", method="spearman")
#Figure 5B---------多因素校正相关性
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

d <- df2 %>%
  select(pt, CuScore, G2M, S, pct_mt) %>%
  drop_na() %>%
  mutate(
    pt_res = resid(lm(pt ~ G2M + S + pct_mt, data = .)),
    cu_res = resid(lm(CuScore ~ G2M + S + pct_mt, data = .))
  )

rho_adj <- cor(d$pt_res, d$cu_res, method = "spearman")

ggplot(d, aes(x = pt_res, y = cu_res, color = pt_res)) +
  geom_point(alpha = 0.55, size = 1.2) +
  geom_smooth(method = "loess", se = TRUE, color = "grey10", linewidth = 1) +
  scale_color_viridis_c(option = "magma", end = 0.95) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  ) +
  labs(
    x = "pt residual (adjusted for G2M + S + pct_mt)",
    y = "CuScore residual (adjusted for G2M + S + pct_mt)",
    color = "pt residual",
    title = sprintf("Adjusted relationship (residual-residual), Spearman rho = %.3f", rho_adj)
  )

#--------------------
#6) Figure 5D/5E：把你已做的 NicheNet targets 变成“CAF-response score”，并接到轨迹上
load("./CUPM/results/V2_202507/5_CuScore_Stromal/ligand_target/full_axis_df.RData")
top_ligands=c("COL1A1", "COL1A2", "CTHRC1", "FN1", "VCAN", "DCN", "ITGB1", "GNAS", "APP", "NEO1")
axis_top_ligand2 <- full_axis_df2 %>%
  filter(ligand %in% top_ligands) %>%
  group_by(ligand) %>%
  arrange(desc(composite_score), .by_group = TRUE)
axis_top_ligand2 <- axis_top_ligand2 %>%
  slice_head(n = 50) %>%   # 每个 ligand 8 条
  ungroup()
#Compute axis activity score (targets → axis)
# keep expressed targets
caf_targets=unique(full_axis_df2$target)
DefaultAssay(tumor)="SCT"
genes_use <- intersect(caf_targets, rownames(tumor[["SCT"]]))
setdiff(caf_targets, genes_use)  # 看看有没有没匹配到的（你说都能匹配到，这里应为空）
expr <- GetAssayData(tumor, assay = "SCT", slot = "data")[genes_use, , drop = FALSE]
dim(expr)  # genes × cells
expr_df <- as.data.frame(t(as.matrix(expr)))
expr_df$cell <- rownames(expr_df)
stopifnot(!anyDuplicated(expr_df$cell))
head(expr_df[, c("cell", genes_use[1])])
# 先看重叠情况
length(intersect(df2$cell, expr_df$cell))
length(setdiff(df2$cell, expr_df$cell))  # df2里有但Seurat里没有的细胞
length(setdiff(expr_df$cell, df2$cell))  # Seurat里有但df2里没有的细胞
# 左连接：以 df2 为主
library(dplyr)
df2 <- df2 %>%
  left_join(expr_df, by = "cell")
#------------------------------AddModuleScore
library(Seurat)
# 确保 seu 的细胞顺序包含 df2 的 cells
seu <- subset(tumor, cells = df2$cell)
seu <- AddModuleScore(
  object = seu,
  features = list(CAF = genes_use),
  assay = "RNA",
  slot  = "data",
  name  = "CAF_AddModule"
)
df2$CAF_AddModule <- FetchData(seu, vars = "CAF_AddModule1")[,1]
cor(df2$pt,df2$CAF_AddModule, method = "spearman")
cor(df2$CuScore,df2$CAF_AddModule, method = "spearman")
#--------------------nested ANOVA / 回归
m0 <- lm(CuScore ~ pt + G2M + S + pct_mt, data = df2)
m1 <- lm(CuScore ~ pt + CAF_AddModule + G2M + S + pct_mt, data = df2)

anova(m0, m1)      # CAF_AddModule 是否显著提升拟合
summary(m1)        # CAF_AddModule 系数方向是否为负
m2 <- lm(CuScore ~ pt * CAF_AddModule + G2M + S + pct_mt, data = df2)
anova(m1, m2)
summary(m2)  # 看 pt:CAF_AddModule
#---------------------------------------把单细胞噪声压下去：pt 分箱 pseudobulk
nbin <- 30
df_bin <- df2 %>%
  mutate(pt_bin = ntile(pt, nbin)) %>%
  group_by(pt_bin) %>%
  summarise(
    pt_mean  = mean(pt, na.rm=TRUE),
    caf_mean = mean(CAF_AddModule, na.rm=TRUE),
    cu_mean  = mean(CuScore, na.rm=TRUE),
    n = n(),
    .groups = "drop"
  ) %>% arrange(pt_bin)

cor.test(df_bin$pt_mean,  df_bin$caf_mean, method="spearman")
cor.test(df_bin$caf_mean, df_bin$cu_mean,  method="spearman")
cor.test(df_bin$pt_mean,  df_bin$cu_mean,  method="spearman")
df_bin <- df_bin %>% mutate(cu_next = lead(cu_mean, 1))
summary(lm(cu_next ~ caf_mean + pt_mean, data = df_bin))
#------------------------plot
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

#----------------------------
# 已知：df2 包含 pt, CuScore, CAF_AddModule, G2M, S, pct_mt
# 已知：df_bin 是你算好的 30-bin pseudobulk（pt_mean/caf_mean/cu_mean/n）
#----------------------------

# 1) 统计量（用你已给出的结果，避免重复计算差异）
rho_pt_caf  <- 0.2204646
rho_cu_caf  <- -0.1568068

# 线性模型（你已拟合）
m1 <- lm(CuScore ~ pt + CAF_AddModule + G2M + S + pct_mt, data = df2)
m2 <- lm(CuScore ~ pt * CAF_AddModule + G2M + S + pct_mt, data = df2)

# 交互项信息（从 summary(m2) 提取）
sm2 <- summary(m2)$coefficients
beta_int <- sm2["pt:CAF_AddModule", "Estimate"]
p_int    <- sm2["pt:CAF_AddModule", "Pr(>|t|)"]

# bin-level spearman（你输出的）
rho_bin_pt_caf <- 0.9337041
p_bin_pt_caf   <- 4.925e-08
rho_bin_caf_cu <- -0.7392659
p_bin_caf_cu   <- 6.355e-06
rho_bin_pt_cu  <- -0.8634038
p_bin_pt_cu    <- 5.831e-07

#----------------------------
# Panel A: single-cell pt vs CAF score (density-aware)
#----------------------------
pA <- ggplot(df2, aes(x = pt, y = CAF_AddModule)) +
  geom_point(size = 0.6, alpha = 0.12) +
  geom_smooth(method = "loess", se = TRUE, span = 0.8,
              color = "black", linewidth = 0.9) +
  labs(
    title = "A  Single-cell association",
    x = "Monocle3 pseudotime (pt)",
    y = "CAF module score (AddModuleScore)"
  ) +
  annotate(
    "text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.2, size = 3.6,
    label = paste0("Spearman ρ = ", round(rho_pt_caf, 3))
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

#----------------------------
# Panel B: single-cell CAF score vs CuScore
#----------------------------
pB <- ggplot(df2, aes(x = CAF_AddModule, y = CuScore)) +
  geom_point(size = 0.6, alpha = 0.12) +
  geom_smooth(method = "loess", se = TRUE, span = 0.8,
              color = "black", linewidth = 0.9) +
  labs(
    title = "B  Single-cell association",
    x = "CAF module score (AddModuleScore)",
    y = "CuScore"
  ) +
  annotate(
    "text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.2, size = 3.6,
    label = paste0("Spearman ρ = ", round(rho_cu_caf, 3))
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

#----------------------------
# Panel C: pseudobulk trend (pt bin means) — CAF vs pt
#----------------------------
pC <- ggplot(df_bin, aes(x = pt_mean, y = caf_mean)) +
  geom_line(linewidth = 1.0, color = "#1f77b4") +
  geom_point(aes(size = n), alpha = 0.85, color = "#1f77b4") +
  scale_size_continuous(range = c(1.5, 6), guide = "none") +
  labs(
    title = "C  Pseudobulk along pseudotime (30 bins)",
    x = "Mean pseudotime per bin",
    y = "Mean CAF module score per bin"
  ) +
  annotate(
    "text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.2, size = 3.6,
    label = paste0("Spearman ρ = ", round(rho_bin_pt_caf, 3),
                   "\nP = ", format.pval(p_bin_pt_caf, digits = 2, eps = 1e-300))
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

#----------------------------
# Panel D: pseudobulk CAF vs CuScore (bin means)
#----------------------------
pD <- ggplot(df_bin, aes(x = caf_mean, y = cu_mean)) +
  geom_path(linewidth = 1.0, color = "#d62728") +
  geom_point(aes(size = n), alpha = 0.85, color = "#d62728") +
  scale_size_continuous(range = c(1.5, 6), guide = "none") +
  labs(
    title = "D  Pseudobulk coupling",
    x = "Mean CAF module score per bin",
    y = "Mean CuScore per bin"
  ) +
  annotate(
    "text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.2, size = 3.6,
    label = paste0("Spearman ρ = ", round(rho_bin_caf_cu, 3),
                   "\nP = ", format.pval(p_bin_caf_cu, digits = 2, eps = 1e-300))
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

#----------------------------
# Panel E: interaction visualization
# 画出在不同 CAF_AddModule 水平下，CuScore 随 pt 的预测斜率差异
#----------------------------
caf_q <- quantile(df2$CAF_AddModule, probs = c(0.2, 0.5, 0.8), na.rm = TRUE)
pt_grid <- seq(min(df2$pt, na.rm=TRUE), max(df2$pt, na.rm=TRUE), length.out = 200)

# 其他协变量固定为中位数（常见投稿做法：展示条件效应）
newdat <- expand.grid(
  pt = pt_grid,
  CAF_AddModule = as.numeric(caf_q),
  G2M = median(df2$G2M, na.rm = TRUE),
  S = median(df2$S, na.rm = TRUE),
  pct_mt = median(df2$pct_mt, na.rm = TRUE)
)

pred <- predict(m2, newdata = newdat, se.fit = TRUE)
newdat$fit <- pred$fit
newdat$se  <- pred$se.fit
newdat$CAF_level <- factor(
  newdat$CAF_AddModule,
  levels = as.numeric(caf_q),
  labels = c("CAF low (20th)", "CAF mid (50th)", "CAF high (80th)")
)

pE <- ggplot(newdat, aes(x = pt, y = fit, color = CAF_level, fill = CAF_level)) +
  geom_line(linewidth = 1.05) +
  geom_ribbon(aes(ymin = fit - 1.96*se, ymax = fit + 1.96*se),
              alpha = 0.15, linewidth = 0) +
  labs(
    title = "E  Stage-dependent effect (interaction model)",
    x = "Monocle3 pseudotime (pt)",
    y = "Predicted CuScore (adjusted)",
    color = "CAF level",
    fill  = "CAF level"
  ) +
  annotate(
    "text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.2, size = 3.6,
    label = paste0("Interaction: β(pt×CAF) = ", round(beta_int, 3),
                   "\nP = ", format.pval(p_int, digits = 2, eps = 1e-300))
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

#----------------------------
# 组合成一个投稿友好的多面板图
#----------------------------
fig <- (pA | pB) / (pC | pD) / pE +
  plot_layout(heights = c(1, 1, 1.1))

fig

# 保存（建议用矢量 PDF）
#ggsave("./CUPM//FigX_CAF_pseudotime_CuScore.pdf", fig, width = 11, height = 12, units = "in")
















