library(Seurat)
library(nichenetr)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
# === 3️⃣ 加载 NicheNet 资源 ===
nichenetr_dir <- "./CUPM/data/Nichenetr/"
# 2️⃣ 加载官方矩阵
ligand_target_matrix <- readRDS(file.path(nichenetr_dir, "ligand_target_matrix_nsga2r_final.rds"))
lr_network           <- readRDS(file.path(nichenetr_dir, "lr_network_human_21122021.rds"))
weighted_networks    <- readRDS(file.path(nichenetr_dir, "weighted_networks_nsga2r_final.rds"))
gr_network           <- readRDS(file.path(nichenetr_dir, "gr_network_human_21122021.rds"))
signaling_network    <- readRDS(file.path(nichenetr_dir, "signaling_network_human_21122021.rds"))
ligand_tf_matrix     <- readRDS("./CUPM/data/Nichenetr/ligand_tf_matrix_nsga2r_final.rds")
#--------------------------------------------
load("./CUPM/results/V2_202507/5_CuScore_Stromal/cellchat_caf/scdata_cellchat.RData")
data1$CuScore_Final=data1@meta.data$CuScore$Final_Score
data1=subset(data1, subset = loc == "metastasis")
#---------------------------------------------
# === 1️⃣ 基础对象 ===
sce_sub <- data1
expr_mat <- sce_sub@assays$RNA@data
meta_data <- sce_sub@meta.data

# 背景基因集
background_genes <- intersect(rownames(expr_mat),rownames(ligand_target_matrix))
Idents(sce_sub) <- sce_sub$interaction_group
# === 2️⃣ sender genes（差异基因）===
de_genes_sender <- FindMarkers(
  sce_sub,
  ident.1 = "ECMCAF",
  ident.2 = c("apCAF", "iCAF"),
  only.pos = TRUE,        # 只取上调基因
  min.pct = 0.3,          # 至少 10% 的细胞表达
  logfc.threshold = 0.5, # logFC 阈值
  test.use = "wilcox"
)
sender_genes <- rownames(de_genes_sender)[de_genes_sender$p_val_adj < 0.05]

# === 2️⃣ Receiver target genes（差异基因）===
other_cells <- setdiff(unique(sce_sub$interaction_group), "Epithelial_cell")
de_genes_target <- FindMarkers(sce_sub, 
    ident.1 = "Epithelial_cell", 
    ident.2 = other_cells,
    only.pos = TRUE, min.pct = 0.3,  
    logfc.threshold = 0.5,
    test.use = "wilcox"  )
target_genes <- rownames(de_genes_target)[de_genes_target$p_val_adj < 0.05]
#===3️⃣ 获取 ECM-CAF 表达的 ligand
# 所有 ligand 在模型中的列名

# === 4️⃣ Sender表达的ligands ===
# 确保与 NicheNet 模型行匹配
target_genes <- intersect(target_genes, rownames(ligand_target_matrix))
length(target_genes)
ligands <- intersect(sender_genes,lr_network$from)
length(ligands)
#===获取 ligand-target links
all_ligand_targets=list()
for(i in 1:length(ligands)){
  tmp <- get_weighted_ligand_target_links(
    ligands[i],
    target_genes,
    ligand_target_matrix
  )
  all_ligand_targets[[i]]=tmp
}
ligand_targets_df=bind_rows(all_ligand_targets)
#===获取 ligand-receptors links
receptors <- lr_network %>% pull(to) %>% unique()
expressed_receptors <- intersect(receptors,background_genes)
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 
#===获取 ligand-TF links
ligand_tf_mat <- ligand_tf_matrix[ligands, , drop = FALSE]
#===获取 receptors-TF links
receptor_tf_mat <- ligand_tf_matrix[expressed_receptors, , drop = FALSE]
#===获取 target-TF links
# 1. 提取 target → TF 权重矩阵
target_tf_mat <- ligand_tf_matrix[target_genes, , drop = FALSE]
# 2. 去掉全 0 的 TF
target_tf_mat <- target_tf_mat[, colSums(target_tf_mat) > 0, drop = FALSE]
tf_target_df=t(target_tf_mat)
#===获取 target-CuScore links
expr_mat <- GetAssayData(sce_sub, slot = "data")
cuscore_vec <- sce_sub$CuScore_Final
target_cuscore_cor <- lapply(target_genes, function(gene) {
  gene_exp <- as.numeric(expr_mat[gene, ])
  ct <- suppressWarnings(
    cor.test(
      gene_exp,
      cuscore_vec,
      method = "spearman"
    )
  )
  data.frame(
    target = gene,
    cor = unname(ct$estimate),
    pvalue = ct$p.value
  )
}) %>%
  bind_rows() %>%
  mutate(
    padj = p.adjust(pvalue, method = "BH")
  )

#--------------------------------------三个mat
#----------ligand_target_mat  # 行 = ligand, 列 = target, 值 = regulatory potential
#----------tf_target_df  # 行 = TF | 列 =target | 值 = tf_target_score
#----------target_cuscore_cor
str(ligand_targets_df)
str(tf_target_df)
str(target_cuscore_cor)
#------------------------------------------------
#-------run
library(dplyr)
library(tidyr)
# -------------------------------
# 1. TF → Target (tf_target_long)
# -------------------------------
tf_target_long <- tf_target_df %>%
  as.data.frame() %>%
  tibble::rownames_to_column("TF") %>%
  pivot_longer(-TF, names_to = "target", values_to = "tf_target_score") %>%
  filter(tf_target_score > 0)

# -------------------------------
# 2. Ligand → Target (ligand_target_df)
# -------------------------------
ligand_target_df <- ligand_targets_df %>%
  mutate(ligand_target_score = as.numeric(weight)) %>%
  dplyr::select(-weight) %>%
  filter(ligand_target_score > 0)

# -------------------------------
# 3. Ligand → TF (ligand_tf_long)
# -------------------------------
ligand_tf_long <- ligand_tf_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ligand") %>%
  pivot_longer(-ligand, names_to = "TF", values_to = "ligand_tf_score") %>%
  filter(ligand_tf_score > 0)

# -------------------------------
# 4. Receptor → TF (receptor_tf_long)
# -------------------------------
receptor_tf_long <- receptor_tf_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("receptor") %>%
  pivot_longer(-receptor, names_to = "TF", values_to = "receptor_tf_score") %>%
  filter(receptor_tf_score > 0)

# -------------------------------
# 5. Ligand → Receptor (ligand_receptor_links_df) 已有
# -------------------------------
ligand_receptor_links_df <- ligand_receptor_links_df %>%
  dplyr::rename(ligand = from, receptor = to)
ligand_receptor_links_df <- ligand_receptor_links_df %>%
  mutate(ligand_receptor_score = as.numeric(weight)) %>%  # 转为普通 numeric 列
  dplyr::select(-weight) %>%
  filter(ligand_receptor_score > 0)
# -------------------------------
# 6. 构建完整 axis: Ligand → Receptor → TF → Target
# -------------------------------
full_axis_df <- ligand_receptor_links_df %>%
  inner_join(receptor_tf_long, by = "receptor") %>%
  inner_join(tf_target_long, by = "TF") %>%
  inner_join(ligand_tf_long, by = c("ligand","TF")) %>%
  inner_join(ligand_target_df, by = c("ligand","target")) %>%
  left_join(target_cuscore_cor %>% dplyr::select(target, cor, padj), by = "target") %>%
  mutate(
    # 标准化各 score（可选）
    ligand_target_norm = (ligand_target_score - min(ligand_target_score)) / (max(ligand_target_score) - min(ligand_target_score)),
    ligand_tf_norm = (ligand_tf_score - min(ligand_tf_score)) / (max(ligand_tf_score) - min(ligand_tf_score)),
    receptor_tf_norm = (receptor_tf_score - min(receptor_tf_score)) / (max(receptor_tf_score) - min(receptor_tf_score)),
    # 综合 axis score
    axis_score = ligand_target_norm * ligand_tf_norm * receptor_tf_norm * tf_target_score * abs(cor)
  )
#--------------------------------
library(dplyr)
library(scales)

full_axis_df2 <- full_axis_df %>%
  mutate(
    # 1. 各层取 rank（比直接用数值稳定）
    lr_rank   = percent_rank(ligand_receptor_score),
    ltf_rank  = percent_rank(ligand_tf_score),
    rtf_rank  = percent_rank(receptor_tf_score),
    tf_rank   = percent_rank(tf_target_score),
    cor_rank  = percent_rank(abs(cor)),

    # 2. 综合 score（可调权重）
    composite_score =
      0.30 * lr_rank +
      0.20 * ltf_rank +
      0.20 * rtf_rank +
      0.20 * tf_rank +
      0.10 * cor_rank
  )
summary(full_axis_df2$composite_score)
sum(full_axis_df2$composite_score>0.9)
# -------------------------------
# 7. 筛选 top axes
# -------------------------------
top_axes_net <- full_axis_df2 %>%
  filter(
    composite_score >= quantile(composite_score, 0.99, na.rm = TRUE)
  ) %>%
  arrange(desc(composite_score))
#统计 ligand频率
top_ligands <- top_axes_net %>%
  count(ligand, sort = TRUE)
#选择 top5 ligand
top5_ligands <- top_ligands %>%
  slice_head(n = 10) %>%
  pull(ligand)
#提取对应axis
top_ligands=c("COL1A1", "COL1A2", "CTHRC1", "FN1", "VCAN", "DCN", "ITGB1", "GNAS", "APP", "NEO1")
top_axes <- top_axes_net %>%
  filter(ligand %in% top_ligands)
#---------在每个 ligand 内再排序
axis_top_ligand <- full_axis_df2 %>%
  filter(ligand %in% top_ligands) %>%
  group_by(ligand) %>%
  arrange(desc(composite_score), .by_group = TRUE)
#---------每个 ligand 只保留 Top N
axis_plot_df <- axis_top_ligand %>%
  slice_head(n = 8) %>%   # 每个 ligand 8 条
  ungroup()
#---------------------------------

#---------------------------------
library(ggalluvial)
library(ggplot2)
ggplot(axis_plot_df,
       aes(axis1 = ligand,
           axis2 = receptor,
           axis3 = TF,
           axis4 = target,
           y = composite_score)) +
  geom_alluvium(aes(fill = ligand),
                alpha = 0.7,
                width = 1/12) +
  geom_stratum(width = 1/12,
               color = "grey30",
               fill = "grey90") +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 3) +
  scale_x_discrete(limits = c("Ligand", "Receptor", "TF", "Target"),
                   expand = c(0.05, 0.05)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  )
#------------ligand-target heatmap
top5_ligands <- top_ligands %>%
  slice_head(n = 10) %>%
  pull(ligand)
axis_top_ligand2 <- full_axis_df2 %>%
  filter(ligand %in% top5_ligands) %>%
  group_by(ligand) %>%
  arrange(desc(composite_score), .by_group = TRUE)
axis_top_ligand2 <- axis_top_ligand2 %>%
  slice_head(n = 50) %>%   # 每个 ligand 8 条
  ungroup()

heatmap_mat <- axis_top_ligand2 %>%
  dplyr::select(ligand, target, ligand_target_score) %>%
  tidyr::pivot_wider(
    names_from  = ligand,
    values_from = ligand_target_score,
    values_fn   = max,  # or mean/sum/etc.
    values_fill = list(ligand_target_score = 0)
  ) %>%
  tibble::column_to_rownames("target") %>%
  as.matrix()
heatmap_mat=t(heatmap_mat)
make_heatmap_ggplot(heatmap_mat, "Prioritized ECM-CAF ligands", "Targets exprssed by epithelial cells",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

#---------------------------------
axis_top_ligand2 <- axis_top_ligand2 %>%
  slice_head(n = 5) %>%   # 每个 ligand 8 条
  ungroup()
ggplot(axis_top_ligand2,
       aes(axis1 = ligand,
           axis2 = receptor,
           axis3 = TF,
           axis4 = target,
           y = composite_score)) +
  geom_alluvium(aes(fill = ligand),
                alpha = 0.7,
                width = 1/12) +
  geom_stratum(width = 1/12,
               color = "grey30",
               fill = "grey90") +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 3) +
  scale_x_discrete(limits = c("Ligand", "Receptor", "TF", "Target"),
                   expand = c(0.05, 0.05)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  )
write.table()
#------------------------------------------------------------------------
#         绘图展示Target和Cuscore的相关性以及axis activity与CuScore的相关性
#------------------------------------------------------------------------
meta_data$CuScore_Final=meta_data$CuScore$Final_Score
#Compute axis activity score (targets → axis)
# keep expressed targets
axis_targets <- unique(axis_top_ligand2$target)
# axis score: mean expression of target genes
axis_score <- colMeans(expr_mat[axis_targets, , drop = FALSE])
# add to metadata
meta_data$AxisScore <- axis_score
#---------------------------------------
#Target–CUScore correlation (gene-level)
target_cor_df <- lapply(axis_targets, function(g) {
  ct <- cor.test(expr_mat[g, ], meta_data$CuScore_Final, method = "spearman")
  data.frame(
    gene = g,
    rho  = ct$estimate,
    pval = ct$p.value
  )
}) |> bind_rows()
# adjust p values
target_cor_df$FDR <- p.adjust(target_cor_df$pval, method = "BH")
#---------------------------
ggplot(target_cor_df, aes(x = reorder(gene, rho), y = rho)) +
  geom_col(fill = "#3C6E71") +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Axis target genes",
    y = "Spearman correlation with CUScore",
    title = "Correlation between axis targets and cuproptosis score"
  ) +
  theme_classic()
#-------------------------------------
ggplot(meta_data, aes(x = AxisScore, y = CuScore_Final)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "#9E2A2B") +
  labs(
    x = "Axis activity score",
    y = "CuScore",
    title = "Higher axis activity corresponds to lower cuproptosis-related scores"
  ) +
  theme_classic()
#--------------
