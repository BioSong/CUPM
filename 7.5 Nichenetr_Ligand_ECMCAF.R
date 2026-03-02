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
gr_network            <- readRDS(file.path(nichenetr_dir, "gr_network_human_21122021.rds"))
signaling_network     <- readRDS(file.path(nichenetr_dir, "signaling_network_human_21122021.rds"))
ligand_tf_matrix <- readRDS("./CUPM/data/Nichenetr/ligand_tf_matrix_nsga2r_final.rds")
dim(ligand_target_matrix)
str(lr_network)
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
#===4️⃣ 预测 ligand 活性
ligand_activities <- predict_ligand_activities(
  geneset = target_genes,
  background_expressed_genes = background_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = ligands
)
summary(ligand_activities$aupr)
#===5️⃣ 筛选 top ligands 并可视化
top_ligands <- ligand_activities %>%
  arrange(desc(aupr)) %>%
  slice(1:30)
# Barplot
library(ggplot2)
ggplot(top_ligands, aes(x = reorder(test_ligand, aupr), y = aupr)) +
  geom_col(fill="#E31A1C") +
  coord_flip() +
  labs(x="Ligands (ECM-CAF)", y="Regulatory activity",
       title="ECM-CAF → Epithelial signaling (NicheNet)") +
  theme_minimal(base_size = 14)
#-----------------------------------------始终为NA，放弃
#===6️⃣ 获取 ligand-target links 并绘制热图
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
#----------------------
df1=ligand_targets_df[ligand_targets_df$weight>0.05,]
topN <- 20  # 20 或 30 都合理
lt_top <- df1 %>%
  group_by(ligand) %>%
  arrange(desc(weight), .by_group = TRUE) %>%
  slice_head(n = topN) %>%
  ungroup()
ligands_use <- unique(lt_top$ligand)
targets_use <- unique(lt_top$target)
length(ligands_use)
length(targets_use)
heatmap_mat <- ligand_targets_df %>%
  filter(
    ligand %in% ligands_use,
    target %in% targets_use
  ) %>%
  select(ligand, target, weight) %>%
  pivot_wider(
    names_from  = ligand,
    values_from = weight,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("target") %>%
  as.matrix()
heatmap_mat=t(heatmap_mat)
make_heatmap_ggplot(heatmap_mat, "Prioritized ECM-CAF ligands", "Targets exprssed by epithelial cells",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

#--------------------------
receptors <- lr_network %>% pull(to) %>% unique()
expressed_receptors <- intersect(receptors,background_genes)
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  ligands_use, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  ligands_use,
  order_hclust = "both") 
p=make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Prioritized ECM-CAF ligands", x_name = "Receptors expressed by epithelial cells",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential")
print(p)

#-----------------------------------------------
#主线 1：Ligand–Receptor–Target 机制轴
## =========================
## Downstream line 1: Target → TF
## =========================

# 1. 提取 target → TF 权重矩阵
target_tf_mat <- ligand_tf_matrix[targets_use, , drop = FALSE]

# 2. 去掉全 0 的 TF
target_tf_mat <- target_tf_mat[, colSums(target_tf_mat) > 0, drop = FALSE]

# 3. 计算每个 TF 的总体调控强度
tf_activity_score <- colSums(target_tf_mat)

tf_activity_df <- data.frame(
  TF = names(tf_activity_score),
  score = tf_activity_score
) %>%
  arrange(desc(score))

# 4. 取 top TF
top_tf_n <- 30
tf_use <- tf_activity_df$TF[1:top_tf_n]

# 5. TF–target heatmap
tf_target_heatmap_mat <- target_tf_mat[, tf_use, drop = FALSE]
tf_target_heatmap_mat <- t(tf_target_heatmap_mat)

p=make_heatmap_ggplot(
  tf_target_heatmap_mat,
  y_name = "Predicted TFs",
  x_name = "Epithelial target genes",
  color  = "darkorange",
  legend_title = "TF–target regulatory potential"
)
p
## =========================
## Downstream line 2: Target → Functional enrichment
## =========================
enrich_tar <- enrichGO(
  gene          = targets_use,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL"
)

# 明确手动过滤（避免“阈值没生效”的误解）
enrich_tar <- enrich_tar@result %>%
  dplyr::filter(p.adjust < 0.05)

head(enrich_tar[, c("Description", "p.adjust")])
res <- enrich_tar@result
library(dplyr)

res2 <- res %>%
  mutate(
    module = case_when(
      Description %in% c(
        "cadherin binding",
        "cadherin binding involved in cell-cell adhesion",
        "cell adhesion mediator activity",
        "cell-cell adhesion mediator activity",
        "integrin binding",
        "alpha-catenin binding",
        "gamma-catenin binding"
      ) ~ "Cell–cell / Cell–ECM adhesion",

      Description %in% c(
        "extracellular matrix binding",
        "proteoglycan binding",
        "heparan sulfate proteoglycan binding",
        "glycosaminoglycan binding"
      ) ~ "ECM interaction",

      Description %in% c(
        "protein tyrosine kinase activator activity",
        "neuropilin binding",
        "virus receptor activity",
        "exogenous protein binding"
      ) ~ "Receptor / signal binding",

      Description %in% c(
        "DNA-binding transcription factor binding",
        "RNA polymerase II-specific DNA-binding transcription factor binding",
        "STAT family protein binding",
        "nuclear receptor binding",
        "nuclear glucocorticoid receptor binding",
        "nuclear retinoid X receptor binding"
      ) ~ "TF / nuclear regulation",

      Description %in% c(
        "actin binding",
        "actin filament binding"
      ) ~ "Cytoskeleton remodeling",

      TRUE ~ "Other"
    )
  )
module_order <- c(
  "Cell–cell / Cell–ECM adhesion",
  "ECM interaction",
  "Receptor / signal binding",
  "TF / nuclear regulation",
  "Cytoskeleton remodeling"
)

res2 <- res2 %>%
  filter(module %in% module_order) %>%
  mutate(
    module = factor(module, levels = module_order),
    Description = factor(Description, levels = rev(Description))
  )
library(ggplot2)

p_main <- ggplot(res2,
                 aes(x = -log10(p.adjust),
                     y = Description,
                     color = module)) +
  geom_point(size = 3) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text.y = element_text(angle = 0, face = "bold"),
    axis.title.y = element_blank()
  ) +
  xlab("-log10(adj. P value)") +
  scale_color_brewer(palette = "Set2")

p_main
#-----------------------------------
#-----绘制ligands-target-tf的桑基图
ligand_list <- colnames(vis_ligand_receptor_network[,colMeans(vis_ligand_receptor_network)>0.065])
target_list1 <- colnames(heatmap_mat[,colMeans(heatmap_mat)>0.005])
target_list2 <- rownames(target_tf_mat[rowMeans(target_tf_mat)>1.0e-06,],)
target_list=intersect(target_list1,target_list2)
# 1. 提取 ligand-target axis
ligand_target_df <- heatmap_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ligand") %>%
  filter(ligand %in% ligand_list) %>%
  tidyr::pivot_longer(-ligand, names_to = "target", values_to = "ligand_target_score") %>%
  filter(ligand_target_score > 0)  # 可根据阈值筛选重要 axis

# 2. 提取 TF-target axis，只保留 ligand-target 中的 targets
tf_target_df <- tf_target_heatmap_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("TF") %>%
  tidyr::pivot_longer(-TF, names_to = "target", values_to = "tf_target_score") %>%
  filter(target %in% ligand_target_df$target) %>%
  filter(tf_target_score > 0)

# 3. 合并成 ligand-target-TF 数据框，用于 alluvial plot
alluvial_df <- dplyr::inner_join(
  ligand_target_df,
  tf_target_df,
  by = "target",
  relationship = "many-to-many"  # 明确表明多对多关系
)
# 选择需要的列
alluvial_df <- dplyr::select(alluvial_df, ligand, target, TF, ligand_target_score, tf_target_score)

# 将 score 标准化到 [0,1] 用于颜色/线宽
alluvial_df <- alluvial_df %>%
  mutate(
    ligand_target_norm = ligand_target_score / max(ligand_target_score),
    tf_target_norm = tf_target_score / max(tf_target_score)
  )

# 转换为 alluvial plot 需要的格式
alluvial_data <- alluvial_df %>%
  dplyr::select(ligand, target, TF) %>%
  mutate(freq = 1)  # 用 freq 表示每条轴流量，如果要用 score 可替换

# 绘制 alluvial plot
ggplot(alluvial_data,
       aes(axis1 = ligand, axis2 = target, axis3 = TF,
           y = freq)) +
  geom_alluvium(aes(fill = ligand), width = 0.2) +
  geom_stratum(width = 0.2, fill = "grey80", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Ligand", "Target", "TF"), expand = c(0.1, 0.05)) +
  theme_minimal() +
  ggtitle("Ligand-Target-TF regulatory axes") +
  theme(axis.title = element_blank(),
        panel.grid = element_blank())




axis_df <- alluvial_df %>%
  mutate(axis_score = ligand_target_score * tf_target_score) %>%
  arrange(desc(axis_score))

head(axis_df, 30)
print(n=30)

