library(ggplot2); library(dplyr); library(tidyr)
library(pheatmap); library(cowplot); library(Matrix); library(reshape2)
library(SingleCellExperiment); library(scater); library(slingshot)
library(edgeR); library(dorothea); library(viper);library(tibble)
library(Seurat)

load("./CUPM/results/V2_202507/FunctionalAnalysis/Seurat_meta_cluster.RData")
seurat_meta$Cu_group=paste("K",seurat_meta$RNA_snn_res.0.1,sep="")

# ---------------------------
# ② cluster 的 cell-type 组成（barplot + fraction heatmap + stats）
# compute counts and fractions
comp_tab <- seurat_meta@meta.data %>%
  group_by(Cu_group, celltype) %>%
  summarise(n = n()) %>%
  group_by(Cu_group) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()


# stacked barplot (fraction)
p_bar <- ggplot(comp_tab, aes(x = Cu_group, y = frac, fill = celltype)) +
  geom_bar(stat="identity", position="stack") +
  theme_minimal() +
  labs(x="Cluster (Cu_group)", y="Fraction of cells", fill="Cell type",
       title="Cell type composition per cluster (metastasis subset)") +
  theme(axis.text.x = element_text(angle=45, hjust=1))
p_bar
# heatmap of fractions (celltype x cluster)
mat_frac <- comp_tab %>% 
  select(Cu_group, celltype, frac) %>%
  pivot_wider(
    names_from = Cu_group,
    values_from = frac,
    values_fill = 0
  ) %>%
  column_to_rownames("celltype") %>%
  as.matrix()
pheatmap(mat_frac, cluster_rows = T, cluster_cols = T,show_rownames=T,show_colnames=T)

# statistical test: for each celltype, test enrichment across clusters (chi-square)
chi_list <- lapply(unique(comp_tab$celltype), function(ct) {
  tab <- table(seurat_meta@meta.data$Cu_group, seurat_meta@meta.data$celltype == ct)
  res <- chisq.test(tab)
  data.frame(celltype = ct, p.value = res$p.value, stringsAsFactors = FALSE)
})
chi_df <- do.call(rbind, chi_list)
chi_df$adj.p <- p.adjust(chi_df$p.value, method="BH")
write.csv(chi_df, file.path(outdir, "celltype_cluster_chisq.csv"), row.names = FALSE)