library(tidyverse)
library(reshape2)
library(ggpubr)        # for stat_compare_means()
library(pheatmap)      # for correlation heatmap
library(FactoMineR)    # PCA
library(factoextra)
#-----------------------------load data
load("./CUPM/results/V2_202507/3_CuScore_CellType/GSE190609_Cell_Proportions_v2.RData")
#est_prop,sce,bulk_expr,bulk_cli
prop_mat <- est_prop$Est.prop.weighted
prop_df  <- as.data.frame(prop_mat) %>% rownames_to_column("Sample")
meta=bulk_cli[,c("Sample_ID","Tissue_type")]
colnames(meta)=c("Sample","Location")
#------------------------------------
bulk_cuscore=read.table("./CUPM/results/V2_202507/2_CuScore_Location_Compare/GSE190609_CuScore.txt",header=T,row.names=1,as.is=T,quote="\"")
cupro=bulk_cuscore %>% rownames_to_column("Sample")
colnames(cupro)[ncol(cupro)]="CuproScore"
cupro=cupro[,c("Sample","CuproScore")]
#----------------------------
df <- prop_df %>%
  left_join(meta, by="Sample") %>%
  left_join(cupro, by="Sample")
#1. Compare Cuproptosis score across locations --------------------------------
df=df[which(df$Location=="primary tumor"|df$Location=="peritoneal metastasis"|df$Location=="liver metastasis"),]
ggplot(df, aes(x=Location, y=CuproScore, fill=Location)) +
  geom_boxplot() +
  stat_compare_means(method="kruskal.test", label="p.format") +
  theme_minimal() +
  labs(title="Cuproptosis Score by CRC Location")

#-------------------------------only include primary and peritoneal metastasis
df=df[which(df$Location=="primary tumor"|df$Location=="peritoneal metastasis"),]
df[which(df$Location=="primary tumor"),"Location"]="Primary"
df[which(df$Location=="peritoneal metastasis"),"Location"]="PeritonealMet"
df$Location=factor(df$Location,levels=c("Primary", "PeritonealMet"))
#Merge all into one long dataframe
df_melt <- df %>%
  gather(key="CellType", value="Proportion", -Sample, -Location, -CuproScore)
#--------------------------------Compare TME composition across locations
# 3.1 Boxplots + Kruskal–Wallis test per cell type
p_list <- df_melt %>%
  group_by(CellType) %>%
  do({
    gg = ggplot(data=., aes(x=Location, y=Proportion, fill=Location)) +
      geom_boxplot(outlier.shape=21, outlier.fill="white") +
      stat_compare_means(method="kruskal.test", label="p.format") +
      theme_classic() +
      labs(title=unique(.$CellType), y="Proportion") +
      theme(legend.position="none")
    tibble(plot=list(gg))
  })

# Arrange into a grid
plot_list <- p_list$plot
n_ct <- length(plot_list)
ncol <- 3
do.call(gridExtra::grid.arrange, c(plot_list, ncol=ncol))


# ──────────────────────────────────────────────────────────────────────────────
# 5. Correlation between TME proportions & CuproScore --------------------------
# ──────────────────────────────────────────────────────────────────────────────

# 5.1 Compute Spearman correlations
cors <- df %>%
  select(-Sample, -Location) %>%
  gather(key="Feature", value="Value", -CuproScore) %>%
  group_by(Feature) %>%
  summarize(
    rho = cor(Value, CuproScore, method="spearman"),
    p   = cor.test(Value, CuproScore, method="spearman")$p.value
  ) %>%
  arrange(desc(abs(rho)))

# 5.2 Barplot of correlations
ggplot(cors, aes(x=reorder(Feature, rho), y=rho, fill=rho)) +
  geom_col() +
  coord_flip() +
  geom_text(aes(label=signif(rho,2)), hjust=ifelse(cors$rho>0, -0.1, 1.1)) +
  theme_minimal() +
  labs(title="Spearman Correlation: Cell Proportions vs. CuproScore",
       x="", y="Spearman rho")

# 5.3 Heatmap of full correlation matrix (all cell types + CuproScore)
cor_mat <- df %>%
  select(-Sample, -Location) %>%
  cor(method="spearman")
sort_cell=c("Epithelial cells","T cells","Mast cells","B cells","Myeloids","Stromal cells","CuproScore")
cor_mat=cor_mat[match(sort_cell,rownames(cor_mat)),match(sort_cell,colnames(cor_mat))]
pheatmap(cor_mat, 
         main="Spearman Correlation Heatmap",
         cluster_rows=F, cluster_cols=F)
# ──────────────────────────────────────────────────────────────────────────────
# 6. PCA on TME composition (optional) -----------------------------------------
# ──────────────────────────────────────────────────────────────────────────────

# 6.1 Perform PCA on proportions only
prop_df <- df %>% 
  select(-Sample, -Location, -CuproScore)
prop_df[] <- lapply(prop_df, function(col) as.numeric(as.character(col)))

# 1.3 Turn into a bare data.frame and set rownames
prop_df <- as.data.frame(prop_df)
rownames(prop_df) <- df$Sample

# ──────────────────────────────────────────────────────────────────────────────
# 2. PCA with FactoMineR -------------------------------------------------------
# ──────────────────────────────────────────────────────────────────────────────
prop_df=prop_df[,c(1,2,3,6)]
pca_fm <- PCA(prop_df, graph = FALSE)

# 2.1 Plot PC1 vs PC2 colored by Location
fviz_pca_ind(
  pca_fm,
  geom.ind = "point",
  habillage = df$Location,
  addEllipses = TRUE,
  title = "PCA of TME Composition by Location"
)

# 2.2 Plot PC1 vs PC2 colored by CuproScore
# Extract coordinates and combine with CuproScore
pca_coords <- as.data.frame(pca_fm$ind$coord) %>%
  rownames_to_column("Sample") %>%
  left_join(df %>% select(Sample, CuproScore), by = "Sample")

ggplot(pca_coords, aes(x = Dim.1, y = Dim.2, color = CuproScore)) +
  geom_point(size = 3) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(
    title = "PCA of TME Composition Colored by Cuproptosis Score",
    x = "PC1", y = "PC2"
  )



ggplot(cell_fraction, aes(x = sample_type, y = fraction, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "", y = "Cell fraction (%)", fill = "Cell type") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )
