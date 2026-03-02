library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

#------------------------------------
load('./CUPM/results/V2_202507/5_CuScore_Stromal/Seurat_Stromal.RData')
stromal$CuScore_Final=stromal@meta.data$CuScore$Final_Score
FeaturePlot(stromal, features = "CuScore_Final") +
  scale_color_gradientn(colours = c("#1B4F72", "#17A589", "#F5B041", "#A04000")) # 或其他颜色方案
#----------------------
data_for_boxplot=data.frame(CuScore=stromal$CuScore_Final,group=stromal$subtype)
ggplot(data_for_boxplot,aes(fill=group,y=CuScore,x=factor(group))) +
  geom_boxplot(width=0.2,alpha=0.2,show.legend =T)+
  theme(axis.text.x =element_text(angle = -35))+
  labs(title="GSE183916_scRNAseq_StromalCell")
#----------------------------------每种subtype比较primary和PM间的CuScore差异
loc=stromal@meta.data$loc
cuScore1=stromal@meta.data$CuScore_Final
celltype=stromal@meta.data$subtype
loc[which(loc=="primary")]="1_primary";loc[which(loc=="metastasis")]="2_metastasis"
data_for_boxplot=data.frame(cuScore=cuScore1,group=loc,cellType=celltype)
data_for_boxplot=data_for_boxplot[!is.na(data_for_boxplot$cellType),]
ggplot(data_for_boxplot,aes(fill=group,y=cuScore,x=factor(cellType))) +
  geom_boxplot(width=0.2,alpha=0.2,show.legend =T)+
  facet_wrap(~cellType, scales = "free")+
  theme(axis.text.x =element_text(angle = -35))+
  labs(title="GSE183916_scRNAseq_celltype")
#----------------------------------每种cluster比较primary和PM间的CuScore差异
cluster=stromal@meta.data$SCT_snn_res.0.8
loc[which(loc=="primary")]="1_primary";loc[which(loc=="metastasis")]="2_metastasis"
data_for_boxplot=data.frame(cuScore=cuScore1,group=loc,cellType=cluster)
data_for_boxplot=data_for_boxplot[!is.na(data_for_boxplot$cluster),]
ggplot(data_for_boxplot,aes(fill=group,y=cuScore,x=factor(cluster))) +
  geom_boxplot(width=0.2,alpha=0.2,show.legend =T)+
  facet_wrap(~cluster, scales = "free")+
  theme(axis.text.x =element_text(angle = -35))+
  labs(title="GSE183916_scRNAseq_cluster")
#-------------------------------量化每种subtype对总体 CuScore 的贡献
# Step 1: summarize by sample and subtype
tbl <- stromal@meta.data %>%
  group_by(orig.ident, loc, subtype) %>% 
  summarise(
    n = n(),
    mean_cucell = mean(CuScore_Final),
    .groups="drop"
  ) %>% 
  group_by(orig.ident) %>% 
  mutate(prop = n / sum(n)) %>% 
  ungroup()

# Step 2: calculate contributions
tbl <- tbl %>%
  mutate(contribution = prop * mean_cucell)

# Step 3: average contribution for primary vs PM
avg_contrib <- tbl %>%
  group_by(loc, subtype) %>%
  summarise(
    mean_contrib = mean(contribution),
    se_contrib = sd(contribution)/sqrt(n()),
    .groups='drop'
  ) %>%
  mutate(subtype = fct_reorder(subtype, mean_contrib))
#Publication-ready 贡献对比图
ggplot(avg_contrib,
       aes(x = subtype,
           y = mean_contrib,
           fill = loc)) +
  geom_col(position="dodge") +
  geom_errorbar(aes(ymin = mean_contrib-se_contrib,
                    ymax = mean_contrib+se_contrib),
                width = 0.2,
                position = position_dodge(.9)) +
  labs(
    x = "Stromal Subtype",
    y = "CuScore Contribution",
    fill = "Tumor Type"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

