
# 加载AUCell包
library(AUCell)
library(Seurat)
library(ggplot2)
#---------scRNAseq annotated with cell type
load("./CUPM/data/scRNAseq/scaled_GSE183916.RData")
exp1=data1@assays$SCT@scale.data
loc=data1@meta.data$loc
cuScore1=data1@meta.data$CuScore
pat_id=data1@meta.data$pat_id
celltype=data1@meta.data$celltype
celltype[celltype=="Endothelial_cell"]="Stromal_cell"
loc[which(loc=="primary")]="1_primary";loc[which(loc=="metastasis")]="2_metastasis"
#-----------------------------每种细胞，比较metastasis and primary
#--------------------------------------
celltype1=c("Epithelial_cell","Myeloid_cell","Stromal_cell","B_cell","T_cell")
tresult=matrix(,5,7);rownames(tresult)=celltype1;colnames(tresult)=colnames(cuScore1)
for(i in 1:5){
        cuScore12=cuScore1[which(celltype==celltype1[i]),]
        loc_1=loc[which(celltype==celltype1[i])]
    for(j in 1:7){
        test1=t.test(cuScore12[which(loc_1=="2_metastasis"),j],cuScore12[which(loc_1=="1_primary"),j])
        tresult[i,j]=test1$p.value
    }
}
tresult
#------------------------------
data_for_boxplot=data.frame(cuScore=cuScore1[,7],group=loc,cellType=celltype)
data_for_boxplot=data_for_boxplot[!is.na(data_for_boxplot$cellType),]
p5=ggplot(data_for_boxplot,aes(fill=group,y=cuScore,x=factor(cellType))) +
  geom_boxplot(width=0.2,alpha=0.2,show.legend =T)+
  facet_wrap(~cellType, scales = "free")+
  theme(axis.text.x =element_text(angle = -35))+
  labs(title="GSE183916_scRNAseq_celltype")
p5
#-----------------------------patient B
cuScore11=cuScore1[which(pat_id=="B"),]
loc_B=loc[which(pat_id=="B")]
celltype_B=celltype[which(pat_id=="B")]
#--------------------------
tresult1=matrix(,6,7);rownames(tresult1)=celltype1;colnames(tresult1)=colnames(cuScore11)
for(i in 1:6){
        cuScore12=cuScore11[which(celltype_B==celltype1[i]),]
        loc_B1=loc_B[which(celltype_B==celltype1[i])]
    for(j in 1:7){
        test1=t.test(cuScore12[which(loc_B1=="2_metastasis"),j],cuScore12[which(loc_B1=="1_primary"),j])
        tresult1[i,j]=test1$p.value
    }
}
tresult1
#------------------------------
data_for_boxplot=data.frame(cuScore=cuScore11[,4],group=loc_B,cellType=celltype_B)
p6=ggplot(data_for_boxplot,aes(fill=group,y=cuScore,x=factor(cellType))) +
  geom_boxplot(width=0.2,alpha=0.2,show.legend =T)+
  facet_wrap(~cellType, scales = "free")+
  theme(axis.text.x =element_text(angle = -35))+
  labs(title="GSE183916_scRNAseq_celltype")
p6
#_----------------------------------展示各种TME Cell在PT和PM之间的分布
library(dplyr)
library(tidyr)
cells_with_B <- rownames(data1@meta.data[data1@meta.data$pat_id == "B", ])
data1_B=subset(data1,cells=cells_with_B)
#-------------------------
DimPlot(data1_B, group.by = "celltype", reduction = "umap", label = TRUE)

# 按样本和细胞类型统计细胞数
cell_counts <- data1_B@meta.data %>%
  group_by(pat_id, loc, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(loc) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# 查看前几行
head(cell_counts)
# 对每种 celltype 构建 2x2 表并做 Fisher 检验
fisher_results <- lapply(unique(cell_counts$celltype), function(ct){
  tmp <- cell_counts %>%
    filter(celltype == ct) %>%
    pivot_wider(names_from = loc, values_from = n, values_fill = 0)

  # 确保两组都存在
  if (all(c("primary", "metastasis") %in% colnames(tmp))) {
    mat <- as.matrix(tmp[, c("primary", "metastasis")])
    # 行数必须 ≥2 才能做 Fisher test
    if (nrow(mat) >= 2) {
      test <- fisher.test(mat)
      data.frame(
        cell_type = ct,
        p_value = test$p.value,
        odds_ratio = test$estimate
      )
    } else {
      data.frame(cell_type = ct, p_value = NA, odds_ratio = NA)
    }
  } else {
    data.frame(cell_type = ct, p_value = NA, odds_ratio = NA)
  }
}) %>% bind_rows()

# 多重检验校正
fisher_results$FDR <- p.adjust(fisher_results$p_value, method = "BH")
fisher_results <- arrange(fisher_results, p_value)

# 查看结果
print(fisher_results)

#----------------------可视化------number
ggplot(cell_counts, aes(x = pat_id, y = n, fill = celltype)) +
  geom_bar(stat = "identity") +
  facet_wrap(~loc, scales = "free_x") +
  theme_bw() +
  labs(y = "Cell Type Number", x = "Sample")

#----------------------可视化------proportion
ggplot(cell_counts, aes(x = pat_id, y = prop, fill = celltype)) +
  geom_bar(stat = "identity") +
  facet_wrap(~loc, scales = "free_x") +
  theme_bw() +
  labs(y = "Cell Type Proportion", x = "Sample")


