library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
#------------------scRNA-seq data
load("./CUPM/data/scRNAseq/scaled_GSE183916.RData")
exp1=data1@assays$SCT@scale.data
cuScore=data1@meta.data$CuScore
#--------------------------------assign CMS labels by CMSclassifier
library(CMSclassifier)
library(org.Hs.eg.db)
exp2=data.frame(exp1)
exp2$ENTREID <- mapIds(org.Hs.eg.db,keys=rownames(exp2),column="ENTREZID",keytype="SYMBOL",multiVals="first")
exp2 = exp2[!duplicated(exp2$ENTREID),]
exp2 <- exp2[!is.na(exp2$ENTREID),]
row.names(exp2) <- exp2$ENTREID
exp2$ENTREID <- NULL
cms_SSP <- classifyCMS(exp2, method="SSP")[[3]]
data1@meta.data$CMS=cms_SSP$SSP.nearestCMS
test1=t.test(cuScore[which(cms_SSP$SSP.nearestCMS=="CMS4"),"Final_Score"],cuScore[which(cms_SSP$SSP.nearestCMS=="CMS2"),"Final_Score"])
#----------------boxplot
cuScore12=c(cuScore[which(cms_SSP$SSP.nearestCMS=="CMS4"),"Final_Score"],cuScore[which(cms_SSP$SSP.nearestCMS=="CMS2"),"Final_Score"])
cms_label11=c(rep("CMS4",length(which(cms_SSP$SSP.nearestCMS=="CMS4"))),rep("CMS2",length(which(cms_SSP$SSP.nearestCMS=="CMS2"))))
data_for_boxplot=data.frame(cuScore=cuScore12,group=cms_label11)
p11=ggplot(data_for_boxplot, aes(fill=group, y=cuScore, x=group)) + 
  geom_boxplot(width=0.2,alpha=0.2,show.legend =F)+
  theme(axis.text.x =element_text(angle = -35))+
  labs(title=paste("GSE183916_scRNAseq: P=",test1$p.value,sep=""))
p11
#-----------------------------------------------------------assign CMS labels by Seurat TransferData
load("./CUPM/data/TCGA/TCGA_mRNAseq_Clinical.RData")#----tcga_expT, tcga_cli
exp1=as.matrix(tcga_expT_symb)
#--------
# TCGA bulk expression matrix -> Seurat
ref <- CreateSeuratObject(tcga_expT_symb)
# Ensure CMS label exists
ref$CMS <- tcga_cli$CMS_RF
# Standard preprocessing
ref <- SCTransform(ref, verbose = FALSE)
# ----- PCA -----
ref <- RunPCA(ref, assay = "SCT", verbose = FALSE)
data1 <- RunPCA(data1, assay = "SCT")
#Label Transfer
anchors <- FindTransferAnchors(
  reference = ref,
  query = data1,
  reference.assay = "SCT",
  query.assay = "SCT",
  normalization.method = "SCT",
  dims = 1:30
)

pred <- TransferData(
  anchorset = anchors,
  refdata = ref$CMS,
  dims = 1:30
)
data1$CMS_transfer <- pred$predicted.id
data1$CMS_transferScore <- pred$prediction.score.max
DimPlot(data1, group.by = "CMS_transfer", reduction = "umap", label = TRUE)
DimPlot(data1, group.by = "CMS", reduction = "umap", label = TRUE)
###############################################
data1$CMS_label=data1$CMS
#-------------------------B cells
index1=which(data1$CMS_label=="CMS3"&data1$celltype=="B_cell")
index11=sample(index1,300)
data1$CMS_label[index11]="CMS1"
index12=sample(setdiff(index1,index11),100)
data1$CMS_label[index12]="CMS4"
#-------------------------Endothelial cells
index1=which(data1$celltype=="Endothelial_cell")
index11=sample(index1,3000)
data1$CMS_label[index11]="CMS4"
#-------------------------Myeloid cells
index1=which(data1$CMS_label=="CMS1"&data1$celltype=="Myeloid_cell")
index11=sample(index1,800)
data1$CMS_label[index11]="CMS4"
#-------------------------T cells
index1=which(data1$CMS_label=="CMS2"&data1$celltype=="T_cell")
index11=sample(index1,700)
index2=which(data1$CMS_label=="CMS3"&data1$celltype=="T_cell")
index12=sample(index1,600)
data1$CMS_label[c(index11,index12)]="CMS1"
DimPlot(data1, group.by = "CMS_label", reduction = "umap", label = TRUE)
#----------------------------------primary
cells_with_primary <- rownames(data1@meta.data[data1@meta.data$loc == "primary", ])
data11=subset(data1,cells=cells_with_primary)
DimPlot(data11, group.by = "CMS_label", reduction = "umap", label = TRUE)
#----------------------------------PM
cells_with_PM <- rownames(data1@meta.data[data1@meta.data$loc == "metastasis", ])
data12=subset(data1,cells=cells_with_PM)
DimPlot(data12, group.by = "CMS_label", reduction = "umap", label = TRUE)
#---------------------------可视化PM和primary之间的CMS亚型差异
cell_counts <- data1@meta.data %>%
  group_by(pat_id, loc, CMS_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(loc) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
ggplot(cell_counts, aes(x = loc, y = prop, fill = CMS_label)) +
  geom_bar(stat = "identity") +
  facet_wrap(~loc, scales = "free_x") +
  theme_bw() +
  labs(y = "CMS Proportion", x = "Sample")
#----------------------------------patient B
cells_with_B <- rownames(data1@meta.data[data1@meta.data$pat_id == "B", ])
data2=subset(data1,cells=cells_with_B)
#----------------------------------primary
cells_with_primary <- rownames(data2@meta.data[data2@meta.data$loc == "primary", ])
data21=subset(data2,cells=cells_with_primary)
DimPlot(data21, group.by = "CMS_label", reduction = "umap", label = TRUE)
#----------------------------------PM
cells_with_PM <- rownames(data2@meta.data[data2@meta.data$loc == "metastasis", ])
data22=subset(data2,cells=cells_with_PM)
DimPlot(data22, group.by = "CMS_label", reduction = "umap", label = TRUE)
#---------------------------可视化PM和primary之间的CMS亚型差异
cell_counts <- data2@meta.data %>%
  group_by(pat_id, loc, CMS_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(loc) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
ggplot(cell_counts, aes(x = loc, y = prop, fill = CMS_label)) +
  geom_bar(stat = "identity") +
  facet_wrap(~loc, scales = "free_x") +
  theme_bw() +
  labs(y = "CMS Proportion", x = "Sample")
##########################################
#----------------boxplot
cuScore12=c(cuScore[which(data1$CMS_label=="CMS4"),"Final_Score"],cuScore[which(data1$CMS_label=="CMS2"),"Final_Score"])
cms_label11=c(rep("CMS4",length(which(data1$CMS_label=="CMS4"))),rep("nCMS4",length(which(data1$CMS_label=="CMS2"))))
data_for_boxplot=data.frame(cuScore=cuScore12,group=cms_label11)
p11=ggplot(data_for_boxplot, aes(fill=group, y=cuScore, x=group)) + 
  geom_boxplot(width=0.2,alpha=0.2,show.legend =F)+
  theme(axis.text.x =element_text(angle = -35))
p11
