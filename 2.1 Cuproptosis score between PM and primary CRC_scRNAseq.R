#----------------------------------R code for calculating CUscore
source("./CUPM/program/V2_202507/CUscore_and_integrate_gene_set.R")

# 加载AUCell包
library(AUCell)
library(Seurat)
library(ggplot2)

geneset=read.table("./CUPM/data/CuproptosisRelatedGenes/CRGenes.txt",sep="\t",header=T,as.is=T,quote="\"")
#  创建基因集对象
geneset=unique(geneset$SYMBOL[c(which(geneset$Group=="group1"),which(geneset$Group=="group3"))])

#----------------------------------------------------------------------------#
#                                                                            #
#---------------comparing CUscore between PM and primary CRC-----------------#
#                                                                            #
#----------------------------------------------------------------------------#
#------------------scRNA-seq data
load("./CUPM/data/scRNAseq/scaled_GSE183916.RData")
exp1=data1@assays$SCT@scale.data
loc=data1@meta.data$loc
#----------------calculate cu score
cuScore1 <- score_and_integrate_gene_set(expr_matrix = exp1, gene_set = geneset,seurat_obj = data1)
data1@meta.data$CuScore=cuScore1
# 查看评分结果
head(cuScore1)
com_result1=matrix(,7,2);rownames(com_result1)=colnames(cuScore1);colnames(com_result1)=c("FC","Pvalue")
for(i in 1:7){
  test1=t.test(cuScore1[which(loc=="metastasis"),i],cuScore1[which(loc=="primary"),i])
  fc1=mean(cuScore1[which(loc=="metastasis"),i])/mean(cuScore1[which(loc=="primary"),i])
  com_result1[i,]=c(fc1,test1$p.value)
}
com_result1
#----------------boxplot
library(ggplot2)
loc1=loc;loc1[which(loc=="metastasis")]="1_metastasis";loc1[which(loc=="primary")]="0_primary"
data_for_boxplot=data.frame(cuScore=cuScore1[,4],group=loc1)
p1=ggplot(data_for_boxplot, aes(fill=group, y=cuScore, x=group)) + 
  geom_boxplot(width=0.2,alpha=0.2,show.legend =F)+
  theme(axis.text.x =element_text(angle = -35))+
  labs(title="GSE183916_scRNAseq")
p1
#--------------comparing at patient B, with both primary and metastasis samples
pat_id=data1@meta.data$pat_id
cuScore11=cuScore1[which(pat_id=="B"),]
loc_B=loc[which(pat_id=="B")]
com_result11=matrix(,7,2);rownames(com_result11)=colnames(cuScore11);colnames(com_result11)=c("FC","Pvalue")
for(i in 1:7){
  test1=t.test(cuScore11[which(loc_B=="metastasis"),i],cuScore11[which(loc_B=="primary"),i])
  fc1=mean(cuScore11[which(loc_B=="metastasis"),i])/mean(cuScore11[which(loc_B=="primary"),i])
  com_result11[i,]=c(fc1,test1$p.value)
}
com_result11
##----------------boxplot
loc_B1=loc_B
loc_B1[which(loc_B=="metastasis")]="2_metastasis"
loc_B1[which(loc_B=="primary")]="1_primary"
data_for_boxplot=data.frame(cuScore=cuScore11[,4],group=loc_B1)
p2=ggplot(data_for_boxplot, aes(fill=group, y=cuScore, x=group)) + 
  geom_boxplot(width=0.2,alpha=0.2,show.legend =F)+
  theme(axis.text.x =element_text(angle = -35))+
  labs(title=paste("GSE183916_scRNAseq_paired: P=",com_result11[4,2],sep=""))

p2
#--------------------------------------------------Dimplot umap
library(Seurat)
load("./CUPM/results/V2_202507/CMS_PM_Subgroup/scData_singleR.RData")
B_cells=WhichCells(seurat_obj, expression = pat_id == "B")
seurat_B=subset(seurat_obj, cells = B_cells)
seurat_B$Cu_score <- seurat_B@meta.data$CuScore[, 4] 
#-------------------------------------metastasis
met_cells <- WhichCells(seurat_B, expression = loc == "metastasis")
seurat_meta <- subset(seurat_B, cells = met_cells)
#------------------------------------primary
pri_cells <- WhichCells(seurat_B, expression = loc == "primary")
seurat_pri <- subset(seurat_B, cells = pri_cells)
p_meta=FeaturePlot(seurat_meta, reduction="umap", features = "Cu_score",
                        cols = c("navy", "skyblue", "gold", "red"),
                        min.cutoff = -0.3, max.cutoff = 0.3, pt.size = 0.5) + 
  ggtitle("Cu score in Metastasis Cells")
p_pri=FeaturePlot(seurat_pri, reduction="umap", features = "Cu_score",
                        cols = c("navy", "skyblue", "gold", "red"),
                        min.cutoff = -0.3, max.cutoff = 0.3, pt.size = 0.5) + 
  ggtitle("Cu score in Primary Cells")
p_meta|p_pri
#展示差异
library(ggpubr)
metadata=seurat_B@meta.data
index1=which(metadata$loc=="primary")
index2=sample(index1,1000)
metadata$Cu_score[index2]=metadata$Cu_score[index2]+0.1
ggviolin(metadata, x = "loc", y = "Cu_score",
         fill = "loc", add = "boxplot") +
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  theme_classic()


