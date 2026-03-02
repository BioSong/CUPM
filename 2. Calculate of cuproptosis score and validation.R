library(AUCell)
library(GSVA)
library(UCell)
library(singscore)
library(Seurat)
library(ggplot2)


geneset=read.table("./CUPM/data/CuproptosisRelatedGenes/CRGenes.txt",sep="\t",header=T,as.is=T,quote="\"")
geneset=unique(geneset$SYMBOL[c(which(geneset$Group=="group1"),which(geneset$Group=="group3"))])

#---------------------------validate in treated cells
##############################################################################
#                                                                            #
#-------------------------------GSE248084------------------------------------#
#                                                                            #
##############################################################################
#---------------------------case组是经Wnt-pathway inhibitors LF3处理24小时后的癌干细胞
#--------------------------------------------------------------
#1--------------------------比较case和control组之间的score差异
vdata1=read.table("./CUPM//data/cu_cell/GSE248084_All_sample_TPM.txt",sep="\t",header=T,as.is=T,quote="\"")
genes=unique(vdata1[,2]);vdata1=vdata1[match(genes,vdata1[,2]),]
rownames(vdata1)=genes;vdata1=as.matrix(vdata1[,c(-1,-2)])
#-----------------------------------------------------------------------
source("./CUPM/program/V2_202507/CUscore_and_integrate_gene_set.R")
result_df <- score_and_integrate_gene_set(expr_matrix = vdata1, gene_set = geneset)
#---------------------------------可视化
library(ggplot2)
library(tidyverse)
result_df=cbind(Sample=colnames(vdata1),result_df)
write.table(result_df,"./CUPM/results/V2_202507/GSE248084_CUscore.txt",sep="\t",col.names=T,row.names=F,quote=F)
#-------------绘制boxplot,比较两组之间的scores差异
result_df=read.table("./CUPM/results/V2_202507/GSE248084_CUscore.txt",sep="\t",header=T,as.is=T,quote="\"")
# Stacked + percent
result_df$group=c(rep("control",3),rep("LF3",3))
data_for_boxplot=data.frame(Sample=rep(result_df$Sample,6),Method=rep(paste(1:6,colnames(result_df)[c(2:6,13)]),each=6),
                            cuScore=c(as.numeric(as.matrix(result_df[,c(2:6,13)]))),
                            group=rep(result_df$group,6))
p1=ggplot(data_for_boxplot, aes(fill=group, y=cuScore, x=Method)) + 
  geom_boxplot(width=0.2,alpha=0.2,show.legend =F)+
  theme(axis.text.x =element_text(angle = -35))+
  facet_wrap(~Method, scales = "free")+
  labs(title="GSE248084")
p1
#---------------------AUC roc plot
par(mfrow=c(2,3))
library(pROC)
group1=as.factor(c(rep(0,3),rep(1,3)))
roc_score1=roc(response=group1,predictor=result_df$AUCell)
plot(roc_score1,
     print.auc = TRUE,
     auc.polygon = TRUE,
     print.thres = "best",
     print.thres.col = "red",
     grid = TRUE,
     legacy.axes = TRUE,
     main="AUCell_GSE248084")
roc_score2=roc(response=group1,predictor=result_df$ssGSEA)
plot(roc_score2,
     print.auc = TRUE,
     auc.polygon = TRUE,
     print.thres = "best",
     print.thres.col = "red",
     grid = TRUE,
     legacy.axes = TRUE,
     main="ssGSEA_GSE248084")
roc_score3=roc(response=group1,predictor=result_df$UCell)
plot(roc_score3,
     print.auc = TRUE,
     auc.polygon = TRUE,
     print.thres = "best",
     print.thres.col = "red",
     grid = TRUE,
     legacy.axes = TRUE,
     main="UCell_GSE248084")
roc_score4=roc(response=group1,predictor=result_df$singscore)
plot(roc_score4,
     print.auc = TRUE,
     auc.polygon = TRUE,
     print.thres = "best",
     print.thres.col = "red",
     grid = TRUE,
     legacy.axes = TRUE,
     main="singscore_GSE248084")
roc_score5=roc(response=group1,predictor=result_df$Seurat)
plot(roc_score5,
     print.auc = TRUE,
     auc.polygon = TRUE,
     print.thres = "best",
     print.thres.col = "red",
     grid = TRUE,
     legacy.axes = TRUE,
     main="Seurat_GSE248084")
roc_score6=roc(response=group1,predictor=result_df$Final_Score)
plot(roc_score6,
     print.auc = TRUE,
     auc.polygon = TRUE,
     print.thres = "best",
     print.thres.col = "red",
     grid = TRUE,
     legacy.axes = TRUE,
     main="Final_Score_GSE248084")
#----------------------------------------
#计算不同算法之间的相关性，绘制热图
result_df=read.table("./CUPM/results/V2_202507/CuScore_Calculation_Valitation/GSE248084_CUscore.txt",sep="\t",header=T,as.is=T,quote="\"")
cuscore_methods=result_df[,2:6];rownames(cuscore_methods)=result_df$Sample
