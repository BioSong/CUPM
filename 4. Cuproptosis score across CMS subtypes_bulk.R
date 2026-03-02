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
#---------------comparing CUscore between CMS4 and CMS2----------------------#
#                                                                            #
#----------------------------------------------------------------------------#

#-------------------------------------------------------------TCGA
load("./CUPM/data/TCGA/TCGA_mRNAseq_Clinical.RData")
exp1=as.matrix(tcga_expT_symb)

cuScore1 <- score_and_integrate_gene_set(expr_matrix = exp1, gene_set = geneset)
# 查看评分结果
head(cuScore1)
#--------------------------------------
comp_result11=matrix(NA,7,2);rownames(comp_result11)=colnames(cuScore1);colnames(comp_result11)=c("FC","Pvalue")
for(i in 1:ncol(cuScore1)){
  test1=t.test(cuScore1[which(tcga_cli$CMS_RF=="CMS4"),i],cuScore1[which(tcga_cli$CMS_RF=="CMS2"),i])
  fc1=mean(cuScore1[which(tcga_cli$CMS_RF=="CMS4"),i])/mean(cuScore1[which(tcga_cli$CMS_RF=="CMS2"),i])
  comp_result11[i,]=c(fc1,test1$p.value)
}
comp_result11
#----------------boxplot
cuScore12=c(cuScore1[which(tcga_cli$CMS_RF=="CMS4"),7],cuScore1[which(tcga_cli$CMS_RF=="CMS2"),7])
cms_label11=c(rep("CMS4",length(which(tcga_cli$CMS_RF=="CMS4"))),rep("CMS2",length(which(tcga_cli$CMS_RF=="CMS2"))))
data_for_boxplot=data.frame(cuScore=cuScore12,group=cms_label11)
p11=ggplot(data_for_boxplot, aes(fill=group, y=cuScore, x=group)) + 
  geom_boxplot(width=0.2,alpha=0.2,show.legend =F)+
  theme(axis.text.x =element_text(angle = -35))+
  labs(title=paste("TCGA_RNAseq: P=",comp_result11[7,2],sep=""))
p11


#---------------------------------------ARGO
load("./CUPM/data/ARGO.RData")
exp2=as.matrix(argo_exp_symbol)
cuScore2 <- score_and_integrate_gene_set(expr_matrix = exp2, gene_set = geneset)
# 查看评分结果
head(cuScore2)
#--------------------------------------
comp_result21=matrix(NA,7,2);rownames(comp_result21)=colnames(cuScore2);colnames(comp_result21)=c("FC","Pvalue")
for(i in 1:ncol(cuScore2)){
  test1=t.test(cuScore2[which(argo_cli$RF.predictedCMS=="CMS4"),i],cuScore2[which(argo_cli$RF.predictedCMS=="CMS2"),i])
  fc1=mean(cuScore2[which(argo_cli$RF.predictedCMS=="CMS4"),i])/mean(cuScore2[which(argo_cli$RF.predictedCMS=="CMS2"),i])
  comp_result21[i,]=c(fc1,test1$p.value)
}
comp_result21
#----------------boxplot
cuScore22=c(cuScore2[which(argo_cli$RF.predictedCMS=="CMS4"),7],cuScore2[which(argo_cli$RF.predictedCMS=="CMS2"),7])
cms_label21=c(rep("CMS4",length(which(argo_cli$RF.predictedCMS=="CMS4"))),rep("CMS2",length(which(argo_cli$RF.predictedCMS=="CMS2"))))
data_for_boxplot=data.frame(cuScore=cuScore22,group=cms_label21)
p21=ggplot(data_for_boxplot, aes(fill=group, y=cuScore, x=group)) + 
  geom_boxplot(width=0.2,alpha=0.2,show.legend =F)+
  theme(axis.text.x =element_text(angle = -35))+
  labs(title=paste("ARGO_RNAseq: P=",comp_result21[7,2],sep=""))
p21
#----------------------------------------------



