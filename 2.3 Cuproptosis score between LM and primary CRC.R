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
#---------------comparing CUscore between LM and primary CRC-----------------#
#                                                                            #
#----------------------------------------------------------------------------#

#----------------------------------------------------------------------------#
#                                                                            #
#---------------             bulk data                      -----------------#
#                                                                            #
#----------------------------------------------------------------------------#
#------GSE190609------35 primary tumor regions and 59 paired PM from 12 patients
exp5=read.table("./CUPM/data/RNAseq/GSE190609/GSE190609_rpkm_SYMBOL.txt",sep="\t",header=T,row.names=1,as.is=T,quote="\"")
cli5=read.table("./CUPM/data/RNAseq/GSE190609/GSE190609_Clinical.txt",sep="\t",header=T,as.is=T,quote="\"")
cli5=cli5[match(colnames(exp5),cli5[,1]),]
exp5=as.matrix(exp5)
cuScore5 <- score_and_integrate_gene_set(expr_matrix = exp5, gene_set = geneset)
#--------------------------------------
comp_result51=matrix(0,7,4);rownames(comp_result51)=colnames(cuScore5);colnames(comp_result51)=c("LM_FC","LM_Pvalue","LN_FC","LN_Pvalue")
for(i in 1:ncol(cuScore5)){
  test1=t.test(cuScore5[which(cli5$Tissue_type=="liver metastasis"),i],cuScore5[which(cli5$Tissue_type=="primary tumor"),i])
  fc1=mean(cuScore5[which(cli5$Tissue_type=="liver metastasis"),i])/mean(cuScore5[which(cli5$Tissue_type=="primary tumor"),i])
  test2=t.test(cuScore5[which(cli5$Tissue_type=="lymph node metastasis"),i],cuScore5[which(cli5$Tissue_type=="primary tumor"),i])
  fc2=mean(cuScore5[which(cli5$Tissue_type=="lymph node metastasis"),i])/mean(cuScore5[which(cli5$Tissue_type=="primary tumor"),i])
  comp_result51[i,]=c(fc1,test1$p.value,fc2,test2$p.value)
}
comp_result51
#---------------------
cuScore52=c(cuScore5[which(cli5$Tissue_type=="liver metastasis"),7],cuScore5[which(cli5$Tissue_type=="primary tumor"),7])
loc_label51=c(rep("2_LM",length(which(cli5$Tissue_type=="liver metastasis"))),rep("1_primary",length(which(cli5$Tissue_type=="primary tumor"))))
data_for_boxplot=data.frame(cuScore=cuScore52,group=loc_label51)
p51=ggplot(data_for_boxplot, aes(fill=group, y=cuScore, x=group)) + 
  geom_boxplot(width=0.2,alpha=0.2,show.legend =F)+
  theme(axis.text.x =element_text(angle = -35))+
  labs(title=paste("GSE190609_bulk: P=",comp_result51[7,2],sep=""))
p51

cuScore52=rbind(cuScore5[which(cli5$Tissue_type=="liver metastasis"),-6],
                  cuScore5[which(cli5$Tissue_type=="primary tumor"),-6])
group=c(rep("1_LM",length(which(cli5$Tissue_type=="liver metastasis"))),rep("0_primary",length(which(cli5$Tissue_type=="primary tumor"))))
data_for_boxplot=data.frame(Sample=rep(rownames(cuScore52),6),Method=rep(paste(1:6,colnames(cuScore52)),each=nrow(cuScore52)),
                            cuScore=c(as.numeric(as.matrix(cuScore52))),
                            group=rep(group,6))
p1=ggplot(data_for_boxplot, aes(fill=group, y=cuScore, x=Method)) + 
  geom_boxplot(width=0.2,alpha=0.2,show.legend =T)+
  theme(axis.text.x =element_text(angle = -35))+
  facet_wrap(~Method, scales = "free")+
  labs(title="GSE190609_LM")
p1

#----------------------------------------------------------------------------#
#                                                                            #
#---------------                scRNAseq data               -----------------#
#                                                                            #
#----------------------------------------------------------------------------#
library(limma)
load("./CUPM/data/scRNAseq/LM/GSE178318/scaled_GSE178318.RData")
#---seu
sampleid=names(seu@active.ident)
sampleid=strsplit2(sampleid,split="_")
seu@meta.data$pat_id=sampleid[,2]
seu@meta.data$loc=sampleid[,3]
#----------------calculate cu score
seu_crc=subset(x=seu,subset=loc=="CRC")
exp1=seu_crc@assays$RNA@scale.data
cuScore1 <- score_and_integrate_gene_set(expr_matrix = exp1, gene_set = geneset,seurat_obj = seu_crc)
seu_lm=subset(x=seu,subset=loc=="LM")
exp2=seu_lm@assays$RNA@scale.data
cuScore2 <- score_and_integrate_gene_set(expr_matrix = exp2, gene_set = geneset,seurat_obj = seu_lm)
seu@meta.data$CuScore=cuScore1
# 查看评分结果
head(cuScore1)
loc=seu@meta.data$loc
com_result1=matrix(,7,2);rownames(com_result1)=colnames(cuScore1);colnames(com_result1)=c("FC","Pvalue")
for(i in 1:7){
  test1=t.test(cuScore1[,i],cuScore2[,i])
  fc1=mean(cuScore1[,i])/mean(cuScore2[,i])
  com_result1[i,]=c(fc1,test1$p.value)
}
com_result1
#----------------boxplot
library(ggplot2)
loc1=loc;
data_for_boxplot=data.frame(cuScore=cuScore1[,4],group=loc1)
p1=ggplot(data_for_boxplot, aes(fill=group, y=cuScore, x=group)) + 
  geom_boxplot(width=0.2,alpha=0.2,show.legend =F)+
  theme(axis.text.x =element_text(angle = -35))+
  labs(title="GSE178318_scRNAseq_LM")
p1


