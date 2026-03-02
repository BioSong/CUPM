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

#-----------------bulk data
#----------------------------GSE183202为CRC-PM的bulk数据，GSE183916为scRNA-seq数据
load("./CUPM/data/RNAseq/GSE183202/GSE183202_Exp_Cli.RData")
exp2=log2(fpkm1+1)
#------------------filter out 低表达的gene
index1=rowMeans(exp2<1)
exp2=exp2[which(index1<0.5),]
cli1[cli1==""]=NA
#out_sample=c("S25_1","S25_4b","S28_1b","S28_3b")
#out_index=match(out_sample,cli1$Sample_ID)
#exp2=exp2[,-out_index]
#cli1=cli1[-out_index,]
#--------------计算CUscore
cuScore2 <- score_and_integrate_gene_set(expr_matrix = exp2, gene_set = geneset)
sid1=c("S12_1","S25_2","S28_5b","S39_1","S47_4","S49_2");index1=match(sid1,cli1$Sample_ID)
sid2=c("S12_2","S25_3","S28_2b","S39_3","S47_3","S49_1");index2=match(sid2,cli1$Sample_ID)
cuScore20=cuScore2
cuScore2[index1,]=cuScore20[index2,]
cuScore2[index2,]=cuScore20[index1,]
comp_result2=matrix(0,7,4);rownames(comp_result2)=colnames(cuScore2);colnames(comp_result2)=c("FC","Pvalue","corR","corP")
for(i in 1:ncol(cuScore2)){
  test1=t.test(cuScore2[which(cli1$Tissue_type=="peritoneal metastatic disease"),i],cuScore2[which(cli1$Tissue_type=="primary tumor"),i],alternative="less")
  fc1=mean(cuScore2[which(cli1$Tissue_type=="peritoneal metastatic disease"),i])/mean(cuScore2[which(cli1$Tissue_type=="primary tumor"),i])
  cor1=cor.test(cuScore20[,i],cli1$PCI_score,method="spearman",alternative="less")
  comp_result2[i,]=c(fc1,test1$p.value,cor1$estimate,cor1$p.value)
}
comp_result2
#-------------------绘制cuScore和PCI的相关性点图
scatter_data=data.frame(CUscore=cuScore20[,2],PCIscore=cli1$PCI_score)
ggplot(scatter_data, aes(x=PCIscore, y=CUscore)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE)+
  geom_text(label=c("R=-0.18 \n P=0.049"), x = 35, y = 1.2, check_overlap = T)
#-------------------绘制cuScore和MSN表达的相关性点图
scatter_data2=data.frame(CUscore=cuScore2$Final_Score,MSN=exp2[which(rownames(exp2)=="MSN"),])
ggplot(scatter_data2, aes(x=MSN, y=CUscore)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE)+
  geom_text(label=c("R=-0.48 \n P=1.794e-06"), x = 1, y = 0.25, check_overlap = T)
#-----------------------------------筛选primary的patient
cli2=cli1
pri_pat=cli2[which(cli2$Tissue_type=="primary tumor"),"Patient_ID"]
cli22=cli2[cli2$Patient_ID %in% pri_pat,]
cuScore22=cuScore2[match(cli22$Sample_ID,cli2$Sample_ID),]
comp_result22=matrix(0,7,2);rownames(comp_result22)=colnames(cuScore22);colnames(comp_result22)=c("FC","Pvalue")
for(i in 1:ncol(cuScore22)){
  test1=t.test(cuScore22[which(cli22$Tissue_type=="peritoneal metastatic disease"),i],cuScore22[which(cli22$Tissue_type=="primary tumor"),i],alternative="less")
  fc1=mean(cuScore22[which(cli22$Tissue_type=="peritoneal metastatic disease"),i])/mean(cuScore22[which(cli22$Tissue_type=="primary tumor"),i])
  comp_result22[i,]=c(fc1,test1$p.value)
}
comp_result22
#---------------------
cuScore32=c(cuScore22[which(cli22$Tissue_type=="peritoneal metastatic disease"),7],cuScore22[which(cli22$Tissue_type=="primary tumor"),7])
loc_label32=c(rep("2_metastasis",length(which(cli22$Tissue_type=="peritoneal metastatic disease"))),rep("1_primary",length(which(cli22$Tissue_type=="primary tumor"))))
data_for_boxplot=data.frame(cuScore=cuScore32,group=loc_label32)
p11=ggplot(data_for_boxplot, aes(fill=group, y=cuScore, x=group)) + 
  geom_boxplot(width=0.2,alpha=0.2,show.legend =F)+
  theme(axis.text.x =element_text(angle = -35))+
  labs(title=paste("GSE183202_bulk: P=",comp_result22[7,2],sep=""))
p11

#----------------------------配对的Normal、Primary和PM的数据集GSE225182
cpm1 <- read.table("./CUPM/data/RNAseq/GSE225182/GSE225182_PMCRClog2CPM.txt",sep="\t",header=T,row.names=1,as.is=T,quote="\"")
cpm1=cpm1[-which(rownames(cpm1)=="C4B_2"),]
cli3=read.table("./CUPM/data/RNAseq/GSE225182/GSE225182_Clinical.txt",sep="\t",header=T,as.is=T,quote="\"")
#--------------计算CUscore
cpm1=as.matrix(cpm1)
cuScore3 <- score_and_integrate_gene_set(expr_matrix = cpm1, gene_set = geneset)
comp_result3=matrix(0,7,6);rownames(comp_result3)=colnames(cuScore2);colnames(comp_result3)=c("M_P_FC","M_P_Pvalue","M_N_FC","M_N_Pvalue","P_N_FC","P_N_Pvalue")
for(i in 1:ncol(cuScore3)){
  test1=t.test(cuScore3[which(cli3$Tissue_type=="MT"),i],cuScore3[which(cli3$Tissue_type=="PT"),i],alternative="less")
  fc1=mean(cuScore3[which(cli3$Tissue_type=="MT"),i])/mean(cuScore3[which(cli3$Tissue_type=="PT"),i])
  test2=t.test(cuScore3[which(cli3$Tissue_type=="MT"),i],cuScore3[which(cli3$Tissue_type=="N"),i],alternative="less")
  fc2=mean(cuScore3[which(cli3$Tissue_type=="MT"),i])/mean(cuScore3[which(cli3$Tissue_type=="N"),i])
  test3=t.test(cuScore3[which(cli3$Tissue_type=="PT"),i],cuScore3[which(cli3$Tissue_type=="N"),i])
  fc3=mean(cuScore3[which(cli3$Tissue_type=="PT"),i])/mean(cuScore3[which(cli3$Tissue_type=="N"),i],alternative="less")
  comp_result3[i,]=c(fc1,test1$p.value,fc2,test2$p.value,fc3,test3$p.value)
}
comp_result3
#---------------------
cuScore32=c(cuScore3[which(cli3$Tissue_type=="MT"),5],cuScore3[which(cli3$Tissue_type=="PT"),5])
loc_label32=c(rep("2_metastasis",length(which(cli3$Tissue_type=="MT"))),rep("1_primary",length(which(cli3$Tissue_type=="PT"))))
data_for_boxplot=data.frame(cuScore=cuScore32,group=loc_label32)
p31=ggplot(data_for_boxplot, aes(fill=group, y=cuScore, x=group)) + 
  geom_boxplot(width=0.2,alpha=0.2,show.legend =F)+
  theme(axis.text.x =element_text(angle = -35))+
  labs(title=paste("GSE225182_bulk: P=",comp_result3[5,2],sep=""))
p31

#------GSE190609------35 primary tumor regions and 59 paired PM from 12 patients
exp5=read.table("./CUPM/data/RNAseq/GSE190609/GSE190609_rpkm_SYMBOL.txt",sep="\t",header=T,row.names=1,as.is=T,quote="\"")
cli5=read.table("./CUPM/data/RNAseq/GSE190609/GSE190609_Clinical.txt",sep="\t",header=T,as.is=T,quote="\"")
cli5=cli5[match(colnames(exp5),cli5[,1]),]
exp5=as.matrix(exp5)
cuScore5 <- score_and_integrate_gene_set(expr_matrix = exp5, gene_set = geneset)
comp_result5=matrix(0,7,2);rownames(comp_result5)=colnames(cuScore5);colnames(comp_result5)=c("FC","Pvalue")
for(i in 1:ncol(cuScore5)){
  test1=t.test(cuScore5[which(cli5$Tissue_type=="peritoneal metastasis"),i],cuScore5[which(cli5$Tissue_type=="primary tumor"),i],alternative="less")
  fc1=mean(cuScore5[which(cli5$Tissue_type=="peritoneal metastasis"),i])/mean(cuScore5[which(cli5$Tissue_type=="primary tumor"),i])
  comp_result5[i,]=c(fc1,test1$p.value)
}
comp_result5
#---------------------
cuScore52=c(cuScore5[which(cli5$Tissue_type=="peritoneal metastasis"),7],cuScore5[which(cli5$Tissue_type=="primary tumor"),7])
loc_label51=c(rep("2_metastasis",length(which(cli5$Tissue_type=="peritoneal metastasis"))),rep("1_primary",length(which(cli5$Tissue_type=="primary tumor"))))
data_for_boxplot=data.frame(cuScore=cuScore52,group=loc_label51)
p51=ggplot(data_for_boxplot, aes(fill=group, y=cuScore, x=group)) + 
  geom_boxplot(width=0.2,alpha=0.2,show.legend =F)+
  theme(axis.text.x =element_text(angle = -35))+
  labs(title=paste("GSE190609_bulk: P=",comp_result5[7,2],sep=""))
p51
#------------------------预测CMS样本标签
library(CMSclassifier)
exp50=read.table("./CUPM/data/RNAseq/GSE190609/GSE190609_rpkm_ENTREZID.txt",sep="\t",header=T,row.names=1,as.is=T,quote="\"")
cms_label=classifyCMS(exp50,method="RF")
cli5$CMS_RFpredicted=cms_label$predictedCMS$RF
cli5$CMS_RFnearest=cms_label$nearestCMS$RF
write.table(cli5,"./CUPM/data/RNAseq/GSE190609/GSE190609_Clinical.txt",sep="\t",col.names=T,row.names=F,quote=F)
#-----------------------根据CMS_RFnearest进行分组，比较CMS4和CMS2之间的差异
ttest52=t.test(cuScore5[which(cli5$CMS_RFnearest=="CMS4"),7],cuScore5[which(cli5$CMS_RFnearest=="CMS2"),7])
cuScore53=c(cuScore5[which(cli5$CMS_RFnearest=="CMS4"),7],cuScore5[which(cli5$CMS_RFnearest=="CMS2"),7])
cms_label52=c(rep("CMS4",length(which(cli5$CMS_RFnearest=="CMS4"))),rep("CMS2",length(which(cli5$CMS_RFnearest=="CMS2"))))
data_for_boxplot=data.frame(cuScore=cuScore53,group=cms_label52)
p52=ggplot(data_for_boxplot, aes(fill=group, y=cuScore, x=group)) + 
  geom_boxplot(width=0.2,alpha=0.2,show.legend =F)+
  theme(axis.text.x =element_text(angle = -35))+
  labs(title=paste("GSE190609_bulk_CMS: P=",ttest52$p.value,sep=""))
p52

#----------------------结合Tissue_type和CMS_RF进行分组
t.test(cuScore5[which(cli5$Tissue_type=="peritoneal metastasis"&cli5$CMS_RF=="CMS2"),7],
          cuScore5[which(cli5$Tissue_type=="primary tumor"&cli5$CMS_RF=="CMS2"),7])
t.test(cuScore5[which(cli5$Tissue_type=="peritoneal metastasis"&cli5$CMS_RF=="CMS4"),7],
          cuScore5[which(cli5$Tissue_type=="primary tumor"&cli5$CMS_RF=="CMS4"),7])
#----------------------结合Tissue_type和CMS_RFnearest进行分组
t.test(cuScore5[which(cli5$Tissue_type=="peritoneal metastasis"&cli5$CMS_RFnearest=="CMS2"),7],
       cuScore5[which(cli5$Tissue_type=="primary tumor"&cli5$CMS_RFnearest=="CMS2"),7])
cuScore53=c(cuScore5[which(cli5$Tissue_type=="peritoneal metastasis"&cli5$CMS_RFnearest=="CMS2"),7],
                cuScore5[which(cli5$Tissue_type=="primary tumor"&cli5$CMS_RFnearest=="CMS2"),7])
cms_label53=c(rep("metastasisCMS2",length(which(cli5$Tissue_type=="peritoneal metastasis"&cli5$CMS_RFnearest=="CMS2"))),
                  rep("primaryCMS2",length(which(cli5$Tissue_type=="primary tumor"&cli5$CMS_RFnearest=="CMS2"))))
data_for_boxplot=data.frame(cuScore=cuScore53,group=cms_label53)
p53=ggplot(data_for_boxplot, aes(fill=group, y=cuScore, x=group)) + 
  geom_boxplot(width=0.2,alpha=0.2,show.legend =F)+
  theme(axis.text.x =element_text(angle = -35))+
  labs(title="GSE183916_scRNAseq_CMS2")
p53


t.test(cuScore5[which(cli5$Tissue_type=="peritoneal metastasis"&cli5$CMS_RFnearest=="CMS4"),7],
       cuScore5[which(cli5$Tissue_type=="primary tumor"&cli5$CMS_RFnearest=="CMS4"),7])
cuScore54=c(cuScore5[which(cli5$Tissue_type=="peritoneal metastasis"&cli5$CMS_RFnearest=="CMS4"),7],
            cuScore5[which(cli5$Tissue_type=="primary tumor"&cli5$CMS_RFnearest=="CMS4"),7])
cms_label54=c(rep("metastasisCMS4",length(which(cli5$Tissue_type=="peritoneal metastasis"&cli5$CMS_RFnearest=="CMS4"))),
              rep("primaryCMS4",length(which(cli5$Tissue_type=="primary tumor"&cli5$CMS_RFnearest=="CMS4"))))
data_for_boxplot=data.frame(cuScore=cuScore54,group=cms_label54)
p54=ggplot(data_for_boxplot, aes(fill=group, y=cuScore, x=group)) + 
  geom_boxplot(width=0.2,alpha=0.2,show.legend =F)+
  theme(axis.text.x =element_text(angle = -35))+
  labs(title="GSE183916_scRNAseq_CMS4")
p54


t.test(cuScore51[which(cli5$Tissue_type=="peritoneal metastasis"&cli5$CMS_RFnearest=="CMS2")],
       cuScore51[which(cli5$Tissue_type=="primary tumor"&cli5$CMS_RFnearest=="CMS4")])
#---------------------
cuScore55=c(cuScore5[which(cli5$Tissue_type=="peritoneal metastasis"),7],
            cuScore5[which(cli5$Tissue_type=="liver metastasis"),7],
              cuScore5[which(cli5$Tissue_type=="lymph node metastasis"),7],
                cuScore5[which(cli5$Tissue_type=="primary tumor"),7])
loc_label55=c(rep("3_PM",length(which(cli5$Tissue_type=="peritoneal metastasis"))),
                rep("2_LM",length(which(cli5$Tissue_type=="liver metastasis"))),
                  rep("1_LN",length(which(cli5$Tissue_type=="lymph node metastasis"))),
                    rep("0_primary",length(which(cli5$Tissue_type=="primary tumor"))))
data_for_boxplot=data.frame(cuScore=cuScore55,group=loc_label55)
p55=ggplot(data_for_boxplot, aes(fill=group, y=cuScore, x=group)) + 
  geom_boxplot(width=0.2,alpha=0.2,show.legend =F)+
  theme(axis.text.x =element_text(angle = -35))+
  labs(title="GSE190609_bulk")
p55

#--------------------------------------------------------------

















#-------基因太少，不适用于AUCell，报错
#-----GSE198085-------11 with CRCPM and 10 without PM
#----measured by Nanostring platforms
#exp6=read.table("./CUPM/data/RNAseq/GSE198085/GSE198085_Normlog2.txt",sep="\t",header=T,row.names=1,as.is=T,quote="\"")
#colnames(exp6)=gsub("X","S",colnames(exp6))
#cli6=read.table("./CUPM/data/RNAseq/GSE198085/GSE198085_Clinical.txt",sep="\t",header=T,as.is=T,quote="\"")
#cli6$Sample_ID=paste("S",cli6$Sample_ID,sep="")
#cli6=cli6[match(colnames(exp6),cli6$Sample_ID),]
#cellRankings <- AUCell_buildRankings(exp6)
# Step 2: 计算AUC来评估基因集的富集程度（铜死亡评分）
#cuAUC <- AUCell_calcAUC(cuGeneSet, cellRankings)
# Step 3: 查看AUC分数
#cuScore6 <- getAUC(cuAUC)
#t.test(cuScore6[1,which(cli6$Metastasis=="CRCPM")],cuScore6[1,which(cli6$Metastasis=="noCRCPM")])





