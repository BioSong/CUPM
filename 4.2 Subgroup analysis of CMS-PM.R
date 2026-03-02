#----------------------------------R code for calculating CUscore
source("./CUPM/program/V2_202507/CUscore_and_integrate_gene_set.R")
# 加载AUCell包
library(AUCell)
library(Seurat)
library(ggplot2)
library(ggpubr)

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
#----------------------------------------------------------
cli11=cli1[!is.na(cli1$groups),]
cuScore21=cuScore2[!is.na(cli1$groups),]
data_for_boxplot=data.frame(cuScore=cuScore21$Final_Score,group=cli11$groups)

kruskal.test(cuScore ~ group, data = data_for_boxplot)
# 如果显著（p < 0.05），进行两两比较（Dunn's test）
library(FSA)
dunnTest(cuScore ~ group, data = data_for_boxplot, method = "bh")
#----------------boxplot
ggplot(data_for_boxplot, aes(x = group, y = cuScore, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  stat_compare_means(method = "kruskal.test", label.y = max(data_for_boxplot$cuScore) * 1.05) +
  theme_minimal() +
  labs(title = "Copper Death Score Across Groups",
       x = "Group", y = "Copper Score")
#--------------------------------------------------------------heatmap between high and low cuScore group
library(ComplexHeatmap)
library(circlize)
library(dplyr)
df=cli11
copper_score=cuScore21$Final_Score
# 将铜死亡 score 分组
df <- df %>%
  mutate(copper_group = ifelse(copper_score >= median(copper_score), "High", "Low"))
#------------------
# 构造 annotation 行，确保行为样本名，列为各变量
# 举例：挑选出以下列
anno_df <- df %>%
  select(APC_mutation, RAS_mutation, TP53_mutation, SMAD4_mutation, PI3K_mutation, groups,
         Tumour_type,HIPEC_chemo,gender, pcsynmeta, tumor_content, Age, PCI_score, copper_group)
rownames(anno_df)=df$Sample_ID
#--------------------
# 颜色映射
col_list <- list(
  copper_group = c("High" = "#E41A1C", "Low" = "#377EB8"),
  gender = c("male" = "#4DAF4A", "female" = "#984EA3"),
  groups = c("k1" = "#66C2A5", "k2" = "#FC8D62", "k3" = "#8DA0CB"),
  pcsynmeta = c("metachronous" = "#FFD92F", "synchronous" = "#A6D854"),
  Tumour_type = c("mucinous_adenocarcinoma" = "#FFD92F", "adenocarcinoma" = "#A6D854", "signet_cell" = "#8DA0CB"),
  HIPEC_chemo = c("mitomycin_c" = "#FFD92F", "none" = "#A6D854"),
  APC_mutation = c("mt" = "#E41A1C", "wt" = "gray90"),
  RAS_mutation = c("mt" = "#E41A1C", "wt" = "gray90"),
  TP53_mutation = c("mt" = "#E41A1C", "wt" = "gray90"),
  SMAD4_mutation = c("mt" = "#E41A1C", "wt" = "gray90"),
  PI3K_mutation = c("mt" = "#E41A1C", "wt" = "gray90"),
  PCI_score = colorRamp2(c(min(df$PCI_score), max(df$PCI_score)), c("white", "black")),
  tumor_content = colorRamp2(c(min(df$tumor_content), max(df$tumor_content)), c("white", "darkgreen")),
  Age = colorRamp2(c(min(df$Age), max(df$Age)), c("white", "purple"))
)
#--------------------------------
