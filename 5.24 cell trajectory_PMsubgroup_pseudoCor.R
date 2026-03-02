######################################################################
library(Seurat)
library(SingleR)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(matrixStats)
library(mclust)
library(cowplot)
library(ggplot2)
library(dplyr)
library(scrapper)

#--------------------------bulk data
load("./CUPM/data/RNAseq/GSE183202/GSE183202_Exp_Cli.RData")
exp1=log2(fpkm1+1)
#------------------filter out 低表达的gene
index1=rowMeans(exp1<1)
exp1=exp1[which(index1<0.5),]
cli1[cli1==""]=NA
ref_exp=exp1[,!is.na(cli1$groups)]
ref_label=cli1[!is.na(cli1$groups),"groups"]

######################################################################
#------------------scRNA-seq data
load("./CUPM/data/scRNAseq/scaled_GSE183916.RData")
sct_exp <- GetAssayData(data1, assay = "SCT", layer = "data")
common_genes <- intersect(rownames(sct_exp), rownames(ref_exp))
test_mat <- sct_exp[common_genes, , drop=FALSE]         # genes x cells
ref_mat  <- ref_exp[common_genes, , drop=FALSE]    # genes x samples
bulk_labels <- ref_label
seurat_obj <- data1


