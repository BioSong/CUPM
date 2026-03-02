library(Seurat)
library(SeuratWrappers)   # as.cell_data_set helper
library(monocle)
library(ggplot2)
library(dplyr)

load("./CUPM/results/V2_202507/FunctionalAnalysis/Seurat_meta_cluster.RData")
seurat_meta$Cu_group=paste("k",seurat_meta$RNA_snn_res.0.1,sep="")

#--------------------Cuproptosis Related Genes
geneset=read.table("./CUPM/data/CuproptosisRelatedGenes/CRGenes.txt",sep="\t",header=T,as.is=T,quote="\"") 
cuGeneSet <- geneset$SYMBOL[c(which(geneset$Group=="group1"),which(geneset$Group=="group3"))]
#Set these genes in the scRNAseq data object, because the next several functions will depend on them.
data2 <- as.CellDataSet(seurat_meta, assay = "SCT", reduction = "umap")
data2 <- estimateSizeFactors(data2)
data2 <- estimateDispersions(data2)
data21 <- setOrderingFilter(data2, cuGeneSet)
plot_ordering_genes(data21)
data21 <- reduceDimension(data21, max_components = 2, method = 'DDRTree')
data21 <- orderCells(data21)
plot_cell_trajectory(data21, color_by = "Cu_group")
