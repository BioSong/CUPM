library(Seurat)
library(dplyr)
load("./CUPM/results/V2_202507/FunctionalAnalysis/Seurat_meta_cluster.RData")
seurat_meta$Cu_group=paste("k",seurat_meta$RNA_snn_res.0.1,sep="")
# 提取表达矩阵
exprMat <- as.matrix(GetAssayData(seurat_meta, assay = "SCT", slot = "data"))
library(Matrix)
# 至少在 1% 细胞中表达的基因
min.cells <- ncol(exprMat) * 0.01  
exprMat_filtered <- exprMat[rowSums(exprMat > 0) >= min.cells, ]
dim(exprMat_filtered)
# 提取元数据
meta <- seurat_meta@meta.data
exprMat_high <- exprMat_filtered[, meta$Cu_group == "k0"]
exprMat_low  <- exprMat_filtered[, meta$Cu_group == "k1"]
library(GENIE3)
# 运行 GENIE3，推荐先对表达矩阵过滤（高变基因）
weightMatrix_high <- GENIE3(exprMat_high)
weightMatrix_low  <- GENIE3(exprMat_low)
library(RcisTarget)
# 加载数据库
dbFile <- "/Users/kai/SCI/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
motifRankings <- importRankings(dbFile)
# TF motif 富集
regulons_high <- cisTarget(weightMatrix_high, motifRankings)
regulons_low  <- cisTarget(weightMatrix_low, motifRankings)
