library(Seurat)
library(monocle)
#---------scRNAseq annotated with cell type
#load("./CUPM/results/1_UmapAnalysis/scaled_GSE183916.RData")#---------data1
load("./CUPM/data/scRNAseq/scaled_GSE183916.RData")
#--------------------Cuproptosis Related Genes
geneset=read.table("./CUPM/data/CuproptosisRelatedGenes/CRGenes.txt",sep="\t",header=T,as.is=T,quote="\"") 
cuGeneSet <- geneset$SYMBOL[c(which(geneset$Group=="group1"),which(geneset$Group=="group3"))]
#Set these genes in the scRNAseq data object, because the next several functions will depend on them.
data2 <- as.CellDataSet(data1, assay = "RNA", reduction = "umap")
data2 <- estimateSizeFactors(data2)
data2 <- estimateDispersions(data2)
data21 <- setOrderingFilter(data2, cuGeneSet)
plot_ordering_genes(data21)
data21 <- reduceDimension(data21, max_components = 2,
    method = 'DDRTree')
data21 <- orderCells(data21)
plot_cell_trajectory(data21, color_by = "loc")
plot_cell_trajectory(data21, color_by = "celltype")
plot_cell_trajectory(data21, color_by = "CMS")
plot_cell_trajectory(data21, color_by = "CuScore$Final_Score")+
  scale_color_viridis_c(option = "inferno")
#----------------------------
pData(data21)$celltype_highlight <- ifelse(
  pData(data21)$celltype %in% c("Epithelial_cell", "Stromal_cell"),
  pData(data21)$celltype,
  "Other"
)
plot_cell_trajectory(data21, color_by = "celltype_highlight") +
  scale_color_manual(
    values = c(
      "Epithelial_cell" = "#1f77b4",     # 蓝色
      "Stromal_cell" = "#ff7f0e",        # 橙色
      "Other" = "lightgray"              # 灰色
    )
  ) +
  ggtitle("Highlight of Epithelial and Stromal Cells in Full Trajectory") +
  theme_minimal()

