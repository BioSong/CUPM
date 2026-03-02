BiocManager::install("TOAST")
devtools::install_github('xuranw/MuSiC')
#-------------------
library(Matrix)
library(MuSiC)
library(SingleCellExperiment)
#Step 1: Load scRNA-seq Data
# Read raw UMI matrix (genes × cells)
umi_mat <- read.table(gzfile("./CUPM/data/scRNAseq/CRC/GSE132465/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"), header = TRUE, row.names = 1, check.names = FALSE)
exprs_mat <- as.matrix(umi_mat)

# Read cell annotation
annot <- read.table(gzfile("./CUPM/data/scRNAseq/CRC/GSE132465/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"), header = TRUE, sep = "\t")
head(annot)

sce <- SingleCellExperiment(
  assays = list(counts = exprs_mat),
  colData = annot
)
table(sce$Cell_type) 

#Step 2: Load Bulk RNA-seq Data
# Format: genes × samples matrix (raw counts or TPM)
bulk_expr <- read.table("./CUPM/data/RNAseq/GSE190609/GSE190609_rpkm_SYMBOL.txt", header = TRUE, row.names = 1, check.names = FALSE)
bulk_expr=as.matrix(bulk_expr)
bulk_expr=log2(bulk_expr+1)
bulk_cli=read.table("./CUPM/data/RNAseq/GSE190609/GSE190609_Clinical.txt",sep="\t",header=T,as.is=T,quote="\"")
bulk_cli=bulk_cli[match(colnames(bulk_expr),bulk_cli[,1]),]

tcell_marker_list <- list("T cell" = cellmarker[which(cellmarker[,1]=="T_cell"),2])
# Run deconvolution
est_prop <- music_prop(
  bulk.mtx = bulk_expr,
  sc.sce = sce,
  clusters = 'Cell_type',
  samples = 'Patient',
  normalize = TRUE
)

# View estimated proportions
head(est_prop$Est.prop.weighted)

save(est_prop,sce,bulk_expr,bulk_cli,file="./CUPM/results/V2_202507/3_CuScore_CellType/GSE190609_Cell_Proportions_v2.RData")

#Visualization
library(ggplot2)
library(reshape2)

df <- melt(est_prop$Est.prop.weighted)
colnames(df) <- c("Sample", "CellType", "Proportion")

ggplot(df, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Estimated Cell Type Proportions in CRC Bulk Samples GSSE190609") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#-----------------------------------------------------------
