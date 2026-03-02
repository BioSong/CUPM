library(clusterProfiler)
library(enrichplot)
load("./CUPM/data/RNAseq/GSE183202/GSE183202_Exp_Cli.RData")
exp2=log2(fpkm1+1)
#------------------filter out 低表达的gene
index1=rowMeans(exp2<1)
vdata1=exp2[which(index1<0.5),]
cli1[cli1==""]=NA
#------------------------------significantly related pathways
cuscore=read.table("./CUPM/results/V2_202507/CuScore_Calculation_Valitation/CuScore_GSE283202.txt",sep="\t",header=T,as.is=T,quote="\"")
#-------GSEA
library(GSVA)
#------------------准备pathways
load("./CUPM/data/MsigDb/MsigDb_genesets.RData")
glist=msigdb_geneid_list$genesets
names(glist)=msigdb_geneid_list$geneset.names
paths=read.table("./CUPM/data/MsigDb/pathway_for_analysis.txt",sep="\t",header=T,as.is=T,quote="\"")
glist=glist[match(paths[,1],names(glist))]
#glist=glist[grepl("WNT|BETA_CATENIN|BETA- CATENIN|WNT_BETA",names(glist))]
#---------------------------将vdata1的Gene ID转化为Gene symbol
library(org.Hs.eg.db)
trans_vdata1 <- select(
  org.Hs.eg.db,
  keys = rownames(vdata1),
  keytype = "SYMBOL",
  columns = "ENTREZID"
)
trans_vdata1=na.omit(trans_vdata1)
gene=intersect(rownames(vdata1),trans_vdata1[,1])
vdata11=vdata1[match(gene,rownames(vdata1)),];
trans_vdata1=trans_vdata1[match(gene,trans_vdata1[,1]),];
rownames(vdata11)=trans_vdata1$ENTREZID
#-----------ssGSEA,correlation analysis
param <- ssgseaParam(
  exprData = vdata11,
  geneSets = glist
)
gsea1=gsva(param)
#-------------------
corresult=matrix(,nrow(gsea1),3)
for(i in 1:nrow(gsea1)){
  cor1=cor.test(cuscore$Final_Score,gsea1[i,])
    corresult[i,1:2]=c(cor1$estimate,cor1$p.value)
}
corresult[,3]=p.adjust(corresult[,2],method="BH")
rownames(corresult)=rownames(gsea1)
colnames(corresult)=c("R","pvalue","FDR")
corresult=corresult[order(corresult[,3]),]
#-------------------------------群体分析
cuscore1=cuscore$Final_Score
group <- ifelse(cuscore1 > median(cuscore1), "High", "Low")
log2FC <- log2((rowMeans(vdata11[, group == "High"]) + 1e-6) /
               (rowMeans(vdata11[, group == "Low"]) + 1e-6))
geneList <- sort(log2FC, decreasing = TRUE)
geneSets <- do.call(rbind, lapply(names(glist), function(gs){
  data.frame(gs_name = gs, gene_symbol = glist[[gs]], stringsAsFactors = FALSE)
}))
gsea_res <- GSEA(geneList, TERM2GENE = geneSets)
gsea_res1=gsea_res@result[,2:6]
# 提取结果
gsea_df <- as.data.frame(gsea_res)
top_pathways <- gsea_df[gsea_df$p.adjust < 0.05, ]
#top_pathways <- gsea_df[gsea_df$pvalue < 0.05, ]
top_pathway1=top_pathways[grepl("KEGG|HALLMARK",rownames(top_pathways)),]
top_pathway1=top_pathway1[1:20,]
ggplot(top_pathway1, aes(x = reorder(Description, NES), y = NES, fill = NES > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal()
#-----------------------------单挑通路
gseaplot2(gsea_res, geneSetID ="KEGG_WNT_SIGNALING_PATHWAY",pvalue_table=T)
gseaplot2(gsea_res, geneSetID ="KEGG_CITRATE_CYCLE_TCA_CYCLE",pvalue_table=T)
gseaplot2(gsea_res, geneSetID =c("KEGG_WNT_SIGNALING_PATHWAY","KEGG_CITRATE_CYCLE_TCA_CYCLE"),pvalue_table=T)
gseaplot2(gsea_res, geneSetID =c("KEGG_WNT_SIGNALING_PATHWAY","HALLMARK_APOPTOSIS"),pvalue_table=T)
#----------------------------计算CUscore与COX17、SCO1、COA6、SOD1、ATP7A等铜稳态基因相关性
cugene1=c("COX17","SCO1","COA6","SOD1","ATP7A")
vdata12=vdata1[na.omit(match(cugene1,rownames(vdata1))),]
cor.test(cuscore$Final_Score,vdata12[1,])
cor.test(cuscore$Final_Score,vdata12[2,])
cor.test(cuscore$Final_Score,vdata12[3,])
cor.test(cuscore$Final_Score,vdata12[4,])
scater_df=data.frame(CuScore=cuscore$Final_Score,COX17=vdata12[1,],SOD1=vdata12[4,])
library(scatterplot3d)
scatterplot3d(vdata12[1,], vdata12[4,],cuscore$Final_Score,
              pch = 16, color = "blue",
              xlab = "COX17", ylab = "SOD1", zlab = "CuScore")
library(GGally)
GGally::ggpairs(scater_df[, c("CuScore", "COX17", "SOD1")],
                upper = list(continuous = wrap("cor", size = 3)),
                lower = list(continuous = wrap("points", alpha = 0.5, size = 1.5)))
