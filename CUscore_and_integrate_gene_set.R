score_and_integrate_gene_set <- function(expr_matrix, gene_set, seurat_obj = NULL) {
  library(AUCell)
  library(GSVA)
  library(UCell)
  library(singscore)
  library(Seurat)
  library(RobustRankAggreg)
  library(SummarizedExperiment)
  
  samples <- colnames(expr_matrix)
  genes <- rownames(expr_matrix)
  
  # Ensure gene set is in data
  gene_set <- intersect(gene_set, genes)
  if (length(gene_set) < 5) stop("Too few genes matched in expression matrix.")
  
  score_df <- data.frame(row.names = samples)
  
  ### 1. AUCell
  cells_rankings <- AUCell_buildRankings(expr_matrix,plotStats=FALSE)
  auc <- AUCell_calcAUC(gene_set, cells_rankings)
  score_df$AUCell <- as.numeric(getAUC(auc)[1, ])
  
  ### 2. ssGSEA (GSVA)
  param <- ssgseaParam(exprData = expr_matrix,geneSets = list(GS = gene_set))
  ss_res=gsva(param)
  score_df$ssGSEA <- ss_res[1,]
    ### 3. UCell
  ucell_scores <- ScoreSignatures_UCell(expr_matrix, features = list(GS = gene_set))
  score_df$UCell <- ucell_scores[,1]

  ### 4. singscore
  ranked <- rankGenes(expr_matrix)
  sing <- simpleScore(ranked, upSet = gene_set,centerScore = TRUE)
  score_df$singscore <- sing[,1]

  ### 5. AddModuleScore
  if (is.null(seurat_obj)) {
    warning("Seurat object not provided; Create it")
    #expr_df=as.data.frame(expr_matrix)
	expr_df=expr_matrix
    seurat_obj <- CreateSeuratObject(counts = expr_df, project = "BulkRNAseq", assay = "RNA")
    seurat_obj <- SetAssayData(object = seurat_obj,layer="data",new.data = expr_df,assay = "RNA")
  }
  seurat_obj <- AddModuleScore(seurat_obj, features = list(gene_set), name = "Module")
  mod_scores <- seurat_obj@meta.data$Module
  score_df$Seurat <- mod_scores

  ### 分位数归一化排名函数
  quantile_rank <- function(x) {
    ecdf(x)(x)  # 经验累积分布函数
  }

  # 3. 应用分位数归一化排名
  score_df$rank_AUCell <- quantile_rank(score_df$AUCell)
  score_df$rank_ssGSEA <- quantile_rank(score_df$ssGSEA)
  score_df$rank_UCell <- quantile_rank(score_df$UCell)
  score_df$rank_singscore <- quantile_rank(score_df$singscore)
  score_df$rank_Seurat <- quantile_rank(score_df$Seurat)

  # 4. 计算综合得分（加权平均）
  # 基于得分组间相关性赋予权重
  cor_matrix <- cor(score_df[, 6:10])
  weight <- apply(cor_matrix, 1, mean)  # 平均相关性作为权重

  score_df$Composite_Score <- (
    weight[1] * score_df$rank_AUCell +
    weight[2] * score_df$rank_ssGSEA +
    weight[3] * score_df$rank_UCell +
    weight[4] * score_df$rank_singscore +
    weight[5] * score_df$rank_Seurat
  ) / sum(weight)

  # 5. 最终归一化到0-1范围（可选）
  score_df$Final_Score <- (score_df$Composite_Score - min(score_df$Composite_Score)) /
                          (max(score_df$Composite_Score) - min(score_df$Composite_Score))
  score_df=score_df[,c("AUCell","ssGSEA","UCell","singscore","Seurat","Composite_Score","Final_Score")]
  return(score_df)
}
