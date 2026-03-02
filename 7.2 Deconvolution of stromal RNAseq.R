library(Seurat)
library(MuSiC)
library(SingleCellExperiment)
library(tidyverse)
#-----------------------------
load('./CUPM/results/V2_202507/5_CuScore_Stromal/Seurat_Stromal.RData')
stromal$CuScore_Final=stromal@meta.data$CuScore$Final_Score
stromal <- subset(stromal, subset = loc == "metastasis")
# 转换为 SCE 对象
sce_ref <- as.SingleCellExperiment(stromal)
# 提取表达矩阵
expr_ref <- assay(sce_ref, "counts")
# 提取 metadata
meta_ref <- data.frame(
  cell_type = stromal$subtype,
  subject = stromal$patient_id,
  stringsAsFactors = FALSE
)
#------------bulk RNAseq
load("./CUPM/results/V2_202507/3_CuScore_CellType/GSE190609_Cell_Proportions_v2.RData")
#est_prop,sce,bulk_expr,bulk_cli
est_prop <- music_prop(
  bulk.mtx = bulk_expr,
  sc.sce = sce_ref,
  clusters = 'subtype',
  samples = 'pat_id2',
  normalize = TRUE
)
prop_mat <- est_prop$Est.prop.weighted
prop_df  <- as.data.frame(prop_mat) %>% rownames_to_column("Sample")
meta=bulk_cli[,c("Sample_ID","Tissue_type")]
colnames(meta)=c("Sample","Location")
#------------------
bulk_cuscore=read.table("./CUPM/results/V2_202507/2_CuScore_Location_Compare/GSE190609_CuScore.txt",header=T,row.names=1,as.is=T,quote="\"")
cupro=bulk_cuscore %>% rownames_to_column("Sample")
colnames(cupro)[ncol(cupro)]="CuproScore"
cupro=cupro[,c("Sample","CuproScore")]
df <- prop_df %>%
  left_join(cupro, by="Sample")
#--------------------
rownames(df)=df$Sample
df_subtype=df[,c(1,8,2:7)]
head(df)
colnames(df_subtype)[2]="CuScore"
#-----------------------------------
library(car)
library(relaimpo)
library(ggplot2)
library(dplyr)
library(broom)
# 选择模型变量
model_vars <- c("others", "MSCs", "TECs", "CAAs", "pericyte")  # 去掉others为baseline
formula <- as.formula(paste("CuScore ~", paste(model_vars, collapse = " + ")))
lm_fit <- lm(formula, data = df_subtype)
relimp <- calc.relimp(lm_fit, type = "lmg", rela = TRUE)
#-----------------------
model_vars <- c("CAF", "MSCs", "TECs", "CAAs", "others")
formula <- as.formula(
  paste("CuScore ~", paste(model_vars, collapse = " + "))
)
# 构建线性模型
lm_fit <- lm(formula, data = df_subtype)
# 提取回归系数
coef_df <- broom::tidy(lm_fit) %>% 
  filter(term != "(Intercept)") %>%
  rename(subtype = term, beta = estimate, pvalue = p.value)
# Semipartial R2 (LMG)
relimp <- calc.relimp(lm_fit, type = "lmg")
lmg_df <- data.frame(
  subtype = names(relimp$lmg),
  lmg = relimp$lmg
)
# 合并结果
result_df <- left_join(coef_df, lmg_df, by = "subtype") %>%
  mutate(
    ContributionRank = rank(-lmg),
    beta_label = sprintf("β = %.3f", beta),
    lmg_label = sprintf("R² = %.3f", lmg)
  )
print(result_df)

#---------------------------------------------
result_df$subtype <- factor(
  result_df$subtype, 
  levels = result_df$subtype[order(result_df$lmg)]
)

ggplot(result_df, aes(x = lmg, y = subtype)) +
  geom_point(size = 4) +
  geom_segment(aes(x = 0, xend = lmg, y = subtype, yend = subtype)) +
  geom_text(aes(label = beta_label), vjust = -0.5, size = 4) +
  theme_classic(base_size = 14) +
  labs(
    x = "Semipartial R² (LMG) – Contribution to CuScore",
    y = "Stromal Subtype",
    title = "Subtype Contribution to Bulk CuScore Variance"
  )
