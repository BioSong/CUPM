#1) CopyKAT：按 pat_id+loc 跑（上皮 + 免疫/内皮/周细胞做参考），并写回 malignant_status
library(Seurat)
library(copykat)
load("./CUPM/results/V2_202507/5_CuScore_Stromal/cellchat_caf/scdata_cellchat.RData")

DefaultAssay(data1) <- "RNA"
data1$pat_id2=data1$pat_id
data1$pat_id2[which(data1$pat_id2!="B")]="C"
ref_types <- c("T_cell","B_cell","Myeloid_cell","TECs","pericyte")

run_copykat_one_sample <- function(seu, pid, loc_value,
                                   ref_types = ref_types,
                                   min_epi = 200, min_ref = 200,
                                   ncores = 4) {
  seu_s <- subset(seu, subset = pat_id2 == pid & loc == loc_value)

  epi_cells <- colnames(seu_s)[seu_s$interaction_group == "Epithelial_cell"]
  ref_cells <- colnames(seu_s)[seu_s$interaction_group %in% ref_types]

  message(pid, "/", loc_value, ": epi=", length(epi_cells), ", ref=", length(ref_cells))

  if (length(epi_cells) < min_epi || length(ref_cells) < min_ref) {
    warning(paste0("Skip CopyKAT for ", pid, "/", loc_value, " (epi<", min_epi, " or ref<", min_ref, ")."))
    out <- seu_s[, epi_cells, drop = FALSE]
    out$copykat_pred <- NA
    return(out)
  }

  use_cells <- c(epi_cells, ref_cells)
  rawmat <- GetAssayData(seu_s[, use_cells], slot = "counts")

  ck <- copykat(
    rawmat = as.matrix(rawmat),
    id.type = "S",  # gene symbol
    ngene.chr = 5,
    win.size = 25,
    KS.cut = 0.1,
    distance = "euclidean",
    norm.cell.names = ref_cells,
    n.cores = ncores
  )

  #pred <- ck$prediction
  #epi_cells_corrected <- gsub("\\-1_1$", ".1_1", epi_cells)
  #pred_epi <- pred[match(epi_cells_corrected, rownames(pred)), ]

  #out <- seu_s[, epi_cells]
  #out$copykat_pred <- pred_epi$copykat.pred
  #return(out)
  return(ck)
}

samples_to_run <- unique(data.frame(pat_id = data1$pat_id2, loc = data1$loc))
samples_to_run <- subset(samples_to_run, loc %in% c("primary","metastasis"))

epi_ck_list <- lapply(seq_len(nrow(samples_to_run)), function(i) {
  run_copykat_one_sample(data1, samples_to_run$pat_id[i], samples_to_run$loc[i],
                                  ref_types,
                                  min_epi = 200, min_ref = 200,
                                  ncores = 6)
})

# 1) 合并所有 prediction
pred_df <- do.call(rbind, lapply(epi_ck_list, function(ck) {  
  p <- ck$prediction  
  p <- p[, c("cell.names", "copykat.pred")]  
  p}))
# 2) 建立映射：names = cell.names, values = copykat.pred
pred_map <- setNames(pred_df$copykat.pred, pred_df$cell.names)
# 4) 回填到 data1
data1$copykat_pred <- NA
# 先按原始名字匹配
hit1 <- intersect(colnames(data1), names(pred_map))
data1$copykat_pred[hit1] <- pred_map[hit1]
# 写回全对象
# 5) 生成 malignant_status
data1$malignant_status <- ifelse(  
  data1$copykat_pred == "aneuploid", "malignant",  
  ifelse(data1$copykat_pred == "diploid", "non_malignant", NA))
# 6) 快速验收
cat("Matched copykat_pred:", sum(!is.na(data1$copykat_pred)), "/", ncol(data1), "\n")
print(table(data1$copykat_pred, useNA = "ifany"))
print(table(data1$malignant_status, useNA = "ifany"))
save(data1,"./CUPM/")
# QC：各病人/部位 malignant 比例
#----------------绘图
library(dplyr)
library(ggplot2)
library(scales)
df_epi <- data1@meta.data |>  
    mutate(pat_id = as.character(pat_id), loc = as.character(loc)) |>  
    filter(!is.na(copykat_pred), copykat_pred %in% c("aneuploid","diploid"))
df_all <- df_epi |>
  mutate(group = dplyr::case_when(
    pat_id == "B" & loc == "primary"     ~ "B_primary",
    pat_id == "B" & loc == "metastasis"  ~ "B_metastasis",
    loc == "metastasis"                  ~ paste0(pat_id, "_metastasis"),
    TRUE ~ NA_character_
  )) |>
  filter(!is.na(group)) |>
  group_by(group) |>
  summarise(
    n_epi = n(),
    frac_aneu = mean(copykat_pred == "aneuploid"),
    .groups = "drop"
  )

ggplot(df_all, aes(x = group, y = frac_aneu)) +
  geom_col(fill = "firebrick", width = 0.75) +
  geom_text(aes(label = paste0("n=", n_epi)), vjust = -0.4, size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = NULL, y = "Aneuploid fraction (epithelial)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#------------------------------only patient B
df_B2 <- df_epi |>
  filter(pat_id == "B") |>
  group_by(loc) |>
  summarise(
    n_epi = n(),
    frac_aneu = mean(copykat_pred == "aneuploid"),
    .groups = "drop"
  )

ggplot(df_B2, aes(x = loc, y = frac_aneu)) +
  geom_col(width = 0.35, fill = "firebrick") +
  geom_text(aes(label = paste0("n=", n_epi)), vjust = -0.4, size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = NULL, y = "Aneuploid fraction (epithelial)") +
  theme_classic()

#-----------------------------------------------
save(data1,file="./CUPM/results/V2_202507/6_CuScore_upstream/scdata_copykat.RData")
