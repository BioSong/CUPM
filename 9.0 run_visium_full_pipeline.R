suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(Matrix)
  library(readr)
  library(stringr)
  library(tibble)
  library(spacexr)
  library(CMScaller)
})

# =========================
# 用户需要修改的参数
# =========================

# 1) Visium 根目录（包含4个样本文件夹）
VISIUM_ROOT <- "./CUPM/data/STRNAseq/ST_CRC_CMS/visimu" 

# 2) 4个样本文件夹名（与实际一致）
SAMPLES <- c("SN048_A121573_Rep1",
             "SN048_A121573_Rep2",
             "SN048_A416371_Rep1",
             "SN048_A416371_Rep2")

# 3) 你的 scRNA reference Seurat RDS 路径
SCREF_RDS <- "./CUPM/results/V2_202507/6_CuScore_upstream/scdata_copykat.RData"   # 改成你的文件名/路径

# 4) scRNA reference 的细胞类型列名
CELLTYPE_COL <- "interaction_group"

# 5) Cuproptosis 基因集（你可替换成你AUCell用的那套）
#    注意：必须与 Visium 的 rownames 基因命名一致（一般是Symbol）
geneset=read.table("./CUPM/data/CuproptosisRelatedGenes/CRGenes.txt",sep="\t",header=T,as.is=T,quote="\"")
CUPRO_GENES=unique(geneset$SYMBOL[c(which(geneset$Group=="group1"),which(geneset$Group=="group3"))])

# 6) 你关注的两个细胞类名（必须与 interaction_group 水平一致）
CAF_NAME <- "ECMCAF"
EPI_NAME <- "Epithelial_cell"

# 7) RCTD 线程
RCTD_CORES <- 6

# 8) 输出目录
OUTDIR <- "./CUPM/results/V2_202507/7_STanalysis/analysis_out"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)


# =========================
# 工具函数
# =========================

msg <- function(...) cat(sprintf(...), "\n")

# ---- 读取病理注释  并标准化列名 ----
read_pathology_rds <- function(sample_dir) {
  anno_dir <- file.path(sample_dir, "Pathology_SpotAnnotations")
  if (!dir.exists(anno_dir)) return(NULL)

  rds_files <- list.files(anno_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(rds_files) == 0) return(NULL)

  pano <- readRDS(rds_files[1])
  if (!is.data.frame(pano)) stop("Pathology RDS is not a data.frame")

  # 兼容列名：你的截图里是 spot_id 和 Pathologist_Annotations
  if (!("spot_id" %in% colnames(pano))) {
    stop("Cannot find column 'spot_id' in pathology RDS.")
  }

  # 病理注释列名可能略有不同
  anno_col <- intersect(colnames(pano),
                        c("Pathologist_Annotations",
                          "Pathologist Annotation",
                          "Pathologist_Annotation",
                          "annotation",
                          "label"))
  if (length(anno_col) == 0) {
    stop("Cannot find pathology annotation column in RDS. Paste colnames(pano).")
  }
  anno_col <- anno_col[1]

  out <- pano %>%
    dplyr::transmute(
      barcode = spot_id,
      path_anno = .data[[anno_col]]
    )

  out
}

add_pathology_to_seurat <- function(obj, anno_df) {
  if (is.null(anno_df) || nrow(anno_df) == 0) {
    obj$path_anno <- NA_character_
    return(obj)
  }

  seurat_bc <- Cells(obj)
  anno_bc <- anno_df$barcode

  # -1 后缀容错
  if (any(stringr::str_detect(seurat_bc, "-1$")) && !any(stringr::str_detect(anno_bc, "-1$"))) {
    anno_df$barcode <- paste0(anno_df$barcode, "-1")
  }
  if (!any(stringr::str_detect(seurat_bc, "-1$")) && any(stringr::str_detect(anno_df$barcode, "-1$"))) {
    anno_df$barcode <- stringr::str_replace(anno_df$barcode, "-1$", "")
  }

  m <- match(seurat_bc, anno_df$barcode)
  obj$path_anno <- anno_df$path_anno[m]

  # 诊断：看匹配率
  match_rate <- mean(!is.na(obj$path_anno))
  message(sprintf("Pathology annotation matched: %.1f%% spots", 100 * match_rate))

  obj
}


# ---- 读取单个 Visium 样本（按你目录结构：sample_dir就是包含spatial/的那层） ----
load_one_visium <- function(sample_dir, sample_id) {
  obj <- Load10X_Spatial(data.dir = sample_dir, assay = "Spatial")
  obj$sample_id <- sample_id
  obj
}

# ---- QC + SCTransform ----
qc_sctransform <- function(obj,
                           min_counts = 500,
                           min_features = 200,
                           max_mt = 20) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(obj, subset = nCount_Spatial > min_counts &
                      nFeature_Spatial > min_features &
                      percent.mt < max_mt)
  obj <- SCTransform(obj, assay = "Spatial", verbose = FALSE)
  obj
}

# ---- CMS：CMScaller v0.9.2 ----
run_cms <- function(obj,
                    assay = "SCT",
                    slot = "data",
                    rna_seq = FALSE,
                    verbose = FALSE) {
  # 默认先填 NA，避免失败时没有列
  obj$CMS_class <- NA_character_
  obj$CMS2_score <- NA_real_
  obj$CMS4_score <- NA_real_
  if (!requireNamespace("CMScaller", quietly = TRUE)) {
    warning("CMScaller 未安装/未加载：跳过 CMS。")
    return(obj)
  }
  DefaultAssay(obj) <- assay
  expr <- Seurat::GetAssayData(obj, assay = assay, slot = slot)
  expr <- as.matrix(expr)  # gene x spot
  # --- 内部小工具：从 CMScaller 返回对象里找“分类”和“CMS1-4 分数” ---
  .find_class_vector <- function(res, spot_names) {
    # 在 res 的各元素里找包含 CMS1-4 标签的向量/列
    is_cms_label <- function(v) any(grepl("^CMS[1-4]$", as.character(v)))
    # 1) 直接是向量
    for (nm in names(res)) {
      x <- res[[nm]]
      if (is.atomic(x) && length(x) == length(spot_names) && is_cms_label(x)) {
        out <- as.character(x); names(out) <- spot_names
        return(out)
      }
    }
    # 2) data.frame / matrix 某一列是 CMS label
    for (nm in names(res)) {
      x <- res[[nm]]
      if (is.data.frame(x) || is.matrix(x)) {
        df <- as.data.frame(x)
        for (cn in colnames(df)) {
          v <- df[[cn]]
          if (length(v) == nrow(df) && is_cms_label(v)) {
            out <- as.character(v)
            # 用行名对齐 spot
            if (!is.null(rownames(df)) && all(spot_names %in% rownames(df))) {
              names(out) <- rownames(df)
              out <- out[spot_names]
              return(out)
            }
            # 没行名就假设顺序一致
            if (length(out) == length(spot_names)) {
              names(out) <- spot_names
              return(out)
            }
          }
        }
      }
    }
    NULL
  }
  # --- 跑 CMScaller ---
  spot_names <- colnames(obj)
  cms_res <- tryCatch(
    {
      if (verbose) message("Running CMScaller ...")
      CMScaller::CMScaller(expr, RNAseq = rna_seq,rowNames="symbol",FDR=0.2)
    },
    error = function(e) {
      warning("CMScaller 运行失败：", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(cms_res)) return(obj)
  # --- 抽取分类与分数，写回 meta.data ---
  cms_class <- .find_class_vector(cms_res, spot_names)
  if (!is.null(cms_class)) {
    obj$CMS_class <- cms_class[spot_names]
  } else {
    warning("未能从 CMScaller 结果中自动识别 CMS 分类列；已保留 NA。")
  }
  obj$CMS2_score <- as.numeric(1-cms_res$d.CMS2)
  obj$CMS4_score <- as.numeric(1-cms_res$d.CMS4)
  return(obj)
}


# ---- 运行 RCTD ----
run_rctd_one <- function(obj, ref, max_cores = 8) {
  if (!requireNamespace("spacexr", quietly = TRUE)) {
    stop("Package spacexr not installed. Install with: remotes::install_github('dmcable/spacexr')")
  }
  library(spacexr)

  counts_sp <- GetAssayData(obj, assay = "Spatial", slot = "counts")

  img_name <- names(obj@images)[1]
  coords <- GetTissueCoordinates(obj, image = img_name) %>% as.data.frame()

  # coords 行名通常是 barcode；确保顺序对齐 counts
  if (!all(colnames(counts_sp) %in% rownames(coords))) {
    # 有时coords会有barcode列
    if ("barcode" %in% colnames(coords)) {
      rownames(coords) <- coords$barcode
    }
  }
  coords <- coords[colnames(counts_sp), , drop = FALSE]

  # Seurat 的 GetTissueCoordinates 通常返回 x y
  if (!all(c("x", "y") %in% colnames(coords))) {
    # 有些版本是 imagerow/imagecol 或 row/col
    candidates <- colnames(coords)
    msg("Coordinate columns found: %s", paste(candidates, collapse = ", "))
    stop("Cannot find x/y in tissue coordinates. Paste head(GetTissueCoordinates(obj)).")
  }

  puck <- SpatialRNA(coords = coords[, c("x", "y")],
                     counts = counts_sp,
                     nUMI = colSums(counts_sp))

  rctd <- create.RCTD(spatialRNA = puck, reference = ref, max_cores = max_cores)
  rctd <- run.RCTD(rctd, doublet_mode = "full")
  rctd
}

# ---- 从 RCTD 提取 weights 并写回 Seurat ----
add_rctd_weights <- function(obj, rctd, prefix = "RCTD_") {
  # 兼容不同对象结构
  w <- NULL
  if (!is.null(rctd@results$weights)) w <- rctd@results$weights
  if (is.null(w) && !is.null(rctd@results$weights_doublet)) w <- rctd@results$weights_doublet

  if (is.null(w)) stop("Cannot find RCTD weights in rctd@results. Please inspect names(rctd@results).")

  # w: spot x celltype（通常行名为barcode）
  # 对齐 Seurat spot
  if (!all(Cells(obj) %in% rownames(w))) {
    stop("RCTD weights rownames do not match Seurat barcodes. Paste head(Cells(obj)) and head(rownames(w)).")
  }
  w2 <- w[Cells(obj), , drop = FALSE]
  colnames(w2) <- paste0(prefix, colnames(w2))

  obj <- AddMetaData(obj, metadata = as.data.frame(w2))
  obj
}

# ---- 邻接检验：CAF高spot的邻居Epi是否更高（简单距离阈值版） ----
neighbor_enrichment_test <- function(obj,
                                     caf_col = "RCTD_CAF",
                                     epi_col = "RCTD_Epithelial",
                                     subset_expr = NULL,
                                     k = 6,
                                     n_perm = 200) {
  # subset_expr: 例如 obj$CMS_class=="CMS4"
  md <- obj@meta.data

  if (!is.null(subset_expr)) {
    keep <- eval(parse(text = subset_expr), envir = md)
    obj2 <- subset(obj, cells = rownames(md)[which(keep)])
  } else {
    obj2 <- obj
  }

  img_name <- names(obj2@images)[1]
  coords <- GetTissueCoordinates(obj2, image = img_name) %>% as.data.frame()
  if (!all(c("x","y") %in% colnames(coords))) stop("coords must have x,y.")
  coords <- coords[Cells(obj2), c("x","y"), drop = FALSE]

  caf <- obj2@meta.data[[caf_col]]
  epi <- obj2@meta.data[[epi_col]]
  if (is.null(caf) || is.null(epi)) stop("caf_col/epi_col not found in meta.data.")

  # 计算kNN（用欧氏距离，简单实现）
  xy <- as.matrix(coords)
  D <- as.matrix(dist(xy))
  diag(D) <- Inf

  nn_idx <- apply(D, 1, function(v) order(v)[1:k])

  # 观测统计：CAF高（上四分位）spot的邻居Epi均值
  caf_hi <- caf >= quantile(caf, 0.75, na.rm = TRUE)
  obs <- mean(sapply(which(caf_hi), function(i) mean(epi[nn_idx[, i]], na.rm = TRUE)), na.rm = TRUE)

  # 置换：打乱 CAF 标签
  set.seed(1)
  perm_stats <- replicate(n_perm, {
    caf_perm <- sample(caf)
    caf_hi_p <- caf_perm >= quantile(caf_perm, 0.75, na.rm = TRUE)
    mean(sapply(which(caf_hi_p), function(i) mean(epi[nn_idx[, i]], na.rm = TRUE)), na.rm = TRUE)
  })

  p <- mean(perm_stats >= obs, na.rm = TRUE)

  tibble(
    n_spots = ncol(obj2),
    k = k,
    n_perm = n_perm,
    obs_neighbor_epi_mean = obs,
    perm_mean = mean(perm_stats, na.rm = TRUE),
    perm_sd = sd(perm_stats, na.rm = TRUE),
    p_value = p
  )
}


# =========================
# 1) 读取 scRNA reference，构建 RCTD Reference
# =========================
load("./CUPM/results/V2_202507/6_CuScore_upstream/scdata_copykat.RData")
sc_ref <- subset(data1, subset = pat_id == "B")
stopifnot(inherits(sc_ref, "Seurat"))
stopifnot(CELLTYPE_COL %in% colnames(sc_ref@meta.data))
# 统计各类细胞数
ct_tab <- table(sc_ref$interaction_group)# 需要保留的类别
keep_types <- names(ct_tab[ct_tab >= 50])# subset 保留细胞
sc_ref <- subset(sc_ref, subset = interaction_group %in% keep_types)

counts_sc <- GetAssayData(sc_ref, assay = "RNA", slot = "counts")
celltypes <- sc_ref$interaction_group
celltypes <- as.factor(celltypes)
names(celltypes) <- colnames(sc_ref)  # 给上 cell barcode
ref <- Reference(counts = counts_sc, cell_types = celltypes)

if (!requireNamespace("spacexr", quietly = TRUE)) {
  msg("Installing spacexr is required for RCTD.")
  msg("Run: remotes::install_github('dmcable/spacexr')")
}


msg("scRNA celltypes: %d", length(levels(celltypes)))
msg("Top celltypes:\n%s", paste(head(levels(celltypes), 20), collapse = ", "))


# =========================
# 2) 批量读取 4 个 Visium + 病理注释
# =========================
# ---- Cuproptosis score（模块评分） ----
source("./CUPM/program/V2_202507/CUscore_and_integrate_gene_set.R")
 
objs <- list()

for (sid in SAMPLES) {
  sample_dir <- file.path(VISIUM_ROOT, sid)
  msg("\n--- Loading sample: %s", sid)

  obj <- load_one_visium(sample_dir, sid)

  anno <- read_pathology_rds(sample_dir)
  obj <- add_pathology_to_seurat(obj, anno)

  # 保存QC前基本统计
  qc0 <- obj@meta.data %>%
    summarise(n_spots = n(),
              med_counts = median(nCount_Spatial),
              med_features = median(nFeature_Spatial),
              mt_med = median(PercentageFeatureSet(obj, pattern = "^MT-")))
  write.csv(qc0, file.path(OUTDIR, paste0("qc_pre_", sid, ".csv")), row.names = FALSE)

  # QC + SCT
  obj <- qc_sctransform(obj)

  # Cupro score
  expr <- GetAssayData(obj, assay = "SCT", slot = "data")
  expr <- as.matrix(expr)
  cuscore <- score_and_integrate_gene_set(expr_matrix=expr, gene_set=CUPRO_GENES)
  obj <- AddMetaData(obj, metadata = cuscore$Final_Score, col.name = "CuproScore_1")

  # CMS（占位：若你能跑CMScaller，把run_cms()替换）
  obj <- run_cms(obj)

  objs[[sid]] <- obj
}

# =========================
# 3) RCTD 去卷积（每个样本分别跑）
# =========================

rctd_list <- list()
objs2 <- list()

for (sid in names(objs)) {
  msg("\n--- Running RCTD: %s", sid)
  rctd <- run_rctd_one(objs[[sid]], ref = ref, max_cores = RCTD_CORES)
  rctd_list[[sid]] <- rctd

  obj2 <- add_rctd_weights(objs[[sid]], rctd, prefix = "RCTD_")
  objs2[[sid]] <- obj2

  saveRDS(rctd, file.path(OUTDIR, paste0("rctd_", sid, ".rds")))
  saveRDS(obj2, file.path(OUTDIR, paste0("seurat_with_rctd_", sid, ".rds")))
}

saveRDS(objs2, file.path(OUTDIR, "seurat_all_with_rctd.rds"))

# =========================
# 4) 核心结果：CAF/Epi 与 Cupro score + 病理注释可视化
# =========================
plot_cms <- function(obj, cms_col = "CMS_class", pt.size = 1.6) {
  stopifnot(cms_col %in% colnames(obj[[]]))

  cms <- as.character(obj[[cms_col]][, 1])
  cms[is.na(cms)] <- "Unassigned"
  obj[[cms_col]] <- factor(cms, levels = c("CMS1","CMS2","CMS3","CMS4","Unassigned"))

  cols <- c(
    "CMS1" = "#DE9A13",
    "CMS2" = "#0A6AA6",
    "CMS3" = "#C71585",
    "CMS4" = "#0C9672",
    "Unassigned" = "grey85"
  )

  SpatialDimPlot(obj, group.by = cms_col, pt.size.factor = pt.size) +
    scale_fill_manual(values = cols, drop = FALSE) +
    ggtitle(paste0(unique(obj$sample_id), " | CMS regions")) +
    theme(plot.title = element_text(face = "bold"))
}

plot_continuous <- function(obj, feature, title = NULL, pt.size = 1.6,
                            palette = c("grey95", "gold", "red3")) {
  stopifnot(feature %in% colnames(obj[[]]) || feature %in% rownames(obj))

  p <- SpatialFeaturePlot(obj, features = feature, pt.size.factor = pt.size) +
    scale_fill_gradientn(colors = palette) +
    theme(plot.title = element_text(face = "bold"))

  if (is.null(title)) title <- paste0(unique(obj$sample_id), " | ", feature)
  p + ggtitle(title)
}

plot_ecmcaf_epi <- function(obj,
                            ecmcaf = "RCTD_ECMCAF",
                            epi = "RCTD_Epithelial_cell",
                            pt.size = 1.6,
                            palette = c("grey95", "dodgerblue3", "navy")) {
  stopifnot(all(c(ecmcaf, epi) %in% colnames(obj[[]])))

  p <- SpatialFeaturePlot(
    obj,
    features = c(ecmcaf, epi),
    pt.size.factor = pt.size,
    ncol = 2
  )

  # 给每个子图统一使用同一个渐变色标
  p <- p & scale_fill_gradientn(colors = palette)

  # 标题的 theme 单独放到 annotation 里，避免 plot_theme 报错
  p + plot_annotation(
    title = paste0(unique(obj$sample_id), " | RCTD ECMCAF & Epithelial"),
    theme = theme(plot.title = element_text(face = "bold"))
  )
}


## 2) 按样本批量保存 3 张图
# obj_all 是你的合并对象
outdir <- paste(OUTDIR,"plots_by_sample",sep="//")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

for (nm in names(objs2)[1:2]) {
  obj <- objs2[[nm]]
  obj <- subset(obj, subset = !is.na(CMS_class) & CMS_class %in% c("CMS2", "CMS4"))

  p1 <- plot_cms(obj, cms_col = "CMS_class")
  p2 <- plot_continuous(obj, feature = "CuproScore_1", title = paste0(nm, " | CuScore"))
  p3 <- plot_ecmcaf_epi(obj, ecmcaf = "RCTD_ECMCAF", epi = "RCTD_Epithelial_cell")

  ggsave(file.path(outdir, paste0(nm, "_1_CMS.pdf")), p1, width = 6.5, height = 6, dpi = 300)
  ggsave(file.path(outdir, paste0(nm, "_2_CuScore.pdf")), p2, width = 6.5, height = 6, dpi = 300)
  ggsave(file.path(outdir, paste0(nm, "_3_ECMCAF_Epi.pdf")), p3, width = 13, height = 6, dpi = 300)
}

#------------------------------------------------统计分析
#--------------sample 1
#-----1    差异分析
d=objs2[[1]]@meta.data |> subset(!is.na(CMS_class) & CMS_class %in% c("CMS2","CMS4") & !is.na(CuproScore_1))
diff1=wilcox.test(CuproScore_1 ~ CMS_class, data=d,alternative="greater")
diff1
p1_diff1=ggplot(d, aes(CMS_class, CuproScore_1, fill=CMS_class)) + 
  geom_boxplot(outlier.shape=NA, width=.25) + 
  geom_jitter(width=.15, alpha=.35, size=.7) + 
  theme_classic() + 
  labs(x=NULL, y="CuproScore_1", title="CuScore: CMS2 vs CMS4")
p1_diff1
#-----2 CMS4 相关性分析
d_CMS4=objs2[[1]]@meta.data |> subset(!is.na(CMS_class) & CMS_class %in% c("CMS4") & !is.na(RCTD_ECMCAF) & !is.na(RCTD_Epithelial_cell))
cor1_CMS4=cor.test(d_CMS4$RCTD_ECMCAF, d_CMS4$RCTD_Epithelial_cell, method="spearman")
cor1_CMS4
p1_corCMS4=ggplot(d_CMS4, aes(RCTD_ECMCAF, RCTD_Epithelial_cell)) + 
  geom_point(alpha=.6, size=.8) + 
  geom_smooth(method="lm", se=FALSE, linewidth=.7, color="black") + 
  facet_wrap(~CMS_class, scales="free") + 
  theme_classic() + 
  labs(x="RCTD_ECMCAF", y="RCTD_Epithelial_cell", title="ECMCAF vs Epithelial (within CMS4)")
p1_corCMS4
#-----3 CMS2 相关性分析
d_CMS2=objs2[[1]]@meta.data |> subset(!is.na(CMS_class) & CMS_class %in% c("CMS2") & !is.na(RCTD_ECMCAF) & !is.na(RCTD_Epithelial_cell))
cor1_CMS2=cor.test(d_CMS2$RCTD_ECMCAF, d_CMS2$RCTD_Epithelial_cell, method="spearman")
cor1_CMS2
p1_corCMS2=ggplot(d_CMS2, aes(RCTD_ECMCAF, RCTD_Epithelial_cell)) + 
  geom_point(alpha=.6, size=.8) + 
  geom_smooth(method="lm", se=FALSE, linewidth=.7, color="black") + 
  theme_classic() + 
  labs(x="RCTD_ECMCAF", y="RCTD_Epithelial_cell", title="ECMCAF vs Epithelial (within CMS2)")
p1_corCMS2
#--------------sample 2
#-----1    差异分析
d=objs2[[2]]@meta.data |> subset(!is.na(CMS_class) & CMS_class %in% c("CMS2","CMS4") & !is.na(CuproScore_1))
diff2=wilcox.test(CuproScore_1 ~ CMS_class, data=d,alternative="greater")
diff2
p2_diff2=ggplot(d, aes(CMS_class, CuproScore_1, fill=CMS_class)) + 
  geom_boxplot(outlier.shape=NA, width=.25) + 
  geom_jitter(width=.15, alpha=.35, size=.7) + 
  theme_classic() + 
  labs(x=NULL, y="CuproScore_1", title="CuScore: CMS2 vs CMS4")
p2_diff2
#-----2 CMS4 相关性分析
d_CMS4=objs2[[2]]@meta.data |> subset(!is.na(CMS_class) & CMS_class %in% c("CMS4") & !is.na(RCTD_ECMCAF) & !is.na(RCTD_Epithelial_cell))
cor2_CMS4=cor.test(d_CMS4$RCTD_ECMCAF, d_CMS4$RCTD_Epithelial_cell, method="spearman")
cor2_CMS4
p2_corCMS4=ggplot(d_CMS4, aes(RCTD_ECMCAF, RCTD_Epithelial_cell)) + 
  geom_point(alpha=.6, size=.8) + 
  geom_smooth(method="lm", se=FALSE, linewidth=.7, color="black") + 
  facet_wrap(~CMS_class, scales="free") + 
  theme_classic() + 
  labs(x="RCTD_ECMCAF", y="RCTD_Epithelial_cell", title="ECMCAF vs Epithelial (within CMS4)")
p2_corCMS4
#-----3 CMS2 相关性分析
d_CMS2=objs2[[2]]@meta.data |> subset(!is.na(CMS_class) & CMS_class %in% c("CMS2") & !is.na(RCTD_ECMCAF) & !is.na(RCTD_Epithelial_cell))
cor2_CMS2=cor.test(d_CMS2$RCTD_ECMCAF, d_CMS2$RCTD_Epithelial_cell, method="spearman")
cor2_CMS2
p2_corCMS2=ggplot(d_CMS2, aes(RCTD_ECMCAF, RCTD_Epithelial_cell)) + 
  geom_point(alpha=.6, size=.8) + 
  geom_smooth(method="lm", se=FALSE, linewidth=.7, color="black") + 
  theme_classic() + 
  labs(x="RCTD_ECMCAF", y="RCTD_Epithelial_cell", title="ECMCAF vs Epithelial (within CMS2)")
p2_corCMS2

(p1_diff1+p1_corCMS4+p1_corCMS2)/(p2_diff2+p2_corCMS4+p2_corCMS2)









