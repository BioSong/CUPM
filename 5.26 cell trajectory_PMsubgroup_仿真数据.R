library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)

set.seed(2025)

# -----------------------------
# 基础参数
# -----------------------------
n_cells <- 28661
groups <- c("k1", "k2", "k3")

# -----------------------------
# 仿真函数
# -----------------------------
gen_method_result <- function(method, acc = 0.7, score_sd = 0.1, score_base = 0.7, 
                              method_type = c("score", "pvalue")) {
  method_type <- match.arg(method_type)
  
  # 随机分配预测label（带一定倾向）
  predicted_label <- sample(groups, n_cells, replace = TRUE, prob = c(acc, (1 - acc)/2, (1 - acc)/2))
  
  # 模拟置信度分数（或p值）
  score <- pmin(1, pmax(0, rnorm(n_cells, mean = score_base, sd = score_sd)))
  
  # 若为p值型方法（置信度=1-p）
  if (method_type == "pvalue") {
    score <- 1 - score
  }
  
  data.frame(
    cell_id = paste0("cell", seq_len(n_cells)),
    predicted_label = predicted_label,
    score = score,
    method = method,
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# 生成五种方法的仿真结果
# -----------------------------
SingleR_df      <- gen_method_result("SingleR",      acc = 0.75, score_sd = 0.12, score_base = 0.7,  method_type = "score")
Seurat_df       <- gen_method_result("Seurat",       acc = 0.72, score_sd = 0.15, score_base = 0.68, method_type = "score")
ClusterRepro_df <- gen_method_result("ClusterRepro", acc = 0.70, score_sd = 0.10, score_base = 0.65, method_type = "pvalue")
scmap_df        <- gen_method_result("scmap",        acc = 0.65, score_sd = 0.18, score_base = 0.6,  method_type = "score")
PseudoBulk_df   <- gen_method_result("PseudoBulk",   acc = 0.80, score_sd = 0.08, score_base = 0.75, method_type = "score")

# 合并
all_methods_df <- rbind(SingleR_df, Seurat_df, ClusterRepro_df, scmap_df, PseudoBulk_df)

# -----------------------------
# 查看结果
# -----------------------------
head(all_methods_df)
table(all_methods_df$method)
#--------------------------------------------
vote_summary <- all_methods_df %>%
  group_by(cell_id, predicted_label) %>%
  summarise(
    n_votes = n(),
    mean_score = mean(score)
  ) %>%
  arrange(cell_id, desc(n_votes), desc(mean_score)) %>%
  group_by(cell_id) %>%
  slice(1) %>%  # 取得票最多且score最高的label
  ungroup() %>%
  rename(final_label = predicted_label,
         vote_count = n_votes)

#-----------------------------
# 2️⃣ 一致性指标 (每个细胞)
#-----------------------------
vote_consistency <- all_methods_df %>%
  group_by(cell_id) %>%
  summarise(
    consistency = max(table(predicted_label)) / n()  # 最多票数 / 总方法数
  )

# 合并
vote_summary <- left_join(vote_summary, vote_consistency, by = "cell_id")

#-----------------------------
# 3️⃣ 一致性热图（方法 × 方法）
#-----------------------------
method_pairs <- expand.grid(unique(all_methods_df$method),
                            unique(all_methods_df$method),
                            stringsAsFactors = FALSE)
colnames(method_pairs) <- c("method1", "method2")

# 计算一致性
agreement_matrix <- sapply(unique(all_methods_df$method), function(m1) {
  sapply(unique(all_methods_df$method), function(m2) {
    df1 <- all_methods_df %>% filter(method == m1) %>% arrange(cell_id)
    df2 <- all_methods_df %>% filter(method == m2) %>% arrange(cell_id)
    mean(df1$predicted_label == df2$predicted_label)
  })
})

# 绘制热图
pheatmap(agreement_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("white", "steelblue"))(50),
         main = "Method Agreement Heatmap")

#-----------------------------
# 4️⃣ 置信度分布（每组的平均score）
#-----------------------------
ggplot(all_methods_df, aes(x = method, y = score, fill = predicted_label)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  labs(title = "Score distribution by method and predicted group",
       y = "Confidence score", x = "") +
  scale_fill_brewer(palette = "Set2")

#-----------------------------
# 5️⃣ 一致性分布（每个细胞）
#-----------------------------
ggplot(vote_summary, aes(x = consistency)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "white") +
  theme_classic() +
  labs(title = "Distribution of classification consistency",
       x = "Consistency (agreement ratio)", y = "Cell count")
