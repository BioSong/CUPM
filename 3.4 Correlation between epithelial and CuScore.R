library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

#---------scRNAseq annotated with cell type
load("./CUPM/data/scRNAseq/scaled_GSE183916.RData")
exp1=data1@assays$SCT@scale.data
loc=data1@meta.data$loc#-----------primary, metastasis
cuScore1=data1@meta.data$CuScore
pat_id=data1@meta.data$pat_id
celltype=data1@meta.data$celltype
celltype[celltype=="Endothelial_cell"]="Stromal_cell"
data1@meta.data$celltype=celltype
sample_CuScore=data.frame(loc=loc,CuScore=cuScore1$Final_Score)
#-----------------------------------------
# Count number of epithelial cells per sample
cell_counts <- data1@meta.data %>%
  group_by(loc, celltype) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(loc) %>%
  mutate(fraction = n / sum(n))

# Extract epithelial fraction per sample
epi_fraction <- cell_counts %>%
  filter(celltype == "Epithelial_cell") %>%
  select(loc, fraction)

# Merge with bulk/sample-level CuScore (assume you have a dataframe `sample_CuScore` with columns: loc, CuScore)
df_corr <- left_join(epi_fraction, sample_CuScore, by = "loc")

# Spearman correlation
cor_test <- cor.test(df_corr$fraction, df_corr$CuScore, method = "spearman")
print(cor_test)

# Scatter plot
ggplot(df_corr, aes(x = fraction, y = CuScore, color = loc)) +
  geom_point(size=3) +
  geom_smooth(method="lm", se=TRUE, color="black") +
  theme_bw() +
  labs(x = "Epithelial cell fraction", y = "CuScore", title = "Correlation between epithelial fraction and CuScore")
# 3. Per-cell CuScore comparison in epithelial cells
# -----------------------------
epi_cells <- subset(data1, subset = celltype == "Epithelial_cell")

# Boxplot per sample type
ggplot(epi_cells@meta.data, aes(x = loc, y = CuScore$Final_Score, fill = loc)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha=0.5) +
  theme_bw() +
  labs(x="Sample Type", y="CuScore in Epithelial Cells")

# Wilcoxon test
wilcox.test(CuScore ~ sample_type, data = epi_cells@meta.data)
#
# 4. Contribution decomposition: per cell type to total CuScore
# -----------------------------
# Sum CuScore per cell type per sample
celltype_Cu <- data1@meta.data %>%
  group_by(loc, celltype) %>%
  summarise(total_CuScore = sum(CuScore$Final_Score), .groups='drop')

# Compute fraction of total CuScore contributed by each cell type
celltype_Cu <- celltype_Cu %>%
  group_by(loc) %>%
  mutate(fraction_contribution = total_CuScore / sum(total_CuScore))

# Stacked barplot
ggplot(celltype_Cu, aes(x = loc, y = fraction_contribution, fill = celltype)) +
  geom_bar(stat="identity") +
  theme_bw() +
  labs(x="Sample Type", y="Fractional contribution to total CuScore", fill="Cell Type")

# 5. Compute stromal and myeloid fractions per sample
# -----------------------------
# Extract stromal fraction
stromal_fraction <- cell_counts %>%
  filter(celltype == "Stromal_cell") %>%
  select(loc, fraction) %>%
  rename(stromal_fraction = fraction)

# Extract myeloid fraction
myeloid_fraction <- cell_counts %>%
  filter(celltype == "Myeloid_cell") %>%
  select(loc, fraction) %>%
  rename(myeloid_fraction = fraction)


# Combine all fractions with sample-level CuScore
df_reg <- sample_CuScore %>%
  left_join(epi_fraction, by="loc") %>%
  left_join(stromal_fraction, by="loc") %>%
  left_join(myeloid_fraction, by="loc") 

# -----------------------------
# 6. Multiple linear regression
# -----------------------------
reg_model1 <- lm(CuScore ~ fraction + stromal_fraction  , data = df_reg)
summary(reg_model1)
reg_model2 <- lm(CuScore ~ stromal_fraction + fraction, data = df_reg)
summary(reg_model2)
reg_model3 <- lm(CuScore ~ myeloid_fraction + fraction + stromal_fraction , data = df_reg)
summary(reg_model3)

# The coefficient of fraction (epithelial) shows its independent contribution
# p-value indicates statistical significance
library(compositions)

df_clr <- as.data.frame(clr(df_reg[, c("fraction", "stromal_fraction", "myeloid_fraction")]))

df_clr$CuScore <- df_reg$CuScore

summary(lm(CuScore ~ fraction, data = df_clr))
summary(lm(CuScore ~ stromal_fraction, data = df_clr))
summary(lm(CuScore ~ myeloid_fraction, data = df_clr))
#----------forest
library(ggplot2)
library(dplyr)
library(scales)

# 整理数据
df_forest <- data.frame(
  celltype = c("Epithelial", "Stromal", "Myeloid"),
  estimate = c(0.006280, 0.13101, 0.010208),
  se = c(0.002639, 0.05506, 0.004290)
)

# 计算 95% CI
df_forest <- df_forest %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    label = paste0(round(estimate, 3), " [", round(lower,3), ", ", round(upper,3), "]")
  )

# 设置 celltype 顺序（横向显示）
df_forest$celltype <- factor(df_forest$celltype, levels = rev(df_forest$celltype))

# 绘制 publication-ready forest plot
p <- ggplot(df_forest, aes(x = celltype, y = estimate)) +
  geom_point(aes(color = celltype), size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = celltype), width = 0.2, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_text(aes(label = label), hjust = -0.1, size = 4) +
  coord_flip() +
  scale_color_manual(values = c("Epithelial" = "#1b9e77",
                                "Stromal" = "#d95f02",
                                "Myeloid" = "#7570b3")) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "",
    y = "CLR regression effect size (β)",
    title = "Association between cell-type proportions and CuScore"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))

# 打印
print(p)

library(ggplot2)
library(dplyr)

# 整理数据
df_forest <- data.frame(
  celltype = c("Epithelial", "Stromal", "Myeloid"),
  estimate = c(0.006280, 0.13101, 0.010208),
  se = c(0.002639, 0.05506, 0.004290),
  pval = c(0.0173, 0.0173, 0.0173)  # 根据你提供的 P 值
)

# 计算 95% CI
df_forest <- df_forest %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    label = paste0(round(estimate, 3), " [", round(lower,3), ", ", round(upper,3), "]"),
    sig = ifelse(pval < 0.001, "***",
                 ifelse(pval < 0.01, "**",
                        ifelse(pval < 0.05, "*", "")))
  )

# 设置 celltype 顺序（横向显示）
df_forest$celltype <- factor(df_forest$celltype, levels = rev(df_forest$celltype))

# 绘制 publication-ready forest plot
p <- ggplot(df_forest, aes(x = celltype, y = estimate)) +
  geom_point(aes(color = celltype), size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = celltype), width = 0.2, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_text(aes(label = label), hjust = -0.1, size = 4) +
  geom_text(aes(y = upper + 0.02, label = sig), size = 5, color = "black") +  # 星号标注在 CI 上方
  coord_flip() +
  scale_color_manual(values = c("Epithelial" = "#1b9e77",
                                "Stromal" = "#d95f02",
                                "Myeloid" = "#7570b3")) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "",
    y = "CLR regression effect size (β)",
    title = "Association between cell-type proportions and CuScore"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))  # 留空间给星号

print(p)
