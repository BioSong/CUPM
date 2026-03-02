# Packages
library(dplyr)
library(ggplot2)

# ---- Example data (replace with yours if needed) ----
df <- tibble::tribble(
  ~cell,  ~protein, ~group,  ~rep, ~fold,
  "SW480","FDX1",   "Blank", 1,    1.53,
  "SW480","FDX1",   "Blank", 2,    0.95,
  "SW480","FDX1",   "Blank", 3,    1.03,
  "SW480","FDX1",   "NC",    1,    1.00,
  "SW480","FDX1",   "NC",    2,    1.00,
  "SW480","FDX1",   "NC",    3,    1.00,
  "SW480","FDX1",   "FOXO1", 1,    2.19,
  "SW480","FDX1",   "FOXO1", 2,    1.40,
  "SW480","FDX1",   "FOXO1", 3,    2.58,
  "SW480","FDX1",   "SOX11", 1,    2.11,
  "SW480","FDX1",   "SOX11", 2,    1.50,
  "SW480","FDX1",   "SOX11", 3,    2.37,

  "SW480","DLAT",   "Blank", 1,    0.80,
  "SW480","DLAT",   "Blank", 2,    0.86,
  "SW480","DLAT",   "Blank", 3,    1.49,
  "SW480","DLAT",   "NC",    1,    1.00,
  "SW480","DLAT",   "NC",    2,    1.00,
  "SW480","DLAT",   "NC",    3,    1.00,
  "SW480","DLAT",   "FOXO1", 1,    1.46,
  "SW480","DLAT",   "FOXO1", 2,    1.30,
  "SW480","DLAT",   "FOXO1", 3,    1.87,
  "SW480","DLAT",   "SOX11", 1,    1.43,
  "SW480","DLAT",   "SOX11", 2,    1.29,
  "SW480","DLAT",   "SOX11", 3,    1.95,

  "LOVO", "FDX1",   "Blank", 1,    0.90,
  "LOVO", "FDX1",   "Blank", 2,    0.80,
  "LOVO", "FDX1",   "Blank", 3,    0.63,
  "LOVO", "FDX1",   "NC",    1,    1.00,
  "LOVO", "FDX1",   "NC",    2,    1.00,
  "LOVO", "FDX1",   "NC",    3,    1.00,
  "LOVO", "FDX1",   "FOXO1", 1,    1.02,
  "LOVO", "FDX1",   "FOXO1", 2,    1.03,
  "LOVO", "FDX1",   "FOXO1", 3,    1.04,
  "LOVO", "FDX1",   "SOX11", 1,    1.03,
  "LOVO", "FDX1",   "SOX11", 2,    1.03,
  "LOVO", "FDX1",   "SOX11", 3,    1.03,

  "LOVO", "DLAT",   "Blank", 1,    0.82,
  "LOVO", "DLAT",   "Blank", 2,    1.13,
  "LOVO", "DLAT",   "Blank", 3,    0.83,
  "LOVO", "DLAT",   "NC",    1,    1.00,
  "LOVO", "DLAT",   "NC",    2,    1.00,
  "LOVO", "DLAT",   "NC",    3,    1.00,
  "LOVO", "DLAT",   "FOXO1", 1,    1.93,
  "LOVO", "DLAT",   "FOXO1", 2,    1.85,
  "LOVO", "DLAT",   "FOXO1", 3,    0.96,
  "LOVO", "DLAT",   "SOX11", 1,    1.70,
  "LOVO", "DLAT",   "SOX11", 2,    1.69,
  "LOVO", "DLAT",   "SOX11", 3,    1.00
)

# factor order
df <- df %>%
  mutate(
    group   = factor(group, levels = c("Blank", "NC", "FOXO1", "SOX11")),
    cell    = factor(cell, levels = c("SW480", "LOVO")),
    protein = factor(protein, levels = c("FDX1", "DLAT"))
  )

# summary for bars/error bars (SEM)
sumdf <- df %>%
  group_by(cell, protein, group) %>%
  summarise(
    mean = mean(fold, na.rm = TRUE),
    sd   = sd(fold, na.rm = TRUE),
    n    = dplyr::n(),
    sem  = sd / sqrt(n),
    .groups = "drop"
  )

# ---- Plot: bar (mean) + errorbar (SEM) + points (replicates) ----
set.seed(1)
p <- ggplot(sumdf, aes(x = group, y = mean, fill = group)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                width = 0.18, linewidth = 0.5) +
  geom_point(data = df,
             aes(x = group, y = fold),
             inherit.aes = FALSE,
             position = position_jitter(width = 0.12, height = 0),
             size = 2.2, shape = 21, stroke = 0.4, fill = "white", color = "black") +
  facet_grid(cell ~ protein) +
  scale_fill_manual(values = c("Blank" = "#BDBDBD", "NC" = "#4D4D4D",
                               "FOXO1" = "#2C7FB8", "SOX11" = "#7FCDBB")) +
  labs(x = NULL, y = "Fold change (NC = 1.0)") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    strip.text = element_text(face = "bold")
  )

p

# Save
ggsave("barplot_foldchange_meanSEM_points.pdf", p, width = 7.5, height = 5.5)
ggsave("barplot_foldchange_meanSEM_points.png", p, width = 7.5, height = 5.5, dpi = 300)
