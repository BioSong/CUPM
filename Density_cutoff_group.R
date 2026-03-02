get_density_cutoff_group <- function(score_group1, score_group2,
                                     name_group1 = "Group1",
                                     name_group2 = "Group2",
                                     plot = TRUE) {
  # 合并数据
  df <- data.frame(
    score = c(score_group1, score_group2),
    source_group = c(rep(name_group1, length(score_group1)),
                     rep(name_group2, length(score_group2)))
  )

  # 计算密度
  dens1 <- density(score_group1)
  dens2 <- density(score_group2)

  # 估算交点（密度差为0的点）
  common_x <- seq(max(min(dens1$x), min(dens2$x)),
                  min(max(dens1$x), max(dens2$x)), length.out = 1000)
  y1_interp <- approx(dens1$x, dens1$y, xout = common_x)$y
  y2_interp <- approx(dens2$x, dens2$y, xout = common_x)$y
  diff_y <- y1_interp - y2_interp
  idx <- which(diff_y[-1] * diff_y[-length(diff_y)] < 0)  # 查找变号点（交点）

  if (length(idx) == 0) {
    warning("No intersection found between densities.")
    cutoff <- median(df$score)
  } else {
    cutoff <- common_x[idx[which.min(abs(diff_y[idx]))]]
  }

  # 分组
  df$cupro_group <- ifelse(df$score >= cutoff, "High", "Low")

  # 画图（可选）
  if (plot) {
    library(ggplot2)
    p <- ggplot(df, aes(x = score, color = source_group)) +
      geom_density(size = 1) +
      geom_vline(xintercept = cutoff, linetype = "dashed", color = "black") +
      annotate("text", x = cutoff, y = 0, label = paste0("Cutoff: ", round(cutoff, 3)),
               hjust = -0.1, vjust = -1, size = 4) +
      labs(title = "Density plot with cutoff", x = "Cuprotosis score", y = "Density") +
      theme_minimal()
    print(p)
  }

  # 返回数据和 cutoff
  return(list(cutoff = cutoff,
              df = df,
              plot = if (plot) p else NULL))
}
