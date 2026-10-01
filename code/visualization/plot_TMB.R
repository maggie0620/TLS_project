
######################################## 計算兩組的 TMB 是否有顯著差異 ##########################################

filtermaf_paad <- read.maf("./data/PAAD/maf/filtermaf_maftools.maf")
filtermaf_chol <- read.maf("./data/CHOL/maf/filtermaf_maftools.maf")
filtermaf_stad <- read.maf("./data/STAD/maf/filtermaf_maftools.maf")

## TMB 沒有顯著區別
TMB_paad <- tmb(maf = filtermaf_paad, captureSize = 35.8, plotType = NA)$total_perMB
TMB_chol <- tmb(maf = filtermaf_chol, captureSize = 35.8, plotType = NA)$total_perMB
TMB_stad <- tmb(maf = filtermaf_stad, captureSize = 35.8, plotType = NA)$total_perMB
TMB_pvalue1 <- wilcox.test(TMB_paad, TMB_chol, alternative = "two.sided", exact = FALSE)$p.value
TMB_pvalue2 <- wilcox.test(TMB_paad, TMB_stad, alternative = "two.sided", exact = FALSE)$p.value
TMB_pvalue3 <- wilcox.test(TMB_chol, TMB_stad, alternative = "two.sided", exact = FALSE)$p.value

# 準備資料
tmb_df <- data.frame(TMB = c(TMB_paad, TMB_chol, TMB_stad),
                     Group = c(rep("PAAD", length(TMB_paad)), rep("CHOL", length(TMB_chol)), rep("STAD", length(TMB_stad))))
tapply(tmb_df$TMB, tmb_df$Group, median)

# 決定 y 軸上限（例如排除最極端的 outlier）
library(ggplot2)
y_limit <- quantile(tmb_df$TMB, 0.99) # 調整這個數字

tiff(filename = "~/TLS_meth/tiff_figure/TMB_cohort.tiff", width = 6, height = 6, res = 300, 
     compression = "lzw", units = "in", bg = "white")
ggplot(tmb_df, aes(x = Group, y = TMB, fill = Group)) +
  stat_boxplot(geom = "errorbar", width = 0.1) + 
  geom_boxplot(outlier.shape = 1, outlier.size = 1.5, width = 0.5) +
  annotate("text",
           x = 1.5,
           y = y_limit * 0.89,             # p value 位置跟著調整
           label = paste0("p = ", format(TMB_pvalue1, digits = 3)),
           size = 4.5) +
  annotate("segment",
           x = 1, xend = 2,
           y = y_limit * 0.85, yend = y_limit* 0.85) +  # 畫橫線連接兩組
  annotate("text",
           x = 2.5,
           y = y_limit * 0.94,             # p value 位置跟著調整
           label = paste0("p = ", format(TMB_pvalue2, digits = 3)),
           size = 4.5) +
  annotate("segment",
           x = 2, xend = 3,
           y = y_limit * 0.9, yend = y_limit* 0.9) +  # 畫橫線連接兩組
  annotate("text",
           x = 2,
           y = y_limit * 1.02,             # p value 位置跟著調整
           label = paste0("p = ", format(TMB_pvalue3, digits = 3)),
           size = 4.5) +
  annotate("segment",
           x = 1, xend = 3,
           y = y_limit * 0.98, yend = y_limit* 0.98) +  # 畫橫線連接兩組
  scale_fill_manual(values = c("PAAD" = "#7E5A92", "CHOL" = "#4FA6A6", "STAD" = "#ffcf33")) +
  coord_cartesian(ylim = c(0, y_limit)) +  # 限制顯示範圍但保留 outlier
  labs(x = NULL, y = "TMB (mutations/Mb)") +
  theme_bw() +                            
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()
