
library(ggplot2)
library(ggpubr)

# 1. 我们重新定义绝对干净的分组和批次标签
batch_info <- c(rep("GSE233973", 36), rep("GSE60427", 32))
group_info <- c(
  rep("Gastritis", 27), rep("Control", 9),
  rep("Control", 4), rep("Gastritis", 12), rep("Control", 4), rep("Gastritis", 12)
)

# 看看前几个名字对不对
print(table(batch_info))

# 2. 画左图（校正前）
pca_before_res <- prcomp(t(combined_mat), scale. = TRUE)
df_before <- data.frame(pca_before_res$x, Batch = batch_info, Group = group_info)
var_before <- round(pca_before_res$sdev^2 / sum(pca_before_res$sdev^2) * 100, 2)

p_before <- ggplot(df_before, aes(x = PC1, y = PC2, color = Batch, shape = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(color = Batch), level = 0.95, linetype = 2, size = 1) +
  scale_color_manual(values = c("GSE233973" = "red", "GSE60427" = "blue")) +
  scale_shape_manual(values = c("Control" = 16, "Gastritis" = 17)) +
  labs(title = "Before ComBat Correction", 
       x = paste0("PC1 (", var_before[1], "%)"), 
       y = paste0("PC2 (", var_before[2], "%)")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 3. 画右图（校正后）
pca_after_res <- prcomp(t(corrected_data), scale. = TRUE)
df_after <- data.frame(pca_after_res$x, Batch = batch_info, Group = group_info)
var_after <- round(pca_after_res$sdev^2 / sum(pca_after_res$sdev^2) * 100, 2)

p_after <- ggplot(df_after, aes(x = PC1, y = PC2, color = Batch, shape = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(color = Batch), level = 0.95, linetype = 2, size = 1) +
  scale_color_manual(values = c("GSE233973" = "red", "GSE60427" = "blue")) +
  scale_shape_manual(values = c("Control" = 16, "Gastritis" = 17)) +
  labs(title = "After ComBat Correction", 
       x = paste0("PC1 (", var_after[1], "%)"), 
       y = paste0("PC2 (", var_after[2], "%)")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 4. 拼图并保存
final_pca <- ggarrange(p_before, p_after, ncol = 2, common.legend = TRUE, legend = "right")

ggsave("PCA_Final_Clean.png", plot = final_pca, width = 10, height = 5, dpi = 300)
cat("✅ 图已保存为 PCA_Final_Clean.png，\n")




