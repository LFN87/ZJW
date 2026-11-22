# 1. 安装并加载包
# install.packages("pheatmap")
library(pheatmap)

# 2. 设置工作目录（假设你还在这个目录下）
setwd("C:")

# 3. 准备画图数据
# 选出 P 值最小的前 50 个基因
top50_genes <- head(rownames(allDiff)[order(allDiff$adj.P.Val)], 50)

# 从去批次后的矩阵中提取这 50 个基因的表达量
heatmap_data <- corrected_data[top50_genes, ]

# 4. 准备样本注释信息（给热图顶部加彩色条）
annotation_col <- data.frame(Group = group_labels)
rownames(annotation_col) <- colnames(heatmap_data)

# 自定义颜色（青色代表Control，粉色代表Gastritis，和你原图保持一致）
ann_colors <- list(Group = c(Control = "#00BFC4", Gastritis = "#F8766D"))

# 5. 绘制热图
pheatmap(
  heatmap_data,
  scale = "row",              # 按行标准化（Z-score），让红蓝对比更明显
  cluster_rows = TRUE,        # 基因聚类
  cluster_cols = TRUE,        # 样本聚类（看看Control和Gastritis是否各自抱团）
  show_rownames = FALSE,      # 隐藏基因名（50个名字太挤了，文章里通常不显示或只显示部分）
  show_colnames = FALSE,      # 隐藏样本名
  annotation_col = annotation_col, # 顶部注释条
  annotation_colors = ann_colors,  # 注释条颜色
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100), # 蓝-白-红渐变
  main = "Top 50 DEGs Heatmap (After Correction)",
  fontsize = 12,
  filename = "Heatmap_Corrected.png", # 直接保存图片
  width = 7,
  height = 8
)

cat("✅ 完看 Heatmap_Corrected.png\n")