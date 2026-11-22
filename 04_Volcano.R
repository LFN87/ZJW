# 1. 加载包
library(ggplot2)
library(ggrepel)
library(dplyr)

# 2. 设置路径并读取数据
setwd("C:\\") 
data_table <- read.table("Corrected_DE.txt", header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

# 3. 定义阈值和显著性标签
filter_logFC <- 1
filter_adjPVal <- 0.05

data_table <- data_table %>%
  mutate(Sig = case_when(
    adj.P.Val < filter_adjPVal & logFC > filter_logFC ~ "Up",
    adj.P.Val < filter_adjPVal & logFC < -filter_logFC ~ "Down",
    TRUE ~ "Not"
  ))

# 4. 准备要添加标签的基因（上调前5，下调前5，外加ICAM1）
top_up <- data_table %>% filter(Sig == "Up") %>% arrange(adj.P.Val) %>% head(5)
top_down <- data_table %>% filter(Sig == "Down") %>% arrange(adj.P.Val) %>% head(5)
icam1_gene <- data_table %>% filter(id == "ICAM1")

# 合并标签数据并去重
label_data <- bind_rows(top_up, top_down, icam1_gene) %>% distinct(id, .keep_all = TRUE)

# 5. 开始画图
p <- ggplot(data_table, aes(x = logFC, y = -log10(adj.P.Val))) +
  # 画散点
  geom_point(aes(color = Sig), alpha = 0.8, size = 1.5) +
  # 设置颜色
  scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5", "Not" = "grey85")) +
  # 添加阈值虚线
  geom_vline(xintercept = c(-filter_logFC, filter_logFC), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(filter_adjPVal), linetype = "dashed", color = "black", alpha = 0.5) +
  # 添加基因标签（带引导线）
  geom_text_repel(data = label_data, aes(label = id), 
                  size = 3.5, 
                  box.padding = 0.5, 
                  point.padding = 0.3, 
                  max.overlaps = Inf,
                  color = "black", fontface = "bold",
                  min.segment.length = 0,    
                  segment.color = "grey30",  
                  segment.size = 0.5) +
  # 【新增：强制设置 X 轴刻度为 -4, -2, 0, 2, 4, 6, 8】
  scale_x_continuous(breaks = seq(-4, 8, by = 2)) +
  # 图表主题美化
  theme_bw() +
  labs(
    title = "Volcano Plot of Corrected DEGs",
    x = expression(log[2](Fold~Change)),
    y = expression(-log[10](adj.~P~value))
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = element_blank(),
    legend.position = "top",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    # 去掉背景细碎网格线，画面更干净
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90")
  )

# 6. 保存图片（PDF和PNG各存一份）
pdf("Volcano_Final_Submit.pdf", width = 7, height = 6)
print(p)
dev.off()

png("Volcano_Final_Submit.png", width = 7, height = 6, units = "in", res = 300)
print(p)
dev.off()

cat("✅ 终查看 Volcano_Final_Submit.png\n")