library(limma)
library(sva)


# ========== 2. 路径设置 ==========
work_dir <- "C:\\" 
out_file <- "merge_corrected.txt"
setwd(work_dir)

# 安全读取函数
safe_read <- function(file) {
  read.table(file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}

# ========== 3. 指定读取顺序并处理数据 ==========
txt_files <- c("GSE233973.txt", "GSE60427.txt") 

data_list <- lapply(txt_files, function(f) {
  cat("正在处理文件:", f, "\n")
  dat <- safe_read(f)
  # 第一列设为行名（基因名）
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  mat <- as.matrix(dat)
  # 处理重复基因（取平均值）
  mat <- avereps(mat, ID = rownames(mat))
  return(mat)
})

# ========== 4. 找共有基因 & 横向合并 ==========
gene_lists <- lapply(data_list, rownames)
common_genes <- Reduce(intersect, gene_lists)
cat("共有基因数量：", length(common_genes), "\n")

# 只保留共有基因
data_list <- lapply(data_list, function(m) m[common_genes, ])
combined_mat <- do.call(cbind, data_list)

# ========== 5. 构建批次和分组标签 ==========
# 批次：GSE233973 为 Batch1，GSE60427 为 Batch2
batch_labels <- c(rep("Batch1", 36), rep("Batch2", 32))
batch_labels <- as.factor(batch_labels)

# 分组
group_labels <- c(
  # GSE233973 (GSM7440xxx)：前27个是treat，后9个是control
  rep("Gastritis", 27), rep("Control", 9),
  
  # GSE60427 (GSM1479xxx)：4个control，12个treat，4个control，12个treat
  rep("Control", 4), rep("Gastritis", 12), rep("Control", 4), rep("Gastritis", 12)
)
group_labels <- factor(group_labels, levels = c("Control", "Gastritis"))


print("检查前 10 个样本名和对应的分组是否一致：")
print(data.frame(Sample = head(colnames(combined_mat), 10), Group = head(group_labels, 10)))


# 检查样本总数是否对齐 (27+9+4+12+4+12 = 68)
if(length(group_labels) != ncol(combined_mat)) {
  stop("❌ 错误：分组标签数量与表达矩阵样本数量不一致，请检查！")
}

# ========== 6. Log2 转换 ==========
# 如果数据最大值大于 50，说明是原始荧光值，自动做 log2(x+1) 转换
if (max(combined_mat, na.rm = TRUE) > 50) {
  cat("检测到数据跨度大，正在执行 log2(x + 1) 转换...\n")
  combined_mat <- log2(combined_mat + 1)
} else {
  cat("数据无需 log2 转换\n")
}

# ========== 7. Quantile 标准化 ==========
combined_mat <- normalizeBetweenArrays(combined_mat)

# ========== 8. ComBat 批次校正（保留生物学分组） ==========
cat("正在执行 ComBat 去批次效应...\n")
# 构建生物学模型矩阵（告诉 ComBat 保留疾病差异）
mod <- model.matrix(~ group_labels)

corrected_data <- ComBat(
  dat = combined_mat,
  batch = batch_labels,
  mod = mod,         # 核心参数：保留真实的 Control vs Gastritis 差异
  par.prior = TRUE
)

# ========== 9. 准备输出 ==========
output_data <- data.frame(
  geneNames = rownames(corrected_data),
  corrected_data,
  check.names = FALSE
)

write.table(
  output_data,
  file = out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("✅ 全部完成！已去除批次效应，结果保存到：", out_file, "\n")


