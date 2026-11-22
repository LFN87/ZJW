library(limma)


# 2. 构建设计矩阵
design <- model.matrix(~0 + group_labels)
colnames(design) <- c("Control", "Gastritis") # 这里的顺序要和 levels 对应
fit <- lmFit(corrected_data, design)

# 3. 构建对比矩阵
cont.matrix <- makeContrasts(Gastritis - Control, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# 4. 提取结果
# ⚠️ 这里调整回原稿的阈值 logFC = 1
logFCfilter <- 1 
adj.P.Val.Filter <- 0.05

allDiff <- topTable(fit2, adjust = 'fdr', number = Inf)

# 输出完整差异结果
allDiffOut <- rbind(id=colnames(allDiff), allDiff)
write.table(allDiffOut, file = "Corrected_DE.txt", sep = "\t", quote = F, col.names = F)

# 输出筛选后的差异基因
diffSig <- allDiff[with(allDiff, (abs(logFC) > logFCfilter & adj.P.Val < adj.P.Val.Filter)), ]
diffSigOut <- rbind(id=colnames(diffSig), diffSig)
write.table(diffSigOut, file = "Corrected_DEGs.txt", sep = "\t", quote = F, col.names = F)

