library(DESeq2)
library(openxlsx)
library(dplyr)

# 1. 读取数据
rawdata <- read.xlsx("1.rawdata.xlsx")

# 2. 检查gene_name列名
colnames(rawdata)[1] <- "gene_name"

# 3. 去掉gene_name缺失
rawdata <- rawdata[!is.na(rawdata$gene_name), ]

# 4. 若有重复gene_name，合并取均值（一般不需要）
rawdata <- rawdata %>%
  group_by(gene_name) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  as.data.frame()

# 5. 提取基因名与计数矩阵
gene_names <- rawdata$gene_name
counts <- rawdata[, -1]
rownames(counts) <- gene_names

# 6. 转为数值型
counts <- as.matrix(counts)
mode(counts) <- "numeric"

# 7. 转整数型（DESeq2需要整数）
counts <- round(counts)
mode(counts) <- "integer"

# 8. 过滤低表达基因（平均表达大于1）
counts <- counts[rowMeans(counts) > 1, ]

# 9. 构建分组信息
sample_table <- data.frame(
  condition = factor(c(rep("WT", 3), rep("KOT", 3)))
)

# 10. 检查行列对应
stopifnot(ncol(counts) == nrow(sample_table))
counts <- counts[complete.cases(counts), ]

# 11. 构建DESeq2对象
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_table, design = ~condition)

# 12. 差异分析
dds <- DESeq(dds)
res <- results(dds)
all_diff <- as.data.frame(res)
all_diff$gene_id <- rownames(all_diff)
all_diff <- all_diff[, c(ncol(all_diff), 1:(ncol(all_diff)-1))]
all_diff <- na.omit(all_diff)

# 13. 重命名
colnames(all_diff)[colnames(all_diff) == "log2FoldChange"] <- "logFC"
colnames(all_diff)[colnames(all_diff) == "pvalue"] <- "P.Value"

# 14. 导出
write.xlsx(all_diff, file = "2_diffALL.xlsx")



# 设置调整后的p值
adjust_p <- 0.05
log_fold_change <- 1
# 筛选差异显著的基因
diff_sig <- all_diff[with(all_diff, (abs(logFC) > log_fold_change &
                                       P.Value < adjust_p)), ]
# 保存差异显著的基因到文件
write.xlsx(diff_sig, "3_diff.xlsx")

# 保存上调基因到文件
diff_up <- all_diff[with(all_diff, (logFC > log_fold_change &
                                      P.Value < adjust_p)), ]
write.xlsx(diff_up, "4_up.xlsx")
# 保存下调基因到文件
diff_down <- all_diff[with(all_diff, (logFC < -log_fold_change &
                                        P.Value < adjust_p)), ]
write.xlsx(diff_down, "5_down.xlsx")





# rm(list = ls())
library(tidyverse)
library(readxl)
library(openxlsx)
# setwd("F:\\postgraduate\\PNH_plan\\2.SmallArticle\\9.scRNA_result\\粒细胞簇差异分析")

# 读取差异分析结果
# data_miRNA <- read.table("13.diffALL.xls",header = T,sep = "\t")
data_miRNA <- read.xlsx("2_diffALL.xlsx")
head(data_miRNA)

#筛选条件
Log2FC <- 1
pvalue <- 0.05

# 替换列的名字
colnames(data_miRNA)[colnames(data_miRNA) == "logFC"] <- "Log2FC"
colnames(data_miRNA)[colnames(data_miRNA) == "P.Value"] <- "pvalue"


data_miRNA$regulated <- "Normal"
data_miRNA$regulated[data_miRNA$Log2FC > Log2FC & data_miRNA$pvalue< pvalue ] <- "Up"
data_miRNA$regulated[data_miRNA$Log2FC < -Log2FC & data_miRNA$pvalue < pvalue ] <- "Down"
table(data_miRNA$regulated)

## 颜色设置：使用的微信截图获取的文章图片配色，当然也有一些包可以根据输入的图片提取其中的颜色
## 这里颜色比较少就用了微信截图获取，方便
# 设置点图的内圈填充色
data_miRNA$color <- ""
data_miRNA$color[data_miRNA$regulated=="Normal"] <- "#ABABA6"
data_miRNA$color[data_miRNA$regulated=="Up"] <- "#FA8260"
data_miRNA$color[data_miRNA$regulated=="Down"] <- "#4D8FD1"
table(data_miRNA$color)

# 设置点的边圈颜色，比内部填充色颜色深一些
data_miRNA$color1 <- ""
data_miRNA$color1[data_miRNA$regulated=="Normal"] <- "#7B7B73"
data_miRNA$color1[data_miRNA$regulated=="Up"] <- "#C85A40"
data_miRNA$color1[data_miRNA$regulated=="Down"] <- "#315E94"

# 计算上调和下调的百分比
downregulated_count <- sum(data_miRNA$color == "#4D8FD1")
upregulated_count <- sum(data_miRNA$color == "#FA8260")
total_count <- nrow(data_miRNA)
downregulated_percentage <- (downregulated_count / total_count) * 100
upregulated_percentage <- (upregulated_count / total_count) * 100

min_value <- min(data_miRNA$pvalue, na.rm = TRUE)
y_max <- floor((-log10(min_value))+1)



# 替换0为很小的值
data_miRNA$pvalue[data_miRNA$pvalue == 0] <- 1e-300
min_value <- min(data_miRNA$pvalue, na.rm = TRUE)
y_max <- floor((-log10(min_value)) + 1)




# 设置步长 y_by
if (y_max <= 10) {
  y_by <- 1
} else if (y_max <= 50) {
  y_by <- 5
} else if (y_max <= 100) {
  y_by <- 10
} else if (y_max <= 500) {
  y_by <- 50
} else {
  y_by <- 100
}



min_logfc <- min(data_miRNA$Log2FC, na.rm = TRUE)
max_logfc <- max(data_miRNA$Log2FC, na.rm = TRUE)

x_min <- ceiling(min_logfc - 1)
x_max <- floor(max_logfc + 1)





pdf(file = "6_volcano.pdf", width = 7.5, height = 6.5)
par(bty="n", mgp = c(1.5,.33,0), mar=c(4.5,3.5,3,3)+.1, las=1, tcl=-.3)

plot(data_miRNA$Log2FC, y = -log10(data_miRNA$pvalue), 
     xlab = "", ylab = "-log P-value", 
     col = data_miRNA$color1, pch = 21, bg = data_miRNA$color,
     cex = 1, lwd = 0.5, yaxt = "n", xaxt = "n", xlim = c(x_min, x_max))
# xlim：x 轴的显示范围

axis(side = 2, at = seq(0, y_max, by = y_by), las = 1, lwd = 2.5)#y轴
axis(side = 1, at = seq(x_min, x_max, by = 1), las = 1, lwd = 2.5)#x轴

mtext(expression(paste(Log[2], " foldchange")), 
      cex = 1.2, col = "black", side = 1, line = 2)

abline(h = -log10(0.05), col = "grey60", lwd = 2.5, lty = 3)
abline(v = c(-1, 1), col = "grey60", lwd = 2.5, lty = 3)
abline(v = 0, col = "black", lwd = 2.5, lty = 1)

# 使用segments函数添加线段
segments(-Log2FC/4, seq(0, y_max, by = y_by), Log2FC/4, seq(0, y_max, by = y_by), 
         col = "black", lwd = 2.5)


# 添加向右的箭头
arrows(x0 = c(Log2FC/2, -Log2FC/2), y0 = y_max*0.9, x1 = c(x_max, x_min), y1 =  y_max*0.9,
       length = 0.1, col = c("#FA8260", "#4D8FD1"), lwd = 4.5)


# 添加文本
text(x = c(Log2FC+2,-Log2FC-2), 
     y = rep(y_max*0.9, 2), 
     labels = c("Higher in PNH", "Higher in control"),
     pos = 3, # 左对齐和右对齐
     col = c("#FA8260", "#4D8FD1"), 
     cex = 1.1, 
     font = 3)

# 百分比
text(c(x_max-0.2,x_min+0.2), rep(y_max*0.9, 1), 
     c(paste(round(upregulated_percentage, 1), "%"), paste(round(downregulated_percentage, 1), "%")),
     pos = 3, col = c("#FA8260", "#4D8FD1"), cex = 1.1, font = 3)

title(main = "mRNA")

dev.off()


