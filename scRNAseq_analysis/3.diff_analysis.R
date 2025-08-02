library(Seurat)
library(dplyr)
library(openxlsx)

# 设置工作目录
setwd("F:\\postgraduate\\PNH_plan\\2.Paper\\12.scRNAseq\\中性粒细胞差异分析")

# 读取Seurat对象
seurat_obj <- qread("../27.Neutrophil_seurat.qs")

# 查看 orig.ident 列中的所有分组
table(seurat_obj$orig.ident)

# 确保 ident 列是 orig.ident
Idents(seurat_obj) <- seurat_obj$orig.ident

# 进行差异表达分析，比较 PNH 相对于 Control 的差异基因
diff_genes <- FindMarkers(seurat_obj, ident.1 = "PNH", ident.2 = "Control")

# 查看差异基因分析结果
head(diff_genes)

# 全部差异结果
significant_genes <- diff_genes %>%
  filter(p_val_adj < 1)

# 输出显著差异基因
head(significant_genes)

# 将行名（基因名）保存到数据框中
significant_genes_with_rownames <- cbind(gene_name = rownames(significant_genes), significant_genes)

# 保存结果为xlsx文件，并包括行名
write.xlsx(significant_genes_with_rownames, file = "1_all_diffgenes.xlsx")





#筛选条件
Log2FC <- 4
pvalue <- 0.05


# up-regulated genes
up_genes <- subset(significant_genes_with_rownames, p_val < pvalue & avg_log2FC > Log2FC)
write.xlsx(up_genes, file = "2_up.xlsx", rowNames = FALSE)

# down-regulated genes
down_genes <- subset(significant_genes_with_rownames, p_val < pvalue & avg_log2FC < -Log2FC)
write.xlsx(down_genes, file = "3_down.xlsx", rowNames = FALSE)









# rm(list = ls())
library(tidyverse)
library(readxl)
library(openxlsx)
data_miRNA <- significant_genes_with_rownames
head(data_miRNA)

#筛选条件
Log2FC <- 4
pvalue <- 0.05
# 替换列的名字
colnames(data_miRNA)[colnames(data_miRNA) == "avg_log2FC"] <- "Log2FC"
colnames(data_miRNA)[colnames(data_miRNA) == "p_val"] <- "pvalue"

head(data_miRNA)

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

pdf(file = "volcano.pdf", width = 7.5, height = 6.5)
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
text(x = c(Log2FC+3,-Log2FC-3), 
     y = rep(y_max*0.9, 2), 
     labels = c("Higher in PNH", "Higher in control"),
     pos = 3, # 左对齐和右对齐
     col = c("#FA8260", "#4D8FD1"), 
     cex = 1.1, 
     font = 3)

# 百分比
text(c(x_max-2,x_min+2), rep(y_max*0.9, 1), 
     c(paste(round(upregulated_percentage, 1), "%"), paste(round(downregulated_percentage, 1), "%")),
     pos = 3, col = c("#FA8260", "#4D8FD1"), cex = 1.1, font = 3)

title(main = "mRNA")

dev.off()
