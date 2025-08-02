# 加载必要的库
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(DOSE)
library(gridExtra)
library(dplyr)
library(tidyr)
library(cols4all)
library(stringr)
library(openxlsx)

# setwd("C:\\Users\\ASUS\\Desktop\\PNH_plan\\试验方案\\0.患者测序\\多组学\\药物预测\\己酮可可碱")

# setwd("F:\\postgraduate\\PNH_plan\\4.自生成\\5.富集分析\\transformed")
# 读取数据
# data <- read.table("C:\\Users\\ASUS\\Desktop\\PNH_plan\\2.final_result\\1.public_data\\1.patient\\1.差异分析\\4.DEseq2_diffALL.xls", header = TRUE, sep = "\t")
data <- read.xlsx("2_diffALL.xlsx")
colnames(data)[1] <- "gene_name"

# 添加 ENTREZID 列并按 log2FoldChange 排序
ids_symbol_and_entrezid <- bitr(geneID = data$gene_name,
                                fromType = "SYMBOL",
                                toType = "ENTREZID",
                                OrgDb = org.Hs.eg.db)

data$ENTREZID <- ids_symbol_and_entrezid[match(data$gene_name,
                                               ids_symbol_and_entrezid$SYMBOL),
                                         2]
data <- na.omit(data)

# 处理重复的 log2FoldChange 值
set.seed(123)
data$logFC <- jitter(data$logFC)
data <- data[order(data$logFC, decreasing = TRUE), ]

# 构建包含 ENTREZID 名称的 log2FoldChange 向量
gene_list <- data$logFC
names(gene_list) <- data$ENTREZID
Sys.setenv(all_proxy="socks5://127.0.0.1:7890")
# 使用 gseKEGG 进行 KEGG 通路富集分析
gsea_result <- gseKEGG(geneList = gene_list,
                       organism = "hsa",
                       minGSSize = 1,
                       maxGSSize = 10000,
                       pvalueCutoff = 1,
                       verbose = FALSE,
                       eps = 0)

# 将 GSEA 结果转换为数据框
gsea_result_df <- as.data.frame(gsea_result)

# 过滤 P 值小于 0.05 的结果
filtered_results <- gsea_result_df[gsea_result_df$pvalue < 0.05, ]
write.table(gsea_result_df, "14_GSEA.xls", sep = "\t", quote = F, row.names = F)
write.table(filtered_results, "15_GSEA_filtered.xls", sep = "\t", quote = F, row.names = F)
# 过滤 NES 列大于 0 的行，保存为 13.GSEA_UP.xlsx
gsea_up <- filtered_results[filtered_results$NES > 0, ]
write.xlsx(gsea_up, file = "16_GSEA_UP.xlsx")

# 过滤 NES 列小于 0 的行，保存为 14.GSEA_DOWN.xlsx
gsea_down <- filtered_results[filtered_results$NES < 0, ]
write.xlsx(gsea_down, file = "17_GSEA_DOWN.xlsx")


# 绘制山脊图并添加横轴标签
ridgeplot(gsea_result) + 
  labs(title = "GSEA KEGG Pathway Enrichment Ridgeplot") + 
  xlab("NES")

# 保存山脊图为 PDF 文件
pdf("15.GSEA_KEGG_ridgeplot.pdf", width = 25/2.54, height = 45/2.54)  # 单位用英寸，1英寸=2.54厘米
ridgeplot(gsea_result, showCategory = 20, core_enrichment = TRUE) + 
  labs(title = "GSEA KEGG Pathway Enrichment Ridgeplot") + 
  xlab("NES")
dev.off()



# 按 NES 值降序排列正的 NES 值并选择前五个
positive_nes <- filtered_results[filtered_results$NES > 0, ]
positive_nes <- positive_nes[order(positive_nes$NES, decreasing = TRUE), ]
top10_positive_nes <- head(positive_nes, nrow(positive_nes))

# 按 NES 值升序排列负的 NES 值并选择前五个
negative_nes <- filtered_results[filtered_results$NES < 0, ]
negative_nes <- negative_nes[order(negative_nes$NES, decreasing = FALSE), ]
top10_negative_nes <- head(negative_nes, nrow(negative_nes))

# 检查是否找到感兴趣的通路
if (nrow(top10_positive_nes) == 0 && nrow(top10_negative_nes) == 0) {
  stop("没有找到符合条件的通路")
}



for (i in 1:nrow(positive_nes)) {
  # 可视化每个通路
  plot <- gseaplot2(gsea_result,
                    geneSetID = top10_positive_nes$ID[i],
                    color = "red",
                    base_size = 11,
                    rel_heights = c(1.5, 0.5, 1),
                    subplots = 1:3,
                    pvalue_table = TRUE,
                    ES_geom = "line")
  # 生成文件名
  file_name <- paste0("18_GSEA_positive_NES_pathway_", i, ".pdf")
  # 保存为 PDF 格式文件
  ggsave(file_name, plot = plot, device = "pdf",
         width = 10, height = 8, units = "in", dpi = 300)
}


# 循环遍历前五个正的 NES 通路
for (i in 1:nrow(negative_nes)) {
  # 可视化每个通路
  plot <- gseaplot2(gsea_result,
                    geneSetID = top10_negative_nes$ID[i],
                    color = "blue",
                    base_size = 11,
                    rel_heights = c(1.5, 0.5, 1),
                    subplots = 1:3,
                    pvalue_table = TRUE,
                    ES_geom = "line")
  # 生成文件名
  file_name <- paste0("19_GSEA_negative_NES_pathway_", i, ".pdf")
  # 保存为 TIFF 格式文件
  ggsave(file_name, plot = plot, device = "pdf",
         width = 12, height = 8, units = "in", dpi = 300)
}


# GSEA结果做基因棒棒糖图
#棒棒糖图：反应显著的通路的正负向调节
# 读取 GSEA 数据文件
dat <- filtered_results

# 根据NES确定通路的方向
dat <- dat %>%
  mutate(Direction = ifelse(NES > 0, "Up", "Down"))

# 分别挑选上调和下调通路的前10个，按Count排序
top_up <- dat %>%
  filter(Direction == "Up") %>%
  arrange(desc(setSize)) %>%
  head(5)

top_down <- dat %>%
  filter(Direction == "Down") %>%
  arrange(desc(setSize)) %>%
  head(5)

# 合并上调和下调的前10个通路
top_paths <- bind_rows(top_up, top_down)

# 为了使得 up 向右延伸，down 向左延伸，我们将 NES 进行调整
top_paths <- top_paths %>%
  mutate(adjusted_NES = ifelse(Direction == "Up", NES, NES))

# 将长标签分成两行
top_paths <- top_paths %>%
  mutate(Description = str_wrap(Description, width = 30)) # 这里假设宽度为30，根据实际情况调整

# 创建图形
p1 <- ggplot(top_paths, aes(x = reorder(Description, adjusted_NES), y = adjusted_NES)) +
  coord_flip() +
  geom_col(aes(fill = Direction), width = 0.05) + # 修改 width 参数调整柱子的宽度
  geom_point(aes(size = setSize, color = pvalue)) + # 添加散点/气泡，颜色表示p值的大小
  scale_size_continuous(range = c(4, 10)) +
  scale_color_gradient2(midpoint = 0.05,
                        low = "#16888b", mid = "#c4d3d2", high = "green",
                        limits = c(0, 0.05), breaks = seq(0, 0.05, 0.01),
                        labels = seq(0, 0.05, 0.01)) + # 手动设置颜色渐变
  scale_fill_manual(values = c("Up" = "#e3c4a5", "Down" = "#97a2dc")) + # 手动设置填充颜色
  geom_hline(yintercept = 0, color = "white",
             size = 0.5, lty = "dashed") +
  geom_vline(xintercept = 0, color = "white",
             size = 0.5, lty = "dashed") + # 添加垂直线分割象限
  labs(x = NULL, y = "NES",
       title = "KEGG Pathway Enrichment",
       fill = "Direction", color = "p.adjust") +
  scale_y_continuous(expand = expansion(add = c(0.1, 0.1)),
                     limits = c(-3, 3),
                     breaks = seq(-3, 3, 1), #根据调整后NES值修改X轴刻度
                     labels = c(-3, -2, -1, 0, 1, 2, 3)) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10, color = "black"), # 确保正确设置 y 轴刻度标签样式
    axis.text.x = element_text(size = 8, color = "black", angle = 0, vjust = 1, hjust = 1) # 调整 x 轴刻度标签
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) + # 将 x 轴标签换行
  guides(size = guide_legend(title = "setSize"),
         fill = guide_legend(title = "Direction"))
# 保存图形为 TIFF 格式
ggsave(filename = "20_BBT_GSEA_result.pdf", plot = p1,
       device = "pdf", dpi = 300, width = 7, height = 8)
