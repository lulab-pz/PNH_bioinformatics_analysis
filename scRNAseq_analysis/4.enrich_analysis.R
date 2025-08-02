# 清空环境
rm(list = ls())
# 加载必要的包
library(Seurat)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(enrichplot)
library(org.Hs.eg.db)
library(tidyverse)

setwd("F:\\postgraduate\\PNH_plan\\2.Paper\\12.scRNAseq\\富集分析")

# 加载单细胞数据
scRNA <- qread("../26.聚类0.1手动注释.qs")
# 查看细胞类型分布
print("细胞类型分布：")
print(table(scRNA$celltype))


celltype_markers <- FindAllMarkers(scRNA, 
                                   group.by = "celltype",
                                   logfc.threshold = 0.5, 
                                   min.pct = 0.25,
                                   only.pos = TRUE)

# 保存标记基因结果，便于后续参考
write.csv(celltype_markers, "1_celltype_markers.csv", row.names = FALSE)

# 创建基因-细胞类型分组数据框
gene_groups <- data.frame(gene = celltype_markers$gene, 
                          celltype = celltype_markers$cluster)


# 将基因符号转换为Entrez ID，用于通路分析
gene_id_mapping <- bitr(celltype_markers$gene, 
                        fromType = "SYMBOL", 
                        toType = "ENTREZID", 
                        OrgDb = "org.Hs.eg.db")

# 输出基因映射统计信息
print(paste0("成功映射 ", nrow(gene_id_mapping), " 个基因，共 ", 
             length(unique(celltype_markers$gene)), " 个"))

# 合并基因ID与细胞类型信息
enrichment_data <- merge(gene_id_mapping, gene_groups, 
                         by.x = 'SYMBOL', 
                         by.y = 'gene')



# 按细胞类型进行KEGG通路富集分析
kegg_results <- compareCluster(ENTREZID ~ celltype,
                               data = enrichment_data,
                               fun = "enrichKEGG",
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05,  # 设置较为宽松的阈值以捕获更多通路
                               qvalueCutoff = 0.1,   
                               organism = "hsa")

# 将结果中的Entrez ID转回基因符号，提高可读性
kegg_results <- setReadable(kegg_results, 
                            OrgDb = "org.Hs.eg.db", 
                            keyType = "ENTREZID")

# 导出完整的富集分析结果
kegg_table <- kegg_results@compareClusterResult
write.csv(kegg_table, "2_celltype_kegg_enrichment.csv", row.names = FALSE)
table(kegg_table$Cluster)

library(dplyr)

# 假设kegg_table已经存在
kegg_table <- kegg_table %>%
  filter(pvalue < 0.05) %>%                             # 筛选p值小于0.05
  filter(!grepl("Human Diseases", category)) %>%        # 去掉Category含有Human Diseases的行
  group_by(Cluster) %>%                                 # 按Cluster分组
  arrange(desc(Count), .by_group = TRUE) %>%            # 每组按Count降序排列
  ungroup()

write.csv(kegg_table, "3_filtered_celltype_kegg_enrichment.csv", row.names = FALSE)





# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(scales)
library(stringr)
# install.packages("ggh4x")  # 如果没有安装ggh4x，请先运行这行代码

# 读取数据
data <- read.csv("3_filtered_celltype_kegg_enrichment.csv", stringsAsFactors = FALSE)

# 查看数据结构
str(data)
head(data)

# 查看所有可用的细胞类型
cat("可用的细胞类型：\n")
print(unique(data$celltype))

# 用户选择细胞类型（修改这里来选择您想要的细胞类型）
# 示例：选择前6个细胞类型，用户可以根据需要修改
selected_celltypes <- c("Neutrophil", "Macrophage", "Megakaryocyte", "Monocyte", "Classical monocyte")

# 如果您想要使用所有细胞类型，可以用下面这行代替上面的代码：
# selected_celltypes <- unique(data$celltype)

# 过滤数据，只保留选定的细胞类型
filtered_data <- data %>%
  filter(celltype %in% selected_celltypes)

# 数据预处理：为每个celltype选择zScore最高的前15条通路
top_pathways <- filtered_data %>%
  group_by(celltype) %>%
  arrange(desc(zScore)) %>%
  slice_head(n = 15) %>%
  ungroup()

# 为了更好的可视化，截断过长的Description
top_pathways$Description_short <- str_wrap(top_pathways$Description, width = 100)

# 创建颜色调色板
unique_celltypes <- unique(top_pathways$celltype)
n_celltypes <- length(unique_celltypes)

# 使用多样化的颜色
colors <- c("#E8A5A5", "#A5E8A5", "#A5A5E8", "#E8E8A5", "#E8A5E8", 
            "#A5E8E8", "#FFB366", "#66B3FF", "#66FFB3", "#B366FF",
            "#FF66B3", "#B3FF66", "#FF6666", "#66FF66", "#6666FF",
            "#FFA366", "#66FFA3", "#A366FF", "#FF66A3", "#A3FF66")

# 如果细胞类型数量超过预设颜色，使用RColorBrewer生成更多颜色
if(n_celltypes > length(colors)) {
  colors <- rainbow(n_celltypes)
}

# 创建颜色映射
color_mapping <- setNames(colors[1:n_celltypes], unique_celltypes)

# 为每个通路创建唯一的y轴位置
plot_data <- top_pathways %>%
  arrange(celltype, desc(zScore)) %>%
  mutate(
    # 为每个通路创建唯一的y轴标识
    pathway_id = paste0(celltype, "_", row_number()),
    # 创建用于排序的数字
    celltype_num = as.numeric(as.factor(celltype)),
    # 在同一细胞类型内按zScore排序
    within_celltype_rank = ave(zScore, celltype, FUN = function(x) rank(-x))
  ) %>%
  # 按细胞类型和zScore排序
  arrange(celltype_num, desc(zScore))

# 使用ggh4x包来实现不同的strip背景色
# 如果没有安装ggh4x，请先运行：install.packages("ggh4x")
library(ggh4x)

# 创建每个细胞类型对应的背景色列表
strip_colors <- color_mapping[unique_celltypes]
names(strip_colors) <- unique_celltypes

# 创建水平条形图 - 使用ggh4x实现不同背景色
p_modified <- ggplot(plot_data, aes(x = zScore + 1, y = reorder(Description_short, zScore))) +  # 将条形图向右偏移1个单位
  geom_col(aes(fill = celltype), alpha = 0.8, width = 0.8) +  # 增加条形宽度，使字体能被包围在条形里
  # 在颜色条块里添加通路名字，垂直居中，左侧偏移1个单位
  geom_text(aes(label = Description_short), 
            hjust = 0, x = 1.1, size = 6, vjust = 0.5, color = "black") +  # 增大字体到6，让文字更清晰
  scale_fill_manual(values = color_mapping, guide = "none") +
  facet_grid2(celltype ~ ., scales = "free_y", space = "free_y", switch = "y",
              strip = strip_themed(background_y = elem_list_rect(fill = strip_colors))) +  # 使用ggh4x设置不同背景色
  theme_minimal() +
  theme(
    # 隐藏y轴文本，因为文本现在在条形图内部
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 16),  # 增大x轴字体
    axis.title = element_text(size = 18),  # 增大轴标题字体，不加粗
    plot.title = element_text(size = 22, hjust = 0.5),  # 增大标题字体，不加粗
    # 调整细胞名显示：左对齐，显示在左侧
    strip.text.y.left = element_text(angle = 0, hjust = 0.5, size = 18, 
                                     margin = margin(l = 8, r = 15, t = 8, b = 8),
                                     color = "black"),  # 增大字体，增加边距
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.6, "lines"),  # 增加面板间距，使每个条形图的上下边距增大
    # 调整strip位置到左侧
    strip.placement = "outside",
    # 增加图的整体边距
    plot.margin = margin(30, 30, 30, 30)
  ) +
  labs(
    title = "Enrichment of celltype marker genes",
    x = "Z-Score",
    y = NULL
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  # 调整y轴的扩展，为每个条形图增加上下边距
  scale_y_discrete(expand = expansion(add = c(0.2, 0.2)))  # 为每个条形图增加上下0.2个单位的边距

# 创建一个更大的绘图设备并保存为 PDF
pdf("celltype_kegg_enrichment_modified.pdf", width = 15, height = 30)  # 增大图像宽度，让元素更清晰
print(p_modified)
dev.off()

# 也可以直接显示
# print(p_modified)


