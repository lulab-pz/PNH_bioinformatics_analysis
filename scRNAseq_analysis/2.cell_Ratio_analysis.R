
library(Seurat)
library(qs)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(SeuratData)
allcolour = c("#DC143C","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#FFFF00","#0000FF",
              "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E",
              "#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22",
              "#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC",
              "#8B008B","#8B4513","#DEB887")
setwd("./细胞比例")
## 读取数据
pbmc3k.final <- qread("../26.聚类0.1手动注释.qs")

head(pbmc3k.final@meta.data)

pbmc3k<- UpdateSeuratObject(pbmc3k.final)
pbmc3k$group <- pbmc3k$orig.ident

pbmc3k$celltype <- pbmc3k@meta.data$celltype
DimPlot(pbmc3k,group.by = 'celltype',split.by = 'group')




metadata <- FetchData(pbmc3k, c("group","celltype"))

cell_counts <- dplyr::count(metadata, group, celltype)



total_cells_by_celltype <- dplyr::group_by(cell_counts, celltype) %>%
  dplyr::summarise(total = sum(n))

total_cells <- cell_counts %>%
  group_by(group) %>%
  summarise(total = sum(n))

cell_counts <- cell_counts %>%
  left_join(total_cells, by = "group")

cell_counts <- cell_counts %>%
  mutate(percentage = n / total*100)

write.csv(cell_counts,'cell_percentage.csv')





# 绘制细胞数量柱状图

p1 <- ggplot(total_cells_by_celltype, aes(x = total, y = celltype, fill = celltype)) +
  geom_bar(stat = "identity") +
  theme_classic() + 
  theme(text = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none") +
  labs(x = '# cells', y = NULL) +
  scale_fill_manual(values = allcolour)+  # 使用 scale_fill_manual 指定颜色+
  geom_text(aes(label = total),         # 添加数字标签
            hjust = 0.3,               # 调整标签位置（右对齐）
            size = 4, 
            color = "black")
p1




## 绘制分组比例柱状图
# 修改细胞比例计算方式（每个细胞类型内各组的比例）
cell_counts <- metadata %>%
  dplyr::count(group, celltype) %>%                           # 按组和细胞类型计数
  dplyr::group_by(celltype) %>%                               # 按细胞类型分组
  dplyr::mutate(
    total_celltype = sum(n),                                  # 计算细胞类型总细胞数
    percent_in_celltype = n / total_celltype * 100            # 计算组内占比
  ) %>%
  dplyr::ungroup()


p2 <- ggplot(cell_counts, aes(x = celltype, y = n, fill = group)) +
  geom_bar(
    stat = "identity",
    position = "fill",          # 自动将 y 轴标准化为 0-1
    width = 0.8,
    size = 0,
    colour = "black"
  ) +
  geom_text(
    aes(label = sprintf("%.1f", percent_in_celltype)),
    position = position_fill(vjust = 0.5),
    size = 4,
    color = "black"
  ) +
  scale_fill_manual(values = allcolour) +
  scale_y_continuous(  # 关键修改：将 y 轴标签转为百分比
    labels = scales::percent  # 使用 scales 包的百分比格式化函数
  ) +
  theme_classic() +
  labs(x = NULL, y = "# percentage") +
  coord_flip() +
  theme(
    axis.ticks.length = unit(0.1, "cm"),
    text = element_text(size = 15, colour = "black"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(colour = "black")
  )

p2




# 将两个图合并
library(patchwork)
p1 + p2 + plot_layout( nrow = 1)
ggsave('barplot_counts_percentage.pdf',width = 20, height = 8)
