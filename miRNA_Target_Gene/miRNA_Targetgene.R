library(multiMiR)
library(openxlsx)

miRNA<-read.xlsx("F:\\postgraduate\\PNH_plan\\2.Paper\\6.miRNA\\9.diff.xlsx")
head(miRNA)

example <- get_multimir(org = "hsa",
                        mirna = miRNA$gene_id,
                        table = "all", 
                        summary = TRUE,
                        predicted.cutoff.type = "p",
                        predicted.cutoff = 10,
                        use.tibble = TRUE)


table(example@data$type)

example1_result <- example@data # 提取data数据，检索结果列表
head(example1_result)

###筛选有双荧光素酶报告实验的
example1_Lucifer <- example1_result[grep("Luciferase", example1_result$experiment), ]
unique(example1_result$experiment)


write.xlsx(example1_result,file = "2_multiMiR_result.xlsx")



# 基于预测结果的下游可视化
library(ggalluvial)
library(readxl)
library(ggplot2)
library(ggsankey)
library(dplyr)
library(viridis)
library(openxlsx)

# 桑基图不适合展示过多的交互节点，这里每个miRNA展示TOP3的靶基因
miRNA_list<-unique(example1_result$mature_mirna_id)
miRNA_list


##根据预测得分降序排列
example1_result <- example1_result[!is.na(example1_result$score), ]
example1_result$score <- as.numeric(example1_result$score)
example1_result<-example1_result[order(-example1_result$score), ]

##提取TOP3合并
sankey_input<-data.frame()
for (i in miRNA_list) {
  sankey_data<-example1_result[example1_result$mature_mirna_id == i,]
  sankey_data<-sankey_data[1:3,]
  sankey_input<-rbind(sankey_input,sankey_data)
}


sankey_input[1:5,]

sankey_input<-sankey_input[,3:4]
colnames(sankey_input)

# 格式化数据
sankey_long <- sankey_input %>%
  make_long(mature_mirna_id,target_symbol)

# 使用ggsankey绘图
p <- ggplot(sankey_long, aes(x = x,
                             next_x = next_x,
                             node = node,
                             next_node = next_node,
                             fill = factor(node),
                             label = node)) +
  geom_sankey(flow.alpha = 0.5, node.color = 1) +
  geom_sankey_label(size = 3, color = 1, fill = "white", hjust = 0) +
  scale_fill_viridis_d() +
  theme_sankey(base_size = 16) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16), # 设置x轴标签的字体大小
        axis.title.x = element_text(size = 16)) + # 设置x轴标题的字体大小
  labs(x = NULL, y = NULL) # 设置x轴标题
# 保存结果为tiff文件
ggsave("4_sankey_plot.tiff", plot = p, width = 18, height = 30, dpi = 300)

