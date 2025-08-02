# Seurat 包：用于处理单细胞数据
library(Seurat)
# CellChat 包：用于细胞通讯分析
library(CellChat)
# future 包：用于并行计算
library(future)
# ggplot2 包：用于绘制图形
library(ggplot2)
# 其他可能用到的包
library(dplyr)  # 用于数据处理
library(tidyr)  # 用于数据整理

setwd("F:\\postgraduate\\PNH_plan\\2.Paper\\12.scRNAseq\\细胞通讯")
# 细胞通讯----
seurat=qread("../26.聚类0.1手动注释.qs")

head(seurat@meta.data)


# 创建cellchat对象，注意需要用Normalized data, 这里用cellchat提供的`normalizeData`函数
cellchat <- createCellChat(normalizeData(seurat[['RNA']]$counts),
                           meta = seurat@meta.data,
                           group.by = 'celltype')

cellchat <- setIdent(cellchat, ident.use="celltype")
groupSize <- as.numeric(table(cellchat@idents))      #每个细胞类型的细胞数目


#导入配体和受体数据库
CellChatDB <- CellChatDB.human       #如果物种是小鼠改成CellChatDB.mouse

# cellchat <- subsetData(cellchat)
CellChatDB.use <- subsetDB(CellChatDB)
cellchat@DB <- CellChatDB.use

#查看配体和受体对的类型
pdf(file="40.COMM01.DatabaseCategory.pdf", width=7, height=5)
showDatabaseCategory(CellChatDB)
dev.off()

#对表达数据进行预处理
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 10) #设置使用的电脑核心数，并行计算节省时间
cellchat <- identifyOverExpressedGenes(cellchat)      #识别每种细胞中过表达的基因
# 增加 global max size 的限制
options(future.globals.maxSize = 9e9)  # 设置为 2GB

cellchat <- identifyOverExpressedInteractions(cellchat)      #识别基因的相互作用
cellchat <- smoothData(cellchat, adj = PPI.human)



#计算细胞通讯的概率
cellchat <- computeCommunProb(cellchat)
#过滤掉小于10个细胞的细胞通讯
cellchat <- filterCommunication(cellchat, min.cells = 10)
#输出细胞间的通讯关系
df.net=subsetCommunication(cellchat)
write.table(file="41.COMM02.Comm.network.xls", df.net, sep="\t", row.names=F, quote=F)


#在信号通路的水平进一步推测胞间的通讯, 推断通路水平的互作网络
cellchat <- computeCommunProbPathway(cellchat)

#对计算结果汇总整合，展示整体细胞通讯状态
cellchat <- aggregateNet(cellchat)
#输出细胞通讯的图形(互作数量的图形)
pdf(file="42.COMM03.cellNetworkCount.pdf", width=10, height=12)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()


#输出细胞通讯的图形(互作强度的图形)
pdf(file="43.COMM04.cellNetworkWeight.pdf", width=10, height=12)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction strength")
dev.off()


#分细胞类型展示(把单个细胞类型提取出来,观察这个细胞与其他细胞的通讯)
pdf(file="44.COMM05.singleCell.pdf", width=18, height=15)
weight_mat <- cellchat@net$weight
par(mfrow = c(2,3), mgp=c(0,0,0), xpd=TRUE)
for (cel in unique(cellchat@idents)){
  cir_mat <- matrix(0, nrow = nrow(weight_mat), ncol = ncol(weight_mat), dimnames = dimnames(weight_mat))
  cir_mat[cel, ] <- weight_mat[cel, ]
  netVisual_circle( cir_mat, vertex.weight= groupSize, weight.scale= T,edge.weight.max = max(weight_mat), vertex.label.cex=1,title.name=cel)
}
dev.off()


