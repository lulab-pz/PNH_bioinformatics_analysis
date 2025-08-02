#加载R包----
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SingleR)
library(CCA)
library(clustree)
library(cowplot)
library(monocle)
library(tidyverse)
library(SCpubr)
library(UCell)
library(irGSEA)
library(GSVA)
library(GSEABase)
library(harmony)
library(plyr)
library(celldex)
library(assertthat)
library(Matrix)
library(stringr)
library(tricycle)   
library(scattermore)
library(scater)
library(randomcoloR)
library(CellChat)
library(future)
library(ggforce)
library(ggsci)
library(BiocParallel)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(velocyto.R)
library(SeuratWrappers)
library(cellcall)
library(limma)
library(NMF)
library(ggalluvial)
library(svglite)
library(ggpubr)
library(RColorBrewer)
library(reshape2)
library(openxlsx)
library(qs)





# 疾病样本单细胞对象创建----
data_dir_pnh <- "./treat/treat_rawdata"
# 获取子文件夹的名字
samples_pnh <- list.files(data_dir_pnh, full.names = TRUE)
# 创建一个空列表来存储每个子文件夹对应的 Seurat 对象
seurat_objects_pnh <- list()
# 循环遍历每个子文件夹，读取数据并创建 Seurat 对象
for (folder in samples_pnh) {
  # 检查文件夹路径是否存在
  if (file.exists(folder)) {
    # 读取单细胞数据
    data_pnh <- Read10X(data.dir = folder)
    
    # 创建 Seurat 对象
    seurat_pnh <- CreateSeuratObject(
      counts = data_pnh,
      project = "PNH",
      min.cells = 3,
      min.features = 200
    )
    # 获取子文件夹名称作为type列
    seurat_pnh$type <- basename(folder)
    
    # 将 Seurat 对象存入列表，使用子文件夹名称作为键
    seurat_objects_pnh[[basename(folder)]] <- seurat_pnh
  } else {
    message(paste("文件夹不存在:", folder))
  }
}

# 将所有 Seurat 对象合并成一个 Control 的 Seurat 对象
seurat_pnh <- Reduce(function(x, y) merge(x, y), seurat_objects_pnh)

# 查看合并后的 Seurat 对象
head(seurat_pnh@meta.data)
dim(seurat_pnh)


#进行数据质控
seurat_pnh[["percent.mt"]] <- PercentageFeatureSet(seurat_pnh, pattern = "^MT-")#线粒体基因
seurat_pnh[["percent.rb"]] <- PercentageFeatureSet(seurat_pnh, pattern = "^RP")#核糖体基因
randomColor <- function() {
  paste0("#",paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE),collapse = ""))
}
randomColors <- replicate(100,randomColor())

pdf("1.pnh_VlnPlot_QC.pdf")
VlnPlot(seurat_pnh, features = c("nCount_RNA","nFeature_RNA", "percent.mt"), ncol = 4)+scale_fill_manual(values =randomColors)
dev.off()

#过滤异常值
mask1 <- seurat_pnh$nCount_RNA<= 30000#细胞内检测到的分子总数
mask2 <- seurat_pnh$nFeature_RNA <= 5000 #每个细胞中检测到的基因数量
mask3 <- seurat_pnh$percent.mt <= 8#结合线粒体基因（percent.mt）去除异常值
# mask4 <- seurat_AML$percent.rb<= 1#核糖体基因（percent.rb）除去异常值
#可除去大多数双峰/死细胞/空液滴
seurat_pnh <- seurat_pnh[, mask1 & mask2 & mask3]

# head(seurat_PNH@meta.data)
dim(seurat_pnh)
# 保存PNH----
qsave(seurat_pnh,"1.seurat_all_pnh.qs")
head(seurat_pnh@meta.data)





# 正常样本单细胞对象创建----
data_dir_con <- "F:\\postgraduate\\PNH_plan\\2.Paper\\12.scRNAseq\\normal\\GSE157344\\GSE157344_RAW"
# 获取子文件夹的名字
samples_con <- list.files(data_dir_con, full.names = TRUE)
# 创建一个空列表来存储每个子文件夹对应的 Seurat 对象
seurat_objects_con <- list()
# 循环遍历每个子文件夹，读取数据并创建 Seurat 对象
for (folder in samples_con) {
  # 检查文件夹路径是否存在
  if (file.exists(folder)) {
    # 读取单细胞数据
    data_con <- Read10X(data.dir = folder)
    
    # 创建 Seurat 对象
    seurat_con <- CreateSeuratObject(
      counts = data_con,
      project = "Control",
      min.cells = 3,
      min.features = 200
    )
    # 获取子文件夹名称作为type列
    seurat_con$type <- basename(folder)
    
    # 将 Seurat 对象存入列表，使用子文件夹名称作为键
    seurat_objects_con[[basename(folder)]] <- seurat_con
  } else {
    message(paste("文件夹不存在:", folder))
  }
}

# 将所有 Seurat 对象合并成一个 Control 的 Seurat 对象
seurat_con <- Reduce(function(x, y) merge(x, y), seurat_objects_con)

# 查看合并后的 Seurat 对象
head(seurat_con@meta.data)
dim(seurat_con)


#进行数据质控
seurat_con[["percent.mt"]] <- PercentageFeatureSet(seurat_con, pattern = "^MT-")#线粒体基因
seurat_con[["percent.rb"]] <- PercentageFeatureSet(seurat_con, pattern = "^RP")#核糖体基因
randomColor <- function() {
  paste0("#",paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE),collapse = ""))
}
randomColors <- replicate(100,randomColor())

pdf("1.con_VlnPlot_QC.pdf")
VlnPlot(seurat_con, features = c("nCount_RNA","nFeature_RNA", "percent.mt"), ncol = 4)+scale_fill_manual(values =randomColors)
dev.off()

#过滤异常值
mask1 <- seurat_con$nCount_RNA<= 10000#细胞内检测到的分子总数
mask2 <- seurat_con$nFeature_RNA <= 2000 #每个细胞中检测到的基因数量
mask3 <- seurat_con$percent.mt <= 15#结合线粒体基因（percent.mt）去除异常值
# mask4 <- seurat_AML$percent.rb<= 1#核糖体基因（percent.rb）除去异常值
#可除去大多数双峰/死细胞/空液滴
seurat_con <- seurat_con[, mask1 & mask2 & mask3]

# head(seurat_PNH@meta.data)
dim(seurat_con)
# 保存PNH----
qsave(seurat_con,"1.seurat_all_con.qs")
head(seurat_con@meta.data)




#合并疾病和正常的seurat样本----
seurat_obj <- merge(seurat_con, seurat_pnh, add.cell.ids = c("Control", "PNH"))

# 查看合并后的 Seurat 对象的信息
head(seurat_obj@meta.data)
dim(seurat_obj)


pdf("2.Con_PNH_VlnPlot_QC.pdf")
VlnPlot(seurat_obj, features = c("nCount_RNA","nFeature_RNA", "percent.mt"), ncol = 4)+scale_fill_manual(values =randomColors)
dev.off()




# V5版本的需要joinLayers
LayerData(seurat_obj, assay = "RNA", layer = "counts")
seurat_obj <- JoinLayers(seurat_obj)
dim(seurat_obj[["RNA"]]$counts )



# 保存合并之后的变量----
qsave(seurat_obj, "2.Con_PNH_Seurat_obj.qs")




#标准化LogNormalize----
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

#鉴定高变基因，定义前2000个为高变基因----
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# 提取前10的高变基因----
top10 <- head(VariableFeatures(seurat_obj), 10)

# 展示高变基因----
plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

pdf(file = "3.展示高变基因.pdf",width =7,height = 8)
plot2
dev.off()





#对所有基因归一化----
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)



# PCA降维----
seurat_obj <- Seurat::RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
seurat_obj <- Seurat::RunUMAP(seurat_obj,dims = 1:20)


pdf(file = "4.降维UMAP图.pdf",width =7.5,height = 5.5)
DimPlot(seurat_obj, reduction = "umap",pt.size = 0.5)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right") #top为图列位置最上方，除此之外还有right,left,bottom(意思同英文)
dev.off()
pdf(file = "5.降维pca图.pdf",width =7.5,height = 5.5)
DimPlot(seurat_obj, reduction = "pca",pt.size = 0.5)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()


colnames(seurat_obj@meta.data)
# 根据 orig.ident 设置 Type 列
seurat_obj$Type <- ifelse(seurat_obj$orig.ident == "Control", "Control",
                          ifelse(seurat_obj$orig.ident == "PNH", "PNH", NA))

# 查看修改后的数据
head(seurat_obj@meta.data)

colaa=distinctColorPalette(100)
# 生成一组具有明显区别的颜色，避免颜色之间过于相似，适用于需要区分多个类别或项的图表
pdf(file = "6.降维UMAP单样品分布图.pdf",width =12,height = 7.5)
do_DimPlot(sample = seurat_obj,
           plot.title = "",
           reduction = "umap",
           legend.position = "bottom",
           dims = c(1,2),split.by = "orig.ident",pt.size =0.5
)
dev.off()
head(seurat_obj@meta.data)






# 去除批次效应----
library(future)
plan("multicore", workers = 10)  # 或 plan("multiprocess", workers = 4)

options(future.globals.maxSize = 60 * 1024^3)  # 允许最大内存占用（这里是100GB）

seurat2 <- RunHarmony(seurat_obj, group.by.vars = "orig.ident")
# ## 鉴定高变基因（由于去除了存在批次的细胞，高变基因可能会发生改变，因此需要重新鉴定高变基因）
seurat2 <- FindVariableFeatures(seurat2, selection.method = "vst", nfeatures = 2000)
# # 提取前10的高变基因
top10 <- head(VariableFeatures(seurat2), 10)
# # 展示高变基因
plot1 <- VariableFeaturePlot(seurat2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf(file = "7.去批次后鉴定高变基因.pdf",width =7,height = 6)
plot2
dev.off()
#3.查看去批次后的PCA,T-SNE及其单样本分布图
pdf(file = "8.harmony去批次pca图.pdf",width =7.5,height = 5.5)
DimPlot(seurat2, reduction = "harmony",pt.size = 0.5)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
seurat2 <- Seurat::RunUMAP(seurat2,dims = 1:20,reduction ='harmony')
pdf(file = "9.去批次后umap图.pdf",width =7.5,height = 5.5)
DimPlot(seurat2, reduction = "umap",pt.size = 0.5)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
collist=c(ggsci::pal_nejm()(8))
names(collist)=names(table(seurat2$orig.ident))
pdf(file = "10.去批次后umap单样本分布图orig.ident.pdf",width =12,height = 7.5)
do_DimPlot(sample = seurat2,
           plot.title = "",
           reduction = "umap",
           legend.position = "bottom",
           dims = c(1,2),split.by = "orig.ident",pt.size =0.5
)
dev.off()
##保存去批次后数据
qsave(seurat2,"11.去批次后seurat.qs")









seurat_obj <- qread("11.去批次后seurat.qs")
head(seurat_obj@meta.data)
dim(seurat_obj)

collist=c(ggsci::pal_nejm()(8))
# pal_nejm：这是ggsci包中的一个函数，返回基于《新英格兰医学杂志》（NEJM）使用的颜色调色板，设计上更适合用于科学图表
# 返回NEJM调色板中的8种颜色
names(collist)=names(table(seurat_obj$orig.ident))

pdf(file = "12.ElbowPlot.pdf",width =5,height = 4)
ElbowPlot(seurat_obj,ndims = 50)
dev.off()
#选择一个节点，其后续的节点都趋于平滑，不出现断层

#热图可视化前20个PC
pdf(file = "13.选20个PC热图.pdf",width =7.5,height = 11)
DimHeatmap(seurat_obj, dims = 1:20, cells = 1000, balanced = TRUE)
dev.off()
#主要看哪些PC中的基因能很好的将其分为两个明显的模块

# seurat_obj <- JackStraw(seurat_obj, num.replicate = 100)
# seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:20)
# pdf(file = "14.Jackstrawplot.pdf",width =7.5,height = 5.5)
# JackStrawPlot(seurat_obj, dims = 1:20)
# dev.off()
#展示每个PC的p值

#选择PC数
seuratPC=30
##对细胞聚类
seurat_obj=FindNeighbors(seurat_obj, dims = 1:seuratPC, reduction = "pca")


for (res in c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) { 
  seurat_obj=FindClusters(seurat_obj, graph.name = "RNA_snn", resolution = res, algorithm = 1)}
apply(seurat_obj@meta.data[,grep("RNA_snn_res",colnames(seurat_obj@meta.data))],2,table)

p2_tree=clustree(seurat_obj@meta.data, prefix = "RNA_snn_res.")
pdf(file = "15.挑选分辨率.pdf",width =25,height =20)
p2_tree
dev.off()
#增加分辨率后每个细胞群还能往下细分为多少个子细胞群




#聚类0.25----
#选择分辨率进行降维
px=0.1
seurat_obj1 <- FindClusters(seurat_obj, resolution = px)



# 查看结果
head(seurat_obj1@meta.data)


DefaultAssay(seurat_obj1)
seurat_obj1 <- JoinLayers(seurat_obj1, assay = "RNA")

# only.pos：只保留上调差异表达的基因
seurat_obj1.markers <- FindAllMarkers(seurat_obj1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(seurat_obj1.markers,file = "16.每个聚类0.1的marker基因.csv")
# 新增 'difference' 列，计算 pct.1 - pct.2
seurat_obj1.markers$difference <- seurat_obj1.markers$pct.1 - seurat_obj1.markers$pct.2

# 按照每个簇内的 'difference' 降序排列
seurat_obj1.markers_sorted <- seurat_obj1.markers %>%
  group_by(cluster) %>%
  arrange(desc(difference)) %>%
  ungroup()

# 保存为xlsx文件
write.xlsx(seurat_obj1.markers_sorted, file = "17.每个聚类0.1的marker基因.xlsx")

# 提取每个细胞簇的前十个标记基因
top_genes_per_cluster <- lapply(unique(seurat_obj1.markers_sorted$cluster), function(cluster) {
  top_genes <- seurat_obj1.markers_sorted[seurat_obj1.markers_sorted$cluster == cluster,]
  top_genes <- head(top_genes$gene, 10)  # 获取前10个基因
  paste(top_genes, collapse = ",")  # 将基因名称拼接成字符串
})

# 为 top_genes_per_cluster 添加簇的名称
names(top_genes_per_cluster) <- unique(seurat_obj1.markers_sorted$cluster)

# 检查 top_genes_per_cluster 和 seurat_obj.markers_sorted 的内容
print(top_genes_per_cluster)
print(head(seurat_obj1.markers_sorted))

# 将每个簇的前十个标记基因连接成字符串并保存在txt文件中
txt_content <- sapply(names(top_genes_per_cluster), function(cluster) {
  paste0("cluster_", cluster, ": ", top_genes_per_cluster[[cluster]])
})

# 强制将 txt_content 转换为字符类型
txt_content <- as.character(txt_content)

# 直接尝试保存到指定路径
writeLines(txt_content, con = "18.每个聚类0.1前十个marker基因.txt")


#选择每个聚类前5各基因绘制热图
top5seurat_obj1.markers <- seurat_obj1.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = difference)
# 对于每个群集，挑选出5个差异表达最显著的基因  

col <- c(ggsci::pal_npg()(9),ggsci::pal_jco()(9),ggsci::pal_jama()(7),ggsci::pal_nejm()(8))

pdf(file = "19.聚类0.1后热图.pdf",width =22,height = 16)
DoHeatmap(seurat_obj1,features = top5seurat_obj1.markers$gene,
          group.colors = col) +
  ggsci::scale_colour_npg() +
  scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',
                       name = 'Z-score')
dev.off()


## 将细胞在低维空间可视化UMAP/tSNE
seurat_obj1 <- RunUMAP(seurat_obj1, dims = 1:seuratPC, reduction = "pca")
# 可视化UMAP/tSNE3
pdf(file = "20.聚类0.1后UMAP.pdf",width =6.5,height = 5.5)
DimPlot(seurat_obj1, reduction = "umap", label = T, label.size = 3.5,pt.size = 0.1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()


head(seurat_obj1@meta.data)
#大数据用UMAP，小数据用TSNE
# UMAP在保留数据的全局结构的同时，更注重保留数据的局部结构，因此相对于TSNE更容易发现群簇的细节信息
# TSNE则在较低维度时能够展示出更好的聚类效果，但随着维度的增加，其聚类效果可能会受到影响

# 只保留需要的列
seurat_obj1@meta.data <- seurat_obj1@meta.data[, c("orig.ident",
                                                   "nCount_RNA",
                                                   "nFeature_RNA",
                                                   "percent.mt",
                                                   "Type",
                                                   "RNA_snn_res.0.1")]#每次运行需要修改一下

head(seurat_obj1@meta.data)

# 保存聚类后数据
qsave(seurat_obj1,"21.聚类0.1后seuratUMAP.qs")




#细胞注释0.5----
seurat=qread("21.聚类0.1后seuratUMAP.qs")#读取数据

marker_cluster <- c(
  "ALAS2","HBD","AHSP",#0.Erythroblast
  "IL7R","IL32","LTB",#1.T cell
  "FCGR3B","G0S2","NAMPT",#2.Neutrophil
  "CXCL8","G0S2","NAMPT",#3.Neutrophil
  "TRAC","IL7R",#4.T cell
  "NKG7","PRF1","CST7",#5.T cell
  "CST3","FCN1","VCAN",#6.Classical monocyte
  "GZMH","NKG7","CCL4",#7.NK cell
  "SEPT7","ATP5L","ATP5E",#8.Megakaryocyte-derived cell
  # 8. 高线粒体活性、结构蛋白丰富、缺乏传统免疫细胞特异性标记
  "FCGR3B","ALPL","PROK2",#9.Neutrophil
  "RETN","MPEG1","MS4A6A",#10.Monocyte
  "CA1","AHSP","SNCA",#11.Erythrocyte
  "MS4A6A","CD36","CYBB",#12.Classical monocyte
  "CD79A","BANK1","MS4A1",#13.B cell
  "PPBP","PF4","NRGN",#14.Megakaryocyte
  "CD79A","TNFRSF13C",#15.B cell
  "HMBS","ANK1","SLC4A1",#16.Erythroid progenitors/precursors cell
  # 16. ANK1、SLC4A1、HMBS，几乎都是红细胞（erythroid cells）/红细胞祖细胞专属的经典marker基因
  "LCN2","CAMP","LTF",#17.Neutrophil progenitor cell
  "MAFB","C10orf54","GNB2L1",#18.Macrophage
  # 18. MAFB是最具指示性的单核-巨噬细胞谱系marker
  # C10orf54强化了这一判断（常见于免疫抑制型单核/巨噬细胞）
  "HDC","IL3RA"#19.Megakaryocyte-erythroid progenitor cell
)

# 可视化 marker 基因的表达情况
dotplot <- DotPlot(seurat, features = unique(marker_cluster), assay = 'RNA') + coord_flip()
ggsave("22.聚类0.1marker基因点图.pdf", plot = dotplot, width = 15, height = 12)
head(seurat@meta.data)



heatmap <- DoHeatmap(seurat, features = marker_cluster) + NoLegend()
ggsave("23.聚类0.1marker基因热图.pdf", plot = heatmap, width = 15, height = 12)


# 根据 marker 基因手动注释细胞类型
new.cluster.ids <- c(
  "Erythroblast",	#0
  "T cell",	#1
  "Neutrophil",	#2
  "Neutrophil",	#3
  "T cell",	#4
  "T cell",	#5
  "Classical monocyte",	#6
  "NK cell",	#7
  "Megakaryocyte-derived cell",	#8
  "Neutrophil",	#9
  "Monocyte",	#10
  "Erythrocyte",	#11
  "Classical monocyte",	#12
  "B cell",	#13
  "Megakaryocyte",	#14
  "B cell",	#15
  "Erythroid progenitors_precursors cell",	#16
  "Neutrophil progenitor cell",	#17
  "Macrophage",	#18
  "Megakaryocyte-erythroid progenitor cell"	#19
)

# 更新 Seurat 对象的聚类标识符
names(new.cluster.ids) <- levels(seurat)
seurat <- RenameIdents(seurat, new.cluster.ids)

# 绘制注释后的 UMAP 图
annotated_umap <- DimPlot(seurat, reduction = "umap", group.by = "ident", pt.size = 0.5, label = TRUE)
ggsave("24.聚类0.1细胞类型注释UMAP图.pdf", plot = annotated_umap, width = 20, height =16)


# 根据 RNA_snn_res.0.1 的值为每个细胞类型分配标签
seurat$celltype <- dplyr::case_when(
  seurat$RNA_snn_res.0.1 == 0 ~ "Erythroblast",
  seurat$RNA_snn_res.0.1 == 1 ~ "T cell",
  seurat$RNA_snn_res.0.1 == 2 ~ "Neutrophil",
  seurat$RNA_snn_res.0.1 == 3 ~ "Neutrophil",
  seurat$RNA_snn_res.0.1 == 4 ~ "T cell",
  seurat$RNA_snn_res.0.1 == 5 ~ "T cell",
  seurat$RNA_snn_res.0.1 == 6 ~ "Classical monocyte",
  seurat$RNA_snn_res.0.1 == 7 ~ "NK cell",
  seurat$RNA_snn_res.0.1 == 8 ~ "Megakaryocyte-derived cell",
  seurat$RNA_snn_res.0.1 == 9 ~ "Neutrophil",
  seurat$RNA_snn_res.0.1 == 10 ~ "Monocyte",
  seurat$RNA_snn_res.0.1 == 11 ~ "Erythrocyte",
  seurat$RNA_snn_res.0.1 == 12 ~ "Classical monocyte",
  seurat$RNA_snn_res.0.1 == 13 ~ "B cell",
  seurat$RNA_snn_res.0.1 == 14 ~ "Megakaryocyte",
  seurat$RNA_snn_res.0.1 == 15 ~ "B cell",
  seurat$RNA_snn_res.0.1 == 16 ~ "Erythroid progenitors_precursors cell",
  seurat$RNA_snn_res.0.1 == 17 ~ "Neutrophil progenitor cell",
  seurat$RNA_snn_res.0.1 == 18 ~ "Macrophage",
  seurat$RNA_snn_res.0.1 == 19 ~ "Megakaryocyte-erythroid progenitor cell",
  TRUE ~ "Unknown"  # 对应没有明确分类的情况
)

# 查看结果
head(seurat@meta.data)

# 保存seurat变量
qsave(seurat,"26.聚类0.1手动注释.qs")




#单独拿出来细胞簇----
# 获取meta数据的前几行
head(seurat@meta.data)

# 获取细胞类型的所有levels
cell_types <- levels(seurat)

# 遍历每个细胞类型，进行子集选择和保存
for (cell_type in cell_types) {
  # 创建子集对象
  cell_type_seurat <- subset(seurat, subset = celltype == cell_type)
  
  # 创建保存文件名
  file_name <- paste0("27.", gsub(" ", "_", cell_type), "_seurat.qs")
  
  # 保存为新的 .rds 文件
  qsave(cell_type_seurat, file_name)
}