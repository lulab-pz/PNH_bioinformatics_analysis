#install.packages(c("cowplot", "igraph", "parallel", "doParallel"))
#install.packages("scCustomize", "patchwork", "WGCNA", "Seurat")

#BiocManager::install(c("UCell", "GeneOverlap"))

#library(devtools)
#devtools::install_github('smorabit/hdWGCNA', ref='dev')


#引用包
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(igraph)
library(parallel)
library(scCustomize)
library(doParallel)

#设置分析的参数
set.seed(12345)
theme_set(theme_cowplot())
options(future.globals.maxSize = 1024^10000)
enableWGCNAThreads(nThreads = 5)

cellName="Fibroblasts"       #细胞的名称(根据自己感兴趣的细胞进行修改)
setwd("C:\\Users\\Administrator\\Desktop\\hdWGCNA\\08.hdWGCNA")       #设置工作目录

#导入单细胞的数据
load("Seurat.Rdata")
seurat_obj=pbmc
seurat_obj$cell_type=Idents(seurat_obj)

#对基因进行过滤
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",    #挑选基因的方法
  fraction = 0.05,            #选择至少在5%的细胞中表达的基因
  wgcna_name = "hdWGCNA"     #WGCNA分析的名称
)

#构造metacells用于hdWGCNA分析
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("cell_type", "Type"),   # specify the columns in seurat_obj@meta.data to group by
  reduction = 'harmony',   # select the dimensionality reduction to perform KNN on
  k = 25,   # nearest-neighbors parameter
  max_shared = 10,   # maximum number of shared cells between two metacells
  ident.group = 'cell_type'   # set the Idents of the metacell seurat object
)

#对Metacells表达矩阵进行矫正
seurat_obj <- NormalizeMetacells(seurat_obj)

#创建表达矩阵用于hdWGCNA分析
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = cellName,   # the name of the group of interest in the group.by column
  group.by='cell_type',   # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA',   # using RNA assay
  layer = 'data'   # using normalized data
)


#选择最佳的power值
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed'   # you can also use "unsigned" or "signed hybrid"
)
#绘制power值的图形
plot_list <- PlotSoftPowers(seurat_obj)
pdf(file="WGCNA01.power.pdf", width=8, height=7)
wrap_plots(plot_list, ncol=2)
dev.off()
#获取power值的表格
power_table <- GetPowerTable(seurat_obj)
head(power_table)


#构建共表达网络
selectPower=10       #填写最佳的power值(需要根据power.pdf图形进行修改)
seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power = selectPower,
  tom_name = cellName    # name of the topoligical overlap matrix written to disk
)
#输出聚类的结果(得到每个基因属于哪个模块)
pdf(file="WGCNA02.cluster.pdf", width=8, height=6)
PlotDendrogram(seurat_obj, main=paste0(cellName, ' hdWGCNA Dendrogram'))
dev.off()


#获取TOM矩阵
TOM <- GetTOM(seurat_obj)
#计算模块基因的kME值
# need to run ScaleData first or else harmony throws an error:
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))
# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
	seurat_obj,
	group.by.vars="Type"
)
# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)
# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']
# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cell_type', group_name = cellName
)

#对模块进行重命名
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = paste0(cellName, "-M")
)

# 对每个模块基因的kME值进行可视化，得到模块核心基因的图形
pdf(file="WGCNA03.kME.pdf",width = 15,height = 15)
PlotKMEs(seurat_obj, ncol=5, n_hubs = 25)
dev.off()


#获取模块基因的表格
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')
# show the first 6 columns:
head(modules[,1:6])
write.csv(modules, file="moduleGene.csv", row.names=F)
#获取每个模块的核心基因
hub_df <- GetHubGenes(seurat_obj, n_hubs = 25)
head(hub_df)
write.csv(hub_df, file="hubGene.csv", row.names=F)
#saveRDS(seurat_obj, file='hdWGCNA_object.rds')


# 根据每个模块的top基因对细胞进行打分
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='Seurat'
)

# 绘制每个模块的散点图(根据hMEs进行FeaturePlot)
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  reduction = "tsne",
  features='hMEs',    # plot the hMEs
  order=TRUE     # order so the points with highest hMEs are on top
)
#输出图形
Ncol=5
Height=ifelse(length(mods)/5*3>6, length(mods)/5*3, 6)
pdf(file="WGCNA04.FeaturePlot.pdf",width = 15,height = Height)
wrap_plots(plot_list, ncol=Ncol)
dev.off()


#模块之间相关性的图形
pdf(file="WGCNA05.ModuleCorrelogram.pdf",width = 8,height = 7)
ModuleCorrelogram(seurat_obj)
dev.off()


#绘制模块和细胞之间相关性的图形
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
p6 <- DotPlot(seurat_obj, features = mods, group.by = "cell_type")
p6 <- p6 + RotatedAxis() +     #坐标轴转置
    scale_color_gradient2(high = "red", mid = "grey95", low = "blue")    #修改图形的颜色
pdf(file="WGCNA06.cellModuleCor.pdf", width = 9, height = 6)
print(p6)
dev.off()


#绘制每个模块的小提琴图
p7 <- VlnPlot(
    seurat_obj,
    features = mods,
    group.by = "cell_type",
    ncol=Ncol,
    sort="increasing",
    pt.size = 0
)
p7 <- p7 + geom_boxplot(width = .25, fill = "white")
p7 <- p7 + xlab("") + ylab("hME") + NoLegend()
pdf(file="WGCNA07.VlnPlot.pdf", width =16, height =Height)
print(p7)
dev.off()


#得到每个模块的网络图
ModuleNetworkPlot(seurat_obj = seurat_obj, mods = "all")


#绘制所有模块核心基因的网络图
pdf(file="WGCNA08.hubNetwork.pdf", width =7, height =7)
HubGeneNetworkPlot(
    seurat_obj,
    n_hubs = 3,     # The number of hub genes to plot for each module
    n_other = 5,     # The number of non-hub genes to sample from each module
    edge_prop = 0.75,    # The proportion of edges in the graph to sample
    mods = "all"
)
dev.off()


#观察基因在UMAP中的空间关系
seurat_obj <- RunModuleUMAP(
    seurat_obj, 
    n_hubs = 10,     #number of hub genes to use in the UMAP computation
    n_neighbors = 15,   #The size of local neighborhood used for manifold approximation
    min_dist = 0.1    #The effective minimum distance between embedded points
)
umap_df <- GetModuleUMAP(seurat_obj)
pdf(file="WGCNA09.geneUmap.pdf", width = 7, height = 6)
ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) + 
    geom_point(
        color = umap_df$color,     # 根据模块定义点的颜色
        size = umap_df$kME * 2     # 根据kME值定义点的大小
    ) + 
    umap_theme()
dev.off()


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


