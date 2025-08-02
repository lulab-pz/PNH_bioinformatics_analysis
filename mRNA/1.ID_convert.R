#1.SPOT ID to genesymbol----
setwd("F:\\PNH_plan\\5.PNH_Pro_Ultra_Max\\6.thp1_result\\mRNA")
library(openxlsx)
gse <- read.table("GSE183484_2v2_matrix.txt",header = T)
gpl <- read.table("GPL23159-184565.txt",header = T)
merge <- merge(data,gse,by="ID")
# 删除gene_names为空值的行
merge <- merge %>%
  filter(!is.na(gene_names) & gene_names != "")
library(org.Hs.eg.db)
library(clusterProfiler)

genes<-merge$SPOT_ID
genes=bitr(genes,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db")
colnames(genes)<-c('SPOT_ID','gene_symbol')

mergedata <- merge(genes,merge,by="SPOT_ID")
mergedata <- mergedata[,-1]
write.table(merge,"rawdata.txt",quote = F,sep = "\t",row.names = F)




#2.entrez to genesymbol----
setwd("F:\\postgraduate\\成果\\大会\\预计投稿\\中华医学会第十九次全国血液学学术会议\\BPDCN\\GSE184656RNA")

library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(biomaRt)
library(rentrez)
library(openxlsx)
library(org.Hs.eg.db)

# rawdata <- read.table("GSE80062.txt", header = TRUE, sep = "\t")
rawdata <- read.xlsx("1.rawdata.xlsx")
colnames(rawdata)[1] <- "gene_name"
first_column <- rawdata[, 1]

# 不同行取平均值
rawdata <- rawdata %>%
  group_by(gene_name) %>%
  summarise(across(everything(), mean))

# 查看可用的BioMart数据库
marts <- listMarts()
print(marts)

# 假设第一列是Entrez ID
entrez_ids <- first_column

# 使用正确的数据库和数据集
human <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# 进行Entrez ID到基因符号的转换
converted_ids <- getBM(attributes = c("entrezgene_id", "hgnc_symbol"),
                       filters = "entrezgene_id",
                       values = entrez_ids,
                       mart = human,
                       uniqueRows = TRUE)

# 合并数据
merged_data <- merge(converted_ids, rawdata,
                     by.x = "entrezgene_id", by.y = "gene_name")
merged_data <- merged_data[, -1]  # 删除Entrez ID列
colnames(merged_data)[1] <- "gene_name"

# 不同行取平均值
merged_data <- merged_data %>%
  group_by(gene_name) %>%
  summarise(across(everything(), mean))

# 写入结果
write.xlsx(merged_data, "2.rawdata_convertID.xlsx", rowNames = FALSE)






#3.ensembl to gene symbol----
setwd("C:\\Users\\ASUS\\Desktop\\SHSY5Y RNAseq\\药物预测\\8.mTOR通路基因")

library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(biomaRt)
library(rentrez)
library(openxlsx)
library(org.Hs.eg.db)

# rawdata <- read.table("GSE80062.txt", header = TRUE, sep = "\t")
rawdata <- read.xlsx("mTOR通路基因.xlsx")
colnames(rawdata)[1] <- "gene_name"
first_column <- rawdata[, 1]
#不同行取平均值
rawdata <- rawdata %>%
  group_by(gene_name) %>%
  summarise(across(everything(), mean))

# 假设第一列是Entrez ID
ensembl_ids <- first_column
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# 进行Entrez ID到基因符号的转换
converted_ids <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                       filters = "ensembl_gene_id",
                       values = ensembl_ids,
                       mart = human,
                       uniqueRows = TRUE)
merged_data <- merge(converted_ids, rawdata,
                     by.x = "ensembl_gene_id", by.y = "gene_name")
merged_data <- merged_data[, -1]
colnames(merged_data)[1] <- "gene_name"
#不同行取平均值
merged_data <- merged_data %>%
  group_by(gene_name) %>%
  summarise(across(everything(), mean))
write.table(merged_data, "1.rawdata.xls", quote = FALSE,
            sep = "\t", row.names = FALSE)





#4.mus to human----
setwd("C:\\Users\\zhaop\\Desktop")
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db) # 小鼠数据库
library(biomaRt)
library(rentrez)
library(openxlsx)
library(org.Hs.eg.db) # 人类数据库

rawdata <- read.table("1.rawdata.xls", header = TRUE, sep = "\t")
colnames(rawdata)[1] <- "gene_name"
# 不同行取平均值
rawdata <- rawdata %>%
  group_by(gene_name) %>%
  summarise(across(everything(), mean))

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                 host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl",
                 host = "https://dec2021.archive.ensembl.org/")
hsa2mus_all <- getLDS(attributes = c("mgi_symbol"),
                      filters = "mgi_symbol",
                      values = rawdata$gene_name,
                      mart = mouse,
                      attributesL = c("hgnc_symbol"),
                      martL = human,
                      uniqueRows = TRUE)
merged_data <- merge(hsa2mus_all, rawdata,
                     by.x = "MGI.symbol", by.y = "gene_name")
merged_data <- merged_data[, -1]
colnames(merged_data)[1] <- "gene_name"
# 不同行取平均值
merged_data <- merged_data %>%
  group_by(gene_name) %>%
  summarise(across(everything(), mean))
write.table(merged_data, "1.rawdata.xls", quote = FALSE,
            sep = "\t", row.names = FALSE)





#5.genesymbol to entrez----
setwd("F:\\postgraduate\\成果\\大会\\预计投稿\\中华医学会第十九次全国血液学学术会议\\BPDCN\\GSE184656RNA")

library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(biomaRt)
library(rentrez)
library(openxlsx)
library(org.Hs.eg.db)

# rawdata <- read.table("GSE80062.txt", header = TRUE, sep = "\t")
rawdata <- read.xlsx("F:\\postgraduate\\PNH_plan\\2.Paper\\8.mRNA\\3_diff.xlsx")
colnames(rawdata)[1] <- "gene_name"
first_column <- rawdata[, 1]

# 不同行取平均值
rawdata <- rawdata %>%
  group_by(gene_name) %>%
  summarise(across(everything(), mean))

# 查看可用的BioMart数据库
marts <- listMarts()
print(marts)

# 假设第一列是Entrez ID
entrez_ids <- first_column

# 使用正确的数据库和数据集
# human <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
human <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
# 进行Entrez ID到基因符号的转换
converted_ids <- getBM(attributes = c("entrezgene_id", "hgnc_symbol"),
                       filters = "hgnc_symbol",
                       values = entrez_ids,
                       mart = human,
                       uniqueRows = TRUE)

# 合并数据
merged_data <- merge(converted_ids, rawdata,
                     by.x = "hgnc_symbol", by.y = "gene_name")
merged_data <- merged_data[, -1]  # 删除Entrez ID列
colnames(merged_data)[1] <- "gene_name"

# 不同行取平均值
merged_data <- merged_data %>%
  group_by(gene_name) %>%
  summarise(across(everything(), mean))

# 写入结果
write.xlsx(merged_data, "F:\\postgraduate\\PNH_plan\\2.Paper\\8.mRNA\\3_entrezID_diff.xlsx", rowNames = FALSE)
