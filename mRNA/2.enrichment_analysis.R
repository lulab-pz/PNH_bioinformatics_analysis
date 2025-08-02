
library(clusterProfiler)
library(stringr)
library(ComplexHeatmap)
library(org.Hs.eg.db)
options(stringsAsFactors = F)
library(ggplot2)
library(openxlsx)

kegg.buble = 1
kegg.chord = 1
# setwd("F:\\postgraduate\\PNH_plan\\3.自生成\\10.生成数据的富集分析")
enrich_ID <- read.xlsx("3_diff.xlsx")
# enrich_ID <- read.table("3.diff.xls",sep = "\t",header = T)
if(str_detect(colnames(enrich_ID)[3],'logFC')){
  colnames(enrich_ID)[3] <- 'log2FoldChange'
}
colnames(enrich_ID)[1] <- "gene_name"



#GO----
gene_symbol <- clusterProfiler::bitr(
  geneID = enrich_ID$gene_name,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db")

enrich_go <- enrichGO(gene = gene_symbol[,2],
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENTREZID",
                      ont = "ALL",
                      pAdjustMethod = "fdr",
                      pvalueCutoff = 1,
                      qvalueCutoff = 1,
                      readable = TRUE)

enrich_go <- mutate(enrich_go, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
GO_ALL <- enrich_go@result
if(median(GO_ALL$richFactor) < 0.1){GO_ALL$richFactor<- GO_ALL$richFactor*10}

GO_ALL <- GO_ALL[GO_ALL$pvalue < 0.05,]
BP <- GO_ALL[GO_ALL$ONTOLOGY=='BP', ]
BP <- BP[order(BP$Count,decreasing = T),]
CC <- GO_ALL[GO_ALL$ONTOLOGY=='CC', ]
CC <- CC[order(CC$Count,decreasing = T),]
MF <- GO_ALL[GO_ALL$ONTOLOGY=='MF', ]
MF <- MF[order(MF$Count,decreasing = T),]
paste0("得到",dim(GO_ALL)[[1]],"个结果，其中",dim(BP)[[1]],"个生物学过程，",dim(CC)[[1]],"个细胞组分，",dim(MF)[[1]],"个分子功能")
# [1]"得到297个结果，其中211个生物学过程，35个细胞组分，51个分子功能"
write.table(GO_ALL,file = "7_GO_ALL.csv",sep = ",",quote = F,row.names = F)
write.table(as.data.frame(BP), '8_GO_BP.csv', sep = ',', row.names = FALSE, quote = FALSE)
write.table(as.data.frame(CC), '9_GO_CC.csv', sep = ',', row.names = FALSE, quote = FALSE)
write.table(as.data.frame(MF), '10_GO_MF.csv', sep = ',', row.names = FALSE, quote = FALSE)

display_number = c(5, 5, 5)
ego_result_BP <- as.data.frame(BP)[1:display_number[1], ]
ego_result_CC <- as.data.frame(CC)[1:display_number[2], ]
ego_result_MF <- as.data.frame(MF)[1:display_number[3], ]
paste0(ego_result_BP$Description,collapse = "（）；")
paste0(ego_result_CC$Description,collapse = "（）；")
paste0(ego_result_MF$Description,collapse = "（）；")
col1 <- c("#E64B35B2", "#DC0000B2",'#E63863','#B2182B', '#D6604D','#EB4B17')
col2 <- c("#4DBBD5B2", "#0073C2B2",'#008280B2','#4393C3', '#23452F','#68A180')
col3 <- c("#E4C755", "#9FA3A8",'#D6E7A3','#F4A371FF','#F9C187FF','#FACF94FF')
Col1 <- sample(col1,1)
Col2 <- sample(col2,1)
Col3 <- sample(col3,1)
#GO柱状图-合并
if(T){
  go_enrich_df <- data.frame(
    ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
    Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
    GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
    type=factor(c(rep("BP:biological process", display_number[1]), 
                  rep("CC:cellular component", display_number[2]),
                  rep("MF:molecular function", display_number[3])), 
                levels=c("BP:biological process", "CC:cellular component","MF:molecular function" )))
  go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
  COLS <- c(Col1, Col2,Col3)
  
  p <- ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + #横纵轴取值
    geom_bar(stat="identity", width=0.8) + #柱状图的宽度，可以自己设置
    scale_fill_manual(values = COLS) + ###颜色
    coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
    xlab("GO term") + 
    ylab("Gene Number") + 
    labs(title = "GO Terms")+
    theme_bw()+theme(
      legend.background = element_rect(fill = "white", color = "black", size = 0.2),
      legend.text = element_text(face="bold",color="black",family = "Times",size=12),
      plot.title = element_text(hjust = 0.5, face = "bold",color = "black",family = "Times",size = 18),
      axis.text.x = element_text(face = "bold",color = "black",size = 14),
      axis.text.y = element_text(face = "bold",color = "black",size = 14),
      axis.title.x = element_text(face = "bold",color = "black",family = "Times",size = 18),
      axis.title.y = element_text(face = "bold",color = "black",family = "Times",size = 18),
      plot.subtitle = element_text(hjust = 0.5,family = "Times", size = 14, face = "italic", colour = "black")
    )+scale_x_discrete(labels=function(x) str_wrap(x, width=40))
  ggsave(filename = "11_GO_barplot.pdf", height = 8, width = 12, p)
} 



# KEGG----
Sys.setenv(all_proxy="socks5://127.0.0.1:7890")
gene_symbol <- clusterProfiler::bitr(
  geneID = enrich_ID$gene_name,  
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db")
KEGG <- enrichKEGG( gene = gene_symbol$ENTREZID,#基因列表 
                    organism = "hsa",  #物种
                    keyType = "kegg", #
                    minGSSize = 1, 
                    maxGSSize = 500,
                    pvalueCutoff = 1,  
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.9
)
KEGG <- mutate(KEGG, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
kk <- setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kegg_result <- kk@result
# if(median(kegg_result$richFactor) < 0.1){kegg_result$richFactor<- kegg_result$richFactor*10}
hh <- as.data.frame(kegg_result)
rownames(hh) <- 1:nrow(hh)
hh <- hh[hh$pvalue <= 0.05 & !grepl("Human Diseases", hh$category), ]
dim(hh)#15
hh <- hh[order(hh$Count,decreasing = T),]
write.xlsx(hh, file = "12_KEGG.xlsx")

rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
kk <- hh[c(1:10),] 
paste0(kk$Description,collapse = "（）；")


if(kegg.buble == 1){
  kk <- hh[c(1:20),] 
  p2 <-  ggplot(kk,aes(y=order,x=Count))+
    geom_point(aes(size=Count,color= pvalue))+# 修改点的大小
    scale_color_gradient(high="#EE0000B2",low = "#008B45B2")+
    labs(color=expression(pvalue,size="Count"), 
         x="Gene Number",y="Pathways",title="KEGG Enrichment")+
    theme_bw()+
    theme(
      legend.background = element_rect(fill = "white", color = "black", size = 0.2),
      legend.text = element_text(face="bold",color="black",family = "Times",size=14),
      plot.title = element_text(hjust = 0.5, face = "bold",color = "black",family = "Times",size = 18),
      axis.text.x = element_text(face = "bold",color = "black",size = 14),
      axis.text.y = element_text(face = "bold",color = "black",size = 14),
      axis.title.x = element_text(face = "bold",color = "black",family = "Times",size = 18),
      axis.title.y = element_text(face = "bold",color = "black",family = "Times",size = 18),
      plot.subtitle = element_text(hjust = 0.5,family = "Times", size = 14, face = "bold", colour = "black")
    )+scale_y_discrete(labels=function(x) str_wrap(x, width=40))
  ggsave(filename = "13_KEGG_dot.pdf", height = 8, width = 9, p2)
}