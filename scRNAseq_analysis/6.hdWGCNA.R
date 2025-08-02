### 初始化 ----
rm(list = ls())
setwd("E:/研究生/课题/MPN/analysis/PB CD34+/GSE53482/miRNA")
if (!dir.exists("miRNA-mRNA")) dir.create("miRNA-mRNA")
setwd("miRNA-mRNA")

### 加载包 ----
library(dplyr)
library(igraph)
library(ggplot2)
library(tidygraph)
library(ggraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(miRBaseConverter)

### Step 1: 读取与处理数据库 ----

# miRTarBase
mirtarbase <- read.csv("data_base/hsa_MTI.csv", stringsAsFactors = FALSE)
mirtarbase_clean <- mirtarbase %>%
  dplyr::filter(`Species..miRNA.` == "hsa", `Species..Target.Gene.` == "hsa") %>%
  dplyr::filter(Support.Type == "Functional MTI") %>%
  dplyr::select(miRNA, target = Target.Gene) %>%
  dplyr::distinct()

# miRDB
mirdb <- read.table("data_base/miRDB_v6.0_prediction_result.txt", sep = "\t", header = FALSE)
colnames(mirdb) <- c("miRNA", "transcript_id", "score")
mirdb <- mirdb %>% filter(score >= 80)

mirdb$gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = mirdb$transcript_id,
  column = "SYMBOL",
  keytype = "REFSEQ",
  multiVals = "first"
)

mirdb_clean <- mirdb %>%
  dplyr::filter(!is.na(gene_symbol)) %>%
  dplyr::select(miRNA = miRNA, target = gene_symbol) %>%
  dplyr::distinct()

# TargetScan
tscan <- read.table("data_base/Summary_Counts.all_predictions.txt", sep = "\t", header = TRUE)
colnames(tscan) <- make.names(colnames(tscan))
tscan_clean <- tscan %>%
  dplyr::filter(Total.context...score <= -0.2) %>%
  dplyr::select(miRNA = Representative.miRNA, target = Gene.Symbol) %>%
  dplyr::distinct() %>%
  na.omit()

### Step 2: 处理差异表达数据 ----

# 读取差异表达 miRNA
miRNA_table <- read.table("../DEG-seq/diff.xls", header = TRUE, sep = "\t")
miRNA_list <- unique(miRNA_table$gene_symbol)

# 区分前体与成熟 miRNA
precursors <- miRNA_list[grepl("^hsa-mir-", miRNA_list)]
matures <- miRNA_list[grepl("^hsa-miR-", miRNA_list)]

precursor_df <- miRNA_NameToAccession(precursors, version = "v22") %>% na.omit()
miRNA_table_all <- getMiRNATable()

filtered_info_precursor <- miRNA_table_all[miRNA_table_all$Precursor_Acc %in% precursor_df$Accession, ]
mature_from_precursor <- unique(na.omit(c(filtered_info_precursor$Mature1, filtered_info_precursor$Mature2)))
final_mature_list <- unique(na.omit(c(mature_from_precursor, matures)))
write.table(final_mature_list, "final_mature_miRNA_list.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# 表达值合并
miRNA_expr <- data.frame(
  miRNA = miRNA_table$gene_symbol,
  log2FC = miRNA_table$logFC,
  p_value = miRNA_table$P.Value
)

precursor_to_mature <- data.frame(
  precursor = rep(filtered_info_precursor$Precursor, 2),
  mature = c(filtered_info_precursor$Mature1, filtered_info_precursor$Mature2)
) %>% na.omit()

mature_expr <- precursor_to_mature %>%
  left_join(miRNA_expr, by = c("precursor" = "miRNA")) %>%
  dplyr::select(mature, log2FC, p_value) %>%
  distinct()

direct_mature_expr <- miRNA_expr %>%
  dplyr::filter(grepl("^hsa-miR-", miRNA)) %>%
  dplyr::rename(mature = miRNA)

full_mature_expr <- bind_rows(mature_expr, direct_mature_expr) %>%
  distinct(mature, .keep_all = TRUE)

colnames(full_mature_expr) <- c("miRNA", "log2FC", "p_value")

# 读取 mRNA DEG
mRNA <- read.table("../../mRNA/DEG-seq/diff.xls", sep = "\t", header = TRUE)
mRNA_deg <- data.frame(
  gene_symbol = mRNA$gene_symbol,
  log2FC = mRNA$logFC,
  p_value = mRNA$P.Value
)

# 分类上下调
up_miRNA   <- full_mature_expr %>% filter(log2FC > 0) %>% pull(miRNA)
down_miRNA <- full_mature_expr %>% filter(log2FC < 0) %>% pull(miRNA)
up_mRNA    <- mRNA_deg %>% filter(log2FC > 0) %>% pull(gene_symbol)
down_mRNA  <- mRNA_deg %>% filter(log2FC < 0) %>% pull(gene_symbol)

### Step 3: 靶基因预测函数 ----

get_reverse_regulation <- function(miRNA_list, target_direction, deg_df,
                                   mirtarbase_df, mirdb_df, tscan_df,
                                   min_db_support = 2) {
  
  all_union <- bind_rows(
    mirtarbase_df %>% dplyr::rename(mature_mirna_id = miRNA, target_symbol = target) %>% mutate(source = "mirtarbase"),
    mirdb_df      %>% dplyr::rename(mature_mirna_id = miRNA, target_symbol = target) %>% mutate(source = "mirdb"),
    tscan_df      %>% dplyr::rename(mature_mirna_id = miRNA, target_symbol = target) %>% mutate(source = "targetscan")
  ) %>% dplyr::filter(mature_mirna_id %in% miRNA_list)
  
  supported_pairs <- all_union %>%
    group_by(mature_mirna_id, target_symbol) %>%
    summarise(db_count = n_distinct(source), .groups = "drop") %>%
    filter(db_count >= min_db_support)
  
  reverse_targets <- intersect(supported_pairs$target_symbol, target_direction)
  
  reverse_interactions <- supported_pairs %>%
    filter(target_symbol %in% reverse_targets) %>%
    left_join(full_mature_expr, by = c("mature_mirna_id" = "miRNA")) %>%
    left_join(deg_df, by = c("target_symbol" = "gene_symbol"))
  
  return(reverse_interactions)
}

### Step 4: 构建网络 ----

upmi_dnmr <- get_reverse_regulation(up_miRNA, down_mRNA, mRNA_deg,
                                    mirtarbase_clean, mirdb_clean, tscan_clean,
                                    min_db_support = 2)

dnmi_upmr <- get_reverse_regulation(down_miRNA, up_mRNA, mRNA_deg,
                                    mirtarbase_clean, mirdb_clean, tscan_clean,
                                    min_db_support = 1)

all_interactions <- bind_rows(
  mutate(upmi_dnmr, network = "Up-miRNA_Down-mRNA"),
  mutate(dnmi_upmr, network = "Down-miRNA_Up-mRNA")
)

nodes <- data.frame(
  id = unique(c(all_interactions$mature_mirna_id, all_interactions$target_symbol))
) %>%
  mutate(type = ifelse(grepl("^hsa-", id), "miRNA", "mRNA")) %>%
  left_join(bind_rows(
    dplyr::select(full_mature_expr, id = miRNA, log2FC_miRNA = log2FC),
    dplyr::select(mRNA_deg, id = gene_symbol, log2FC_mRNA = log2FC)
  ), by = "id")

edges <- all_interactions %>% dplyr::select(from = mature_mirna_id, to = target_symbol, network)
net <- graph_from_data_frame(edges, vertices = nodes, directed = TRUE)

### Step 5: 可视化网络 ----

V(net)$color <- case_when(
  V(net)$type == "miRNA" & V(net)$log2FC_miRNA > 0 ~ "#FF6666",
  V(net)$type == "miRNA" & V(net)$log2FC_miRNA < 0 ~ "#6666FF",
  V(net)$type == "mRNA"  & V(net)$log2FC_mRNA  > 0 ~ "#FF9999",
  V(net)$type == "mRNA"  & V(net)$log2FC_mRNA  < 0 ~ "#9999FF"
)
V(net)$size <- ifelse(V(net)$type == "miRNA", abs(V(net)$log2FC_miRNA)*5 + 10, abs(V(net)$log2FC_mRNA)*5 + 8)
V(net)$label.cex <- 0.8
V(net)$label.color <- "black"
E(net)$color <- ifelse(E(net)$network == "Up-miRNA_Down-mRNA", "#FF6666", "#6666FF")
E(net)$arrow.size <- 0.4

set.seed(123)
plot(net, layout = layout_with_fr, main = "miRNA-mRNA Reverse Regulatory Network", vertex.frame.color = "white")

### 可视化 tidygraph 格式 ----

tg_net <- as_tbl_graph(net)

ggraph(tg_net, layout = "fr") +
  geom_edge_link(aes(color = network),
                 arrow = arrow(length = unit(3, 'mm')),
                 end_cap = circle(2, 'mm'), alpha = 0.6) +
  geom_node_point(aes(color = type, size = abs(coalesce(log2FC_miRNA, 0) + coalesce(log2FC_mRNA, 0)))) +
  geom_node_text(aes(label = name), repel = TRUE, size = 2.5) +
  scale_color_manual(values = c(miRNA = "#FF6666", mRNA = "#9999FF")) +
  theme_void()

### Step 6: GO 富集分析 ----

all_targets <- unique(all_interactions$target_symbol)

ego <- enrichGO(
  gene = all_targets,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

if (nrow(ego) > 0) {
  dotplot(ego, showCategory = 15, title = "GO Enrichment: Reverse-Regulated Targets")
}

### Step 7: 保存结果 ----

write.csv(all_interactions, "reverse_regulation_interactions.csv", row.names = FALSE)
write.csv(nodes, "network_nodes.csv", row.names = FALSE)
write_graph(net, "reverse_regulation_network.graphml", format = "graphml")
if (nrow(ego) > 0) write.csv(ego@result, "GO_enrichment_results.csv", row.names = FALSE)

cat("✅ 全流程分析完成！")
