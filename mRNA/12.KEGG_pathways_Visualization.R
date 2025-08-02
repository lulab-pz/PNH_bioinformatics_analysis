library(ggkegg)
library(readxl)
library(dplyr)
library(ggfx)
library(ggraph)
library(interp)
library(akima)
library(geomtextpath)
library(sysfonts)
library(pathview)
library(tidyr)

# 加载差异表达数据
diff_df <- read_excel("3_entrezID_diff.xlsx")
colnames(diff_df) <- tolower(colnames(diff_df))
colnames(diff_df)[colnames(diff_df) == "logfc"] <- "log2FoldChange"

# 加载通路注释数据
kegg_df <- read_excel("12_KEGG.xlsx")

head(diff_df)
print(colnames(kegg_df))



# 假设 diff_df$gene_id 是 ENTREZID
fc_vector <- diff_df$log2FoldChange
names(fc_vector) <- as.character(diff_df$gene_name)  # 必须是命名向量




selected_ids <- c("hsa04613", "hsa04610", "hsa04640")

# 创建路径 ID 到清洗后的文件名的映射（不使用 deframe）
id_to_name <- kegg_df %>%
  filter(ID %in% selected_ids) %>%
  mutate(name_clean = gsub(" ", "_", Description)) %>%
  { setNames(.$name_clean, .$ID) }

# 批量保存 PDF 图
for (pathway_id in selected_ids) {
  output_name <- id_to_name[[pathway_id]]
  output_path <- paste0("F:/postgraduate/PNH_plan/2.Paper/8.mRNA/", output_name, ".pdf")
  
  pathview(
    gene.data = fc_vector,
    pathway.id = pathway_id,
    species = "hsa",
    out.suffix = NULL,
    kegg.native = TRUE,
    print.image = FALSE,
    pdf.width = 9,
    pdf.height = 7,
    output.file = output_path
  )
}


# 自动转换 png 为 pdf 并改名
library(magick)

# 定义映射：KEGG ID 到描述
id_to_name <- c(
  "hsa04613" = "Neutrophil_extracellular_trap_formation",
  "hsa04610" = "Complement_and_coagulation_cascades",
  "hsa04640" = "Hematopoietic_cell_lineage"
)

# 批量处理
for (id in names(id_to_name)) {
  png_path <- paste0("hsa", substr(id, 4, 9), "..png")
  pdf_path <- paste0("F:/postgraduate/PNH_plan/2.Paper/8.mRNA/12.", id_to_name[[id]], ".pdf")
  
  img <- image_read(png_path)
  image_write(img, path = pdf_path, format = "pdf")
}

files_to_delete <- c("hsa04640..png", "hsa04640.png", "hsa04640.xml",
                     "hsa04610..png", "hsa04610.png", "hsa04610.xml",
                     "hsa04613..png", "hsa04613.png", "hsa04613.xml")
file.remove(files_to_delete)
