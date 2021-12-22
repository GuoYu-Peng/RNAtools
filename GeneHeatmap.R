suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(viridis))


what_plot <- "从基因表达矩阵生成热图，默认第一列为基因 ID"
parser <- ArgumentParser(description = what_plot, add_help = TRUE)
parser$add_argument("--expression", dest = "EXPR", help = "csv 格式的表达数据", required = TRUE)
parser$add_argument("--output", dest = "OUTPUT", help = "图片输出路径，要求为 pdf 格式。默认：GeneHeatmap.pdf", 
                    default = "GeneHeatmap.pdf")
parser$add_argument("--normalize", dest = "NORMALIZE", help = "数据 Normalize 方法，注意选择 log2 默认进行 log2(M + 1) 处理。默认：no", 
                    choices = c("no", "log2", "z-score"), default = "no")
parser$add_argument("--list", dest = "LIST", help = "可选参数，显示指定的基因，文件里每个基因一行", default = NULL)
parser$add_argument("--legend_name", dest = "LEGEND", help = "图例名称。默认：Normalized Expression", default = "Normalized Expression")
parser$add_argument("--column_title", dest = "COLUMN_TITLE", help = "列标题。默认：Sample", default = "Sample")
parser$add_argument("--cluster_row", dest = "CLUSTER_ROW", help = "是否行聚类。默认：true", 
                    choices = c("true", "false"), default = "true")
parser$add_argument("--show_row_dend", dest = "SHOW_ROW_DEND", help = "是否显示行聚类树。默认：true", 
                    choices = c("true", "false"), default = "true")
parser$add_argument("--cluster_column", dest = "CLUSTER_COLUMN", help = "是否列聚类。默认：false", 
                    choices = c("true", "false"), default = "false")
parser$add_argument("--show_column_dend", dest = "SHOW_COLUMN_DEND", help = "是否显示列聚类树。默认：false",
                    choices = c("true", "false"), default = "false")
parser$add_argument("--show_row_names", dest = "SHOW_ROW_NAMES", help = "是否显示行名。默认：true", 
                    choices = c("true", "false"), default = "true")
parser$add_argument("--row_names_site", dest = "ROW_NAMES_SITE", help = "行名位置。默认：left", 
                    choices = c("left", "right"), default = "left")
parser$add_argument("--show_column_names", dest = "SHOW_COLUMN_NAMES", help = "是否显示列名。默认：true", 
                    choices = c("true", "false"), default = "true")
parser$add_argument("--column_names_site", dest = "COLUMN_NAMES_SITE", help = "行名位置。默认：top", 
                    choices = c("top", "bottom"), default = "top")
parser$add_argument("--pdf_width", dest = "PDF_WIDTH", help = "pdf 宽度，默认：7", default = 7)
parser$add_argument("--pdf_height", dest = "PDF_HEIGHT", help = "pdf 高度，默认：7", default = 7)


argvs <- parser$parse_args()
expr_path <- file.path(argvs$EXPR)
output_path <- file.path(argvs$OUTPUT)
normalize <- argvs$NORMALIZE
gene_list <- argvs$LIST
column_title <- argvs$COLUMN_TITLE
legend_name <- argvs$LEGEND
cluster_row <- argvs$CLUSTER_ROW
present_row_dend <- argvs$SHOW_ROW_DEND
cluster_col <- argvs$CLUSTER_COLUMN
present_col_dend <- argvs$SHOW_COLUMN_DEND
present_row_names <- argvs$SHOW_ROW_NAMES
row_name_position <- argvs$ROW_NAMES_SITE
present_col_names <- argvs$SHOW_COLUMN_NAMES
col_name_position <- argvs$COLUMN_NAMES_SITE
pdf_width <- as.integer(argvs$PDF_WIDTH)
pdf_height <- as.integer(argvs$PDF_HEIGHT)

expr_data <- read.csv(expr_path, header = TRUE, row.names = 1) %>% as.matrix()
if (normalize == "no") {
  print("原始输入数据画图")
} else if (normalize == "log2") {
  print("输入数据进行 log2(M + 1) 转换后画图")
  expr_data <- log2(expr_data + 1)
} else if (normalize == "z-score") {
  print("输入数据进行 Z-Score 转换后画图")
  expr_data <- t(scale(t(expr_data)))
} else {
  print("--normalize 参数错误，请检查命令")
  q(save = "no")
}

# 如果指定了基因列表
if (!is.null(gene_list)) {
  gene_path <- file.path(gene_list)
  gene_list <- scan(file = gene_path, what = "character", sep = "\n", quiet = TRUE)
  expr_data <- expr_data[which(rownames(expr_data) %in% gene_list),]
}

hm_color <- viridis(200, direction = -1)
hm_plot <- Heatmap(matrix = expr_data, col = hm_color, name = legend_name, column_title = column_title, 
                   column_title_side = "top", cluster_rows = switch(cluster_row, "true" = TRUE, "false" = FALSE), 
                   show_row_dend = switch(present_row_dend, "true" = TRUE, "false" = FALSE), 
                   cluster_columns = switch(cluster_col, "true" = TRUE, "false" = FALSE), 
                   show_column_dend = switch(present_col_dend, "true" = TRUE, "false" = FALSE), 
                   show_row_names = switch(present_row_names, "true" = TRUE, "false" = FALSE), 
                   show_column_names = switch (present_col_names, "true" = TRUE, "false" = FALSE), 
                   row_names_side = row_name_position, column_names_side = col_name_position)

pdf(output_path, width = pdf_width, height = pdf_height)
draw(hm_plot)
dev.off()

q(save = "no")
