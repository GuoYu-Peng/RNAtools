suppressPackageStartupMessages(library(this.path))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(viridis))


src_dir <- this.path::this.dir()
source(file.path(src_dir, "shareobj.R"))

x <- "Create heatmap from expression matrix"
parser <- ArgumentParser(description = x, add_help = TRUE)
parser$add_argument("--expression", dest = "EXPR", required = TRUE, 
                    help = "expression matrix in csv format, the first column should be gene ID")
parser$add_argument("--output", dest = "OUTPUT", default = "GeneHeatmap.pdf", 
                    help = "heatmap output, default: GeneHeatmap.pdf")
parser$add_argument("--normalize", dest = "NORMALIZE", choices = c("no", "log2", "z-score"), default = "no", 
                    help = "normalize data before plot, default: no")
parser$add_argument("--list", dest = "LIST", default = NULL, 
                    help = "select genes from file, per gene one line, default: false")
parser$add_argument("--legend_name", dest = "LEGEND", default = "Normalized Expression", 
                    help = "legend title, default: Normalized Expression")
parser$add_argument("--column_title", dest = "COLUMN_TITLE", default = "Sample", 
                    help = "column title, default: Sample")
parser$add_argument("--cluster_row", dest = "CLUSTER_ROW", choices = c("true", "false"), 
                    default = "true", help = "cluster row, default: true")
parser$add_argument("--show_row_dend", dest = "SHOW_ROW_DEND", choices = c("true", "false"), default = "true", 
                    help = "show row dendrogram, default: true")
parser$add_argument("--cluster_column", dest = "CLUSTER_COLUMN", choices = c("true", "false"), default = "false", 
                    help = "clust column, default: false")
parser$add_argument("--show_column_dend", dest = "SHOW_COLUMN_DEND", choices = c("true", "false"), default = "false", 
                    help = "show column dendrogram, default: false")
parser$add_argument("--show_row_names", dest = "SHOW_ROW_NAMES", choices = c("true", "false"), default = "true", 
                    help = "show row name, default: true")
parser$add_argument("--row_names_site", dest = "ROW_NAMES_SITE", choices = c("left", "right"), default = "left", 
                    help = "row name side, default: left")
parser$add_argument("--show_column_names", dest = "SHOW_COLUMN_NAMES", choices = c("true", "false"), default = "true", 
                    help = "show column name, default: true")
parser$add_argument("--column_names_site", dest = "COLUMN_NAMES_SITE", choices = c("top", "bottom"), default = "top", 
                    help = "column name side, default: top")
parser$add_argument("--pdf_width", dest = "PDF_WIDTH", default = 7, 
                    help = "pdf width, default: 7")
parser$add_argument("--pdf_height", dest = "PDF_HEIGHT", default = 7, 
                    help = "pdf height, default: 7")


argvs <- parser$parse_args()
expr_path <- file.path(argvs$EXPR)
output_path <- file.path(argvs$OUTPUT)
normalize <- argvs$NORMALIZE
list_setting <- argvs$LIST
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
  log_msg("INFO", "raw input plot")
} else if (normalize == "log2") {
  log_msg("INFO", "log2 normalize plot")
  expr_data <- log2(expr_data + 1)
} else if (normalize == "z-score") {
  log_msg("INFO", "z-score normalize plot")
  expr_data <- t(scale(t(expr_data)))
} else {
  log_msg("ERROR", "invalid --normalize:", normalize)
  q(status = 1)
}

# 如果指定了基因列表
if (!is.null(list_setting)) {
  gene_path <- file.path(list_setting)
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
