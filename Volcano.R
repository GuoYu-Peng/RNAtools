suppressPackageStartupMessages(library(this.path))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))


src_dir <- this.path::this.dir()
source(file.path(src_dir, "shareobj.R"))

x <- "DESeq2 volcano plot"
parser <- ArgumentParser(description = x, add_help = TRUE)
parser$add_argument("--DEGs", dest = "DEGs", required = TRUE, 
                    help = "DESeq2 DEGs results in CSV format")
parser$add_argument("--output_dir", dest = "OUTDIR", default = "./", 
                    help = "output directory, default: ./")
parser$add_argument("--prefix", dest = "PREFIX", default = "NoName", 
                    help = "output prefix, default: NoName")
parser$add_argument("--annot_list", dest = "ANNOT_LIST", default = NULL, 
                    help = "annotate gene list file, one gene per line")
parser$add_argument("--annot_gene", dest = "ANNOT_GENE", default = NULL, 
                    help = "annotate gene list separated by comma")
parser$add_argument("--label_column", dest = "LABEL_COLUMN", default = "hgnc_symbol", 
                    help = "gene ID column of annotation, default: hgnc_symbol")
parser$add_argument("--plot_title", dest = "TITLE", default = "Differential Expression Genes", 
                    help = "plot title, default: Differential Expression Genes")
parser$add_argument("--pvalue_cutoff", dest = "PVAL", default = 0.05, 
                    help = "p-value cutoff, default: 0.05")
parser$add_argument("--log2FC_cutoff", dest = "LOG2FC", default = 1, 
                    help = "abs(log2FC) cutoff, default: 1")
parser$add_argument("--max_y", dest = "MAXY", default = 12, 
                  help = "max y-axis, default: 12")
parser$add_argument("--max_x", dest = "MAXX", default = 5, 
                    help = "max x-axis, default: 5")
parser$add_argument("--color1", dest = "COL1", default = "#56B949", 
                    help = "color of filtered DEGs, default: #f56B949")
parser$add_argument("--color2", dest = "COL2", default = "#F0A32F", 
                    help = "color of significant DEGs, default: #F0A32F")
parser$add_argument("--color3", dest = "COL3", default = "#30499B", 
                    help = "color of non DEGs, default: #30499B")
parser$add_argument("--color4", dest = "COL4", default = "#EE4035", 
                    help = "annotation color, default: #EE4035")
parser$add_argument("--plot_width", dest = "PLOT_WIDTH", default = 180, 
                    help = "plot width (mm), default: 180")
parser$add_argument("--plot_height", dest = "PLOT_HEIGHT", default = 200, 
                    help = "plot height (mm), default: 200")


argvs <- parser$parse_args()
input_path <- file.path(argvs$DEGs)
annot_list_pm <- argvs$ANNOT_LIST
annot_gene_pm <- argvs$ANNOT_GENE
label_column <- argvs$LABEL_COLUMN
output_dir <- file.path(argvs$OUTDIR)
prefix <- argvs$PREFIX
plot_title <- argvs$TITLE
max_y <- as.double(argvs$MAXY)
max_x <- as.double(argvs$MAXX)
pval_cutoff <- as.double(argvs$PVAL)
fc_cutoff <- abs(as.double(argvs$LOG2FC))
color1 <- argvs$COL1
color2 <- argvs$COL2
color3 <- argvs$COL3
color4 <- argvs$COL4
plot_width <- as.integer(argvs$PLOT_WIDTH)
plot_height <- as.integer(argvs$PLOT_HEIGHT)




# 被压缩的点用三角形，以示区别
shape_value <- function(p_value, log2_foldchange, max_x, max_y) {
  if (- log10(p_value) <= max_y & (log2_foldchange >= - max_x & log2_foldchange <= max_x)) {
    shape <- "circle"
  } else {
    shape <- "triangle"
  }
  return(shape)
}

# 调整 Y 值
adjust_y <- function(p_value, max_y) {
  y <- - log10(p_value)
  if (y > max_y) {
    y <- max_y
  }
  return(y)
}

# 调整 X 值
adjust_x <- function(log2_foldchange, max_x) {
    if (log2_foldchange > max_x) {
    new_log2 <- max_x
  } else if (log2_foldchange < (- max_x)) {
    new_log2 <- - max_x
  } else {
    new_log2 <- log2_foldchange
  }
  return(new_log2)
}

# 颜色
get_color <- function(log2fc, padj) {
  if (abs(log2fc) >= fc_cutoff & padj < pval_cutoff) {
    data_col <- "col1"
  } else if (abs(log2fc) < fc_cutoff & padj < pval_cutoff) {
    data_col <- "col2"
  } else {
    data_col <- "col3"
  }
  return(data_col)
}

# 火山图添加注释基因
annot_plot <- function(raw_plot_data, gene_list, raw_plot) {
  label_data <- filter(raw_plot_data, .data[[label_column]] %in% gene_list)
  if (nrow(label_data) > 0) {
    annot_plot <- raw_plot +
      ggrepel::geom_text_repel(data = label_data, 
                               mapping = aes(x, y, label = .data[[label_column]]), 
                               show.legend = FALSE, 
                               max.time = 1) +
      geom_point(data = label_data, 
                 mapping = aes(x, y), color = color4)
  } else {
    log_msg("WARN", "no annotation genes in DEGs data")
    annot_plot <- raw_plot
  }
  return(annot_plot)
}

degs_summary <- function(x) {
  up_num <- filter(x, padj < pval_cutoff, log2FoldChange > 0) %>% 
    nrow()
  down_num <- filter(x, padj < pval_cutoff, log2FoldChange < 0) %>% 
    nrow()
  log_msg("INFO", "up regulate significant DEGs number:", up_num)
  log_msg("INFO", "down regulate significant DEGs number:", down_num)
  
  up_num <- filter(x, padj < pval_cutoff, log2FoldChange >= fc_cutoff) %>% 
    nrow()
  down_num <- filter(x, padj < pval_cutoff, log2FoldChange <= (- fc_cutoff)) %>% 
    nrow()
  log_msg("INFO", "up regulate filtered DEGs number:", up_num)
  log_msg("INFO", "down regulate filtered DEGs number:", down_num)
}


plot_data <- read_csv(input_path) %>% 
  filter(!is.na(padj)) %>% 
  mutate(x = map2_dbl(log2FoldChange, max_x, adjust_x), y = map2_dbl(padj, max_y, adjust_y), 
         data_shape = pmap_chr(list(padj, log2FoldChange, max_x, max_y), shape_value), 
         data_color = map2_chr(log2FoldChange, padj, get_color))
degs_summary(plot_data)

# 不显示 Legend
# 设置 expand 让图像框不覆盖点
volcano_plot <- ggplot(plot_data, aes(x, y)) + 
    geom_point(aes(colour = data_color, shape = data_shape), show.legend = FALSE) + 
	  scale_colour_manual(values = c("col1" = color1, "col2" = color2, "col3" = color3)) + 
    scale_shape_manual(values = c("circle" = "circle", "triangle" = "triangle")) +
	  geom_vline(xintercept = c(- fc_cutoff, fc_cutoff), linetype = "dashed", size = 0.5) + 
	  geom_hline(yintercept = - log10(pval_cutoff), linetype = "dashed", size = 0.5) + 
	  labs(x = "Fold Change(log2)", y = "P value(-log10)", title = plot_title) + 
    scale_x_continuous(limits = c(- max_x, max_x), expand = expansion(mult = c(0.005))) +
    scale_y_continuous(limits = c(0, max_y), expand = expansion(c(0, 0.005))) +
    theme(panel.grid = element_blank(), 
          panel.background = element_rect(fill = "white", colour = "black", linetype = "solid", size = 1.2))

# 注释
if (!is.null(annot_list_pm)) {
  gene_list_path <- file.path(annot_list_pm)
  check_path(gene_list_path, TRUE, FALSE)
  genes <- scan(gene_list, what = character(), sep = "\n")
  vplot <- annot_plot(plot_data, genes, volcano_plot)
} else if (!is.null(annot_gene_pm)) {
  genes <- str_split(annot_gene_pm, pattern = ",") %>% 
    unlist()
  vplot <- annot_plot(plot_data, genes, volcano_plot)
} else {
  vplot <- volcano_plot
}


png_file <- paste(prefix, "Volcano.png", sep = "_")
pdf_file <- paste(prefix, "Volcano.pdf", sep = "_")
ggsave(filename = png_file, plot = vplot, dpi = 600, device = "png", path = output_dir, 
       width = plot_width, height = plot_height, units = "mm")
ggsave(filename = pdf_file, plot = vplot, device = "pdf", path = output_dir, 
       width = plot_width, height = plot_height, units = "mm")
