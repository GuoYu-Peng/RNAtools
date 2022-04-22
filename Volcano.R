suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))

# 传参
{
what_plot <- "DESeq2 差异基因火山图"
parser <- ArgumentParser(description = what_plot, add_help = TRUE)
parser$add_argument("--DEGs", dest = "DEGs", help = "csv 格式 DESeq2 差异基因分析结果", required = TRUE)
parser$add_argument("--output_dir", dest = "OUTDIR", help = "火山图保存目录。默认：当前目录", default = ".")
parser$add_argument("--prefix", dest = "PREFIX", help = "输出文件名前缀。默认：NoName", default = "NoName")
parser$add_argument("--annot_list", dest = "ANNOT_LIST", help = "每个基因一行的注释基因列表文件", 
                    default = NULL)
parser$add_argument("--annot_gene", dest = "ANNOT_GENE", help = "注释基因，多个基因用逗号 \",\" 分隔", 
                    default = NULL)
parser$add_argument("--label_column", dest = "LABEL_COLUMN", help = "注释基因列名，默认 hgnc_symbol", 
                    default = "hgnc_symbol")
parser$add_argument("--plot_title", dest = "TITLE", help = "火山图标题。默认：Differential Expression Genes", 
                    default = "Differential Expression Genes")
parser$add_argument("--pvalue_cutoff", dest = "PVAL", help = "差异基因 P 值阈值。默认：0.05", default = 0.05)
parser$add_argument("--log2FC_cutoff", dest = "LOG2FC", help = "差异基因 log2 差异倍数绝对值阈值。默认：1", default = 1)
parser$add_argument("--max_y", dest = "MAXY", help = "Y 轴最大值，超过的点将被调整到此值。默认：12", default = 12)
parser$add_argument("--max_x", dest = "MAXX", help = "X 轴绝对值最大值，超过的点将被调整到此值。默认：5", default = 5)
parser$add_argument("--color1", dest = "COL1", help = "差异基因的颜色。默认：#f56B949", default = "#56B949")
parser$add_argument("--color2", dest = "COL2", help = "P 值显著但差异倍数小基因的颜色。默认：#F0A32F", default = "#F0A32F")
parser$add_argument("--color3", dest = "COL3", help = "非差异基因的颜色。默认：#30499B", default = "#30499B")
parser$add_argument("--color4", dest = "COL4", help = "注释基因的颜色。默认：#EE4035", default = "#EE4035")
parser$add_argument("--plot_width", dest = "PLOT_WIDTH", help = "图片宽度（mm），默认：180", default = 180)
parser$add_argument("--plot_height", dest = "PLOT_HEIGHT", help = "图片高度（mm），默认：200", default = 200)


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
}

# 函数
{
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
    cat("画图数据没有要注释基因，注意检查数据\n")
    annot_plot <- raw_plot
  }
  return(annot_plot)
}

degs_summary <- function(x) {
  up_num <- filter(x, padj < pval_cutoff, log2FoldChange > 0) %>% 
    nrow()
  down_num <- filter(x, padj < pval_cutoff, log2FoldChange < 0) %>% 
    nrow()
  cat("\n差异基因总结：\n")
  cat("P 值显著上调基因数 ")
  cat(up_num)
  cat("\n")
  cat("P 值显著下调基因数 ")
  cat(down_num)
  cat("\n")
  
  up_num <- filter(x, padj < pval_cutoff, log2FoldChange >= fc_cutoff) %>% 
    nrow()
  down_num <- filter(x, padj < pval_cutoff, log2FoldChange <= (- fc_cutoff)) %>% 
    nrow()
  cat("P 值和差异倍数符合条件上调基因数 ")
  cat(up_num)
  cat("\n")
  cat("P 值和差异倍数符合条件下调基因数 ")
  cat(down_num)
  cat("\n\n")
}
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
  stopifnot("基因列表路径错误" = file.exists(gene_list_path))
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
