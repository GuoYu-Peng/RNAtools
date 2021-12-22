suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))

what_plot <- "DESeq2 差异基因火山图"
parser <- ArgumentParser(description = what_plot, add_help = TRUE)
parser$add_argument("--DEGs", dest = "DEGs", help = "csv 格式 DESeq2 差异基因分析结果", required = TRUE)
parser$add_argument("--output_dir", dest = "OUTDIR", help = "火山图输出目录，脚本保存 pdf 和 png 格式图片。默认：当前目录", default = ".")
parser$add_argument("--prefix", dest = "PREFIX", help = "输出文件名前缀。默认：NoName", default = "NoName")
parser$add_argument("--plot_title", dest = "TITLE", help = "火山图标题。默认：Differential Expression Genes", 
                    default = "Differential Expression Genes")
parser$add_argument("--pvalue_cutoff", dest = "PVAL", help = "差异基因 P 值阈值。默认：0.05", default = 0.05)
parser$add_argument("--log2FC_cutoff", dest = "LOG2FC", help = "差异基因 log2 差异倍数绝对值阈值。默认：1", default = 1)
parser$add_argument("--max_y", dest = "MAXY", help = "Y 轴最大值，超过的点将被调整到此值。默认：12", default = 12)
parser$add_argument("--max_x", dest = "MAXX", help = "X 轴绝对值最大值，超过的点将被调整到此值。默认：5", default = 5)
parser$add_argument("--color1", dest = "COL1", help = "差异基因的颜色。默认：#f2300e", default = "#f2300e")
parser$add_argument("--color2", dest = "COL2", help = "P 值显著但差异倍数小基因的颜色。默认：#0c775e", default = "#0c775e")
parser$add_argument("--color3", dest = "COL3", help = "非差异基因的颜色。默认：#352749", default = "#352749")
parser$add_argument("--plot_width", dest = "PLOT_WIDTH", help = "图片宽度（mm），默认：180", default = 180)
parser$add_argument("--plot_height", dest = "PLOT_HEIGHT", help = "图片高度（mm），默认：200", default = 200)


argvs <- parser$parse_args()
input_path <- file.path(argvs$DEGs)
output_dir <- file.path(argvs$OUTDIR)
prefix <- argvs$PREFIX
plot_title <- argvs$TITLE
max_y <- as.double(argvs$MAXY)
max_x <- as.double(argvs$MAXX)
pval_cutoff <- as.double(argvs$PVAL)
fc_cutoff <- as.double(argvs$LOG2FC)
color1 <- argvs$COL1
color2 <- argvs$COL2
color3 <- argvs$COL3
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

plot_data <- read_csv(input_path) %>% 
  filter(!is.na(padj)) %>% 
  mutate(x = map2_dbl(log2FoldChange, max_x, adjust_x), y = map2_dbl(padj, max_y, adjust_y), 
         data_shape = pmap_chr(list(padj, log2FoldChange, max_x, max_y), shape_value), 
         data_color = map2_chr(log2FoldChange, padj, get_color))

# 不显示 Legend
# 设置 expand 让图像框不覆盖点
volcano_plot <- ggplot(plot_data, aes(x, y)) + 
    geom_point(aes(colour = data_color, shape = data_shape), show.legend = FALSE) + 
	  scale_colour_manual(values = c("col1" = color1, "col2" = color2, "col3" = color3)) + 
    scale_shape_manual(values = c("circle" = "circle", "triangle" = "triangle")) +
	  geom_vline(xintercept = c(- fc_cutoff, fc_cutoff), linetype = "dashed", size = 0.7) + 
	  geom_hline(yintercept = - log10(pval_cutoff), linetype = "dashed", size = 0.7) + 
	  labs(x = "Fold Change(log2)", y = "P value(-log10)", title = plot_title) + 
    scale_x_continuous(limits = c(- max_x, max_x), expand = expansion(mult = c(0.005))) +
    scale_y_continuous(limits = c(0, max_y), expand = expansion(c(0, 0.005))) +
    theme(panel.grid = element_blank(), 
          panel.background = element_rect(fill = "white", colour = "black", linetype = "solid", size = 1.2))


png_file <- paste(prefix, "Volcano.png", sep = "_")
pdf_file <- paste(prefix, "Volcano.pdf", sep = "_")
ggsave(filename = png_file, plot = volcano_plot, dpi = 600, device = "png", path = output_dir, 
       width = plot_width, height = plot_height, units = "mm")
ggsave(filename = pdf_file, plot = volcano_plot, device = "pdf", path = output_dir, 
       width = plot_width, height = plot_height, units = "mm")

writeLines("\no(^▽^)o")
