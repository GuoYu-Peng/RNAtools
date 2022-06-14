suppressPackageStartupMessages(library(this.path))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))

src_dir <- this.path::this.dir()
source(file.path(src_dir, "shareobj.R"))

what_plot <- "clusterProfiler GSEA 分析结果柱状图"
parser <- ArgumentParser(description = what_plot, add_help = TRUE)
parser$add_argument("--gsea_result", dest = "GSEA", help = "csv 格式 clusterProfiler GSEA 结果", required = TRUE)
parser$add_argument("--output", dest = "OUTPUT", help = "图片输出路径，要求为 pdf 格式。默认：GSEA.pdf", 
                    default = "GSEA.pdf")
parser$add_argument("--plot_title", dest = "TITLE", help = "图片标题。默认：Pathway Enrichment", 
                    default = "Pathway Enrichment")
parser$add_argument("--list", dest="LIST", help="只展示文件内通路，每个一行，可选参数", default = NULL)
parser$add_argument("--show_number", dest="NUMBER", 
                    help="展示的通路数目，如果指定了 --list 参数，那么此参数无效。默认 25", default=25)
parser$add_argument("--row_name_width", dest="ROWWIDTH", help = "通路名宽度，超过此宽度将换行显示，默认 36", default = 36)
parser$add_argument("--p_min", dest = "P_MIN", help = "图例 P 值刻度（小）。默认：0", default = 0)
parser$add_argument("--p_max", dest = "P_MAX", help = "图例 P 值刻度（大）。默认：0.1", default = 0.1)
parser$add_argument("--p_min_color", dest = "P_MIN_COLOR", help = "图例 P 值刻度（小）颜色。默认：#404040", 
                    default = "#404040")
parser$add_argument("--p_max_color", dest = "P_MAX_COLOR", help = "图例 P 值刻度（大）颜色。默认：#bfbfbf", 
                    default = "#bfbfbf")
parser$add_argument("--bar_width", dest = "BAR_WIDTH", help = "方条宽度，默认：0.8", default = 0.8)
parser$add_argument("--plot_width", dest = "PLOT_WIDTH", help = "图片宽度（mm），默认：240", default = 240)
parser$add_argument("--height_scale", dest = "HEIGHT_SCALE", help = "图片高度缩放比例，默认：1", default = 1)

argvs <- parser$parse_args()
gsea_path <- file.path(argvs$GSEA)
output_path <- file.path(argvs$OUTPUT)
plot_title <- argvs$TITLE
list_para <- argvs$LIST
show_number <- as.integer(argvs$NUMBER)
row_name_width <- as.integer(argvs$ROWWIDTH)
p_min <- as.numeric(argvs$P_MIN)
p_max <- as.numeric(argvs$P_MAX)
p_min_color <- argvs$P_MIN_COLOR
p_max_color <- argvs$P_MAX_COLOR
bar_width <- as.numeric(argvs$BAR_WIDTH)
plot_width <- as.integer(argvs$PLOT_WIDTH)
height_scale <- as.numeric(argvs$HEIGHT_SCALE)

# 选取需要的条目
# 按照 P 值排序 排序
gsea_result <- read_csv(gsea_path)%>% arrange(`p.adjust`)
if (is.null(list_para)) {
  gsea_result <- slice(gsea_result, 1:show_number)
} else {
  keep_pathway <- scan(file = file.path(list_para), what = character(), sep = "\n")
  gsea_result <- filter(gsea_result, ID %in% keep_pathway)
}

# 作图用 NES 数据排序
gsea_result <- arrange(gsea_result, NES) %>% 
  mutate(Pathway=factor(Description, levels=Description))
gsea_plot <- ggplot(gsea_result, aes(Pathway, NES)) +
  geom_bar(aes(fill = `p.adjust`), stat = "identity", width = bar_width, color = "#151515") +
  # geom_hline(yintercept = 0, color = "grey") +
  scale_fill_gradient(high = p_max_color, low = p_min_color, limits = c(p_min, p_max)) +
  labs(y = "Normalized Enrichment Score", title = plot_title, x = "Pathway", fill = "P value") +
  scale_x_discrete(labels = scales::wrap_format(row_name_width)) +
  theme(panel.background = element_rect(fill = "#FFFFFF"), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.text = element_text(size = 13), 
        axis.title = element_text(size = 14), axis.line.x.bottom = element_line(color = "grey"), 
        axis.line.y.left = element_line(colour = "grey"), axis.ticks = element_line(colour = "grey")) +
  coord_flip()

plot_height <- 10 * nrow(gsea_result) * height_scale
ggsave(filename = output_path, plot = gsea_plot, device = "pdf", width = plot_width, 
       height = plot_height, units = "mm")
