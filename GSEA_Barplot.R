suppressPackageStartupMessages(library(this.path))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))

src_dir <- this.path::this.dir()
source(file.path(src_dir, "shareobj.R"))

x <- "barplot of clusterProfiler GSEA results"
parser <- ArgumentParser(description = x, add_help = TRUE)
parser$add_argument("--gsea_result", dest = "GSEA", required = TRUE, 
                    help = "clusterProfiler GSEA results in CSV format")
parser$add_argument("--output", dest = "OUTPUT", default = "GSEA.pdf", 
                    help = "output, default: GSEA.pdf")
parser$add_argument("--plot_title", dest = "TITLE", default = "Pathway Enrichment", 
                    help = "plot title, default: Pathway Enrichment")
parser$add_argument("--show_number", dest="NUMBER", default=25, 
                    help="display top NUMBER gene sets, default: 25")
parser$add_argument("--list", dest="LIST", default = NULL, 
                    help="only display genesets from the file, each geneset one line, this will overwrite --show_number")
parser$add_argument("--row_name_width", dest="ROWWIDTH", default = 40, 
                    help = "wrap row name more than this character width, default: 40")
parser$add_argument("--p_min", dest = "P_MIN", default = 0, 
                    help = "minimal p-value on plot legend, default: 0")
parser$add_argument("--p_max", dest = "P_MAX", default = 0.1, 
                    help = "maximal p-value on plot legend, default: 0.1")
parser$add_argument("--p_min_color", dest = "P_MIN_COLOR", default = "#404040", 
                    help = "color of min p value, default: #404040")
parser$add_argument("--p_max_color", dest = "P_MAX_COLOR", default = "#bfbfbf", 
                    help = "color of max p value, default: #bfbfbf")
parser$add_argument("--bar_width", dest = "BAR_WIDTH", default = 0.8, 
                    help = "bar width, default: 0.8")
parser$add_argument("--plot_width", dest = "PLOT_WIDTH", default = 240, 
                    help = "plot width (mm), default: 240")
parser$add_argument("--height_scale", dest = "HEIGHT_SCALE", default = 8, 
                    help = "plot height (mm) = (geneset numbers) * HEIGHT_SCALE, default: 8")

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

plot_height <- nrow(gsea_result) * height_scale
ggsave(filename = output_path, plot = gsea_plot, device = "pdf", width = plot_width, 
       height = plot_height, units = "mm")
