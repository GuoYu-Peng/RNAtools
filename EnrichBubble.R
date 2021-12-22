suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))


what_plot <- "clusterProfiler 通路富集分析结果泡泡图"
parser <- ArgumentParser(description = what_plot, add_help = TRUE)
parser$add_argument("--enrich_result", dest = "ENRICHMENT", help = "csv 格式 clusterProfiler 通路富集结果", required = TRUE)
parser$add_argument("--output", dest = "OUTPUT", help = "图片输出路径，要求为 pdf 格式。默认：Bubble.pdf", 
                    default = "Bubble.pdf")
parser$add_argument("--plot_title", dest = "TITLE", help = "图片标题。默认：Pathway", default = "pathway")
parser$add_argument("--list", dest="LIST", help="只展示文件内通路，每个一行，可选参数", default = NULL)
parser$add_argument("--show_number", dest="NUMBER", 
                    help="展示的通路数目，如果指定了 --list 参数，那么此参数无效。默认 25", default = 25)
parser$add_argument("--row_name_width", dest="ROWWIDTH", help = "通路名宽度，超过此宽度将换行显示，默认 40", default = 40)
parser$add_argument("--p_min", dest = "P_MIN", help = "图例 P 值刻度（小）。默认：0", default = 0)
parser$add_argument("--p_max", dest = "P_MAX", help = "图例 P 值刻度（大）。默认：0.1", default = 0.1)
parser$add_argument("--p_min_color", dest = "P_MIN_COLOR", help = "图例 P 值刻度（小）颜色。默认：#890011", default = "#890011")
parser$add_argument("--p_max_color", dest = "P_MAX_COLOR", help = "图例 P 值刻度（大）颜色。默认：#f87c5e", default = "#f87c5e")
parser$add_argument("--plot_width", dest = "PLOT_WIDTH", help = "图片宽度（mm），默认：180", default = 200)
parser$add_argument("--height_scale", dest = "HEIGHT_SCALE", help = "图片高度缩放比例，默认：1", default = 1)

argvs <- parser$parse_args()
enrich_path <- file.path(argvs$ENRICHMENT)
output_path <- file.path(argvs$OUTPUT)
plot_title <- argvs$TITLE
list_path <- argvs$LIST
show_number <- as.integer(argvs$NUMBER)
row_name_width <- as.integer(argvs$ROWWIDTH)
p_min <- as.numeric(argvs$P_MIN)
p_max <- as.numeric(argvs$P_MAX)
p_min_color <- argvs$P_MIN_COLOR
p_max_color <- argvs$P_MAX_COLOR
plot_width <- as.integer(argvs$PLOT_WIDTH)
height_scale <- as.numeric(argvs$HEIGHT_SCALE)


# 保留需要做图的数据
# 如果给定了通路列表，那么限定通路数目不起作用
enrich_result <- read_csv(enrich_path) %>% 
  arrange(`p.adjust`) %>% 
  separate(GeneRatio, into=c("k", "n"), sep="/") %>% 
  separate(BgRatio, into=c("M", "N"), sep="/") %>% 
  mutate(RichFactor=as.numeric(k)/as.numeric(M))

if (is.null(list_path)) {
  enrich_result <- slice(enrich_result, 1:show_number)
} else {
  keep_pathway <- scan(file = file.path(list_path), what = character(), sep = "\n")
  enrich_result <- filter(enrich_result, ID %in% keep_pathway)
}

# 按照 RichFactor 排列
enrich_result <- arrange(enrich_result, RichFactor) %>% 
  mutate(Description=factor(Description, levels = Description))
glimpse(enrich_result)

bubble_plot <- ggplot(enrich_result, aes(x = Description, y = RichFactor)) +
  geom_point(aes(size = Count, colour = `p.adjust`), show.legend = TRUE) +
  scale_colour_gradient(low = p_min_color, high = p_max_color, limits = c(p_min, p_max)) +
  scale_x_discrete(labels=scales::wrap_format(row_name_width)) +
  labs(x = "Pathway", y = "RichFactor", title = plot_title, size = "Gene count", colour = "P value") +
  theme(axis.title.y = element_blank(), panel.background = element_rect(fill = "#FFFFFF", color = "#000000", linetype = "solid"), 
        panel.grid.major = element_line(color = "#DCDCDC", linetype = "solid"), 
        axis.text.x = element_text(size = 14)) +
  coord_flip()

plot_height = 8 * nrow(enrich_result) * height_scale
ggsave(filename = output_path, plot = bubble_plot, device = "pdf", width = plot_width, 
       height = plot_height, units = "mm")
