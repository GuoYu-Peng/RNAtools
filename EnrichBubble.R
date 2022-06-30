suppressPackageStartupMessages(library(this.path))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))


src_dir <- this.path::this.dir()
source(file.path(src_dir, "shareobj.R"))

x <- "Scatter plot of geneset enrich analysis by clusterProfiler"
parser <- ArgumentParser(description = what_plot, add_help = TRUE)
parser$add_argument("--enrich_result", dest = "ENRICHMENT", required = TRUE, 
                    help = "enrichment results in csv format")
parser$add_argument("--output", dest = "OUTPUT", default = "Bubble.pdf", 
                    help = "output path, default: Bubble.pdf")
parser$add_argument("--plot_title", dest = "TITLE", default = "Pathway Enrichment", 
                    help = "plot title, default: Pathway Enrichment")
parser$add_argument("--show_number", dest="NUMBER", 
                    help="display top NUMBER genesets, default: 25", 
                    default = 25)
parser$add_argument("--list", dest="LIST", default = NULL, 
                    help="only display genesets from the file, each geneset one line, this will overwrite --show_number")
parser$add_argument("--row_name_width", dest="ROWWIDTH", default = 40, 
                    help = "wrap row name more than this character width, default: 40")
parser$add_argument("--p_min", dest = "P_MIN", default = 0, help = "minimal p-value on plot legend, default: 0")
parser$add_argument("--p_max", dest = "P_MAX", default = 0.1, help = "maximal p-value on plot legend, default: 0.1")
parser$add_argument("--p_min_color", dest = "P_MIN_COLOR", default = "#890011",
                    help = "color of min p value, default：#890011")
parser$add_argument("--p_max_color", dest = "P_MAX_COLOR", default = "#f87c5e", 
                    help = "color of max p-value, default：#f87c5e")
parser$add_argument("--plot_width", dest = "PLOT_WIDTH", default = 200, 
                    help = "plot width in mm, default: 200")
parser$add_argument("--height_scale", dest = "HEIGHT_SCALE", default = 8, 
                    help = "plot height (mm) = (geneset numbers) * HEIGHT_SCALE, default: 8")

argvs <- parser$parse_args()
enrich_path <- file.path(argvs$ENRICHMENT)
output_path <- file.path(argvs$OUTPUT)
plot_title <- argvs$TITLE
list_setting <- argvs$LIST
show_number <- as.integer(argvs$NUMBER)
row_name_width <- as.integer(argvs$ROWWIDTH)
p_min <- as.numeric(argvs$P_MIN)
p_max <- as.numeric(argvs$P_MAX)
p_min_color <- argvs$P_MIN_COLOR
p_max_color <- argvs$P_MAX_COLOR
plot_width <- as.integer(argvs$PLOT_WIDTH)
height_scale <- as.numeric(argvs$HEIGHT_SCALE)

check_path(enrich_path, TRUE, FALSE)
enrich_result <- read_csv(enrich_path) %>% 
  arrange(`p.adjust`) %>% 
  separate(GeneRatio, into=c("k", "n"), sep="/") %>% 
  separate(BgRatio, into=c("M", "N"), sep="/") %>% 
  mutate(RichFactor=as.numeric(k)/as.numeric(M))

# keep genesets
if (is.null(list_setting)) {
  enrich_result <- slice(enrich_result, 1:show_number)
} else {
  keep_pathway <- scan(file = file.path(list_para), what = character(), sep = "\n")
  enrich_result <- filter(enrich_result, ID %in% keep_pathway)
}

# RichFactor
enrich_result <- arrange(enrich_result, RichFactor) %>% 
  mutate(Description=factor(Description, levels = Description))
glimpse(enrich_result)

bubble_plot <- ggplot(enrich_result, aes(x = Description, y = RichFactor)) +
  geom_point(aes(size = Count, colour = `p.adjust`), show.legend = TRUE) +
  scale_colour_gradient(low = p_min_color, high = p_max_color, limits = c(p_min, p_max)) +
  scale_x_discrete(labels=scales::wrap_format(row_name_width)) +
  labs(x = "Pathway", y = "RichFactor", title = plot_title, 
       size = "Gene count", colour = "P value") +
  theme(axis.title.y = element_blank(), 
        panel.background = element_rect(fill = "#FFFFFF", color = "#000000", linetype = "solid"), 
        panel.grid.major = element_line(color = "#DCDCDC", linetype = "solid"), 
        axis.text.x = element_text(size = 14)) +
  coord_flip()

plot_height = nrow(enrich_result) * height_scale
ggsave(filename = output_path, plot = bubble_plot, device = "pdf", width = plot_width, 
       height = plot_height, units = "mm")
