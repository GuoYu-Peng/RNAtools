suppressPackageStartupMessages(library(this.path))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ggrepel))

src_dir <- this.path::this.dir()
source(file.path(src_dir, "shareobj.R"))

x <- "sample distance of RNAseq"
parser <- ArgumentParser(description = x, add_help = TRUE)
parser$add_argument("--expression", dest = "EXPR", required = TRUE, 
                    help = "expression data in CSV format")
parser$add_argument("--group", dest = "GROUP", required = TRUE, 
                    help = "group information in CSV format, first column is sample name.")
parser$add_argument("--DEGs", dest = "DEGs", required = TRUE, 
                    help = "filtered DESeq2 DEGs results in CSV format")
parser$add_argument("--output", dest = "OUTPUT", default = "SampleDIST.pdf", 
                    help = "output, default: SampleDIST.pdf")
parser$add_argument("--column", dest = "COLUMN", default = "Group", 
                    help = "plot grouping column, default: Group")
parser$add_argument("--slice", dest = "SLICE", default = NULL, 
                    help = "select sample groups, use commas to separate multiple groups")
parser$add_argument("--dist_method", dest = "DIST_METHOD", default = "euclidean", 
                    help = "distance method, default: euclidean")


argvs <- parser$parse_args()
expr_path <- file.path(argvs$EXPR)
group_path <- file.path(argvs$GROUP)
output_path <- file.path(argvs$OUTPUT)
deg_path <- file.path(argvs$DEGs)
need_column <- argvs$COLUMN
need_value <- argvs$SLICE
dist_method <- argvs$DIST_METHOD

# ------ 样本分组信息和表达数据 ------
# 保留需要的列和行
sample_group <- read.csv(group_path, header = TRUE, row.names = 1) %>% 
  subset(select = need_column)

if (!is.null(need_value)) {
  keep_value <- strsplit(need_value, split = ",", fixed = TRUE)
  keep_row <- sample_group[, 1] %in% keep_value
  sample_group <- sample_group[keep_row, , drop = FALSE]
}
sample_group <- as_tibble(sample_group, rownames = "SampleName")

# 表达数据只保留需要的基因和样本
deg_data <- read_csv(deg_path)
expr_data <- read_csv(expr_path) %>% 
  dplyr::select(ensembl_gene_id, all_of(sample_group$SampleName)) %>% 
  dplyr::filter(ensembl_gene_id %in% deg_data$ensembl_gene_id) %>% 
  dplyr::distinct(ensembl_gene_id, .keep_all = TRUE) %>% 
  as.data.frame()
rownames(expr_data) <- expr_data$ensembl_gene_id
expr_data$ensembl_gene_id <- NULL

# ------ PCA ------
data_pca <- prcomp(t(expr_data))
pca_plot_data <- data_pca$x %>% 
  as_tibble(rownames = "SampleName") %>% 
  dplyr::select(SampleName, PC1, PC2) %>% 
  left_join(sample_group, by = "SampleName")

# 取得 PC1,PC2 解释占比
pca_summary <- summary(data_pca)
summary_tb <- pca_summary$importance
pc1_prot <- summary_tb[[2, 1]] %>% round(3)
pc2_prot <- summary_tb[[2, 2]] %>% round(3)
x_lab <- str_glue("PC1({pc1_prot})")
y_lab <- str_glue("PC2({pc2_prot})")

pca_plot <- ggplot(pca_plot_data, aes(x = PC1, y = PC2)) +
           geom_point(aes_string(shape = need_column)) +
           ggrepel::geom_text_repel(aes(label = SampleName)) +
           labs(title = "PCA of DEGs Expression", x = x_lab, y = y_lab) +
           theme_bw()


# ------ clust + heatmap ------
data_dist <- dist(t(expr_data), method = dist_method, diag = TRUE, upper = TRUE)
hm_data <- as.matrix(data_dist)

data_clust <- hclust(data_dist, method = "ward.D2")
hm_name <- toupper(str_glue("{dist_method}"))
hm_color <- inferno(200)
hm_plot <- Heatmap(hm_data, name = hm_name, col = hm_color, cluster_rows = TRUE, 
              cluster_columns = TRUE, show_row_names = TRUE, row_names_side = "right", 
              show_row_dend = FALSE, show_column_dend = FALSE, column_title = "Sample Distance", 
              column_title_side = "top")


pdf(output_path, width = 9, height = 9)
print(pca_plot)
plot(data_clust, xlab = "Sample", main = "Sample Clusters")
draw(hm_plot)
dev.off()
