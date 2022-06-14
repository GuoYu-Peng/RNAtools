suppressPackageStartupMessages(library(this.path))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(KEGGREST))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))


src_dir <- this.path::this.dir()
source(file.path(src_dir, "shareobj.R"))

tool_desc <- "将多个通路按重叠基因数目输出到表格适合生成网络图，可添加 GSEA 结果。"
parser <- ArgumentParser(description = tool_desc, add_help = TRUE)
parser$add_argument("--list", dest = "LIST", help = "通路列表文件，每个通路一行，限定为人种 KEGG 或 GO 通路", required = TRUE)
parser$add_argument("--GSEA", dest = "GSEA", help = "csv 格式 GSEA 结果，如果提供将多输出一个表，把 GSEA 数据映射到节点", default = NULL)
parser$add_argument("--output_dir", dest = "OUTDIR", help = "输出目录，默认当前目录", default = ".")
parser$add_argument("--basename", dest = "BASENAME", help = "输出文件前缀，默认 pathway", default = "pathway")
parser$add_argument("--gmt", dest = "GMT", help = "通路 GMT 文件，不提供时脚本会自动获取通路数据", default = NULL)
parser$add_argument("--type", dest = "TYPE", help = "KEGG 或 GO 通路？默认 KEGG", default = "KEGG")
parser$add_argument("--min_percent", dest = "PERCENT", help = "重叠基因数占较小通路最低比例要求，默认 0.25", default = 0.25)

argvs <- parser$parse_args()
list_path <- file.path(argvs$LIST)
stopifnot(file.exists(list_path))
output_dir <- file.path(argvs$OUTDIR)
name_prefix <- argvs$BASENAME
dataset_type <- argvs$TYPE
min_percent <- as.numeric(argvs$PERCENT)

kegg_entrez <- function(x) {
  pathway <- keggGet(x)[[1]]
  pathway_gene <- pathway$GENE
  if (length(pathway_gene) >= 2) {
    pathway_entrez <- pathway_gene[seq(from = 1, to = length(pathway_gene), by = 2)]
  } else {
    pathway_entrez <- ""
  }
  return(pathway_entrez)
}

read_gmt <- function(x) {
  gmt_line <- readLines(x)
  gmt_list <- str_split(gmt_line, pattern = "\t")
  names(gmt_list) <- vapply(gmt_list, function(y) y[1], character(1))
  gmt_list <- lapply(gmt_list, "[", -c(1:2))
  return(gmt_list)
}

pathway_list <- scan(list_path, what = character(), sep = "\n")
pathway_gene <- list()

gmt_status <- argvs$GMT
if (is.null(gmt_status)) {
  if (identical(dataset_type, "KEGG")) {
    for (pathway in pathway_list) {
      entrez_list <- kegg_entrez(pathway)
      pathway_gene[[pathway]] <- entrez_list
    }
    } else {
      for (pathway in pathway_list) {
        entrez_list <- mapIds(x = org.Hs.eg.db, keys = pathway, column = "ENTREZID", 
                              keytype = "GOALL", multiVals = "list") %>% unlist()
        pathway_gene[[pathway]] <- entrez_list
      }
    }
} else {
  # 考虑列表的通路不在 gmt 情况
  gmt_path <- file.path(gmt_status)
  stopifnot(file.exists(gmt_path))
  gmt <- read_gmt(gmt_path)
  gmt_sets <- names(gmt)
  need_sets <- pathway_list[pathway_list %in% gmt_sets]
  pathway_gene <- gmt[need_sets]
}

stopifnot(length(pathway_gene) >= 3)
pathway_list2 <- names(pathway_gene)
nums <- length(pathway_list2)
overlap_count <- c()
overlap_percent <- c()
source_node <- c()
target_node <- c()
# 在出图时作为无向网络处理
for (i in 1:(nums - 1)) {
  i_gene <- pathway_gene[[i]]
  i_name <- pathway_list2[i]
  for (j in (i + 1):nums) {
    j_gene <- pathway_gene[[j]]
    j_name <- pathway_list2[j]
    # 重叠基因数目
    i_j_overlap <- length(intersect(i_gene, j_gene))
    i_j_percent <- i_j_overlap / min(c(length(i_gene), length(j_gene)))
    overlap_percent <- c(overlap_percent, i_j_percent)
    source_node <- c(source_node, i_name)
    target_node <- c(target_node, j_name)
    overlap_count <- c(overlap_count, i_j_overlap)
  }
}
network_raw <- tibble(source = source_node, target = target_node, overlap = overlap_count, 
                      percent_min = overlap_percent)
network_pass <- filter(network_raw, percent_min >= min_percent)
network_fail <- filter(network_raw, percent_min < min_percent)
glimpse(network_pass)
if (nrow(network_fail) >= 1) {
  cat("\n")
  cat("移除以下占比过低条目：\n")
  print(network_fail)
}


node_attrs <- tibble("shared name" = pathway_list2, gene_count = sapply(pathway_gene, length)) %>% 
  filter(`shared name` %in% unique(x = c(network_pass$source, network_pass$target)))

gsea_status <- argvs$GSEA
if (!is.null(gsea_status)) {
  gsea_path <- file.path(gsea_status)
  stopifnot(file.exists(gsea_path))
  gsea <- read_csv(gsea_path) %>% 
    filter(ID %in% pathway_list2)
  node_attrs <- left_join(node_attrs, gsea, by = c("shared name" = "ID")) %>% 
    select(`shared name`, gene_count, Description, enrichmentScore, NES, 
           pvalue, `p.adjust`)
}
glimpse(node_attrs)

out_file1 <- paste(name_prefix, "network.csv", sep = "_")
out_path1 <- file.path(output_dir, out_file1)
out_file2 <- paste(name_prefix, "attributes.csv", sep = "_")
out_path2 <- file.path(output_dir, out_file2)

write_csv(network_pass, out_path1)
write_csv(node_attrs, out_path2)
