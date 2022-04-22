suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(org.Hs.eg.db))

what_plot <- "GSEA 过程 ES 图"
parser <- ArgumentParser(description = what_plot, add_help = TRUE)
parser$add_argument("--gsea", dest = "GSEA", help = "clusterProfiler GSEA 分析结果（rds）", required = TRUE)
parser$add_argument("--output_dir", dest = "OUTPUT_DIR", help = "输出目录，默认当前目录", default = ".")
parser$add_argument("--list_path", dest = "LIST_PATH", help = "画图的通路列表，每个通路 ID 一行", 
                    default = NULL)
parser$add_argument("--geneset", dest = "GENESET", help = "直接指定通路 ID. 如果多条用逗号 \",\" 分隔", 
                    default = NULL)
# 默认颜色取自 BuenColors 的 brewer 系列
parser$add_argument("--line_color", dest = "LINE_COLOR", help = "ES 线颜色，默认：#74c476", default = "#74c476")
parser$add_argument("--point_color", dest = "POINT_COLOR", help = "核心富集基因点颜色，默认：#fd8d3c", default = "#fd8d3c")
parser$add_argument("--segment_color1", dest = "SEG_COL1", help = "核心富集基因线段颜色，默认：#fb6a4a", default = "#fb6a4a")
parser$add_argument("--segment_color2", dest = "SEG_COL2", help = "非核心富集基因线段颜色，默认：#6baed6", default = "#6baed6")
parser$add_argument("--line_width", dest = "LINE_WIDTH", help = "ES 线宽度，默认：0.6", default = 0.6)
parser$add_argument("--segment_width1", dest = "SEG_WIDTH1", help = "核心富集基因线段宽度，默认：0.4", default = 0.4)
parser$add_argument("--segment_width2", dest = "SEG_WIDTH2", help = "非核心富集基因线段宽度，默认：0.4", default = 0.4)
parser$add_argument("--plot_height", dest = "HEIGHT", help = "图片高度（mm），默认：90", default = 90)
parser$add_argument("--plot_width", dest = "WIDTH", help = "图片宽度（mm），默认：180", default = 180)

argvs <- parser$parse_args()
gsea_path <- file.path(argvs$GSEA)
output_dir <- file.path(argvs$OUTPUT_DIR)
list_para <- argvs$LIST_PATH
geneset_para <- argvs$GENESET
line_color <- argvs$LINE_COLOR
point_color <- argvs$POINT_COLOR
seg_col1 <- argvs$SEG_COL1
seg_col2 <- argvs$SEG_COL2
line_width <- as.numeric(argvs$LINE_WIDTH)
seg_width1 <- as.numeric(argvs$SEG_WIDTH1)
seg_width2 <- as.numeric(argvs$SEG_WIDTH2)
plot_height <- as.numeric(argvs$HEIGHT)
plot_width <- as.numeric(argvs$WIDTH)

gseaScores <- function (geneList, geneSet, exponent = 1, fortify = FALSE) {
  geneSet <- intersect(geneSet, names(geneList))
  N <- length(geneList)
  Nh <- length(geneSet)
  Phit <- Pmiss <- numeric(N)
  hits <- names(geneList) %in% geneSet
  Phit[hits] <- abs(geneList[hits])^exponent
  # 让每次 runing 加分和减分的总分都为 1
  NR <- sum(Phit)
  Phit <- cumsum(Phit/NR)
  Pmiss[!hits] <- 1/(N - Nh)
  Pmiss <- cumsum(Pmiss)
  runningES <- Phit - Pmiss
  max.ES <- max(runningES)
  min.ES <- min(runningES)
  if (abs(max.ES) > abs(min.ES)) {
    ES <- max.ES
  } else {
    ES <- min.ES
  }
  df <- data.frame(x = seq_along(runningES), runningScore = runningES,
                   position = as.integer(hits))
  if (fortify == TRUE) {
    return(df)
  }
  df$gene = names(geneList)
  res <- list(ES = ES, runningES = df)
  return(res)
}

gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList
  
  if (is.numeric(geneSetID)) {
    geneSetID <- object@result[geneSetID, "ID"]
  }
    
  
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df_raw <- gseaScores(geneList, geneSet, exponent, fortify=FALSE)
  ES <- df_raw$ES
  df <- df_raw$runningES
  
  df$ymin <- 0
  df$ymax <- 0
  # hit position
  pos <- df$position == 1
  # h 是总 ES 分差的 1/20
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  
  df$Description <- object@result[geneSetID, "Description"]
  res <- list(ES = ES, runningES = df)
  return(res)
}

p_text <- function(padj) {
  if (padj < 0.0001) {
    ptext <- "P.adjust < 0.0001"
  } else {
    p <- round(padj, digits = 4)
    ptext <- paste("P.adjust", p, sep = " = ")
  }
  return(ptext)
}

need_sets <- function(gsea_result, limit1, limit2) {
  gsea_sets <- gsea_result$ID
  if (!is.null(limit1)) {
    list_path <- file.path(limit1)
    limit_sets <- scan(file = list_path, sep = "\n", what = character()) %>% 
      unique()
  } else if (!is.null(limit2)) {
    limit_sets <- limit2 %>% 
      strsplit(split = ",", fixed = TRUE) %>% 
      unlist()
  } else {
    limit_sets <- gsea_sets
  }
  need_sets <- intersect(gsea_sets, limit_sets)
  return(need_sets)
}

gsea_run <- readRDS(gsea_path)
gsea_result <- gsea_run@result
needsets <- need_sets(gsea_result, list_para, geneset_para)
cat("以下通路在富集结果：\n")
print(needsets)
stopifnot("无可画图基因集" = length(needsets) > 0)


for (each_set in needsets) {
  set_data <- gsInfo(gsea_run, each_set)
  set_result <- filter(gsea_result, ID == each_set)
  leading_genes <- pull(set_result, core_enrichment) %>% 
    strsplit(split = "/", fixed = TRUE) %>% 
    unlist()
  set_desc <- pull(set_result, Description)
  
  # 是个矩阵
  running_es <- set_data$runningES
  
  leading_data <- filter(running_es, position == 1, gene %in% leading_genes) %>% 
    mutate(hgnc_symbol = mapIds(x = org.Hs.eg.db, keys = gene, column = "SYMBOL", keytype = "ENTREZID"))
  rest_data <- filter(running_es, position == 1, !(gene %in% leading_genes))
  last_gene_data <- tail(leading_data, n = 1)
  seg_min <- min(running_es$ymin)
  seg_max <- max(running_es$ymax)
  
  ES <- set_data$ES
  mid_pos <- median(running_es$x)
  set_padj <- pull(set_result, `p.adjust`) %>% 
    p_text()

  es_plot <- ggplot() +
    geom_line(data = running_es, mapping = aes(x = x, y = runningScore), 
              color = line_color, size = line_width) +
    geom_point(data = leading_data, mapping = aes(x = x, y = runningScore), color = seg_col1, 
               size = (seg_width1 - 0.1), shape = 8) +
    geom_text(mapping = aes(x = mid_pos, y = ES, label = set_padj), size = 2.6) +
    geom_segment(data = leading_data, mapping = aes(x = x, y = seg_min, xend = x, yend = seg_max), 
                 color = seg_col1, size = seg_width1) +
    geom_segment(data = rest_data, mapping = aes(x = x, y = seg_min, xend = x, yend = seg_max), 
                 color = seg_col2, size = seg_width2) +
    geom_segment(data = last_gene_data, mapping = aes(x = x, y = 0, xend = x, yend = runningScore), 
                 color = seg_col1, size = seg_width1, linetype = "dashed") +
    geom_hline(yintercept = 0, size = 0.5, color = "#000000", alpha = 0.7) +
    ggrepel::geom_text_repel(data = leading_data, mapping = aes(x = x, y = runningScore, label = hgnc_symbol), 
                             max.overlaps = 20, segment.size = 0.2, size = 1.6, segment.alpha = 0.6) +
    labs(x = "Rank in Ordered Dataset", y = "Running Enrichment Score", title = set_desc) +
    theme_bw()
  
  set_name <- gsub(pattern = "[-:/]", replacement = "_", x = each_set)
  file_name <- paste(set_name, "pdf", sep = ".")
  ggsave(filename = file_name, plot = es_plot, device = "pdf", path = output_dir, 
        width = plot_width, height = plot_height, units = "mm")
}
