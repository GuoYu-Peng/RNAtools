suppressPackageStartupMessages(library(this.path))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(tidyverse))


src_dir <- this.path::this.dir()
source(file.path(src_dir, "shareobj.R"))

x <- "Geneset Over-Representation Analysis"
parser <- ArgumentParser(description = x, add_help = TRUE)
parser$add_argument("--DEGs", dest = "DEGS", default = "DEGs.csv", 
                    help = "DEGs result from DESeq2 in csv format. Default: ./DEGs.csv")
parser$add_argument("--outdir", dest = "OUTDIR", default = "geneset_ORA", 
                    help = "output directory. Default: ./geneset_ORA")
parser$add_argument("--prefix", dest = "PREFIX", default = "DEGs", 
                    help = "output prefix. Default: DEGs")
parser$add_argument("--target", dest = "TARGET", default = "KEGG,GOBP", 
                    help = "target gene sets, comma separate if multiple genesets. Default: KEGG,GOBP. Optional: KEGG,GOBP,GOMF,GOCC")
parser$add_argument("--test_all", default = NULL, action="store_true", 
                    help = "test all gene sets, this option will overwrite --target")
parser$add_argument("--genome", dest = "GENOME", default = "GRCh38", 
                    help = "genome, GRCh38 or GRCm39. Default: GRCh38")
argvs <- parser$parse_args()

DEGs_path <- file.path(argvs$DEGS)
output_dir <- file.path(argvs$OUTDIR)
prefix <- argvs$PREFIX
targets <- argvs$TARGET
genome <- argvs$GENOME

check_path(DEGs_path, quit = TRUE)
check_path(output_dir, quit = FALSE, mkdir = TRUE)

if (!is.null(argvs$test_all)) {
  targets <- "KEGG,GOBP,GOMF,GOCC"
}
target_list <- strsplit(targets, split = ",", fixed = TRUE) %>% 
  unlist()
full_target <- c("KEGG", "GOBP", "GOMF", "GOCC")
keep_target <- c()
rm_target <- c()
for (each_target in target_list) {
  if (each_target %in% full_target) {
    keep_target <- c(keep_target, each_target)
  } else {
    rm_target <- c(rm_target, each_target)
  }
}
# check gene sets
if (length(keep_target) == 0) {
  writeLines("\n[ERROR] all gene sets illegal:")
  print(rm_target)
  quit(status = 1)
}
if (length(rm_target) != 0) {
  writeLines("\n[WARN] unknown gene sets:")
  print(rm_target)
}
writeLines("\n[INFO] gene set list:")
print(keep_target)

check_DEGs <- function(x) {
  n <- nrow(x)
  if (n == 0) {
    writeLines("\n[ERROR] no DEGs under \"padj < 0.05 &abs(log2FoldChange) >= 1\" condiction")
    quit(status = 1)
  } else {
    t <- paste("\n[INFO] DEGs after filter: ", n, sep = "")
    writeLines(t)
  }
}

run_KEGG <- function(entrez_list, entrez_bg, output_dir, ...) {
  writeLines("\n[INFO] testing KEGG...")
  csv_path <- make_path(dest_dir = output_dir, suffix = "csv", ...)
  rds_path <- make_path(dest_dir = output_dir, suffix = "rds", ...)
  
  ora_kegg <- enrichKEGG(gene = entrez_list, organism = "hsa", keyType = "ncbi-geneid", 
                         universe = entrez_bg, minGSSize = 15, maxGSSize = 300, pvalueCutoff = 0.1)
  
  if (is.null(ora_kegg)) {
    writeLines("\n[WARN] no KEGG enrichment results")
  } else {
    kegg_tb <- as_tibble(ora_kegg@result)
    write_csv(kegg_tb, csv_path)
    saveRDS(ora_kegg, rds_path)
  }
}

run_GO <- function(entrez_list, entrez_bg, output_dir, ont = "BP", ...) {
  info_msg <- paste("\n[INFO] testing GO", ont, "...", sep = "")
  writeLines(info_msg)
  csv_path <- make_path(dest_dir = output_dir, suffix = "csv", ...)
  rds_path <- make_path(dest_dir = output_dir, suffix = "rds", ...)
  
  ora_go1 <- enrichGO(gene = entrez_list, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", 
                      ont = ont, universe = entrez_bg, minGSSize = 15, maxGSSize = 250, 
                      readable = TRUE, pvalueCutoff = 0.1)
  
  if (is.null(ora_go1)) {
    writeLines("\n[WARN] no GO enrichment results")
  } else {
    ora_go2 <- clusterProfiler::simplify(ora_go1)
    saveRDS(ora_go2, rds_path)
    go_tb <- as_tibble(ora_go2@result)
    write_csv(go_tb, csv_path)
  }
}

if (genome == "GRCh38") {
  DEGs_all <- read_csv(DEGs_path) %>% 
    dplyr::filter(!is.na(entrezgene_id))
  DEGs_filter <- dplyr::filter(DEGs_all, !is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1)
  check_DEGs(DEGs_filter)
  
  DEGs_list <- DEGs_filter$entrezgene_id %>% 
    unique()
  bg_list <- DEGs_all$entrezgene_id %>% 
    unique() %>% 
    as.character()    # background list must be chr
} else if (genome == "GRCm39") {
  DEGs_all <- read_csv(DEGs_path) %>% 
    dplyr::filter(!is.na(hsapiens_homolog_entrez_gene))
  DEGs_filter <- dplyr::filter(DEGs_all, !is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1)
  check_DEGs(DEGs_filter)
  
  DEGs_list <- DEGs_filter$hsapiens_homolog_entrez_gene %>% 
    unique()
  # background list must be chr
  bg_list <- DEGs_all$hsapiens_homolog_entrez_gene %>% 
    unique() %>% 
    as.character()    
} else {
  writeLines("\n[ERROR] unsupported genome:")
  print(genome)
  quit(status = 1)
}

for (each_target in keep_target) {
  if (each_target == "KEGG") {
    run_KEGG(DEGs_list, bg_list, output_dir, prefix, "ORA", "KEGG")
  } else if (each_target == "GOBP") {
    run_GO(DEGs_list, bg_list, output_dir, ont = "BP", prefix, "ORA", "GOBP")
  } else if (each_target == "GOMF") {
    run_GO(DEGs_list, bg_list, output_dir, ont = "MF", prefix, "ORA", "GOMF")
  } else if (each_target == "GOCC") {
    run_GO(DEGs_list, bg_list, output_dir, ont = "CC", prefix, "ORA", "GOCC")
  } else {
    writeLines("\n[WARN] unknown gene set:")
    print(each_target)
  }
}
