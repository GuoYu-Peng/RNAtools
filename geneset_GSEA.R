suppressPackageStartupMessages(library(this.path))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(tidyverse))


src_dir <- this.path::this.dir()
source(file.path(src_dir, "shareobj.R"))

x <- "GSEA of DEGs from DESeq2"
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
                    help = "run all gene sets, this option will overwrite --target")
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


run_KEGG <- function(rank_list, output_dir, ...) {
  writeLines(log_msg(level = "INFO", "testing KEGG..."))
  csv_path <- make_path(dest_dir = output_dir, suffix = "csv", ...)
  rds_path <- make_path(dest_dir = output_dir, suffix = "rds", ...)
  
  gsea_kegg <- gseKEGG(geneList = rank_list, organism = "hsa", keyType = "ncbi-geneid", 
                       minGSSize = 15, maxGSSize = 300, eps = 0, pvalueCutoff = 0.1)
  
  if (is.null(gsea_kegg)) {
    writeLines("\n[WARN] no KEGG GSEA results")
  } else {
    kegg_tb <- as_tibble(gsea_kegg@result)
    write_csv(kegg_tb, csv_path)
    saveRDS(gsea_kegg, rds_path)
  }
}

run_GO <- function(rank_list, output_dir, ont = "BP", ...) {
  writeLines(log_msg(level = "INFO", "testing GO", ont, "..."))
  csv_path <- make_path(dest_dir = output_dir, suffix = "csv", ...)
  rds_path <- make_path(dest_dir = output_dir, suffix = "rds", ...)
  
  gsea_go1 <- gseGO(geneList = rank_list, ont = ont, OrgDb = org.Hs.eg.db, 
                    keyType = "ENTREZID", minGSSize = 15, maxGSSize = 300, 
                    eps = 0, pvalueCutoff = 0.1)
  if (is.null(gsea_go1)) {
    warn_msg <- paste("\n[WARN] no GO", ont, " GSEA results", sep = "")
    writeLines(warn_msg)
  } else {
    gsea_go2 <- clusterProfiler::simplify(gsea_go1)
    go_tb <- as_tibble(gsea_go2@result)
    write_csv(go_tb, csv_path)
    saveRDS(gsea_go2, rds_path)
  }
}

if (genome == "GRCh38") {
  DEGs <- read_csv(DEGs_path) %>% 
    dplyr::filter(!is.na(entrezgene_id)) %>% 
    dplyr::arrange(desc(hgnc_symbol)) %>% 
    dplyr::distinct(entrezgene_id, .keep_all = TRUE) %>%
    arrange(desc(log2FoldChange))
  fc_rank <- DEGs$log2FoldChange
  names(fc_rank) <- DEGs$entrezgene_id
} else if (genome == "GRCm39") {
  DEGs <- read_csv(DEGs_path) %>% 
    dplyr::filter(!is.na(hsapiens_homolog_entrez_gene)) %>% 
    dplyr::arrange(desc(hsapiens_homolog_hgnc_symbol)) %>% 
    dplyr::distinct(hsapiens_homolog_entrez_gene, .keep_all = TRUE) %>%
    arrange(desc(log2FoldChange))
  fc_rank <- DEGs$log2FoldChange
  names(fc_rank) <- DEGs$hsapiens_homolog_entrez_gene
} else {
  err_msg <- paste("\n[ERROR] unsupported genome: ", genome, sep = "")
  writeLines(err_msg)
}
writeLines("\n[INFO] ranked list head and tail:")
print(head(fc_rank))
print(tail(fc_rank))

for (each_target in keep_target) {
  if (each_target == "KEGG") {
    run_KEGG(fc_rank, output_dir, prefix, "GSEA", "KEGG")
  } else if (each_target == "GOBP") {
    run_GO(fc_rank,output_dir, "BP", prefix, "GSEA", "GOBP")
  } else if (each_target == "GOMF") {
    run_GO(fc_rank,output_dir, "MF", prefix, "GSEA", "GOMF")
  } else if (each_target == "GOCC") {
    run_GO(fc_rank,output_dir, "CC", prefix, "GSEA", "GOCC")
  } else {
    writeLines(log_msg(level = "WARN", "unsupported gene set: ", each_target))
  }
}
