suppressPackageStartupMessages(library(this.path))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))


src_dir <- this.path::this.dir()
source(file.path(src_dir, "shareobj.R"))

# remove gene version
short_id <- function(x) {
  id_version <- strsplit(x, split = ".", fixed = TRUE) %>% unlist()
  ensembl_id <- id_version[1]
  return(ensembl_id)
}

x <- "run DESeq2"
parser <- ArgumentParser(description = x, add_help = TRUE)
parser$add_argument("--output_dir", dest = "OUTDIR", default = "./", 
                    help = "output directory, Default: ./")
parser$add_argument("--salmon", dest = "SALMON", default = "salmon", 
                    help = "salmon results directory. Default: ./salmon")
parser$add_argument("--group", dest = "GROUP", default = "sample_group.csv", 
                    help = "group information in CSV format, first column is sample name, second column is group. Default: ./sample_group.csv")
parser$add_argument("--control", dest = "CONTROL", default = NULL, 
                    help = "specific which group is control.")
parser$add_argument("--counts", dest = "COUNTS", default = NULL, 
                    help = "read counts matrix in CSV format. This will overwrite --salmon")
parser$add_argument("--genome", dest = "GENOME", default = "GRCh38", 
                    help = "genome, GRCh38 or GRCm39. Default: GRCh38")
parser$add_argument("--tx2gene", dest = "TX2GENE", required = TRUE, 
                    help = "mapping of transcripts and genes in csv format, first column is transcript, second column is gene.  Create by MakeTx2GeneFromGTF.R from GTF/GFF3")
argvs <- parser$parse_args()

output_dir <- file.path(argvs$OUTDIR)
salmon_dir <- file.path(argvs$SALMON)
group_path <- file.path(argvs$GROUP)
tx2gene_path <- file.path(argvs$TX2GENE)
ctr_setting <- argvs$CONTROL
genome <- argvs$GENOME

check_path(group_path, quit = TRUE)
check_path(output_dir, mkdir = TRUE)

# no header
group_info <- read.csv(group_path, header = TRUE, stringsAsFactors = TRUE)
colnames(group_info) <- c("sample", "treatment")
if (is.null(ctr_setting)) {
  grp_lev <- levels(group_info$treatment)
  ctr_group <- grp_lev[1]
  test_group <- grp_lev[2]
} else {
  ctr_group <- ctr_setting
  group_info$treatment <- relevel(group_info$treatment, ref = ctr_group)
  test_group <- levels(group_info$treatment)[2]
}

x <- table(group_info$treatment)
if (any(x < 3)) {
  log_msg("ERROR", "DESeq2 requires at least 3 repeat samples for each group.")
  quit(status = 1)
}

count_input <- argvs$COUNTS
if (is.null(count_input)) {
  check_path(salmon_dir, quit = TRUE)
  salmon_paths <- file.path(salmon_dir, group_info$sample, "quant.sf")
  names(salmon_paths) <- group_info$sample
  
  for (each_path in salmon_paths) {
    check_path(each_path, quit = TRUE)
  }
  check_path(tx2gene_path, quit = TRUE)
  tx_gene <- read.csv(tx2gene_path, header = TRUE, quote = "", stringsAsFactors = FALSE)
  
  txp <- tximport(salmon_paths, type = "salmon", tx2gene = tx_gene)
  tpm <- txp$abundance
  tpm <- tpm[rowSums(tpm) > 0, ]
  
  dds <- DESeqDataSetFromTximport(txp, colData = group_info, design = ~ treatment)
} else {
  count_path <- file.path(count_input)
  check_path(count_path, quit = TRUE)
  
  cm <- read.csv(count_path, header = TRUE, row.names = 1, stringsAsFactors = FALSE)[, group_info$sample]
  dds <- DESeqDataSetFromMatrix(cm, colData = group_info, design = ~ treatment)
}

# filter out 0 reads genes
keep_genes <- rowSums(counts(dds)) > 0
dds <- dds[keep_genes,]
dds <- DESeq(dds)

# output paths
output_paths <- file.path(output_dir, c("read_counts.csv", "normalized_counts.csv", 
                                      "RLog_blind.csv", "RLog.csv", "DEGs_raw_all.csv", 
                                      "DEGs_shrink_all.csv", "DEGs_shrink_filter.csv", 
                                      "MAplot.pdf"))
names(output_paths) <- c("read_counts", "normalized_counts", "RLog_blind", "RLog", 
                         "DEGs_raw_all", "DEGs_shrink_all", "DEGs_shrink_filter", 
                         "MAplot")

gene_list <- rownames(dds) %>% 
  sapply(short_id, USE.NAMES = FALSE)
if (genome == "GRCh38") {
  ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  id_map <- getBM(mart = ensembl, values = gene_list, filters = "ensembl_gene_id", 
                  attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol")) %>% 
    as_tibble() %>% 
    arrange(entrezgene_id, desc(hgnc_symbol)) %>% 
    distinct(ensembl_gene_id, .keep_all = TRUE)
  
  if (is.null(count_input)) {
    tpm <- as_tibble(tpm, rownames = "gene_id") %>% 
      separate(col = "gene_id", into = c("ensembl_gene_id", "ensembl_gene_version"), 
                                         extra = "drop", fill = "right", sep = "\\.") %>% 
      select(- ensembl_gene_version) %>% 
      left_join(id_map, by = "ensembl_gene_id") %>% 
      select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())
    tpm_path <- file.path(output_dir, "TPM.csv")
    write_csv(tpm, tpm_path)
  }
  
  raw_counts <- counts(dds, normalized = FALSE) %>% 
    as_tibble(rownames = "gene_id") %>% 
    separate(col = "gene_id", into = c("ensembl_gene_id", "ensembl_gene_version"), 
                                       extra = "drop", fill = "right", sep = "\\.") %>% 
    select(- ensembl_gene_version) %>% 
    left_join(id_map, by = "ensembl_gene_id") %>% 
    select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())
  write_csv(raw_counts, output_paths["read_counts"])
  nor_counts <- counts(dds, normalized = TRUE) %>% 
    as_tibble(rownames = "gene_id") %>% 
    separate(col = "gene_id", into = c("ensembl_gene_id", "ensembl_gene_version"), 
                                       extra = "drop", fill = "right", sep = "\\.") %>% 
    select(- ensembl_gene_version) %>% 
    left_join(id_map, by = "ensembl_gene_id") %>% 
    select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())
  write_csv(nor_counts, output_paths["normalized_counts"])
  rlog1 <- rlog(dds, blind = TRUE) %>% 
    assay() %>% 
    as_tibble(rownames = "gene_id") %>% 
    separate(col = "gene_id", into = c("ensembl_gene_id", "ensembl_gene_version"), 
                                       extra = "drop", fill = "right", sep = "\\.") %>% 
    select(- ensembl_gene_version) %>% 
    left_join(id_map, by = "ensembl_gene_id") %>% 
    select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())
  write_csv(rlog1, output_paths["RLog_blind"])
  rlog2 <- rlog(dds, blind = FALSE) %>% 
    assay() %>% 
    as_tibble(rownames = "gene_id") %>% 
    separate(col = "gene_id", into = c("ensembl_gene_id", "ensembl_gene_version"), 
                                       extra = "drop", fill = "right", sep = "\\.") %>% 
    select(- ensembl_gene_version) %>% 
    left_join(id_map, by = "ensembl_gene_id") %>% 
    select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())
  write_csv(rlog2, output_paths["RLog"])
  
  raw_res <- results(dds, contrast = c("treatment", test_group, ctr_group))
  snk_res <- lfcShrink(dds, coef = str_glue("treatment_{test_group}_vs_{ctr_group}"), 
                       type = "apeglm", res = raw_res)
  
  pdf(output_paths["MAplot"], width = 9)
  plotMA(raw_res, ylim = c(-5, 5), main = "LFC raw")
  plotMA(snk_res, ylim = c(-5, 5), main = "LFC shrink")
  dev.off()
  
  raw_degs <- as_tibble(raw_res, rownames = "gene_id") %>% 
    separate(col = "gene_id", into = c("ensembl_gene_id", "ensembl_gene_version"), 
                                       extra = "drop", fill = "right", sep = "\\.") %>% 
    select(- ensembl_gene_version) %>% 
    left_join(id_map, by = "ensembl_gene_id") %>% 
    select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())
  write_csv(raw_degs, output_paths["DEGs_raw_all"])
  snk_degs1 <- as_tibble(snk_res, rownames = "gene_id") %>% 
    separate(col = "gene_id", into = c("ensembl_gene_id", "ensembl_gene_version"), 
                                       extra = "drop", fill = "right", sep = "\\.") %>% 
    select(- ensembl_gene_version) %>% 
    left_join(id_map, by = "ensembl_gene_id") %>% 
    select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())
  write_csv(snk_degs1, output_paths["DEGs_shrink_all"])
  snk_degs2 <- as_tibble(snk_res, rownames = "gene_id") %>% 
    separate(col = "gene_id", into = c("ensembl_gene_id", "ensembl_gene_version"), 
                                       extra = "drop", fill = "right", sep = "\\.") %>% 
    select(- ensembl_gene_version) %>% 
    left_join(id_map, by = "ensembl_gene_id") %>% 
    select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything()) %>% 
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1)
  write_csv(snk_degs2, output_paths["DEGs_shrink_filter"])
  
} else if (genome == "GRCm39") {
  ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
  id_map1 <- getBM(mart = ensembl, values = gene_list, filters = "ensembl_gene_id", 
                  attributes = c("ensembl_gene_id", "mgi_id", "mgi_symbol", "entrezgene_id")) %>% 
    as_tibble() %>% 
    arrange(desc(mgi_id), desc(mgi_symbol)) %>% 
    distinct(ensembl_gene_id, .keep_all = TRUE)
  
  id_map2 <- getBM(mart = ensembl, attributes = c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene"), 
                     values = id_map1$ensembl_gene_id, filters = "ensembl_gene_id") %>% 
    as_tibble() %>% 
    filter(hsapiens_homolog_ensembl_gene != "") %>% 
    distinct(ensembl_gene_id, .keep_all = TRUE)
  
  ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  id_map3 <- getBM(mart = ensembl, attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"), 
                     values = id_map2$hsapiens_homolog_ensembl_gene, 
                     filters = "ensembl_gene_id") %>% 
    as_tibble() %>% 
    arrange(entrezgene_id, desc(hgnc_symbol)) %>% 
    distinct(ensembl_gene_id, .keep_all = TRUE) %>% 
    rename("hsapiens_homolog_ensembl_gene" = ensembl_gene_id, 
                  "hsapiens_homolog_entrez_gene" = entrezgene_id, 
                  "hsapiens_homolog_hgnc_symbol" = hgnc_symbol)
  id_map <- left_join(id_map1, id_map2, by = "ensembl_gene_id") %>% 
    left_join(id_map3, by = "hsapiens_homolog_ensembl_gene")
  
  if (is.null(count_input)) {
    tpm <- as_tibble(tpm, rownames = "gene_id") %>% 
      separate(col = "gene_id", into = c("ensembl_gene_id", "ensembl_gene_version"), 
                                         extra = "drop", fill = "right", sep = "\\.") %>% 
      select(- ensembl_gene_version) %>% 
      left_join(id_map, by = "ensembl_gene_id") %>% 
      select(ensembl_gene_id, mgi_id, mgi_symbol, entrezgene_id, hsapiens_homolog_ensembl_gene, 
             hsapiens_homolog_entrez_gene, hsapiens_homolog_hgnc_symbol, everything())
    tpm_path <- file.path(output_dir, "TPM.csv")
    write_csv(tpm, tpm_path)
  }
  
  raw_counts <- counts(dds, normalized = FALSE) %>% 
    as_tibble(rownames = "gene_id") %>% 
    separate(col = "gene_id", into = c("ensembl_gene_id", "ensembl_gene_version"), 
                                       extra = "drop", fill = "right", sep = "\\.") %>% 
    select(- ensembl_gene_version) %>% 
    left_join(id_map, by = "ensembl_gene_id") %>% 
    select(ensembl_gene_id, mgi_id, mgi_symbol, entrezgene_id, hsapiens_homolog_ensembl_gene, 
           hsapiens_homolog_entrez_gene, hsapiens_homolog_hgnc_symbol, everything())
  write_csv(raw_counts, output_paths["read_counts"])
  nor_counts <- counts(dds, normalized = TRUE) %>% 
    as_tibble(rownames = "gene_id") %>% 
    separate(col = "gene_id", into = c("ensembl_gene_id", "ensembl_gene_version"), 
                                       extra = "drop", fill = "right", sep = "\\.") %>% 
    select(- ensembl_gene_version) %>% 
    left_join(id_map, by = "ensembl_gene_id") %>% 
    select(ensembl_gene_id, mgi_id, mgi_symbol, entrezgene_id, hsapiens_homolog_ensembl_gene, 
           hsapiens_homolog_entrez_gene, hsapiens_homolog_hgnc_symbol, everything())
  write_csv(nor_counts, output_paths["normalized_counts"])
  rlog1 <- rlog(dds, blind = TRUE) %>% 
    assay() %>% 
    as_tibble(rownames = "gene_id") %>% 
    separate(col = "gene_id", into = c("ensembl_gene_id", "ensembl_gene_version"), 
                                       extra = "drop", fill = "right", sep = "\\.") %>% 
    select(- ensembl_gene_version) %>% 
    left_join(id_map, by = "ensembl_gene_id") %>% 
    select(ensembl_gene_id, mgi_id, mgi_symbol, entrezgene_id, hsapiens_homolog_ensembl_gene, 
           hsapiens_homolog_entrez_gene, hsapiens_homolog_hgnc_symbol, everything())
  write_csv(rlog1, output_paths["RLog_blind"])
  rlog2 <- rlog(dds, blind = FALSE) %>% 
    assay() %>% 
    as_tibble(rownames = "gene_id") %>% 
    separate(col = "gene_id", into = c("ensembl_gene_id", "ensembl_gene_version"), 
                                       extra = "drop", fill = "right", sep = "\\.") %>% 
    select(- ensembl_gene_version) %>% 
    left_join(id_map, by = "ensembl_gene_id") %>% 
    select(ensembl_gene_id, mgi_id, mgi_symbol, entrezgene_id, hsapiens_homolog_ensembl_gene, 
           hsapiens_homolog_entrez_gene, hsapiens_homolog_hgnc_symbol, everything())
  write_csv(rlog2, output_paths["RLog"])
  
  raw_res <- results(dds, contrast = c("treatment", test_group, ctr_group))
  snk_res <- lfcShrink(dds, coef = str_glue("treatment_{test_group}_vs_{ctr_group}"), 
                       type = "apeglm", res = raw_res)
  pdf(output_paths["MAplot"], width = 9)
  plotMA(raw_res, ylim = c(-5, 5), main = "LFC raw")
  plotMA(snk_res, ylim = c(-5, 5), main = "LFC shrink")
  dev.off()
  
  raw_degs <- as_tibble(raw_res, rownames = "gene_id") %>%
    separate(col = "gene_id", into = c("ensembl_gene_id", "ensembl_gene_version"), 
                                       extra = "drop", fill = "right", sep = "\\.") %>% 
    select(- ensembl_gene_version) %>% 
    left_join(id_map, by = "ensembl_gene_id") %>% 
    select(ensembl_gene_id, mgi_id, mgi_symbol, entrezgene_id, hsapiens_homolog_ensembl_gene, 
           hsapiens_homolog_entrez_gene, hsapiens_homolog_hgnc_symbol, everything())
  write_csv(raw_degs, output_paths["DEGs_raw_all"])
  snk_degs1 <- as_tibble(snk_res, rownames = "gene_id") %>%
    separate(col = "gene_id", into = c("ensembl_gene_id", "ensembl_gene_version"), 
                                       extra = "drop", fill = "right", sep = "\\.") %>% 
    select(- ensembl_gene_version) %>% 
    left_join(id_map, by = "ensembl_gene_id") %>% 
    select(ensembl_gene_id, mgi_id, mgi_symbol, entrezgene_id, hsapiens_homolog_ensembl_gene, 
           hsapiens_homolog_entrez_gene, hsapiens_homolog_hgnc_symbol, everything())
  write_csv(snk_degs1, output_paths["DEGs_shrink_all"])
  snk_degs2 <- as_tibble(snk_res, rownames = "gene_id") %>%
    separate(col = "gene_id", into = c("ensembl_gene_id", "ensembl_gene_version"), 
                                       extra = "drop", fill = "right", sep = "\\.") %>% 
    select(- ensembl_gene_version) %>% 
    left_join(id_map, by = "ensembl_gene_id") %>% 
    select(ensembl_gene_id, mgi_id, mgi_symbol, entrezgene_id, hsapiens_homolog_ensembl_gene, 
           hsapiens_homolog_entrez_gene, hsapiens_homolog_hgnc_symbol, everything()) %>% 
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1)
  write_csv(snk_degs2, output_paths["DEGs_shrink_filter"])
} else {
  log_msg("ERROR", "genome must be GRCh38 or GRCm39!")
  quit(status = 1)
}
