# 利用 clusterProfiler 取得通路 gmt 文件
# 分别取 KEGG 和 GO 通路
# 同时保存一个合并的，包含 KEGG, GO 全部通路

library(org.Hs.eg.db)
library(clusterProfiler)

genes <- sample(keys(org.Hs.eg.db, keytype = "ENTREZID"), size = 1200)
kegg_enrich <- enrichKEGG(gene = genes, organism = "hsa", keyType = "ncbi-geneid", 
                          minGSSize = 1, maxGSSize = 2000)
go_enrich <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, ont = "ALL", 
                      minGSSize = 1, maxGSSize = 2000)

kegg_sets <- kegg_enrich@geneSets
go_sets <- go_enrich@geneSets

write_gmt <- function(gene_sets, output_path, check_exists = TRUE) {
  if (check_exists) {
    stopifnot("GMT 文件已存在" = !file.exists(output_path))
  }
  
  n <- length(gene_sets)
  set_id <- names(gene_sets)
  for (i in 1:n) {
    id <- set_id[i]
    genes <- gene_sets[[i]]
    entry <- paste(paste(c(id, "NA", genes), sep = "\t", collapse = "\t"), "\n", sep = "")
    cat(entry, file = output_path, append = TRUE)
  }
}

prefix <- gsub("-", "", Sys.Date(), fixed = TRUE)
kegg_path <- paste(prefix, "KEGG.gmt", sep = "_")
go_path <- paste(prefix, "GO.gmt", sep = "_")
merged_path <- paste(prefix, "Merged.gmt", sep = "_")

write_gmt(kegg_sets, kegg_path)
write_gmt(go_sets, go_path)
write_gmt(kegg_sets, merged_path)
write_gmt(go_sets, merged_path, check_exists = FALSE)