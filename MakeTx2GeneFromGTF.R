suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(GenomicFeatures))


x <- "Get mapping between genes and transcripts then save to csv format file."
parser <- ArgumentParser(description = x, add_help = TRUE)
parser$add_argument("--input", dest = "INPUT", help = "genome annotation file, gtf/gff3 format", 
                    required = TRUE)
parser$add_argument("--output", dest = "OUTPUT", help = "output path", required = TRUE)
argvs <- parser$parse_args()

input_path <- file.path(argvs$INPUT)
output_path <- file.path(argvs$OUTPUT)
stopifnot("Input file does not exist!" = file.exists(input_path))

txdb <- makeTxDbFromGFF(input_path, format = "auto")
tx2gene <- AnnotationDbi::select(txdb, keys=keys(txdb, "TXNAME"), 
                                 keytype="TXNAME", columns=c("TXNAME", "GENEID"))
write.csv(tx2gene, file = output_path, quote = FALSE, row.names = FALSE)
