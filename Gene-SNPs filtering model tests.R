library(ranger)
library(SNPRelate)
library(gdsfmt)
library(SeqArray)
library(ggplot2)
library(rrBLUP)
library(BGLR)
library(dplyr)
library(caret)
library(RhpcBLASctl)


gds_file <- "filtered_snps.gds"
gds_file2 <- 'filtered_file.gds'

# Open the GDS file
genofile <- snpgdsOpen(gds_file)


snp.id <- read.gdsn(index.gdsn(genofile, "snp.id"))
chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
pos <- read.gdsn(index.gdsn(genofile, "snp.position"))

snp_gr <- GRanges(seqnames = chr, ranges = IRanges(start = pos, end = pos))

seqlevels(genomic_rangesr_output) <- sub("^Chr", "", seqlevels(genomic_rangesr_output))


hits <- findOverlaps(snp_gr, genomic_rangesr_output)
filtered_snp_ids <- snp.id[queryHits(hits)]

cat("Number of SNPs to retain:", length(filtered_snp_ids), "\n")

snpgdsClose(genofile)

input_gds <- "filtered_snps.gds"
output_gds <- "in_gene_snps.gds"

cat("Total SNPs in original GDS:", length(snp.id), "\n")
cat("Filtered SNPs:", length(filtered_snp_ids), "\n")
cat("Duplicates in filtered:", sum(duplicated(filtered_snp_ids)), "\n")
cat("Missing in filtered vs original:", length(setdiff(filtered_snp_ids, snp.id)), "\n")

filtered_snp_ids <- unique(filtered_snp_ids)

snpgdsCreateGenoSet(
  src.fn       = input_gds,
  dest.fn      = output_gds,
  sample.id    = NULL,
  snp.id       = filtered_snp_ids,
  snpfirstdim  = NULL,
  compress.annotation = "ZIP_RA.max",
  compress.geno        = "",
  verbose = TRUE
)




