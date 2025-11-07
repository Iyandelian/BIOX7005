library(SNPRelate)
library(SeqArray)

vcffile <- "/home/blake/Documents/BIOX7005/new variants/variants/Arabidopsis_2029_Maf005_Filter95.vcf.bgz"

seqVCF2GDS(vcffile, "new_SNPs.gds", verbose = TRUE, parallel = 5)

gdssetfile <- "new_SNPs.gds"

openbook <- seqOpen(gdssetfile)

filtered <- seqSetFilterCond(openbook, maf=0.05, missing.rate=0.1, parallel=5, .progress=TRUE, verbose=TRUE)

seqExport(openbook, "filtered_file.gds")

seqClose(openbook)

gdsfilteredfile <- "filtered_file.gds"

openbook <- seqOpen(gdsfile)

seqGDS2SNP(gdssetfile, "2prt1.gds", dosage=FALSE, optimize=TRUE, verbose=TRUE)

genofile <- snpgdsOpen("2prt1.gds")

LDP_SNPs <- snpgdsLDpruning(genofile, autosome.only=FALSE, remove.monosnp=FALSE, maf=NaN, missing.rate=NaN, method="corr", slide.max.bp=20000, slide.max.n=NaN, ld.threshold=0.8, start.pos="random", num.thread=5, autosave=NULL, verbose=TRUE)


genofile <- snpgdsOpen("2prt1.gds", readonly = FALSE)
pruned_ids <- unlist(LDP_SNPs, use.names = FALSE)
add.gdsn(genofile, "ld_pruned_snp.id", pruned_ids)
snpgdsClose(genofile)

as  r2 = 0.10, a window size of 200 kb, and a step size of 20 kb