library(SNPRelate)
library(data.table)
library(GMMAT)  # add GMMAT library

snpgds_file <- "filtered_snps.gds"
phenotype_files <- c(
  "/home/blake/Documents/BIOX7005/phenotypes/DTF1/DTFvalues.csv",
  "/home/blake/Documents/BIOX7005/phenotypes/CL/CLvalues.csv",
  "/home/blake/Documents/BIOX7005/phenotypes/RL/RLvalues.csv",
  "/home/blake/Documents/BIOX7005/phenotypes/RBN/RBNvalues.csv",
  "/home/blake/Documents/BIOX7005/phenotypes/FRI/FRIvalues.csv",
  "/home/blake/Documents/BIOX7005/phenotypes/RosetteDM/Rosettevalues.csv"
)

output_dir <- "/home/blake/Documents/BIOX7005/gds_filtered_subsets"

for (pheno_path in phenotype_files) {
  
  # --- Load and clean phenotype data ---
  pheno <- fread(pheno_path)
  pheno <- pheno[!is.na(phenotype_value)]
  
  # --- Load GDS sample IDs in true stored order ---
  gds <- snpgdsOpen(snpgds_file)
  all_samples <- read.gdsn(index.gdsn(gds, "sample.id"))
  snpgdsClose(gds)
  
  # --- Restrict to samples with this phenotype only ---
  common_samples <- all_samples[all_samples %in% pheno$accession_id]
  
  # --- Reorder phenotype table to match GDS sample order ---
  pheno_sub <- pheno[match(common_samples, pheno$accession_id), ]
  
  # --- Split THIS phenotype subset into train/test ---
  n_train <- floor(0.7 * length(common_samples))
  set.seed(123)  # Seed here for reproducibility of train/test split
  train_samples <- sample(common_samples, n_train)
  test_samples  <- setdiff(common_samples, train_samples)
  
  # --- Create training GDS containing only phenotype-specific train samples ---
  train_gds <- tempfile(fileext = ".gds")
  snpgdsCreateGenoSet(
    snpgds_file,
    train_gds,
    sample.id = train_samples,
    snpfirstdim = TRUE,
    compress.annotation = "ZIP_RA.max",
    compress.geno = "ZIP_RA",
    verbose = TRUE
  )
  
  # --- Compute GRM from training data ---
  gds_open <- snpgdsOpen(train_gds)
  
  ibs <- snpgdsIBS(
    gds_open,
    autosome.only = FALSE,
    remove.monosnp = FALSE,
    maf = NaN,
    missing.rate = NaN,
    num.thread = 5,
    verbose = TRUE
  )
  
  # --- Extract phenotype for training samples ---
  train_samples_actual <- read.gdsn(index.gdsn(gds_open, "sample.id"))
  pheno_train <- pheno_sub[match(train_samples_actual, pheno_sub$accession_id), ]
  stopifnot(all(train_samples_actual == pheno_train$accession_id))
  
  # --- Extract genotype matrix for training samples ---
  geno_mat <- snpgdsGetGeno(gds_open, sample.id = train_samples_actual, with.id = TRUE)
  snpgdsClose(gds_open)
  
  # geno_mat$genotype is matrix SNP x Sample; transpose for GMMAT: rows = samples, columns = SNPs
  geno_matrix <- t(geno_mat$genotype)
  colnames(geno_matrix) <- geno_mat$snp.id
  rownames(geno_matrix) <- train_samples_actual
  
  # --- Fit null model with GMMAT ---
  pheno_df <- data.frame(id = pheno_train$accession_id,
                         phenotype = pheno_train$phenotype_value)
  
  kinship_mat <- ibs$ibs
  rownames(kinship_mat) <- ibs$sample.id
  colnames(kinship_mat) <- ibs$sample.id
  kinship_mat <- kinship_mat[as.character(pheno_train$accession_id),
                             as.character(pheno_train$accession_id)]
  
  nullmod <- glmmkin(phenotype ~ 1, data = pheno_df, kins = kinship_mat, id = "id", family = gaussian())
  
  # --- Run association test with GMMAT ---
  stopifnot(all(rownames(geno_matrix) == pheno_df$id))
  
  geno_matrix_transposed <- t(geno_matrix)
  
  # Use phenotype-specific unique file names for input/output
  pheno_name <- tools::file_path_sans_ext(basename(pheno_path))
  infile <- paste0("transposed_genotype_matrix_", pheno_name, ".txt")
  outfile <- paste0("output_glmmr_", pheno_name, ".txt")
  
  write.table(geno_matrix_transposed, file = infile, sep = "\t", row.names = TRUE, col.names = NA)
  
  # Create/overwrite the output file
  if (file.exists(outfile)) file.remove(outfile)
  file.create(outfile)
  
  assoc <- glmm.score(nullmod,
                      infile = infile,
                      outfile = outfile,
                      infile.ncol.skip = 1,
                      infile.nrow.skip = 1,
                      verbose = TRUE)
  
  assoc <- read.delim(outfile, header = TRUE, stringsAsFactors = FALSE)
  
  # --- Debug info ---
  cat("Phenotype:", pheno_name, "\n")
  cat("Number of SNPs tested:", nrow(assoc), "\n")
  cat("Number of SNPs with NA p-values:", sum(is.na(assoc$PVAL)), "\n")
  print(summary(assoc$PVAL))
  
  # --- Filter top SNPs ---
  suggestive_threshold <- 1e-2
  top_snps <- assoc$SNP[assoc$PVAL < suggestive_threshold & !is.na(assoc$PVAL)]
  cat("Number of top SNPs (p <", suggestive_threshold, "):", length(top_snps), "\n")
  
  # Calculate lambda genomic inflation factor
  p <- as.numeric(assoc$PVAL)
  p <- p[!is.na(p) & p > 0 & p <= 1]
  lambda <- median(qchisq(1 - p, 1)) / qchisq(0.5, 1)
  cat("λ =", round(lambda, 3), "\n")
  
  # --- Create filtered test GDS using training-selected SNPs ---
  final_gds_path <- file.path(output_dir, paste0(pheno_name, "_filtered_test.gds"))
  
  if (length(top_snps) > 0) {
    snpgdsCreateGenoSet(
      snpgds_file,
      final_gds_path,
      sample.id = test_samples,
      snp.id = as.numeric(top_snps),
      snpfirstdim = TRUE,
      compress.annotation = "ZIP_RA.max",
      compress.geno = "ZIP_RA",
      verbose = TRUE
    )
  } else {
    cat("Warning: No SNPs passed threshold for", pheno_name, "\n")
  }
  
  cat("Completed GWAS training and saved filtered test GDS for:",
      pheno_name, "\nPath:", final_gds_path, "\n\n")
  
  unlink(train_gds)
  
  # Optional: remove temporary input/output files if desired
  # file.remove(infile)
  # file.remove(outfile)
}
