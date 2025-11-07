library(BGLR)
library(SNPRelate)
library(gdsfmt)
library(SeqArray)
library(dplyr)

# make sure the GDS SNP and CSV phenotype files are in the right folders, accessible from the directory containing the R script

# get GDS file
gds_files <- list.files("/Delivered_GDSs", pattern = "\\.gds$", full.names = TRUE)

# Get phenotype file
phenotype_files <- list.files("/Delivered_phenotypes", pattern = "values\\.csv$", full.names = TRUE)


# Extract phenotype names directly from the file names, assuming they are in the format: 'phenotype_values.csv'
phenotype_names <- gsub("_values\\.csv$", "", basename(phenotype_files))


# Loop through each GDS file
for (gds_file in gds_files) {
  gds_name <- basename(gds_file)  # Get the GDS file name
  
  # Open the GDS file
  cat("\nProcessing GDS file:", gds_name, "\n")
  genomedata <- snpgdsOpen(gds_file)
  GRMdata <- snpgdsGRM(
    genomedata,
    autosome.only = FALSE,
    remove.monosnp = FALSE,
    maf = NaN,
    missing.rate = NaN,
    method = c("GCTA"),
    num.thread = 5,
    with.id = TRUE,
    verbose = TRUE
  )
  
  # Extract GRM matrix
  GRM_Matrix <- GRMdata$grm
  rownames(GRM_Matrix) <- GRMdata$sample.id
  colnames(GRM_Matrix) <- GRMdata$sample.id
  
  # Determine phenotypes to process based on the GDS file name (partial matching)
  if (grepl("full", gds_name)) {
    phenotypes_to_process <- phenotype_names  # Process all phenotypes
  } else {
    phenotypes_to_process <- grep(gds_name, phenotype_names, value = TRUE)  # Process phenotypes that match part of the GDS file name
  }
  
  # Loop through each phenotype to process
  for (phenotype_name in phenotypes_to_process) {
    cat("\nProcessing phenotype:", phenotype_name, "for GDS file:", gds_name, "\n")
    
    # Find the correct phenotype file based on the phenotype name (partial matching)
    phenotype_file <- grep(phenotype_name, phenotype_files, value = TRUE)
    if (length(phenotype_file) == 0) {
      cat("No phenotype file found for", phenotype_name, "in", gds_name, "\n")
      next
    }
    
    # Read phenotype data
    phenotype_data <- read.csv(phenotype_file)
    
    # Accessions for this phenotype
    accessions <- phenotype_data$accession_id
    common_accessions <- accessions %in% GRMdata$sample.id
    accessions_common <- accessions[common_accessions]
    
    # Subset GRM matrix for this phenotype
    GRM_subset <- GRM_Matrix[as.character(accessions_common), as.character(accessions_common)]
    
    # Set seed for reproducibility
    set.seed(400)
    
    # Phenotype values (dependent variable)
    y <- phenotype_data$phenotype_value[common_accessions]
    n <- length(y)  # number of samples
    folds <- sample(rep(1:10, length.out = n))  # assign 10 CV folds
    
    # Prepare storage for results
    y_pred <- rep(NA, n)               # to hold predictions
    fold_correlations <- numeric(10)   # to hold correlation per fold
    
    # Perform 10-fold cross-validation
    for (fold in 1:10) {
      cat("\nRunning fold", fold, "for phenotype", phenotype_name, "in GDS file", gds_name, "...\n")
      
      # Define training/test sets
      tst <- which(folds == fold)       # test indices
      y_cv <- y                         # copy phenotype vector
      y_cv[tst] <- NA                   # mask test phenotypes
      
      # Define ETA for GBLUP (linear kernel)
      ETA <- list(list(K = GRM_subset, model = "RKHS"))  # select phenotype
      
      # Fit the model
      fit <- BGLR(
        y = y_cv,
        ETA = ETA,
        nIter = 5000,
        burnIn = 1000,
        verbose = FALSE
      )
      
      # Store predictions for test set
      y_pred[tst] <- fit$yHat[tst]
      
      # Calculate and store fold-wise correlation
      fold_cor <- cor(y[tst], fit$yHat[tst], use = "complete.obs")
      fold_correlations[fold] <- fold_cor
      cat("Fold", fold, "correlation:", round(fold_cor, 4), "\n")
    }
    
    # Compile results
    cv_results <- data.frame(
      name = phenotype_data$accession_name[common_accessions],  # select phenotype
      y = y,
      yhat = y_pred,
      fold = folds
    )
    
    # Evaluation metrics
    overall_cor <- cor(cv_results$y, cv_results$yhat, use = "complete.obs")
    rmse <- (mean((cv_results$y - cv_results$yhat)^2))
    
    cat("Overall correlation (10-fold CV) for phenotype", phenotype_name, "in GDS file", gds_name, ":", round(overall_cor, 4), "\n")
    cat("Overall MSE for phenotype", phenotype_name, "in GDS file", gds_name, ":", rmse, "\n")
    
    # Plot observed vs predicted for this phenotype
    plot(cv_results$y, cv_results$yhat,
         xlab = "Observed", ylab = "Predicted",
         main = paste("10-Fold GBLUP: Observed vs Predicted for", phenotype_name, "in", gds_name))
    abline(0, 1, col = "red", lty = 2)
    
    # Summary of fold-wise correlations
    fold_summary <- data.frame(
      Fold = 1:10,
      Correlation = round(fold_correlations, 4)
    )
    print(fold_summary)
    
    # Optionally: save results
    write.csv(cv_results, paste0("GBLUP_", gds_name, "_", phenotype_name, "_results.csv"), row.names = FALSE)
  }
  
  # Close the GDS file after processing
  snpgdsClose(genomedata)
}
