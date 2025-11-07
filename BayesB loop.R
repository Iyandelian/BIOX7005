library(BGLR)
library(SNPRelate)
library(gdsfmt)
library(SeqArray)
library(dplyr)
library(foreach)
library(doParallel)


#make sure the GDS snp and csv phenotype files are in the right folders, accessible from the directory containing the R script

# get GDS file
gds_files <- list.files("Delivered_GDSs", pattern = "\\.gds$", full.names = TRUE, recursive = TRUE)

# Get phenotype file
phenotype_files <- list.files("Delivered_phenotypes", pattern = "values\\.csv$", full.names = TRUE, recursive = TRUE)


print(gds_files)
print(phenotype_files)

# Extract phenotype names directly from the file names, assuming they are in the format: 'phenotype_values.csv'
phenotype_names <- gsub("_values\\.csv$", "", basename(phenotype_files))


# Loop through each GDS file
for (gds_file in gds_files) {
  gds_name <- basename(gds_file)  # Get the GDS file name
  print(gds_name)
  # Open the GDS file
  cat("\nProcessing GDS file:", gds_name, "\n")
  genomedata <- snpgdsOpen(gds_file)
  snp_ids <- read.gdsn(index.gdsn(genomedata, "snp.id"))
  sample_ids <- read.gdsn(index.gdsn(genomedata, "sample.id"))
  GMdata <- snpgdsGetGeno(genomedata, snp.id = snp_ids, verbose = TRUE)
  GMdata <- t(GMdata)
  rownames(GMdata) <- sample_ids
  colnames(GMdata) <- snp_ids
  snpgdsClose(genomedata)
  
  
  if (grepl("full", gds_name)) {
    phenotypes_to_process <- phenotype_names  # Process all phenotypes
  } else {
    matches_logical <- lapply(phenotype_names, function(p) grepl(p, gds_name))
    phenotypes_to_process <- phenotype_names[unlist(matches_logical)]
  }
  print(phenotypes_to_process)
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
    
    # Match accessions between phenotype and genotype data
    accessions <- phenotype_data$accession_id
    common_accessions <- accessions %in% rownames(GMdata)
    accessions <- accessions[common_accessions]
    
    if (length(accessions) == 0) next  # Skip if no common accessions
    
    
    # Subset GM matrix for this phenotype
    #GM_pheno <- GMdata[as.character(accessions), ]
    
    # Set seed for reproducibility
    set.seed(400)
    
    # Phenotype values (dependent variable)
    # --- Alignment check and fix ---
    common <- intersect(rownames(GMdata), phenotype_data$accession_id)
    GM_pheno <- GMdata[common, , drop = FALSE]
    phenotype_data <- phenotype_data[match(common, phenotype_data$accession_id), ]
    if(!all(rownames(GM_pheno) == phenotype_data$accession_id)) stop("Alignment failed")

    y <- phenotype_data$phenotype_value
    # --- End alignment block ---
    
    n <- length(y)  # number of samples
    
    for (j in 1:ncol(GM_pheno)) {
      na_idx <- which(is.na(GM_pheno[, j]))
      if (length(na_idx) > 0) {
        GM_pheno[na_idx, j] <- mean(GM_pheno[, j], na.rm = TRUE)
      }
    }
    
    folds <- sample(rep(1:10, length.out = n))  # assign 10 CV folds
    
    # Prepare storage for results
    y_pred <- rep(NA, n)               # to hold predictions
    fold_correlations <- numeric(10)   # to hold correlation per fold
    
    # Detect the number of available cores and use all but one
    #num_cores <- detectCores() - 1
    num_cores <- 10
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    
    # Perform 10-fold cross-validation
    foreach (fold = 1:10, .packages = c("BGLR"), .combine = c) %dopar% {
      cat("\nRunning fold", fold, "for phenotype", phenotype_name, "in GDS file", gds_name, "...\n")
      
      # Define training/test sets
      tst <- which(folds == fold)       # test indices
      y_cv <- y                         # copy phenotype vector
      y_cv[tst] <- NA                   # mask test phenotypes
      
      # Define ETA for GBLUP (linear kernel)
      ETA <- list(list(X = GM_pheno, model = "BayesB"))  # select phenotype
      
      # Fit the model
      fit <- BGLR(
        y = y_cv,
        ETA = ETA,
        nIter = 5000,
        burnIn = 1000,
        verbose = FALSE
      )
      
      # The return value of each parallel process is a list
      list(y_pred_fold = fit$yHat[tst], tst = tst)
    } -> predictions_and_cors
    
    
    # The `.combine='c'` argument combines the list of lists into a single list
    # where every other element is the y_pred_fold and the other is the tst indices.
    for (i in 1:10) {
      fold_pred <- predictions_and_cors[[2*i - 1]] # Odd-indexed elements are y_pred_fold
      fold_tst <- predictions_and_cors[[2*i]]      # Even-indexed elements are tst
      y_pred[fold_tst] <- fold_pred
      
      # Calculate and store fold-wise correlation
      fold_correlations[i] <- cor(y[fold_tst], y_pred[fold_tst], use = "complete.obs")
      cat("Fold", i, "correlation:", round(fold_correlations[i], 4), "\n")
    }
    stopCluster(cl)
    
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
    write.csv(cv_results, paste0("BayesB", gds_name, "_", phenotype_name, "_results.csv"), row.names = FALSE)
  }
  
  # Close the GDS file after processing
  #snpgdsClose(genomedata)
}
