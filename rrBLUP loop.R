library(rrBLUP)
library(SNPRelate)
library(gdsfmt)
library(SeqArray)
library(dplyr)
library(caret)
library(RhpcBLASctl)
library(ggplot2)

# make sure the GDS SNP and CSV phenotype files are in the right folders, accessible from the directory containing the R script
# and if needed, set the directory holding the file to the working directory

# get GDS file
gds_files <- list.files("Delivered_GDSs", pattern = "\\.gds$", full.names = TRUE, recursive = TRUE)

# Get phenotype file
phenotype_files <- list.files("Delivered_phenotypes", pattern = "values\\.csv$", full.names = TRUE, recursive = TRUE)


# Extract phenotype names directly from the file names, assuming they are in the format: 'phenotype_values.csv'
phenotype_names <- gsub("_values\\.csv$", "", basename(phenotype_files))

print(gds_files)
print(phenotype_files)
print(phenotype_names)

# Loop through each GDS file
for (gds_file in gds_files) {
  gds_name <- basename(gds_file)  # Get the GDS file name
  print(gds_name)
  # Open the GDS file
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
  
  GMdata_full <- GMdata
  
  for (phenotype_name in phenotypes_to_process) {
    cat("\nProcessing phenotype:", phenotype_name, "for GDS file:", gds_name, "\n")
    
    GMdata <- GMdata_full  # Reset GMdata to full data each phenotype iteration
    
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
    
    

    # Set seed for reproducibility
    set.seed(400)
    
    # --- Alignment check and fix ---
    common <- intersect(rownames(GMdata), phenotype_data$accession_id)
    GMdata <- GMdata[common, , drop = FALSE]
    phenotype_data <- phenotype_data[match(common, phenotype_data$accession_id), ]
    if(!all(rownames(GMdata) == phenotype_data$accession_id)) stop("Alignment failed")
    Z <- GMdata
    y <- phenotype_data$phenotype_value
    # --- End alignment block ---
    
    
    # Handle missing values in genotype matrix (replace with column mean)
    for (j in 1:ncol(Z)) {
      na_idx <- which(is.na(Z[, j]))
      if (length(na_idx) > 0) {
        Z[na_idx, j] <- mean(Z[, j], na.rm = TRUE)
      }
    }
    
    # Ensure dimensions match
    if (length(y) != nrow(Z)) stop("Mismatch between number of phenotypes and genotypes!")
    
    # Set up k-fold cross-validation
    set.seed(400)
    k <- 10
    folds <- createFolds(seq_along(y), k = min(k, length(y)), list = TRUE)

      
    
    # Initialize storage for predictions
    predicted <- rep(NA, length(y))
    per_fold_accuracy <- numeric(k)
    per_fold_mse <- numeric(k)
    
    # Perform k-fold cross-validation
    for (i in seq_along(folds)) {
      test_idx <- folds[[i]]
      y_train <- y
      y_train[test_idx] <- NA  # Mask test phenotypes
      
      # Fit RR-BLUP model
      model <- mixed.solve(y = y_train, Z = Z)
      
      # Predict for all individuals
      GEBVs <- rep(model$beta, nrow(Z)) + Z %*% model$u
      
      # Store predictions
      predicted[test_idx] <- GEBVs[test_idx]
      
      # Evaluate fold performance
      true_vals <- y[test_idx]
      pred_vals <- GEBVs[test_idx]
      
      per_fold_accuracy[i] <- cor(pred_vals, true_vals)
      per_fold_mse[i] <- mean((pred_vals - true_vals)^2)
      
      cat(sprintf("Fold %d: Accuracy = %.3f, MSE = %.3f\n", i, per_fold_accuracy[i], per_fold_mse[i]))
    }
    
    # Calculate overall performance
    overall_accuracy <- cor(predicted, y)
    overall_mse <- mean((predicted - y)^2)
    
    cat("\n===== Cross-Validation Summary for", phenotype_file, "=====\n")
    cat(sprintf("Overall prediction accuracy (r): %.3f\n", overall_accuracy))
    cat(sprintf("Overall MSE: %.3f\n", overall_mse))
    cat(sprintf("Mean fold accuracy: %.3f (SD = %.3f)\n", 
                mean(per_fold_accuracy), sd(per_fold_accuracy)))
    cat(sprintf("Mean fold MSE: %.3f (SD = %.3f)\n", 
                mean(per_fold_mse), sd(per_fold_mse)))
    
    
    # Assign fold numbers for each sample
    n <- length(y)
    fold_assignment <- rep(NA, n)
    
    for (fold_num in seq_along(folds)) {
      test_indices <- folds[[fold_num]]
      fold_assignment[test_indices] <- fold_num
    }
    
    # Build dataframe for CSV output
    pred_df <- data.frame(
      name = phenotype_data$accession_name,
      observed = y,
      predicted = predicted,
      fold = fold_assignment
    )
    
    r2 <- round(cor(pred_df$observed, pred_df$predicted, use = "complete.obs")^2, 3)
    
    ggplot(pred_df, aes(x = observed, y = predicted)) + 
      geom_point(alpha = 0.6, color = "steelblue") + 
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      annotate("text", x = min(pred_df$observed, na.rm = TRUE),
               y = max(pred_df$predicted, na.rm = TRUE),
               label = paste("r2 =", r2), hjust = 0, size = 5) +
      labs(title = paste(phenotype_file, "observed vs predicted"),
           x = "Observed", y = "Predicted") + theme_minimal()
    
    write.csv(pred_df,
              paste0("rrBLUP", gds_name, "_", phenotype_name, "_results.csv"),
              row.names = FALSE)
    
    
    
    

    
  
  }  
}  