library(ranger)
library(SNPRelate)
library(gdsfmt)
library(SeqArray)
library(ggplot2)

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
    
    
    
    # Set seed for reproducibility
    set.seed(400)
    
    # --- Alignment check and fix ---
    common <- intersect(rownames(GMdata), phenotype_data$accession_id)
    GM_pheno_data <- GMdata[common, , drop = FALSE]
    phenotype_data <- phenotype_data[match(common, phenotype_data$accession_id), ]
    if(!all(rownames(GM_pheno_data) == phenotype_data$accession_id)) stop("Alignment failed")
    x <- GM_pheno_data
    y <- phenotype_data$phenotype_value
    # --- End alignment block ---
    
    # Handle missing values in the genotype matrix (imputation)
    for (j in 1:ncol(x)) {
      na_idx <- which(is.na(x[, j]))
      if (length(na_idx) > 0) {
        x[na_idx, j] <- round(mean(x[, j], na.rm = TRUE), digits = 0)
      }
    }
    
    # ---- 2. Settings ----
    set.seed(42)
    k <- 10                            # CV folds
    ntree <- 100                     # Trees in forest
    n_cores <- 5
    top_n <- 10                      # Top N features to plot
    
    # ---- 3. Scale Genotype Matrix ----
    X_scaled <- x

    # ---- 4. Cross-validation ----
    n <- nrow(X_scaled)
    folds <- sample(rep(1:k, length.out = n))
    cv_pred <- rep(NA, n)
    importance_list <- list()
    
    
    for (i in 1:k) {
      cat("Fold", i, "...\n")
      
      test_idx <- which(folds == i)
      train_idx <- setdiff(1:n, test_idx)
      
      X_train <- X_scaled[train_idx, , drop = FALSE]
      y_train <- y[train_idx]
      X_test  <- X_scaled[test_idx, , drop = FALSE]
      
      rf_cv_model <- ranger(
        x = X_train,
        y = y_train,
        num.trees = ntree,
        importance = "none",
        num.threads = n_cores,
        min.node.size = 50,
        max.depth = 15
      )
      
      # Save importance
      importance_list[[i]] <- rf_cv_model$variable.importance
      
      # Predict
      cv_pred[test_idx] <- predict(rf_cv_model, data = X_test)$predictions
    }
    
    # ---- 5. Results ----
    overall_accuracy <- cor(cv_pred, y)
    overall_mse <- mean((cv_pred - y)^2)
    
    cat("\n===== Cross-Validation Summary =====\n")
    cat(sprintf("Overall prediction accuracy (r): %.3f\n", overall_accuracy))
    cat(sprintf("Overall MSE: %.3f\n", overall_mse))
    
    n <- length(y)
    fold_assignment <- rep(NA, n)
    predicted_full <- rep(NA, n)
    
    counter_fold <- 1
    
    pred_df <- data.frame(
      name = phenotype_data$accession_name[common_accessions],  # select phenotype
      observed = y,
      predicted = cv_pred,
      fold = folds
    )
    
    r2 <- round(cor(pred_df$observed, pred_df$predicted)^2, 3)
    
    ggplot(pred_df, aes(x = observed, y = predicted)) + geom_point(alpha = 0.6, color = "steelblue") + 
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      annotate("text", x = min(pred_df$observed), y = max(pred_df$predicted), label = paste("r2 =", r2), hjust = 0, size = 5) +
      labs(title = "observed vs predicted",
           x = "observer", y = "predicted") + theme_minimal()
    
    # Optionally: save results
    write.csv(pred_df, paste0("Random_Forest", gds_name, "_", phenotype_name, "_results.csv"), row.names = FALSE)
  }
  
  # Clear memory after processing the GDS file
  gc()
  

    
  } 
  
  


