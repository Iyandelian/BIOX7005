library(SNPRelate)
library(gdsfmt)
library(SeqArray)
library(ggplot2)

gds_file <- "filtered_snps.gds"
gds_file2 <- 'filtered_file.gds'

# Open the GDS files
genofile <- snpgdsOpen(gds_file)
genofile2 <- seqOpen(gds_file2)

# accession info
reg1001_data <- read.csv("/home/blake/Documents/BIOX7005/new variants/variants/Arabidopsis_2029_accessions_imputation.csv")
sampleIDlist <- read.gdsn(index.gdsn(genofile, "sample.id"))
countrylist <- list()

# --- Split samples into training (70%) and test (30%) sets ---
set.seed(123)
all_samples <- read.gdsn(index.gdsn(genofile, "sample.id"))
n_train <- floor(0.7 * length(all_samples))
train_samples <- sample(all_samples, n_train)
test_samples  <- setdiff(all_samples, train_samples)

for (x in sampleIDlist) {
  for(i in 1:nrow(reg1001_data)) {
    row <- reg1001_data[i,]
    if (x == row$FID) {
      countrylist <- c(countrylist, row$country)
    }
  }
}

print(length(countrylist))

# --- Run PCA on training samples only ---
pca <- snpgdsPCA(genofile, sample.id = train_samples, num.thread = 4)

# print eigenvector values
pc.percent <- pca$varprop * 100
pc.percent.sig <- pc.percent[pc.percent >= 1]
num_pca <- length(pc.percent.sig)

loadings <- snpgdsPCASNPLoading(pca, genofile)
snploading <- loadings$snploading
snp_ids <- loadings$snp.id

num_pca <- min(32, nrow(snploading))
combined_loading <- colSums(snploading[1:num_pca, ]^2)
length(combined_loading)

top_idx <- order(combined_loading, decreasing = TRUE)[1:2500]
top_snps <- snp_ids[top_idx]

# eigenvectors for table
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],
                  EV2 = pca$eigenvect[,2],
                  EV3 = pca$eigenvect[,3],
                  stringsAsFactors = FALSE)
head(tab)

# Draw PCA
plot(tab$EV3, tab$EV2, xlab = "eigenvector 3", ylab = "eigenvector 2")

# PCA with country codes
pop_code <- unlist(countrylist)
head(cbind(sampleIDlist, pop_code))
length(unique(pop_code))
unique(pop_code)

tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(pop_code)[match(pca$sample.id, sampleIDlist)],
                  EV1 = pca$eigenvect[,1],
                  EV2 = pca$eigenvect[,2],
                  EV3 = pca$eigenvect[,3],
                  stringsAsFactors = FALSE)
head(tab)

ggplot(data = tab, aes(x = EV3, y = EV2, label = pop, colour = pop, shape = pop)) +
  geom_point() +
  geom_text() +
  labs(title = "PCA plot by country", x = "eigenvector 2", y = "eigenvector 1") +
  theme_minimal()

all_snps <- read.gdsn(index.gdsn(genofile, "snp.id"))
sum(top_snps %in% all_snps)

# --- Create filtered test GDS using top SNPs from training PCA ---
snpgdsCreateGenoSet(
  gds_file,
  "filtered_test_PCA_snps.gds",
  sample.id = test_samples,
  snp.id = top_snps,
  snpfirstdim = TRUE,
  compress.annotation = "ZIP_RA.max",
  compress.geno = "ZIP_RA",
  verbose = TRUE
)

# Sanity checks
cat("Training samples:", length(train_samples), "\n")
cat("Test samples:", length(test_samples), "\n")
cat("Top SNPs selected from training:", length(top_snps), "\n")

# Close GDS files at end
snpgdsClose(genofile)
seqClose(genofile2)

