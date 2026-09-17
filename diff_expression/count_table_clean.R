# Load required packages
library(readr)
library(stringr)

# Read featureCounts output
cnt <- read_delim("data\all_countMatrix.csv", skip = 1, delim = "\t")

# Keep gene ID and sample count columns only
cnt_clean <- cnt[, c(1,7:ncol(cnt))]

# Simplify sample column names: extract filename and remove bam suffix
old_names <- colnames(cnt_clean)
new_names <- str_extract(old_names, "[^/]+$") %>% str_remove("\\.bam")
colnames(cnt_clean) <- new_names

# Set gene ID as row names, remove Geneid column
cnt_clean <- as.data.frame(cnt_clean)
rownames(cnt_clean) <- cnt_clean$Geneid
cnt_clean$Geneid <- NULL

# Remove genes with zero reads across all samples
cnt_clean <- cnt_clean[rowSums(cnt_clean) > 0, ]

# Export cleaned count matrix
write.csv(cnt_clean, file = "all_countMatrix_clean.csv", row.names = TRUE)
