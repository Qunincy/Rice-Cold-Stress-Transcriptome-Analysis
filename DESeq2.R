# Load required packages
library(readr)
library(stringr)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(patchwork)
library(biomaRt)
library(httr2)
library(AnnotationDbi)

# ========== User configurable parameters ==========
base_dir <- getwd()
out_dir_volcano <- file.path(base_dir, "output/volcano")
out_dir_heatmap <- file.path(base_dir, "output/heatmap")
out_dir_deg <- file.path(base_dir, "output/deg_tables")
count_file <- file.path(base_dir, "input/all_countMatrix.csv")
id_map_file <- file.path(base_dir, "input/RAP-MSU_2026-02-05.txt.gz")

# Auto create output folders
if (!dir.exists(out_dir_volcano)) dir.create(out_dir_volcano, recursive = TRUE)
if (!dir.exists(out_dir_heatmap)) dir.create(out_dir_heatmap, recursive = TRUE)
if (!dir.exists(out_dir_deg)) dir.create(out_dir_deg, recursive = TRUE)

# ---------------------- 2. 读入计数矩阵 + 样本名清洗 ----------------------
cnt <- read_delim(count_file, skip = 1, delim = "\t")
cnt_clean <- cnt[,c(1,7:ncol(cnt))]
old_names <- colnames(cnt_clean)
new_names <- str_extract(old_names, "[^/]+$") %>% str_remove("\\.bam")
colnames(cnt_clean) <- new_names
cnt_clean <- as.data.frame(cnt_clean)
rownames(cnt_clean) <- cnt_clean$Geneid
cnt_clean$Geneid <- NULL
cnt_clean <- cnt_clean[rowSums(cnt_clean) > 0, ]

# ---------------------- 3. 拆分样本分组信息 ----------------------
sample_names <- colnames(cnt_clean)
group_names  <- str_extract(sample_names, "^[^_]+")                    # Genotype
time_names   <- str_extract(sample_names, "(?<=_)[0-9]+(?=_[0-9]$)")    # Time point: 0/24/72
rep_names    <- str_extract(sample_names, "\\d$")                       # Biological replicate ID
sample_info <- data.frame(
  sample = sample_names,
  genotype = group_names,
  time = time_names,
  rep = rep_names
)
rownames(sample_info) <- sample_info$sample
sample_info$sample <- NULL
# 转为因子适配DESeq2
sample_info$genotype <- factor(sample_info$genotype)
sample_info$time     <- factor(sample_info$time, levels = c("0","24","72"))
sample_info$rep      <- factor(sample_info$rep)

# ---------------------- 4. DESeq2 differential expression with interaction model ----------------------
dds <- DESeqDataSetFromMatrix(
  countData = cnt_clean,
  colData   = sample_info,
  design    = ~ genotype + time + genotype:time
)
dds <- DESeq(dds)
resultsNames(dds)

# Extract Y vs C contrast
res <- results(dds, contrast = c("genotype","Y","C"), alpha = 0.05, lfcThreshold = 1)
res_sig <- subset(res, padj < 0.05 & abs(log2FoldChange) >= 1)

# Save core DESeq2 objects
save(dds, res, res_sig, file = file.path(out_dir_deg, "Y_vs_C_deseq2_result.RData"))

# ---------------------- 5. Format results, filter NA, label up/down genes ----------------------
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_filter <- res_df[!is.na(res_df$padj), ]

# Label significance (separate from res object to avoid sync issues)
res_filter$sig <- "Not significant"
res_filter$sig[res_filter$padj < 0.05 & res_filter$log2FoldChange >= 1]  <- "Up"
res_filter$sig[res_filter$padj < 0.05 & res_filter$log2FoldChange <= -1] <- "Down"

# Export up/down gene ID lists and full results table
up_genes   <- res_filter$gene[res_filter$sig == "Up"]
down_genes <- res_filter$gene[res_filter$sig == "Down"]

write.table(up_genes,   file.path(out_dir_deg, "Y_vs_C_up_geneID.txt"),   quote = F, row.names = F, col.names = F)
write.table(down_genes, file.path(out_dir_deg, "Y_vs_C_down_geneID.txt"), quote = F, row.names = F, col.names = F)
write.csv(res_df, file.path(out_dir_deg, "Y_vs_C_all_res.csv"), row.names = T)

