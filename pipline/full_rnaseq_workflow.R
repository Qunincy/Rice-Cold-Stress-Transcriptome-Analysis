# ==============================================================================
# full_rnaseq_workflow.R
# Rice leaf time-series RNA-seq full analysis pipeline
# Treatment: CK, salt, PEG; Time gradient: 0h,30m,1h,3h,6h,10h,24h
#
# Workflow modules:
# 1. Raw count matrix cleaning & filtering
# 2. DESeq2 differential expression analysis
# 3. TPM normalization & log2(TPM FoldChange) calculation
# 4. Heatmap visualization of all DEGs
# 5. K-means clustering of DEG expression profiles (k=6)
# 6. Gene annotation & cluster-wise GO enrichment
# 7. Marker stress-response gene trend plotting
#
# Input: raw_data/leaf_countMatrix.csv (featureCounts output)
# Output: results/ folder (DEG tables, TPM matrix, clustering, GO tables, figures)
#
# Dependencies: DESeq2, clusterProfiler, org.Osativa.eg.db, pheatmap, factoextra
# ==============================================================================

# Load libraries
library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(cluster)
library(factoextra)
library(ggplot2)
library(clusterProfiler)
library(AnnotationDbi)
library(enrichplot)

# Load rice gene annotation database (compatible for HPC & local machine)
if(requireNamespace("org.Osativa.eg.db", quietly = TRUE)){
  library(org.Osativa.eg.db)
}else{
  # HPC custom library path
  library(org.Osativa.eg.db, lib = "/share/home/yuqun/MJ20260724165-R-260724054/go_enrich")
}

# ====================== 全局路径设置（仓库根目录为基准）======================
input_dir  <- "raw_data"
out_deseq  <- "results/deseq2"
out_cluster<- "results/cluster"
out_go     <- "results/go_enrich"
out_fig    <- "results/figures"

# Create output folders automatically
dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_deseq, showWarnings = FALSE, recursive = TRUE)
dir.create(out_cluster, showWarnings = FALSE, recursive = TRUE)
dir.create(out_go, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

# Check input file exists
input_file <- file.path(input_dir, "leaf_countMatrix.csv")
if(!file.exists(input_file)){
  stop(paste0("ERROR: Input count matrix NOT found at: ", input_file, "\nPlease put featureCounts result into raw_data/ folder"))
}

## -------- Step1: clean_data --------
cnt <- read_delim(input_file, skip = 1, delim = "\t")
cnt_length <- cnt[, c(1,6:ncol(cnt))]
old_names <- colnames(cnt_length)
new_names <- str_extract(old_names,"[^/]+$") %>%
  str_remove("^L1MLG\\d+-") %>%
  str_remove("\\.sorted\\.bam$")
colnames(cnt_length) <- new_names
cnt_length <- as.data.frame(cnt_length)
rownames(cnt_length) <- cnt_length$Geneid
cnt_length$Geneid <- NULL
count_cols <- colnames(cnt_length)[colnames(cnt_length) != "Length"]

# Filter low-expression genes: at least 3 samples with count >=10
keep <- rowSums(cnt_length[, count_cols] >= 10) >= 3
cnt_length <- cnt_length[keep,]
cat("After filtering, gene count = ", nrow(cnt_length), "\n")

write.csv(cnt_length,file = file.path(input_dir, "leaf_countMatrix_clean_length.csv"),row.names = TRUE)
cnt_clean <- cnt_length[,count_cols]
write.csv(cnt_clean,file = file.path(input_dir, "leaf_countMatrix_clean1.csv"),row.names = TRUE)

## -------- Step2: DESeq2 --------
sample_names <- colnames(cnt_clean)
split_mat <- str_split_fixed(sample_names, "_", n = 4)
sample_info <- tibble(
  sample = sample_names,
  tissue = split_mat[,1],
  time = split_mat[,2],
  treatment = split_mat[,3],
  rep = split_mat[,4]
) %>%
  mutate(
    rep = as.character(rep),
    condition = str_c(treatment, "_", time)
  )
stopifnot( sample_info$sample == colnames(cnt_clean) )

# Factor level setup
sample_info$time <- factor(sample_info$time,levels = c("0h","30m","1h","3h","6h","10h","24h"))
sample_info$condition <- factor(sample_info$condition)
sample_info$condition <- relevel(sample_info$condition, ref = "ck_0h")

dds <- DESeqDataSetFromMatrix(
  countData = cnt_clean,
  colData   = sample_info,
  design    = ~ condition
)
dds <- DESeq(dds)
cat("DESeq2 model finished. Available contrast names:\n")
print(resultsNames(dds))

## -------- Step3: Extract DEG --------
# All comparisons vs ck_0h
conditions <- levels(sample_info$condition)
conditions <- conditions[conditions != "ck_0h"]
cat("Comparisons to run: ", conditions, "\n")

deg_list <- list()
all_de_genes <- c()
for (cond in conditions) {
  res <- results(dds,contrast = c("condition",cond,"ck_0h"))
  res <- lfcShrink(dds, contrast = c("condition",cond,"ck_0h"), res = res, type = "ashr")
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  deg_list[[cond]] <- res_df
  # DEG threshold: padj <0.05, |log2FC|>1
  de_genes <- res_df %>% filter(padj <0.05,abs(log2FoldChange)>1) %>% pull(gene)
  all_de_genes <- union(all_de_genes,de_genes)
}
cat("Total unique DEGs across all comparisons: ", length(all_de_genes), "\n")

# Save DEG results
for (cond in names(deg_list)) {
  write.csv(deg_list[[cond]],file = file.path(out_deseq, paste0(cond, "_vs_ck_0h_DEG_all.csv")),row.names = FALSE)
}
write.table(data.frame(geneID = all_de_genes),file = file.path(out_deseq,"all_DEG_genes.txt"),sep = "\t", row.names = FALSE)

## -------- Step4: TPM calculation --------
gene_length <- cnt_length$Length
names(gene_length) <- rownames(cnt_length)
stopifnot(rownames(cnt_clean)==names(gene_length))

count2tpm <- function(counts, geneLength){
  rpk <- counts/(geneLength/1000)
  tpm <- t(t(rpk)/colSums(rpk, na.rm = TRUE)*1e6)
  return(tpm)
}
tpm_matrix <- count2tpm(cnt_clean,gene_length)

write.csv(tpm_matrix,file = file.path(out_deseq,"tpm_matrix.csv"),row.names = TRUE)
# Fix bug: variable name all_de_genes
deg_tpm <- tpm_matrix[all_de_genes, ]
write.csv(deg_tpm,file = file.path(out_deseq,"all_deg_gene_tpm.csv"),row.names = TRUE)

## -------- Step5: log2(TPM FoldChange) --------
ck0h_samples <- sample_info %>% filter(condition == "ck_0h") %>% pull(sample)
ck0h_tpm_mean <- rowMeans(tpm_matrix[,ck0h_samples])
write.csv(ck0h_tpm_mean, file = file.path(out_deseq,"ck0h_tpm_mean.csv"),row.names = TRUE)

# +1 to avoid log2(0)
log2_tpm_fc <- log2((tpm_matrix + 1) /(ck0h_tpm_mean + 1))
log2_tpm_fc_de <- log2_tpm_fc[all_de_genes, ]
write.csv(log2_tpm_fc,file = file.path(out_deseq,"log2_tpm_fc.csv"),row.names = TRUE)

# Reorder samples for heatmap
sample_order <- sample_info %>%
  mutate(treatment=factor(treatment, levels=c("ck","salt","PEG"))) %>%
  arrange(treatment, time, rep) %>%
  pull(sample)
log2_tpm_fc_de <- log2_tpm_fc_de[,sample_order]
write.csv(log2_tpm_fc_de,file = file.path(out_deseq,"log2_tpm_fc_de.csv"),row.names = TRUE)

## -------- Step6: Heatmap for all DEGs --------
anno_col <- sample_info %>% select(sample, time, treatment) %>% column_to_rownames("sample")
pheatmap(
  log2_tpm_fc_de,
  annotation_col = anno_col, 
  scale = "none",
  show_rownames = FALSE,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_row = 20,
  treeheight_col = 20,
  color = colorRampPalette(c("#2166ac","white","#b2182b"))(100),
  breaks = seq(-4,4,length.out=101),
  main = "log2(TPM FoldChange vs CK_0h) | All DEGs",
  filename = file.path(out_fig, "DEG_heatmap_log2TPM_FC.pdf")
)

## -------- Step7: K-means clustering of DEG expression pattern --------
log2_tpm <- log2(tpm_matrix + 1)
cluster_matrix <- log2_tpm[all_de_genes,]
# Z-score by row(gene): normalize expression trend
cluster_matrix_z <- t(scale(t(cluster_matrix)))
cat("NA number after Z-score: ", sum(is.na(cluster_matrix_z)), "\n")

# Estimate best k
p_wss <- fviz_nbclust(cluster_matrix_z, kmeans, method = "wss")
p_sil <- fviz_nbclust(cluster_matrix_z, kmeans, method = "silhouette")
ggsave(file.path(out_fig,"kmeans_wss_plot.pdf"), p_wss, width = 8, height =5)
ggsave(file.path(out_fig,"kmeans_silhouette_plot.pdf"), p_sil, width =8, height=5)

set.seed(123)
km <- kmeans(
  cluster_matrix_z,
  centers = 6,
  nstart = 50
)
cluster_gene_6 <- data.frame(
  geneID = names(km$cluster),
  cluster = km$cluster
)
write.csv(
  cluster_gene_6,
  file.path(out_cluster,"cluster_gene_6.csv"),
  row.names = FALSE
)
cat("Gene count per cluster:\n")
print(table(cluster_gene_6$cluster))

cluster_list_6 <- split(
  cluster_gene_6$geneID,
  cluster_gene_6$cluster
)

# Cluster trend plot
center_df <- as.data.frame(km$centers)
center_df <- center_df[,sample_order]
center_df$cluster <- paste0("Cluster ", 1:nrow(center_df))
center_long <- center_df %>%
  pivot_longer(cols=-cluster, names_to="sample", values_to="expression") %>%
  left_join(sample_info %>% dplyr::select(sample,time,treatment), by="sample")
center_summary <- center_long %>%
  group_by(cluster, time, treatment) %>%
  summarise(expression = mean(expression), .groups="drop")

p_cluster_trend <- ggplot(center_summary, aes(x=time, y=expression, group=treatment, color=treatment))+
  geom_line(linewidth=1)+
  facet_wrap(~cluster, ncol=3)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90))
ggsave(file.path(out_fig,"cluster_expression_trend.pdf"), p_cluster_trend, width =12, height=7)

## -------- Step8: Gene annotation & GO mapping statistics --------
gene_anno <- bitr(cluster_gene_6$geneID,
                  fromType = "GID",
                  toType = c("Gene_Name","GO"),
                  OrgDb = org.Osativa.eg.db)

# ==== Fix original bugs: gene_anno_6 -> gene_anno; cluster_gene -> cluster_gene_6 ====
cat("GO mapping stats:\n")
print(table(cluster_gene_6$geneID %in% gene_anno$GID))
cluster_anno <- left_join(cluster_gene_6, gene_anno, by = c("geneID"="GID"))

gene_anno <- bitr(
  cluster_gene_6$geneID,
  fromType = "GID",
  toType = "GO",
  OrgDb = org.Osativa.eg.db
)
mapped_go_genes <- unique(gene_anno$GID)
cat("Mapped GO gene count: ", length(mapped_go_genes), "\n")
cat("Mapping rate: ", length(mapped_go_genes) / nrow(cluster_gene_6), "\n")

gene_anno_goall <- bitr(
  cluster_gene_6$geneID,
  fromType = "GID",
  toType = "GOALL",
  OrgDb = org.Osativa.eg.db
)
cat("GOALL mapping rate: ", length(unique(gene_anno_goall$GID)) / nrow(cluster_gene_6), "\n")

cluster_map_stat <- cluster_gene_6 %>%
  mutate(GO_mapped = geneID %in% mapped_go_genes) %>%
  group_by(cluster) %>%
  summarise(
    total = n(),
    mapped = sum(GO_mapped),
    unmapped = total - mapped,
    mapping_rate = mapped / total
  )
print(cluster_map_stat)

go_background <- unique(gene_anno$GID)
cluster_list_go <- lapply(
  cluster_list_6,
  function(x) intersect(x, go_background)
)
cat("Genes available for GO enrichment per cluster:\n")
print(sapply(cluster_list_go, length))

## -------- Step9: GO enrichment for each cluster --------
go_res_list <- list()
go_res_simple_list <- list()
for (cid in names(cluster_list_go)) {
  gene_set <- cluster_list_go[[cid]]
  if (length(gene_set) < 10) {
    message(sprintf("Cluster %s: gene number <10, skip enrichment", cid))
    next
  }

  ego <- enrichGO(
    gene          = gene_set,
    OrgDb         = org.Osativa.eg.db,
    keyType       = "GID",
    ont           = "ALL",
    pAdjustMethod = "fdr",
    universe      = go_background,
    qvalueCutoff  = 0.05
  )
  go_res_list[[cid]] <- ego

  ego_simple <- simplify(ego, cutoff = 0.7)
  go_res_simple_list[[cid]] <- ego_simple

  write.csv(ego@result,
            file = file.path(out_go, paste0("GO_Cluster_", cid, ".csv")),
            row.names = FALSE)
  write.csv(ego_simple@result,
            file = file.path(out_go, paste0("GO_simplified_Cluster_", cid, ".csv")),
            row.names = FALSE)
}

# Plot GO dotplot
for (cid in names(go_res_simple_list)) {
  p <- dotplot(go_res_simple_list[[cid]],
               showCategory = 10,
               split = "ONTOLOGY",
               color = "p.adjust",
               size = "Count",
               label_format = 40,
               title = paste0("Cluster ",cid," GO enrichment")) +
    facet_grid(ONTOLOGY ~ ., scales = "free_y") +
    scale_size_continuous(range = c(2, 8)) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9),
      axis.title = element_text(size = 11),
      strip.text = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10),
      strip.text.y = element_text(angle = 0)
    )
  ggsave(file.path(out_go, paste0("GO_dotplot_Cluster", cid, ".png")),
         plot = p, width = 10, height =9, dpi=300)
}

## -------- Step10: Marker stress gene expression trend --------
gene_name_map <- c(
  "LOC_Os12g44360" = "OsSOS1",
  "LOC_Os06g40370" = "OsSOS2/CIPK24",
  "LOC_Os05g45810" = "OsSOS3/CBL4",
  "LOC_Os01g20160" = "OsHKT1;5",
  "LOC_Os04g51820" = "OsHKT1;1",
  "LOC_Os01g07120" = "OsDREB2A",
  "LOC_Os03g44380" = "OsNCED3",
  "LOC_Os01g64000" = "OsABI5",
  "LOC_Os09g35030" = "OsDREB1A / OsCBF3"
)
marker_df <- log2_tpm %>%
  as.data.frame() %>%
  rownames_to_column("geneID") %>%
  filter(geneID %in% gene_name_map) %>%
  pivot_longer(-geneID, names_to = "sample", values_to = "log2TPM") %>%
  left_join(sample_info, by = "sample") %>%
  mutate(time = factor(time, levels=c("0h","30m","1h","3h","6h","10h","24h")))

marker_summary <- marker_df %>%
  group_by(geneID, time, treatment) %>%
  summarise(
    mean_log2TPM = mean(log2TPM),
    sd_log2TPM = sd(log2TPM),
    .groups = "drop"
  ) %>%
  mutate(time = factor(time, levels=c("0h","30m","1h","3h","6h","10h","24h")))

marker_df <- marker_df %>%
  mutate(geneName = gene_name_map[geneID])
marker_summary <- marker_summary %>%
  mutate(geneName = gene_name_map[geneID])

p_marker <- ggplot() +
  geom_point(data = marker_df, aes(x=time, y=log2TPM, color=treatment), alpha=0.3, size=1.5) +
  geom_line(data = marker_summary, aes(x=time, y=mean_log2TPM, color=treatment, group=treatment), linewidth=1) +
  geom_errorbar(data = marker_summary, aes(x=time, ymin=mean_log2TPM-sd_log2TPM, ymax=mean_log2TPM+sd_log2TPM, color=treatment), width=0.2) +
  facet_wrap(~geneID, scales="free_y", ncol=3,
             labeller = labeller(geneID = function(x) paste0(gene_name_map[x], "\n(", x, ")"))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        strip.text = element_text(size=10, face="bold")) +
  labs(x="Time", y="log2(TPM+1)", color="Treatment")

ggsave(file.path(out_fig,"marker_gene_mean_trend.png"), p_marker, width=14, height=9, dpi=300)
ggsave(file.path(out_fig,"marker_gene_mean_trend.pdf"), p_marker, width=14, height=9)

cat("===== Analysis finished! All outputs saved under results/ =====\n")
