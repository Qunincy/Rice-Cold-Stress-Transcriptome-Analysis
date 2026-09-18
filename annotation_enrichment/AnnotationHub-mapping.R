library(AnnotationHub)

# Initialize AnnotationHub
hub <- AnnotationHub()
# Query rice OrgDb resources
rice_orgdb_query <- query(hub, c("OrgDb", "Oryza"))
rice_orgdb_query

# Retrieve japonica rice OrgDb
rice_orgdb <- hub[['AH117673']]

# Check available ID types
keytypes(rice_orgdb)

# Preview gene keys (MSU LOC IDs)
head(keys(rice_orgdb, keytype = "GID"), 20)

# Test if MSU7 LOC ID exists in alias mapping
"LOC_Os01g01010" %in% keys(rice_orgdb, keytype = "ALIAS")

# ==============================================================================
# Rice RNA-seq: MSU7 gene ID annotation & multi-source ID mapping
# Resolve the lack of native MSU7 OrgDb package for rice
# Two strategies for RAP -> ENTREZ ID conversion:
# 1. biomaRt package (preferred, simple)
# 2. Biomart REST API via httr2 (fallback, when biomaRt fails / connection blocked)
# ==============================================================================
library(tidyverse)
library(AnnotationHub)
library(httr2)
library(AnnotationDbi)

# ---------------------- Load DE results ----------------------
load("Y_vs_C_deseq2_result.RData")

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_filter <- res_df[!is.na(res_df$padj), ]

# Define up/down regulated genes
up_genes <- rownames(res_filter)[res_filter$padj < 0.05 & res_filter$log2FoldChange >= 1]
down_genes <- rownames(res_filter)[res_filter$padj < 0.05 & res_filter$log2FoldChange <= -1]

up_genes_df <- data.frame(gene = up_genes)
down_genes_df <- data.frame(gene = down_genes)

res_up <- res_filter[res_filter$padj <0.05 & res_filter$log2FoldChange >=1, ]
res_down <- res_filter[res_filter$padj <0.05 & res_filter$log2FoldChange <= -1, ]

# ---------------------- Strategy1: Local MSU-RAP mapping table ----------------------
rap_msu <- read_tsv("RAP-MSU_2026-02-05.txt.gz",
                    col_names = c("RAP", "LOC"),
                    col_types = cols(RAP = col_character(), LOC = col_character()))

rap_msu <- rap_msu %>%
  separate_rows(LOC, sep = ",") %>%
  mutate(LOC = trimws(LOC),
         LOC = ifelse(LOC == "None", NA, LOC),
         RAP = ifelse(RAP == "None", NA, RAP)) %>%
  filter(!is.na(LOC), !is.na(RAP)) %>%
  mutate(LOC_gene = sub("\\..*$", "", LOC)) %>%
  distinct(LOC_gene, .keep_all = TRUE) %>%
  dplyr::select(RAP, LOC_gene)

write.csv(rap_msu, "RAP-MSU_clean.csv", row.names = F)

# Map DEG MSU LOC ID to RAP ID
up_LOC_MSU <- up_genes_df %>% left_join(rap_msu, by = c("gene" = "LOC_gene")) %>% filter(!is.na(RAP))
down_LOC_MSU <- down_genes_df %>% left_join(rap_msu, by = c("gene" = "LOC_gene")) %>% filter(!is.na(RAP))

# ---------------------- Strategy2: RAP to ENTREZ ID, two interchangeable methods ----------------------
gene_ids <- unique(na.omit(up_LOC_MSU$RAP))
rap_ncbi <- NULL

## Method A: biomaRt (primary, recommended)
tryCatch({
  library(biomaRt)
  mart <- useMart(biomart = "plants_mart", host = "https://plants.ensembl.org", dataset = "osativa_eg_gene")
  rap_ncbi <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                    filters = "ensembl_gene_id",
                    values = gene_ids,
                    mart = mart)
  colnames(rap_ncbi) <- c("RAP", "ENTREZID")
  rap_ncbi <- rap_ncbi %>% filter(!is.na(ENTREZID))
  message("Success: ID mapping using biomaRt")
}, error = function(e){
  message("biomaRt failed, switch to Biomart REST API (httr2 fallback)")
  
  ## Method B: httr2 REST API (backup, when biomaRt connection error)
  get_entrez <- function(ids) {
    query_xml <- paste0(
      '<Query virtualSchemaName="plants_mart" formatter="TSV" header="1">',
      '<Dataset name="osativa_eg_gene" interface="default">',
      '<Filter name="ensembl_gene_id" value="',
      paste(ids, collapse = ","),
      '"/>',
      '<Attribute name="ensembl_gene_id"/>',
      '<Attribute name="entrezgene_id"/>',
      '</Dataset>',
      '</Query>'
    )
    resp <- request("https://plants.ensembl.org/biomart/martservice") |>
      req_url_query(query = query_xml) |>
      req_retry(max_tries = 3) |>
      req_perform()
    read.delim(text = resp_body_string(resp), sep = "\t", check.names = FALSE)
  }
  
  groups <- split(gene_ids, ceiling(seq_along(gene_ids)/50))
  result_list <- lapply(groups, get_entrez)
  rap_ncbi <<- bind_rows(result_list)
  colnames(rap_ncbi) <- c("RAP", "ENTREZID")
  rap_ncbi <- rap_ncbi %>% filter(!is.na(ENTREZID))
})

# Merge ENTREZ ID back to DEG table
up_RAP_ENTREZID <- up_LOC_MSU %>%
  left_join(rap_ncbi, by = "RAP", relationship = "many-to-many") %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(gene, .keep_all = TRUE)

down_RAP_ENTREZID <- down_LOC_MSU %>%
  left_join(rap_ncbi, by = "RAP", relationship = "many-to-many") %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(gene, .keep_all = TRUE)

# Save workspace
save(dds,res,res_sig,res_df,res_filter,res_down,res_up,
     up_genes_df,down_genes,up_genes,down_genes_df,
     up_LOC_MSU,down_LOC_MSU,rap_msu,rap_ncbi,
     up_RAP_ENTREZID,down_RAP_ENTREZID,
     file = "up_and_down.RData")

write.csv(rap_ncbi, "RAP-ENTREZID_clean.csv", row.names = F)
