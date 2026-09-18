# ==============================================================================
# Rice RNA-seq: Batch GO enrichment
# Auto ID mapping: MSU-LOC → RAP → ENTREZID
# Batch enrichment, GO simplify, export tables and dotplot figures
# ==============================================================================
library(DESeq2)
library(dplyr)
library(tidyverse)
library(clusterProfiler)
library(openxlsx)

#' GO enrichment wrapper
#' @param up_gene_vec vector of up-regulated MSU LOC IDs
#' @param down_gene_vec vector of down-regulated MSU LOC IDs
#' @param rap_msu MSU-LOC_gene to RAP ID table
#' @param rap_ncbi RAP to ENTREZID table
#' @param rice_orgdb rice OrgDb annotation database
#' @return list(ego_up,ego_down)
run_general_go <- function(up_gene_vec,
                           down_gene_vec,
                           rap_msu,
                           rap_ncbi,
                           rice_orgdb){
  # Step1: MSU LOC to RAP
  up_genes_df <- data.frame(gene = up_gene_vec)
  down_genes_df <- data.frame(gene = down_gene_vec)
  
  up_LOC_MSU <- up_genes_df %>%
    left_join(rap_msu, by = c("gene" = "LOC_gene")) %>%
    filter(!is.na(RAP))
  down_LOC_MSU <- down_genes_df %>%
    left_join(rap_msu, by = c("gene" = "LOC_gene")) %>%
    filter(!is.na(RAP))
  
  # Step2: RAP to ENTREZID
  up_RAP_ENTREZID <- up_LOC_MSU %>%
    left_join(rap_ncbi, by = "RAP", relationship = "many-to-many") %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(gene, .keep_all = TRUE)
  down_RAP_ENTREZID <- down_LOC_MSU %>%
    left_join(rap_ncbi, by = "RAP", relationship = "many-to-many") %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(gene, .keep_all = TRUE)
  
  # GO enrichment for up-regulated genes
  entrez_up <- unique(as.character(up_RAP_ENTREZID$ENTREZID))
  ego_up <- NULL
  if(length(entrez_up) > 0){
    ego_up <- enrichGO(
      gene          = entrez_up,
      OrgDb         = rice_orgdb,
      keyType       = "ENTREZID",
      ont           = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05,
      readable      = TRUE
    )
  }
  
  # GO enrichment for down-regulated genes
  entrez_down <- unique(as.character(down_RAP_ENTREZID$ENTREZID))
  ego_down <- NULL
  if(length(entrez_down) > 0){
    ego_down <- enrichGO(
      gene          = entrez_down,
      OrgDb         = rice_orgdb,
      keyType       = "ENTREZID",
      ont           = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05,
      readable      = TRUE
    )
  }
  return(list(ego_up = ego_up, ego_down = ego_down))
}

#' Extract DE gene vectors from DESeq2 contrast
get_deg_vec <- function(dds,contrast_vec){
  res <- results(dds,contrast = contrast_vec)
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_filter <- res_df[!is.na(res_df$padj), ]
  up_vec <- rownames(res_filter)[res_filter$padj <0.05 & res_filter$log2FoldChange >=1]
  down_vec <- rownames(res_filter)[res_filter$padj <0.05 & res_filter$log2FoldChange <=-1]
  return(list(up = up_vec, down = down_vec))
}

# List of all contrast groups
contrast_list <- list(
  K_24vs0 = list(c("time_24_vs_0","genotypeK.time24")),
  K_72vs0 = list(c("time_72_vs_0","genotypeK.time72")),
  K_24vs72 = list(c("time_24_vs_72","genotypeK.time24")),
  K0_vs_C0 = list(c("genotypeK.0h","genotypeC.0h")),
  K24_vs_C24 = list(c("genotypeK.24h","genotypeC.24h"))
  # Add remaining contrasts here
)

all_go_allcomparison <- list()
all_go_simplify <- list()

# Run GO enrichment for all contrasts
for(nm in names(contrast_list)){
  cat("\n===== Running: ",nm," =====\n")
  ct <- contrast_list[[nm]]
  deg_tmp <- get_deg_vec(dds,ct)
  go_tmp <- run_general_go(
    up_gene_vec = deg_tmp$up,
    down_gene_vec = deg_tmp$down,
    rap_msu = rap_msu,
    rap_ncbi = rap_ncbi,
    rice_orgdb = rice_orgdb
  )
  all_go_allcomparison[[nm]] <- go_tmp
}

# Test plot for K_72vs0
ego <- all_go_allcomparison[["K_72vs0"]]$ego_up
dotplot(ego, showCategory =15)
p_up <- barplot(all_go_allcomparison[["K_72vs0"]]$ego_up, showCategory =15, title = paste0("K_72vs0","_Up_GO"))
ego_bp <- all_go_allcomparison[["K_72vs0"]]$ego_up %>% filter(ONTOLOGY=="BP")
p_up <- dotplot(ego_bp, showCategory =15, title = paste0("K_72vs0","_Up_GO_BP")) 
show(p_up)

# GO simplify for all groups
for(nm in names(all_go_allcomparison)){
  cat("\n==== Simplify GO: ", nm, " ====\n")
  go_res <- all_go_allcomparison[[nm]]
  
  if(!is.null(go_res$ego_up)){
    go_up_simple <- simplify(go_res$ego_up, cutoff = 0.7)
  }else{
    go_up_simple <- NULL
  }
  
  if(!is.null(go_res$ego_down)){
    go_down_simple <- simplify(go_res$ego_down, cutoff = 0.7)
  }else{
    go_down_simple <- NULL
  }
  all_go_simplify[[nm]] <- list(ego_up = go_up_simple, ego_down = go_down_simple)
}

# Export result tables
out_dir <- "./go_all_comparison_result"
if(!dir.exists(out_dir)) dir.create(out_dir,recursive = TRUE)

# Export csv
for(nm in names(all_go_allcomparison)){
  item <- all_go_allcomparison[[nm]]
  if(!is.null(item$ego_up)){
    write_csv(as.data.frame(item$ego_up),file.path(out_dir,paste0(nm,"_up_GO.csv")))
  }
  if(!is.null(item$ego_down)){
    write_csv(as.data.frame(item$ego_down),file.path(out_dir,paste0(nm,"_down_GO.csv")))
  }
}

# Export xlsx
for (i in names(all_go_allcomparison)) {
  if(!is.null(all_go_allcomparison[[i]]$ego_up)){
    write.xlsx(as.data.frame(all_go_allcomparison[[i]]$ego_up), file.path(out_dir, sprintf("%s_up.xlsx", i)))
  }
  if(!is.null(all_go_allcomparison[[i]]$ego_down)){
    write.xlsx(as.data.frame(all_go_allcomparison[[i]]$ego_down), file.path(out_dir, sprintf("%s_down.xlsx", i)))
  }
}

# Batch draw dotplots for simplified GO
for(nm in names(all_go_simplify)){
  res_item <- all_go_simplify[[nm]]
  
  if(!is.null(res_item$ego_up)){
    p_up <- dotplot(res_item$ego_up, showCategory =15, title = paste0(nm,"_Up_GO"))
    ggsave(filename = file.path(out_dir,paste0(nm,"_up_dotplot.pdf")),
           plot = p_up, width =10, height =7)
  }else{
    cat(nm," No significant GO terms for up-regulated genes, skip plot\n")
  }
  
  if(!is.null(res_item$ego_down)){
    p_down <- dotplot(res_item$ego_down, showCategory =15, title = paste0(nm,"_Down_GO"))
    ggsave(filename = file.path(out_dir,paste0(nm,"_down_dotplot.pdf")),
           plot = p_down, width =10, height =7)
  }else{
    cat(nm," No significant GO terms for down-regulated genes, skip plot\n")
  }
}

# Facet dotplot, split by GO ontology(BP/CC/MF)
for(nm in names(all_go_simplify)){
  res_item <- all_go_simplify[[nm]]
  
  if(!is.null(res_item$ego_up)){
    p_up <- dotplot(res_item$ego_up, 
                    showCategory = 10, 
                    split = "ONTOLOGY", 
                    color = "p.adjust", 
                    size = "Count") +
      facet_grid(ONTOLOGY~., scales = "free") +
      ggtitle(paste0(nm, "_Up_GO"))
    
    ggsave(filename = file.path(out_dir, paste0(nm, "_up_go_facet.pdf")),
           plot = p_up, width = 10, height = 9)
  }else{
    cat(nm, " No significant GO terms for up-regulated genes, skip facet plot\n")
  }
  
  if(!is.null(res_item$ego_down)){
    p_down <- dotplot(res_item$ego_down, 
                      showCategory = 10, 
                      split = "ONTOLOGY", 
                      color = "p.adjust", 
                      size = "Count") +
      facet_grid(ONTOLOGY~., scales = "free") +
      ggtitle(paste0(nm, "_Down_GO"))
    
    ggsave(filename = file.path(out_dir, paste0(nm, "_down_go_facet.pdf")),
           plot = p_down, width = 10, height = 9)
  }else{
    cat(nm, " No significant GO terms for down-regulated genes, skip facet plot\n")
  }
}
