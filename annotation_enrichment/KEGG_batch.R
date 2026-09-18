library(clusterProfiler)
library(tidyverse)
library(openxlsx)

# Download rice KEGG gson object (run once, save locally for offline usage)
download_KEGG("osa")
osa_gson <- gson_KEGG(
  species = "osa",
  KEGG_Type = "KEGG",
  keyType = "ncbi-geneid"
)

# Inspect gson object
exists("osa_gson")
class(osa_gson)
osa_gson@species
osa_gson@keytype
osa_gson@version
osa_gson@accessed_date

# Save offline KEGG database
saveRDS(
  osa_gson,
  file = "osa_KEGG_ncbi_geneid_gson.rds"
)

#' General KEGG enrichment function
#' Input MSU‑LOC ID vector, perform two-step ID mapping and KEGG enrichment
#' @param up_gene_vec Character vector of up-regulated MSU‑LOC IDs
#' @param down_gene_vec Character vector of down-regulated MSU‑LOC IDs
#' @param rap_msu MSU‑LOC_gene → RAP mapping table
#' @param rap_ncbi RAP → ENTREZID mapping table
#' @param osa_gson KEGG gson object for rice (osa)
#' @return list(kegg_up, kegg_down), NULL if no valid genes
run_general_kegg <- function(up_gene_vec,
                             down_gene_vec,
                             rap_msu,
                             rap_ncbi,
                             osa_gson){
  # Step1: build gene dataframe
  up_genes_df <- data.frame(gene = up_gene_vec)
  down_genes_df <- data.frame(gene = down_gene_vec)
  
  # Mapping 1: MSU‑LOC_gene → RAP ID
  up_LOC_MSU <- up_genes_df %>%
    left_join(rap_msu, by = c("gene" = "LOC_gene")) %>%
    filter(!is.na(RAP))
  
  down_LOC_MSU <- down_genes_df %>%
    left_join(rap_msu, by = c("gene" = "LOC_gene")) %>%
    filter(!is.na(RAP))
  
  # Mapping 2: RAP → ENTREZID
  up_RAP_ENTREZID <- up_LOC_MSU %>%
    left_join(rap_ncbi, by = "RAP", relationship = "many-to-many") %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(gene, .keep_all = TRUE)
  
  down_RAP_ENTREZID <- down_LOC_MSU %>%
    left_join(rap_ncbi, by = "RAP", relationship = "many-to-many") %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(gene, .keep_all = TRUE)
  
  # KEGG enrichment for up-regulated genes
  entrez_up <- unique(as.character(up_RAP_ENTREZID$ENTREZID))
  ekegg_up <- NULL
  if(length(entrez_up) > 0){
    ekegg_up <- enrichKEGG(
      gene          = entrez_up,
      organism      = osa_gson,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05
    )
  }
  
  # KEGG enrichment for down-regulated genes
  entrez_down <- unique(as.character(down_RAP_ENTREZID$ENTREZID))
  ekegg_down <- NULL
  if(length(entrez_down) > 0){
    ekegg_down <- enrichKEGG(
      gene          = entrez_down,
      organism      = osa_gson,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05
    )
  }
  return(list(kegg_up = ekegg_up, kegg_down = ekegg_down))
}

# ---------------------- Main workflow ----------------------
# Define all comparison groups
contrast_list <- list(
  TB_24vs0 = list(c("time_24_vs_0","genotypeTB.time24")),
  TB_72vs0 = list(c("time_72_vs_0","genotypeTB.time72")),
  NJ_24vs0 = list(c("time_24_vs_0","genotypeNJ.time24")),
  NJ_72vs0 = list(c("time_72_vs_0","genotypeNJ.time72")),
  Y_24vs0 = list(c("time_24_vs_0","genotypeY.time24")),
  Y_72vs0 = list(c("time_72_vs_0","genotypeY.time72")),
  M_24vs0 = list(c("time_24_vs_0","genotypeM.time24")),
  M_72vs0 = list(c("time_72_vs_0","genotypeM.time72")),
  N_24vs0 = list(c("time_24_vs_0","genotypeN.time24")),
  N_72vs0 = list(c("time_72_vs_0","genotypeN.time72")),
  T_24vs0 = list(c("time_24_vs_0","genotypeT.time24")),
  T_72vs0 = list(c("time_72_vs_0","genotypeT.time72")),
  KO_24vs0 = list(c("time_24_vs_0","genotypeKO.time24")),
  KO_72vs0 = list(c("time_72_vs_0","genotypeKO.time72"))
)

all_kegg_allcomparison <- list()
# Run KEGG enrichment for all contrasts
for(nm in names(contrast_list)){
  cat("\n===== Running contrast: ",nm," =====\n")
  deg_tmp <- get_deg_vec(dds,contrast_list[[nm]])
  kegg_tmp <- run_general_kegg(
    up_gene_vec = deg_tmp$up,
    down_gene_vec = deg_tmp$down,
    rap_msu = rap_msu,
    rap_ncbi = rap_ncbi,
    osa_gson = osa_gson
  )
  all_kegg_allcomparison[[nm]] <- kegg_tmp
}

# Save workspace
save(rap_msu,rap_ncbi,contrast_list,all_go_allcomparison,all_kegg_allcomparison,all_go_simplify,groups,mart,get_deg_vec,get_entrez,run_general_go,
     file = "go_kegg_batch.RData")

# Output folder
out_dir <- "./go_kegg_comparison_result"
if(!dir.exists(out_dir)) dir.create(out_dir,recursive = TRUE)

# Export KEGG to xlsx
for (i in names(all_kegg_allcomparison)) {
  if(!is.null(all_kegg_allcomparison[[i]]$kegg_up)){
    write.xlsx(as.data.frame(all_kegg_allcomparison[[i]]$kegg_up), file.path(out_dir, sprintf("%s_up.xlsx", i)))
  }
  if(!is.null(all_kegg_allcomparison[[i]]$kegg_down)){
    write.xlsx(as.data.frame(all_kegg_allcomparison[[i]]$kegg_down), file.path(out_dir, sprintf("%s_down.xlsx", i)))
  }
}

# Export KEGG to csv
for(nm in names(all_kegg_allcomparison)){
  item <- all_kegg_allcomparison[[nm]]
  if(!is.null(item$kegg_up)){
    write_csv(as.data.frame(item$kegg_up),file.path(out_dir,paste0(nm,"_up_KEGG.csv")))
  }
  if(!is.null(item$kegg_down)){
    write_csv(as.data.frame(item$kegg_down),file.path(out_dir,paste0(nm,"_down_KEGG.csv")))
  }
}

# Single plot test
dotplot(all_kegg_allcomparison[["KO_72vs0"]]$kegg_up, showCategory=10, color="p.adjust", size="Count") + theme_bw()

# Batch plot KEGG dotplot
out_dir_kegg <- "./kegg_results"
dir.create(out_dir_kegg, recursive = TRUE, showWarnings = FALSE)
for(nm in names(all_kegg_allcomparison)){
  kegg_item <- all_kegg_allcomparison[[nm]]
  # Up-regulated
  if(!is.null(kegg_item$kegg_up)){
    p_kegg_up <- dotplot(kegg_item$kegg_up, showCategory = 10, color="p.adjust", size="Count") +
      theme_bw() + ggtitle(paste0(nm,"_Up_KEGG"))
    ggsave(file.path(out_dir_kegg,paste0(nm,"_up_kegg.pdf")),p_kegg_up,width=9,height=6)
  }
  # Down-regulated
  if(!is.null(kegg_item$kegg_down)){
    p_kegg_down <- dotplot(kegg_item$kegg_down, showCategory = 10, color="p.adjust", size="Count") +
      theme_bw() + ggtitle(paste0(nm,"_Down_KEGG"))
    ggsave(file.path(out_dir_kegg,paste0(nm,"_down_kegg.pdf")),p_kegg_down,width=9,height=6)
  }
}
