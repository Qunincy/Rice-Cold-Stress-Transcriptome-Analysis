# ---------------------- 8. Heatmap plotting ----------------------
library(pheatmap)

# Variance stabilization transformation for heatmap (avoid raw read counts)
vst_mat <- vst(dds, blind = FALSE)

# Extract vst expression matrix of selected DE genes
expr_sig_genes <- assay(vst_mat)[sig_genes, ]

# Row-wise Z-score normalization for relative expression
expr_z <- t(scale(t(expr_sig_genes)))

pheatmap(
  expr_z,
  show_rownames = FALSE,   # Hide gene names for large DE gene sets; set TRUE for candidate genes
  show_colnames = TRUE,
  scale = "none",          # Manual Z-score applied already, skip built-in scaling
  cluster_rows = TRUE,     # Cluster genes
  cluster_cols = TRUE,     # Cluster samples
  treeheight_row = 20,
  treeheight_col = 20,
  color = colorRampPalette(c("#2166ac","white","#b2182b"))(100),
  main = "Y_vs_C DE genes heatmap"
)

# Subset samples only for Y and C genotype
anno_col <- sample_info[,c("genotype","time")]
keep_samples <- anno_col$genotype %in% c("Y","C")
expr_z_sub <- expr_z[, keep_samples]

pheatmap(
  expr_z_sub,
  show_rownames = FALSE,
  show_colnames = TRUE,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_row = 20,
  treeheight_col = 20,
  color = colorRampPalette(c("#2166ac","white","#b2182b"))(100),
  main = "Y_vs_C DE genes heatmap"
)

# Heatmap generation function
plot_heatmap_de <- function(dds, contrast_name, plot_title, anno_col, vst_mat, padj_cut = 0.05, lfc_cut  = 1.0){
  # 1. Extract DESeq2 results
  res <- results(dds, name = contrast_name, alpha = 0.05, lfcThreshold = 1)
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_filter <- res_df[!is.na(res_df$padj), ]
  
  # 2. Filter significant DE genes
  sig_genes <- res_filter[res_filter$padj < padj_cut & abs(res_filter$log2FoldChange) >= lfc_cut, "gene", drop = TRUE]
  expr_sig <- assay(vst_mat)[sig_genes, ]
  
  # Skip heatmap if fewer than 2 significant genes
  n_gene <- length(sig_genes)
  if(n_gene < 2){
    message(paste0(contrast_name,": ",n_gene," significant genes found, skip heatmap"))
    return(NULL)
  }
  
  # Parse target genotype from contrast name (genotype_Y_vs_C -> Y)
  curr_geno <- str_extract(contrast_name, "(?<=genotype_).*(?=_vs_C)")
  
  # 3. Keep only target genotype and control C samples
  keep_samples <- anno_col$genotype %in% c(curr_geno,"C")
  expr_sub <- expr_sig[, keep_samples]
  
  # 4. Row Z-score normalization, replace NaN generated during scaling
  expr_z_sub <- t(scale(t(expr_sub)))
  expr_z_sub[is.na(expr_z_sub)] <- 0
  
  # 5. Plot heatmap
  p_heat <- pheatmap(
    expr_z_sub,
    show_rownames = FALSE,
    show_colnames = TRUE,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    treeheight_row = 20,
    treeheight_col = 20,
    color = colorRampPalette(c("#2166ac","white","#b2182b"))(100),
    main = plot_title
  )
  return(p_heat)
}

# Prepare sample annotation and vst matrix (run once globally)
anno_col <- sample_info[,c("genotype","time")]
vst_mat <- vst(dds, blind = FALSE)

# Generate heatmaps for multiple genotype comparisons
h_KvsC <- plot_heatmap_de(dds,"genotype_K_vs_C","genotype_K_vs_C DE genes heatmap", anno_col, vst_mat)
h_RvsC <- plot_heatmap_de(dds,"genotype_R_vs_C","genotype_R_vs_C DE genes heatmap", anno_col, vst_mat)
h_LvsC <- plot_heatmap_de(dds,"genotype_L_vs_C","genotype_L_vs_C DE genes heatmap", anno_col, vst_mat)
h_TBvsC <- plot_heatmap_de(dds,"genotype_TB_vs_C","genotype_TB_vs_C DE genes heatmap", anno_col, vst_mat)
h_NJvsC <- plot_heatmap_de(dds,"genotype_NJ_vs_C","genotype_NJ_vs_C DE genes heatmap", anno_col, vst_mat)
h_YvsC <- plot_heatmap_de(dds,"genotype_Y_vs_C","genotype_Y_vs_C DE genes heatmap", anno_col, vst_mat)
h_MvsC <- plot_heatmap_de(dds,"genotype_M_vs_C","genotype_M_vs_C DE genes heatmap", anno_col, vst_mat)
h_NvsC <- plot_heatmap_de(dds,"genotype_N_vs_C","genotype_N_vs_C DE genes heatmap", anno_col, vst_mat)
h_TvsC <- plot_heatmap_de(dds,"genotype_T_vs_C","genotype_T_vs_C DE genes heatmap", anno_col, vst_mat)
h_KOvsC <- plot_heatmap_de(dds,"genotype_KO_vs_C","genotype_KO_vs_C DE genes heatmap", anno_col, vst_mat)

# Print heatmaps
print(h_KvsC)
print(h_RvsC)
print(h_LvsC)
print(h_TBvsC)
print(h_NJvsC)
print(h_YvsC)
print(h_MvsC)
print(h_NvsC)
print(h_TvsC)
print(h_KOvsC)

# Output directory setup
out_dir <- "heatmap_output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
plot_w <-1810
plot_h <-1062
plot_dpi <-144

heat_list <- list(
  "K_vs_C" = h_KvsC,
  "R_vs_C" = h_RvsC,
  "L_vs_C" = h_LvsC,
  "TB_vs_C" = h_TBvsC,
  "NJ_vs_C" = h_NJvsC,
  "Y_vs_C" = h_YvsC,
  "M_vs_C" = h_MvsC,
  "N_vs_C" = h_NvsC,
  "T_vs_C" = h_TvsC,
  "KO_vs_C" = h_KOvsC
)

# Export individual heatmaps
for (plot_name in names(heat_list)) {
  ggsave(
    filename = file.path(out_dir, paste0(plot_name, "_DE_heatmap.png")),
    plot = heat_list[[plot_name]],
    width = plot_w,
    height = plot_h,
    dpi = plot_dpi,
    units = "px",
    device = "png"
  )
}

# Combine all heatmaps into one figure
library(patchwork)
library(grid)
heat_grob_list <- lapply(heat_list, function(x){
  x$gtable
})
big_heat <- wrap_plots(heat_grob_list, nrow =5, ncol =2)
ggsave(
  filename = file.path(out_dir,"all_heatmap_one_big.png"),
  plot = big_heat,
  width = 1810*2,
  height = 1062*5,
  dpi = 144,
  units = "px",
  device = "png"
)

# Save R workspace
save(list = ls(), file = "full_rna_seq_workspace.RData")
