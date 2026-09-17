# Load required packages
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(dplyr)

# ---------------------- 6. Select top 5 most significant up/down-regulated genes for labeling ----------------------
top_up    <- res_filter %>% filter(sig == "Up")   %>% arrange(padj) %>% head(5)
top_down  <- res_filter %>% filter(sig == "Down") %>% arrange(padj) %>% head(5)
top_genes <- rbind(top_up, top_down)

# Count up/down regulated genes for legend annotation
tb     <- table(res_filter$sig)
down_n <- tb["Down"]
up_n   <- tb["Up"]

# ---------------------- 7. Volcano plot function (remove NA, enable gene label annotation) ----------------------
plot_volcano_de <- function(dds, contrast_name, plot_title){
  # Extract differential expression results
  res <- results(dds, name = contrast_name, alpha = 0.05, lfcThreshold = 1)
  
  # Mark significantly DE genes (padj < 0.05 & |log2FoldChange| >=1)
  res$sig <- "Not significant"
  res$sig[res$padj < 0.05 & res$log2FoldChange >= 1] <- "Up"
  res$sig[res$padj < 0.05 & res$log2FoldChange <= -1] <- "Down"
  
  tb <- table(res$sig)
  down_n <- tb["Down"]
  up_n <- tb["Up"]
  
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  # Select top 5 significant genes
  top_up <- res_df %>%
    filter(sig == "Up") %>%
    arrange(padj) %>%
    head(5)
  
  top_down <- res_df %>%
    filter(sig == "Down") %>%
    arrange(padj) %>%
    head(5)
  
  top_genes <- rbind(top_up, top_down)
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), colour = sig)) +
    geom_point(alpha = 1, size = 3.5) +
    geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 20, show.legend = FALSE)+
    ylab("-log10(P-adjust)") +
    xlab("log2FoldChange") +
    scale_color_manual(
      values = c("Down"="blue","Not significant"="grey","Up"="red"),
      labels = c(paste0('Down(',down_n,")"),"Not significant",paste0("Up(",up_n,")"))
    ) +
    geom_vline(xintercept = c(-1, 1), lty = 4, col = "black", lwd = 0.8) +
    geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.8) +
    ggtitle(plot_title)+
    theme_bw()
  return(p)
}

# Generate volcano plots for multiple genotype comparisons
p_KvsC <- plot_volcano_de(dds,"genotype_K_vs_C","genotype_K_vs_C")
p_RvsC <- plot_volcano_de(dds,"genotype_R_vs_C","genotype_R_vs_C")
p_LvsC <- plot_volcano_de(dds,"genotype_L_vs_C","genotype_L_vs_C")
p_TBvsC <- plot_volcano_de(dds,"genotype_TB_vs_C","genotype_TB_vs_C")
p_NJvsC <- plot_volcano_de(dds,"genotype_NJ_vs_C","genotype_NJ_vs_C")
p_YvsC <- plot_volcano_de(dds,"genotype_Y_vs_C","genotype_Y_vs_C")
p_MvsC <- plot_volcano_de(dds,"genotype_M_vs_C","genotype_M_vs_C")
p_NvsC <- plot_volcano_de(dds,"genotype_N_vs_C","genotype_N_vs_C")
p_TvsC <- plot_volcano_de(dds,"genotype_T_vs_C","genotype_T_vs_C")
p_KOvsC <- plot_volcano_de(dds,"genotype_KO_vs_C","genotype_KO_vs_C")

# Print plots to console
print(p_KvsC)
print(p_RvsC)
print(p_LvsC)
print(p_TBvsC)
print(p_NJvsC)
print(p_YvsC)
print(p_MvsC)
print(p_NvsC)
print(p_TvsC)
print(p_KOvsC)

# Plot export parameters
out_dir <- "volcano_output"
plot_w <- 1245
plot_h <- 1042
plot_dpi <- 144

# Store plot objects in a named list
plot_list <- list(
  "genotype_K_vs_C" = p_KvsC,
  "genotype_KO_vs_C" = p_KOvsC,
  "genotype_L_vs_C"  = p_LvsC,
  "genotype_M_vs_C"  = p_MvsC,
  "genotype_N_vs_C"  = p_NvsC,
  "genotype_NJ_vs_C" = p_NJvsC,
  "genotype_R_vs_C"  = p_RvsC,
  "genotype_T_vs_C"  = p_TvsC,
  "genotype_TB_vs_C" = p_TBvsC,
  "genotype_Y_vs_C"  = p_YvsC
)

# Export individual volcano plots
for (plot_name in names(plot_list)) {
  ggsave(
    filename = file.path(out_dir, paste0(plot_name, ".png")),
    plot = plot_list[[plot_name]],
    width = plot_w,
    height = plot_h,
    dpi = plot_dpi,
    units = "px",
    device = "png"
  )
}

# Combine all subplots into one figure (2 rows × 5 columns)
big_plot <- wrap_plots(plot_list, nrow = 2, ncol = 5)

# Export combined figure
ggsave(
  filename = file.path(out_dir,"all_volcano_one_big.png"),
  plot = big_plot,
  width = plot_w *5,
  height = plot_h *2,
  dpi = plot_dpi,
  units = "px",
  device = "png"
)
