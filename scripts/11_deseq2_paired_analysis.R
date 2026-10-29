suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

counts <- readRDS("data/processed/TCGA_LUAD_unstranded_counts.rds")
metadata <- read.csv(
  "results/tables/TCGA_LUAD_selected_samples.csv",
  stringsAsFactors = FALSE
)

# Keep only patients with one selected tumor and one selected normal sample.
paired_ids <- names(which(table(metadata$patient_id) == 2))
paired_meta <- metadata[metadata$patient_id %in% paired_ids, ]

# Confirm that every included patient has both tissue groups.
group_check <- table(paired_meta$patient_id, paired_meta$group)
if (!all(c("Primary Tumor", "Solid Tissue Normal") %in% colnames(group_check)) ||
    !all(group_check[, c("Primary Tumor", "Solid Tissue Normal")] == 1)) {
  stop("Pairing check failed. Inspect the selected sample metadata.")
}

# Align counts exactly to the selected metadata order.
paired_counts <- counts[, paired_meta$sample_id, drop = FALSE]

stopifnot(identical(colnames(paired_counts), paired_meta$sample_id))
if (any(paired_counts < 0) ||
    any(abs(paired_counts - round(paired_counts)) > 0)) {
  stop("Counts are not non-negative integers.")
}

# Set normal tissue as the reference group.
paired_meta$patient_id <- factor(paired_meta$patient_id)
paired_meta$group <- factor(
  paired_meta$group,
  levels = c("Solid Tissue Normal", "Primary Tumor")
)
rownames(paired_meta) <- paired_meta$sample_id

# Remove genes with very low counts.
keep <- rowSums(paired_counts >= 10) >= 10
paired_counts <- round(paired_counts[keep, , drop = FALSE])

cat("Paired samples:", ncol(paired_counts), "\n")
cat("Matched patients:", nlevels(paired_meta$patient_id), "\n")
cat("Genes retained:", nrow(paired_counts), "\n")

# Paired model: patient effects are controlled by the design.
dds <- DESeqDataSetFromMatrix(
  countData = paired_counts,
  colData = paired_meta,
  design = ~ patient_id + group
)

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("group", "Primary Tumor", "Solid Tissue Normal"),
  alpha = 0.05
)

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/logs", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

write.csv(
  res_df,
  "results/tables/DESeq2_paired_all_genes.csv",
  row.names = FALSE
)

sig <- subset(res_df, !is.na(padj) & padj < 0.05)
write.csv(
  sig,
  "results/tables/DESeq2_paired_significant_genes.csv",
  row.names = FALSE
)

write.csv(
  counts(dds, normalized = TRUE),
  "results/tables/DESeq2_paired_normalized_counts.csv"
)

# Volcano plot.
volcano_df <- res_df
volcano_df$neglog10_padj <- -log10(pmax(volcano_df$padj, 1e-300))
volcano_df$status <- ifelse(
  !is.na(volcano_df$padj) & volcano_df$padj < 0.05 &
    volcano_df$log2FoldChange > 1,
  "Up in tumor",
  ifelse(
    !is.na(volcano_df$padj) & volcano_df$padj < 0.05 &
      volcano_df$log2FoldChange < -1,
    "Down in tumor",
    "Not significant"
  )
)

p_volcano <- ggplot(
  volcano_df[!is.na(volcano_df$padj), ],
  aes(x = log2FoldChange, y = neglog10_padj, color = status)
) +
  geom_point(alpha = 0.55, size = 1) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "TCGA-LUAD paired differential expression",
    x = "Log2 fold change (tumor vs normal)",
    y = "-Log10 adjusted p-value",
    color = "Result"
  ) +
  theme_bw()

ggsave(
  "figures/DESeq2_paired_volcano.png",
  p_volcano, width = 9, height = 7, dpi = 300
)

# PCA from variance-stabilized counts.
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = group)) +
  geom_point(size = 2.5, alpha = 0.8) +
  labs(
    title = "TCGA-LUAD paired cohort PCA",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = "Tissue"
  ) +
  theme_bw()

ggsave("figures/DESeq2_paired_PCA.png", p_pca,
       width = 8, height = 6, dpi = 300)

# Heatmap of the 30 most variable retained genes.
mat <- assay(vsd)
gene_vars <- apply(mat, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[
  seq_len(min(30, length(gene_vars)))
]

annotation <- data.frame(group = paired_meta$group)
rownames(annotation) <- rownames(paired_meta)

png("figures/DESeq2_paired_top_variable_heatmap.png",
    width = 1800, height = 1800, res = 180)
pheatmap(
  mat[top_genes, , drop = FALSE],
  scale = "row",
  annotation_col = annotation,
  show_colnames = FALSE,
  fontsize_row = 6,
  main = "Top variable genes: matched TCGA-LUAD samples"
)
dev.off()

writeLines(
  capture.output(sessionInfo()),
  "results/logs/DESeq2_paired_sessionInfo.txt"
)

cat("\nAnalysis complete.\n")
cat("Matched patients:", nlevels(paired_meta$patient_id), "\n")
cat("Samples analysed:", ncol(paired_counts), "\n")
cat("Genes tested:", nrow(res_df), "\n")
cat("Genes with adjusted p < 0.05:", sum(!is.na(res_df$padj) & res_df$padj < 0.05), "\n")
cat("Upregulated (padj < 0.05, log2FC > 1):",
    sum(!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange > 1), "\n")
cat("Downregulated (padj < 0.05, log2FC < -1):",
    sum(!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < -1), "\n")
