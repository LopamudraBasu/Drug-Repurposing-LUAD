suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
})

input_file <- "results/tables/DESeq2_paired_all_genes.csv"
res <- read.csv(input_file, stringsAsFactors = FALSE)

# Remove Ensembl version suffixes, e.g. ENSG00000000003.15
res$ensembl_id <- sub("\\..*$", "", res$gene_id)

# Keep valid, unique Ensembl IDs for enrichment.
res <- res[!is.na(res$ensembl_id) & !duplicated(res$ensembl_id), ]

tested_ids <- unique(res$ensembl_id)
sig <- res[!is.na(res$padj) & res$padj < 0.05, ]
up <- sig[!is.na(sig$log2FoldChange) & sig$log2FoldChange > 1, ]
down <- sig[!is.na(sig$log2FoldChange) & sig$log2FoldChange < -1, ]

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

run_go <- function(gene_ids, label) {
  gene_ids <- unique(gene_ids)
  cat("\n", label, "- input genes:", length(gene_ids), "\n")

  if (length(gene_ids) < 10) {
    cat("Too few genes; skipping enrichment.\n")
    return(invisible(NULL))
  }

  ego <- enrichGO(
    gene = gene_ids,
    universe = tested_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )

  out_csv <- paste0("results/tables/GO_BP_", label, ".csv")
  write.csv(as.data.frame(ego), out_csv, row.names = FALSE)

  if (nrow(as.data.frame(ego)) > 0) {
    p <- dotplot(ego, showCategory = 20) +
      ggtitle(paste("GO Biological Process:", label))

    ggsave(
      paste0("figures/GO_BP_", label, "_dotplot.png"),
      p, width = 11, height = 8, dpi = 300
    )

    cat("Enriched terms:", nrow(as.data.frame(ego)), "\n")
    print(head(as.data.frame(ego)[,
      c("ID", "Description", "GeneRatio", "p.adjust", "Count")], 10))
  } else {
    cat("No terms passed the selected significance cutoffs.\n")
  }

  invisible(ego)
}

cat("Genes tested:", length(tested_ids), "\n")
cat("Significant genes:", nrow(sig), "\n")
cat("Upregulated genes (padj < 0.05, log2FC > 1):", nrow(up), "\n")
cat("Downregulated genes (padj < 0.05, log2FC < -1):", nrow(down), "\n")

run_go(up$ensembl_id, "upregulated")
run_go(down$ensembl_id, "downregulated")

cat("\nPathway enrichment finished.\n")
