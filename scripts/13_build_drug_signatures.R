suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

res <- read.csv(
  "results/tables/DESeq2_paired_all_genes.csv",
  stringsAsFactors = FALSE
)

# Remove Ensembl version suffixes.
res$ensembl_id <- sub("\\..*$", "", res$gene_id)
res <- res[!is.na(res$ensembl_id) & !duplicated(res$ensembl_id), ]

# Map IDs to gene symbols and Entrez IDs.
symbols <- mapIds(
  org.Hs.eg.db,
  keys = res$ensembl_id,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

entrez <- mapIds(
  org.Hs.eg.db,
  keys = res$ensembl_id,
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)

res$symbol <- unname(symbols[res$ensembl_id])
res$entrez_id <- unname(entrez[res$ensembl_id])

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

write.csv(
  res,
  "results/tables/DESeq2_results_with_gene_symbols.csv",
  row.names = FALSE
)

sig_up <- res[
  !is.na(res$padj) & res$padj < 0.05 &
    !is.na(res$log2FoldChange) & res$log2FoldChange > 1 &
    !is.na(res$symbol),
]
sig_up <- sig_up[order(-sig_up$log2FoldChange), ]

sig_down <- res[
  !is.na(res$padj) & res$padj < 0.05 &
    !is.na(res$log2FoldChange) & res$log2FoldChange < -1 &
    !is.na(res$symbol),
]
sig_down <- sig_down[order(sig_down$log2FoldChange), ]

write.csv(
  head(sig_up, 150),
  "results/tables/Drug_signature_up_top150.csv",
  row.names = FALSE
)
write.csv(
  head(sig_down, 150),
  "results/tables/Drug_signature_down_top150.csv",
  row.names = FALSE
)

# GSEA uses the Wald statistic across the ranked tested genes.
rank_data <- res[
  !is.na(res$stat) & !is.na(res$entrez_id),
  c("entrez_id", "stat")
]
rank_data <- rank_data[order(abs(rank_data$stat), decreasing = TRUE), ]
rank_data <- rank_data[!duplicated(rank_data$entrez_id), ]

gene_list <- rank_data$stat
names(gene_list) <- rank_data$entrez_id
gene_list <- sort(gene_list, decreasing = TRUE)

gsea <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pAdjustMethod = "BH",
  verbose = FALSE
)

gsea_df <- as.data.frame(gsea)
write.csv(
  gsea_df,
  "results/tables/GSEA_GO_BP_ranked.csv",
  row.names = FALSE
)

cat("Total tested genes:", nrow(res), "\n")
cat("Mapped gene symbols:", sum(!is.na(res$symbol)), "\n")
cat("Mapped Entrez IDs:", sum(!is.na(res$entrez_id)), "\n")
cat("Upregulated signature genes with symbols:", nrow(sig_up), "\n")
cat("Downregulated signature genes with symbols:", nrow(sig_down), "\n")
cat("Genes in ranked GSEA:", length(gene_list), "\n")
cat("GSEA terms returned:", nrow(gsea_df), "\n")

if (nrow(gsea_df) > 0) {
  print(head(
    gsea_df[order(gsea_df$p.adjust),
            intersect(c("ID", "Description", "NES", "p.adjust", "setSize"),
                      names(gsea_df))],
    10
  ))
}

cat("\nDrug signatures and ranked GSEA files saved.\n")
