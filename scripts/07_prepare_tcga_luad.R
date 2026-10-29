suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

out_dir <- "data/raw/TCGA-LUAD"
processed_dir <- "data/processed"
results_dir <- "results/tables"
logs_dir <- "results/logs"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

message("Step 1: Querying TCGA-LUAD STAR counts")

query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

manifest <- getResults(query)

write.csv(
  manifest,
  file.path(results_dir, "TCGA_LUAD_GDC_query_manifest.csv"),
  row.names = FALSE
)

message("Files available: ", nrow(manifest))

message("Step 2: Downloading files in batches of 20")

GDCdownload(
  query,
  directory = out_dir,
  method = "api",
  files.per.chunk = 20
)

message("Step 3: Preparing expression data")

se <- GDCprepare(
  query,
  directory = out_dir,
  summarizedExperiment = TRUE
)

message("Available assays:")
print(assayNames(se))

if (!"unstranded" %in% assayNames(se)) {
  stop(
    "The unstranded count assay was not found. ",
    "Inspect assayNames(se) before proceeding."
  )
}

counts <- assay(se, "unstranded")
metadata <- as.data.frame(colData(se))
gene_info <- as.data.frame(rowData(se))

stopifnot(ncol(counts) == nrow(metadata))

saveRDS(
  counts,
  file.path(processed_dir, "TCGA_LUAD_unstranded_counts.rds")
)

write.csv(
  metadata,
  file.path(results_dir, "TCGA_LUAD_sample_metadata.csv"),
  row.names = TRUE
)

write.csv(
  gene_info,
  file.path(results_dir, "TCGA_LUAD_gene_annotation.csv"),
  row.names = TRUE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(logs_dir, "TCGA_LUAD_preparation_sessionInfo.txt")
)

cat("\n=== Expression matrix dimensions ===\n")
print(dim(counts))

cat("\n=== Sample types ===\n")
if ("sample_type" %in% colnames(metadata)) {
  print(table(metadata$sample_type, useNA = "ifany"))
} else {
  print(names(metadata))
}

cat("\nPrepared count matrix and metadata successfully.\n")
