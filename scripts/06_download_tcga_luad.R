suppressPackageStartupMessages({
  library(TCGAbiolinks)
})

out_dir <- "data/raw/TCGA-LUAD"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Querying TCGA-LUAD STAR counts...")

query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

sample_info <- getResults(query)

write.csv(
  sample_info,
  file = "results/tables/TCGA_LUAD_GDC_query_manifest.csv",
  row.names = FALSE
)

cat("Files in query:", nrow(sample_info), "\n")
cat("\nSample types:\n")
print(table(sample_info$sample_type, useNA = "ifany"))

cat("\nQuery manifest saved.\n")
