counts <- readRDS("data/processed/TCGA_LUAD_unstranded_counts.rds")

sample_ids <- colnames(counts)

# Extract the TCGA sample-type code from each barcode
sample_code <- sub(
  "^TCGA-[^-]+-[^-]+-([0-9]{2}).*$",
  "\\1",
  sample_ids
)

group <- ifelse(
  sample_code == "01", "Primary Tumor",
  ifelse(sample_code == "11", "Solid Tissue Normal", NA_character_)
)

metadata <- data.frame(
  sample_id = sample_ids,
  patient_id = sub(
    "^((TCGA-[^-]+-[^-]+)).*$",
    "\\1",
    sample_ids
  ),
  sample_code = sample_code,
  group = group,
  stringsAsFactors = FALSE
)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

write.csv(
  metadata,
  "results/tables/TCGA_LUAD_sample_metadata.csv",
  row.names = FALSE
)

cat("Count matrix dimensions:", dim(counts), "\n\n")
cat("Sample-type codes:\n")
print(table(metadata$sample_code, useNA = "ifany"))

cat("\nGroup counts:\n")
print(table(metadata$group, useNA = "ifany"))

cat("\nUnique patients:", length(unique(metadata$patient_id)), "\n")
cat("Duplicated patient IDs:", sum(duplicated(metadata$patient_id)), "\n")
cat("Unclassified samples:", sum(is.na(metadata$group)), "\n")

cat("\nMetadata saved to results/tables/TCGA_LUAD_sample_metadata.csv\n")
