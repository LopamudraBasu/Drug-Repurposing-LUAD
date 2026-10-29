metadata <- read.csv(
  "results/tables/TCGA_LUAD_sample_metadata.csv",
  stringsAsFactors = FALSE
)

# Correctly extract the two-digit TCGA sample-type code.
metadata$sample_code <- substr(
  sub("^TCGA-[^-]+-[^-]+-", "", metadata$sample_id),
  1, 2
)

metadata$group <- ifelse(
  metadata$sample_code == "01", "Primary Tumor",
  ifelse(metadata$sample_code == "11", "Solid Tissue Normal", NA)
)

# Extract the vial letter (A, B, etc.) from the barcode.
metadata$vial <- sub(
  "^TCGA-[^-]+-[^-]+-[0-9]{2}([A-Z]).*$",
  "\\1",
  metadata$sample_id
)

# Prefer vial A, then sort by full barcode for reproducibility.
metadata <- metadata[order(
  metadata$patient_id,
  metadata$group,
  metadata$vial != "A",
  metadata$sample_id
), ]

selected <- metadata[
  !is.na(metadata$group) &
    !duplicated(paste(metadata$patient_id, metadata$group)),
]

write.csv(
  selected,
  "results/tables/TCGA_LUAD_selected_samples.csv",
  row.names = FALSE
)

cat("Selected samples by group:\n")
print(table(selected$group))

cat("\nUnique selected patients:", length(unique(selected$patient_id)), "\n")
cat("Duplicate patient-group combinations:",
    sum(duplicated(paste(selected$patient_id, selected$group))), "\n")
cat("\nSaved: results/tables/TCGA_LUAD_selected_samples.csv\n")
