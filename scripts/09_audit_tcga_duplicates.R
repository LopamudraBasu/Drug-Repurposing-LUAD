metadata <- read.csv(
  "results/tables/TCGA_LUAD_sample_metadata.csv",
  stringsAsFactors = FALSE
)

cat("Samples per patient:\n")
print(sort(table(metadata$patient_id), decreasing = TRUE)[1:15])

cat("\nSample groups per patient:\n")
patient_groups <- aggregate(
  group ~ patient_id,
  data = metadata,
  FUN = function(x) paste(sort(unique(x)), collapse = " + ")
)
print(table(patient_groups$group))

cat("\nPatients with both tumor and normal samples:\n")
print(head(patient_groups[grepl("\\+", patient_groups$group), ], 15))

cat("\nFull barcode examples for patients with multiple samples:\n")
duplicate_patients <- names(which(table(metadata$patient_id) > 1))
print(head(metadata[metadata$patient_id %in% duplicate_patients,
                    c("sample_id", "patient_id", "sample_code", "group")], 20))
