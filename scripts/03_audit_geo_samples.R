meta <- read.csv(
  "data/processed/GSE140343_metadata_1.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Total samples:", nrow(meta), "\n\n")

fields <- c(
  "disease state:ch1",
  "tissue:ch1",
  "source_name_ch1",
  "supplementary_file_1"
)

for (field in fields) {
  cat("\n============================\n")
  cat(field, "\n")
  cat("============================\n")

  if (field %in% colnames(meta)) {
    print(table(meta[[field]], useNA = "ifany"))
  } else {
    cat("Field not found\n")
  }
}

cat("\nSample titles:\n")
print(meta[, c("geo_accession", "title")][1:min(15, nrow(meta)), ])

cat("\nSupplementary file examples:\n")
if ("supplementary_file_1" %in% colnames(meta)) {
  print(head(unique(meta$supplementary_file_1), 10))
}

write.csv(
  meta,
  "results/tables/GSE140343_sample_audit.csv",
  row.names = FALSE
)

cat("\nAudit saved to results/tables/GSE140343_sample_audit.csv\n")
