meta <- read.csv(
  "data/processed/GSE140343_metadata_1.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sample_sheet <- data.frame(
  sample_id = meta$geo_accession,
  title = meta$title,
  source = meta$source_name_ch1,
  disease_state = meta[["disease state:ch1"]],
  tissue = meta[["tissue:ch1"]],
  supplementary_file = meta$supplementary_file_1,
  stringsAsFactors = FALSE
)

sample_sheet$group <- ifelse(
  grepl("clinical normal", sample_sheet$source, ignore.case = TRUE),
  "Normal",
  ifelse(
    grepl("LUAD", sample_sheet$source, ignore.case = TRUE),
    "Tumor",
    NA_character_
  )
)

sample_sheet$patient_id <- sub(
  "^.*_([0-9]+)[TN]\\.gtf\\.gz$",
  "\\1",
  sample_sheet$supplementary_file
)

cat("Group counts:\n")
print(table(sample_sheet$group, useNA = "ifany"))

cat("\nPatient pairing check:\n")
print(table(sample_sheet$patient_id, sample_sheet$group, useNA = "ifany"))

write.csv(
  sample_sheet,
  "data/processed/GSE140343_sample_sheet.csv",
  row.names = FALSE
)

cat("\nSaved sample sheet.\n")
