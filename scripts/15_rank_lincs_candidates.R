up_file <- "results/tables/LINCS_reversal_LUAD_up_vs_drug_down.csv"
down_file <- "results/tables/LINCS_reversal_LUAD_down_vs_drug_up.csv"

up <- read.csv(up_file, stringsAsFactors = FALSE)
down <- read.csv(down_file, stringsAsFactors = FALSE)

# Extract compound names from terms such as:
# "LJP006 HT29 24H-palbociclib-0.37"
extract_compound <- function(term) {
  # Remove plate and cell-line fields plus the perturbation duration.
  # Anchor to the first three whitespace-separated fields so an H in
  # a compound name (e.g. XMD11-85H-10) is not mistaken for duration.
  x <- sub("^[^[:space:]]+[[:space:]]+[^[:space:]]+[[:space:]]+[0-9]+H-", "", term)
  # Remove the final dose suffix.
  x <- sub("-[0-9.]+$", "", x)
  trimws(x)
}

summarize_compounds <- function(dat, direction_label) {
  dat$compound <- vapply(dat$term, extract_compound, character(1))
  dat <- dat[
    !is.na(dat$adjusted_p_value) &
      dat$adjusted_p_value < 0.05 &
      !is.na(dat$compound) & nzchar(dat$compound),
  ]

  compounds <- split(dat, dat$compound)

  out <- lapply(names(compounds), function(nm) {
    x <- compounds[[nm]]
    best <- x[which.min(x$adjusted_p_value), ]

    data.frame(
      compound = nm,
      best_adjusted_p = best$adjusted_p_value,
      best_combined_score = max(x$combined_score, na.rm = TRUE),
      number_of_signatures = nrow(x),
      best_signature = best$term,
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, out)
  names(result)[-1] <- paste0(names(result)[-1], "_", direction_label)
  result
}

up_compounds <- summarize_compounds(up, "up_reversal")
down_compounds <- summarize_compounds(down, "down_reversal")

# Keep compounds present in both opposing-direction comparisons.
candidates <- merge(up_compounds, down_compounds, by = "compound")

if (nrow(candidates) > 0) {
  candidates$combined_evidence_score <-
    -log10(pmax(candidates$best_adjusted_p_up_reversal, 1e-300)) +
    -log10(pmax(candidates$best_adjusted_p_down_reversal, 1e-300))

  candidates <- candidates[
    order(-candidates$combined_evidence_score),
  ]
}

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

write.csv(
  up_compounds,
  "results/tables/LINCS_compounds_reversing_up_signature.csv",
  row.names = FALSE
)
write.csv(
  down_compounds,
  "results/tables/LINCS_compounds_reversing_down_signature.csv",
  row.names = FALSE
)
write.csv(
  candidates,
  "results/tables/LINCS_compounds_in_both_screens.csv",
  row.names = FALSE
)

cat("Unique compounds in LUAD-up vs drug-down screen:",
    nrow(up_compounds), "\n")
cat("Unique compounds in LUAD-down vs drug-up screen:",
    nrow(down_compounds), "\n")
cat("Compounds appearing in both screens:", nrow(candidates), "\n")

if (nrow(candidates) > 0) {
  cat("\nTop 20 compounds appearing in both screens:\n")
  print(head(candidates[, c(
    "compound",
    "best_adjusted_p_up_reversal",
    "best_adjusted_p_down_reversal",
    "number_of_signatures_up_reversal",
    "number_of_signatures_down_reversal",
    "combined_evidence_score"
  )], 20))
} else {
  cat("No compound names overlapped after parsing. Inspect the term formats.\n")
}
