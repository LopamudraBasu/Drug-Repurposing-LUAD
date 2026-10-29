suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
})

base_url <- "https://maayanlab.cloud/Enrichr"

input_file <- "results/tables/DESeq2_results_with_gene_symbols.csv"
if (!file.exists(input_file)) {
  stop("Missing input file: ", input_file,
       ". Run script 13 first.")
}

res <- read.csv(input_file, stringsAsFactors = FALSE)

get_signature <- function(direction) {
  x <- res[
    !is.na(res$padj) & res$padj < 0.05 &
      !is.na(res$log2FoldChange) &
      !is.na(res$symbol) & nzchar(res$symbol),
  ]

  if (direction == "up") {
    x <- x[x$log2FoldChange > 1, ]
    x <- x[order(-x$log2FoldChange), ]
  } else {
    x <- x[x$log2FoldChange < -1, ]
    x <- x[order(x$log2FoldChange), ]
  }

  # Remove duplicated symbols; keep the strongest-ranked occurrence.
  x <- x[!duplicated(x$symbol), ]
  unique(x$symbol)
}

submit_list <- function(genes, description) {
  response <- POST(
    paste0(base_url, "/addList"),
    body = list(
      list = paste(genes, collapse = "\n"),
      description = description
    ),
    encode = "multipart",
    timeout(60)
  )

  stop_for_status(response)
  parsed <- fromJSON(content(response, "text", encoding = "UTF-8"))

  if (is.null(parsed$userListId)) {
    stop("Enrichr did not return a userListId for: ", description)
  }

  parsed$userListId
}

query_library <- function(list_id, library_name) {
  response <- GET(
    paste0(base_url, "/enrich"),
    query = list(
      userListId = list_id,
      backgroundType = library_name
    ),
    timeout(120)
  )

  stop_for_status(response)
  parsed <- fromJSON(
    content(response, "text", encoding = "UTF-8"),
    simplifyVector = FALSE
  )

  entries <- parsed[[library_name]]
  if (is.null(entries) || length(entries) == 0) {
    return(data.frame())
  }

  rows <- lapply(entries, function(z) {
    data.frame(
      rank = as.integer(z[[1]]),
      term = as.character(z[[2]]),
      pvalue = as.numeric(z[[3]]),
      zscore = as.numeric(z[[4]]),
      combined_score = as.numeric(z[[5]]),
      overlapping_genes = paste(unlist(z[[6]]), collapse = ";"),
      adjusted_p_value = as.numeric(z[[7]]),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out[order(out$adjusted_p_value, -out$combined_score), ]
}

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/logs", recursive = TRUE, showWarnings = FALSE)

luad_up <- get_signature("up")
luad_down <- get_signature("down")

cat("LUAD up signature genes:", length(luad_up), "\n")
cat("LUAD down signature genes:", length(luad_down), "\n")

if (length(luad_up) < 10 || length(luad_down) < 10) {
  stop("Too few mapped genes for a useful screen.")
}

# Disease-up genes are compared with genes decreased by compounds.
cat("\nSubmitting LUAD-up signature...\n")
up_id <- submit_list(luad_up, "TCGA-LUAD tumor-up genes")
Sys.sleep(2)

cat("Querying LINCS_L1000_Chem_Pert_down...\n")
up_reversed <- query_library(up_id, "LINCS_L1000_Chem_Pert_down")
write.csv(
  up_reversed,
  "results/tables/LINCS_reversal_LUAD_up_vs_drug_down.csv",
  row.names = FALSE
)

Sys.sleep(2)

# Disease-down genes are compared with genes increased by compounds.
cat("\nSubmitting LUAD-down signature...\n")
down_id <- submit_list(luad_down, "TCGA-LUAD tumor-down genes")
Sys.sleep(2)

cat("Querying LINCS_L1000_Chem_Pert_up...\n")
down_reversed <- query_library(down_id, "LINCS_L1000_Chem_Pert_up")
write.csv(
  down_reversed,
  "results/tables/LINCS_reversal_LUAD_down_vs_drug_up.csv",
  row.names = FALSE
)

cat("\nResults returned:\n")
cat("Drug-down signatures enriched against LUAD-up:",
    nrow(up_reversed), "\n")
cat("Drug-up signatures enriched against LUAD-down:",
    nrow(down_reversed), "\n")

cat("\nTop 10 results: LUAD-up vs drug-down\n")
if (nrow(up_reversed) > 0) {
  print(head(up_reversed[, c(
    "term", "adjusted_p_value", "combined_score",
    "overlapping_genes"
  )], 10))
}

cat("\nTop 10 results: LUAD-down vs drug-up\n")
if (nrow(down_reversed) > 0) {
  print(head(down_reversed[, c(
    "term", "adjusted_p_value", "combined_score",
    "overlapping_genes"
  )], 10))
}

writeLines(
  capture.output(sessionInfo()),
  "results/logs/14_lincs_reversal_sessionInfo.txt"
)

cat("\nLINCS first-pass screen completed.\n")
