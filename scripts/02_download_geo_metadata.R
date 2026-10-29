suppressPackageStartupMessages({
  library(GEOquery)
})

geo_id <- "GSE140343"

message("Downloading GEO metadata for ", geo_id)

gse <- getGEO(
  geo_id,
  GSEMatrix = TRUE,
  getGPL = FALSE
)

message("Number of expression series objects: ", length(gse))

for (i in seq_along(gse)) {
  eset <- gse[[i]]
  meta <- pData(eset)

  out <- sprintf(
    "data/processed/GSE140343_metadata_%d.csv",
    i
  )

  write.csv(meta, out, row.names = TRUE)

  cat("\nSeries:", i, "\n")
  cat("Samples:", nrow(meta), "\n")
  cat("Metadata columns:\n")
  print(colnames(meta))
  cat("\nSample titles:\n")
  print(head(meta$title, 10))
  cat("\nSaved:", out, "\n")
}
