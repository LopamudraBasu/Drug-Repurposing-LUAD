packages <- c(
  "DESeq2",
  "TCGAbiolinks",
  "GEOquery",
  "limma",
  "clusterProfiler",
  "org.Hs.eg.db",
  "ggplot2",
  "pheatmap"
)

audit <- data.frame(
  package = packages,
  version = vapply(
    packages,
    function(pkg) {
      if (requireNamespace(pkg, quietly = TRUE)) {
        as.character(packageVersion(pkg))
      } else {
        NA_character_
      }
    },
    character(1)
  )
)

print(audit, row.names = FALSE)

write.csv(
  audit,
  "results/tables/software_versions.csv",
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  "results/logs/R_sessionInfo.txt"
)
