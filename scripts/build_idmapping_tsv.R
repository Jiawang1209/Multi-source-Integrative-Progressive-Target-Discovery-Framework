#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  parsed <- list(
    output = "src/miptd/resources/idmapping.tsv"
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    value <- if (i + 1 <= length(args)) args[[i + 1]] else NULL
    if (key == "--output") parsed$output <- value
    i <- i + 2
  }

  parsed
}

opt <- parse_args(args)
out_path <- opt$output
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = keys(org.Hs.eg.db, keytype = "UNIPROT"),
  columns = c("SYMBOL"),
  keytype = "UNIPROT"
)

mapping <- unique(mapping[, c("UNIPROT", "SYMBOL")])
names(mapping) <- c("From", "To")
mapping <- mapping[!is.na(mapping$From) & mapping$From != "" & !is.na(mapping$To) & mapping$To != "", ]

write.table(
  mapping,
  file = out_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("Wrote: ", normalizePath(out_path, winslash = "/", mustWork = FALSE))
message("Rows: ", nrow(mapping))
