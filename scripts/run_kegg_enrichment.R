#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(readr)
  library(writexl)
})

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  parsed <- list(case_dir = NULL)
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    value <- if (i + 1 <= length(args)) args[[i + 1]] else NULL
    if (key == "--case-dir") parsed$case_dir <- value
    i <- i + 2
  }
  if (is.null(parsed$case_dir)) {
    stop("Missing --case-dir")
  }
  parsed
}

opt <- parse_args(args)
case_dir <- normalizePath(opt$case_dir, winslash = "/", mustWork = TRUE)
venn_path <- file.path(case_dir, "venn.rds")
out_path <- file.path(case_dir, "2.KEGG_result_2.xlsx")

venn <- readRDS(venn_path)
union_id <- unique(unlist(venn))

geneid_trans <- bitr(
  geneID = union_id,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db"
)

if (nrow(geneid_trans) == 0) {
  stop("No SYMBOL to ENTREZID mapping succeeded.")
}

kegg_result <- enrichKEGG(
  gene = geneid_trans$ENTREZID,
  organism = "hsa"
)

kegg_result_2 <- setReadable(
  kegg_result,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

write_xlsx(as.data.frame(kegg_result_2), out_path)
message("Wrote: ", out_path)
