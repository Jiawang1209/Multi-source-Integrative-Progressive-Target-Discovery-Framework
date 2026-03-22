#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
})

read_table_auto <- function(path, ..., skip = 0, delim = NULL) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    return(read_xlsx(path, skip = skip, ...))
  }
  if (is.null(delim)) {
    delim <- if (ext %in% c("tsv", "txt")) "\t" else ","
  }
  read_delim(path, delim = delim, show_col_types = FALSE, skip = skip, ...)
}

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  parsed <- list(
    case_dir = NULL,
    swiss = NULL,
    sea = NULL,
    chembl = NULL,
    idmapping = NULL,
    ppb2 = NULL,
    kegg = NULL,
    disease_regex = "non-alcoholic fatty liver disease|lipid|fatty|insulin|oxidative|apoptosis|chemokine|MAPK|PI3K|mTOR|AMPK|PPAR|TNF|NF-kappa|FoxO|adipocytokine",
    output_prefix = NULL,
    top_n = 10
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    value <- if (i + 1 <= length(args)) args[[i + 1]] else NULL
    if (key == "--case-dir") parsed$case_dir <- value
    if (key == "--swiss") parsed$swiss <- value
    if (key == "--sea") parsed$sea <- value
    if (key == "--chembl") parsed$chembl <- value
    if (key == "--idmapping") parsed$idmapping <- value
    if (key == "--ppb2") parsed$ppb2 <- value
    if (key == "--kegg") parsed$kegg <- value
    if (key == "--disease-regex") parsed$disease_regex <- value
    if (key == "--output-prefix") parsed$output_prefix <- value
    if (key == "--top-n") parsed$top_n <- as.integer(value)
    i <- i + 2
  }

  if (is.null(parsed$case_dir) || is.null(parsed$swiss) || is.null(parsed$sea) ||
      is.null(parsed$chembl) || is.null(parsed$idmapping) || is.null(parsed$output_prefix)) {
    stop("Missing required arguments.")
  }

  parsed
}

normalize_path <- function(path) {
  if (is.null(path)) return(NULL)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

read_swiss <- function(path) {
  ext <- tolower(tools::file_ext(path))
  swiss_skip <- if (ext %in% c("xlsx", "xls")) 1 else 0
  read_table_auto(path, skip = swiss_skip) %>%
    transmute(
      target = str_to_upper(`Common name`),
      swiss_prob = as.numeric(`Probability*`)
    ) %>%
    filter(!is.na(target), target != "", !is.na(swiss_prob), swiss_prob > 0) %>%
    group_by(target) %>%
    summarise(swiss_prob = max(swiss_prob, na.rm = TRUE), .groups = "drop")
}

read_ppb2 <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    return(tibble(target = character(), ppb2 = integer()))
  }

  read_table_auto(path) %>%
    transmute(target = str_to_upper(Symbol), ppb2 = 1L) %>%
    filter(!is.na(target), target != "") %>%
    distinct()
}

read_sea <- function(path) {
  read_delim(path, delim = ",", show_col_types = FALSE) %>%
    transmute(
      target = str_to_upper(Name),
      sea_p = `P-Value`,
      sea_z = `Z-Score`,
      sea_tc = `Max Tc`
    ) %>%
    filter(!is.na(target), target != "") %>%
    group_by(target) %>%
    summarise(
      sea_p = min(sea_p, na.rm = TRUE),
      sea_z = max(sea_z, na.rm = TRUE),
      sea_tc = max(sea_tc, na.rm = TRUE),
      .groups = "drop"
    )
}

read_chembl <- function(chembl_path, idmapping_path) {
  id_map <- read_tsv(idmapping_path, show_col_types = FALSE) %>%
    distinct(From, To)

  read_table_auto(chembl_path) %>%
    filter(Organism == "Homo sapiens", !is.na(Accessions), Accessions != "") %>%
    select(Accessions) %>%
    separate_rows(Accessions, sep = "\\|") %>%
    left_join(id_map, by = c("Accessions" = "From"), relationship = "many-to-many") %>%
    filter(!is.na(To), To != "") %>%
    transmute(target = str_to_upper(To), chembl = 1L) %>%
    distinct()
}

read_kegg <- function(path, disease_regex) {
  if (is.null(path) || !file.exists(path)) {
    return(list(
      pathways = tibble(),
      genes = tibble(),
      summary = tibble()
    ))
  }

  pathways <- read_table_auto(path) %>%
    filter(grepl(disease_regex, Description, ignore.case = TRUE)) %>%
    select(ID, Description, p.adjust, geneID, Count)

  genes <- pathways %>%
    separate_rows(geneID, sep = "/") %>%
    transmute(
      target = str_to_upper(geneID),
      pathway = Description,
      pathway_id = ID,
      pathway_padj = p.adjust
    )

  summary <- genes %>%
    group_by(target) %>%
    summarise(
      pathway_n = n(),
      best_pathway_padj = min(pathway_padj, na.rm = TRUE),
      pathways = str_c(unique(pathway), collapse = " | "),
      .groups = "drop"
    )

  list(pathways = pathways, genes = genes, summary = summary)
}

write_markdown <- function(path, all_targets, key_targets, top_targets, pathway_table) {
  lines <- c(
    "# Target Prediction Summary",
    "",
    "## Disease-relevant KEGG pathways",
    ""
  )

  if (nrow(pathway_table) == 0) {
    lines <- c(lines, "No KEGG pathway file was provided.")
  } else {
    lines <- c(lines, "| Pathway | Count | Adjusted p-value |", "| --- | ---: | ---: |")
    for (i in seq_len(nrow(pathway_table))) {
      row <- pathway_table[i, ]
      lines <- c(
        lines,
        sprintf("| %s | %s | %.3g |", row$Description, row$Count, row$p.adjust)
      )
    }
  }

  lines <- c(lines, "", "## Top predicted targets", "", "| Target | Score | Platform vote | Pathway count | Swiss | SEA p-value | Sources |", "| --- | ---: | ---: | ---: | ---: | ---: | --- |")

  for (i in seq_len(nrow(top_targets))) {
    row <- top_targets[i, ]
    source_text <- row$sources
    lines <- c(
      lines,
      sprintf(
        "| %s | %.3f | %d | %d | %s | %s | %s |",
        row$target,
        row$consensus_score,
        row$platform_vote,
        row$pathway_n,
        ifelse(is.na(row$swiss_prob), "-", sprintf("%.3f", row$swiss_prob)),
        ifelse(is.na(row$sea_p), "-", sprintf("%.3g", row$sea_p)),
        source_text
      )
    )
  }

  lines <- c(
    lines,
    "",
    sprintf("Total merged targets: %d", nrow(all_targets)),
    sprintf("Key targets after disease filtering: %d", nrow(key_targets))
  )

  writeLines(lines, con = path)
}

opt <- parse_args(args)

case_dir <- normalize_path(opt$case_dir)
swiss_path <- normalize_path(file.path(case_dir, opt$swiss))
sea_path <- normalize_path(file.path(case_dir, opt$sea))
chembl_path <- normalize_path(file.path(case_dir, opt$chembl))
idmapping_path <- normalize_path(file.path(case_dir, opt$idmapping))
ppb2_path <- if (!is.null(opt$ppb2)) normalize_path(file.path(case_dir, opt$ppb2)) else NULL
kegg_path <- if (!is.null(opt$kegg)) normalize_path(file.path(case_dir, opt$kegg)) else NULL
output_prefix <- normalize_path(opt$output_prefix)

dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

swiss <- read_swiss(swiss_path)
ppb2 <- read_ppb2(ppb2_path)
sea <- read_sea(sea_path)
chembl <- read_chembl(chembl_path, idmapping_path)
kegg <- read_kegg(kegg_path, opt$disease_regex)

all_targets <- full_join(swiss, ppb2, by = "target") %>%
  full_join(sea, by = "target") %>%
  full_join(chembl, by = "target") %>%
  full_join(kegg$summary, by = "target") %>%
  mutate(
    source_swiss = !is.na(swiss_prob),
    source_ppb2 = !is.na(ppb2),
    source_sea = !is.na(sea_p),
    source_chembl = !is.na(chembl),
    platform_vote = rowSums(cbind(source_swiss, source_ppb2, source_sea, source_chembl)),
    sea_score = pmin(ifelse(is.na(sea_p), 0, -log10(sea_p) / 20), 5),
    consensus_score = platform_vote + coalesce(swiss_prob, 0) + sea_score + coalesce(pathway_n, 0) / 2,
    sources = pmap_chr(
      list(source_swiss, source_ppb2, source_sea, source_chembl),
      function(a, b, c, d) {
        out <- c()
        if (a) out <- c(out, "Swiss")
        if (b) out <- c(out, "PPB2")
        if (c) out <- c(out, "SEA")
        if (d) out <- c(out, "ChEMBL")
        str_c(out, collapse = ", ")
      }
    )
  ) %>%
  arrange(desc(consensus_score), desc(platform_vote), desc(pathway_n), sea_p, desc(swiss_prob))

key_targets <- all_targets %>%
  filter(!is.na(pathway_n), platform_vote >= 2) %>%
  arrange(desc(consensus_score), desc(pathway_n), desc(platform_vote), sea_p, desc(swiss_prob))

top_targets <- key_targets %>%
  slice_head(n = opt$top_n)

write_csv(all_targets, paste0(output_prefix, "_all_targets.csv"))
write_csv(key_targets, paste0(output_prefix, "_key_targets.csv"))
if (nrow(kegg$pathways) > 0) {
  write_csv(kegg$pathways, paste0(output_prefix, "_disease_pathways.csv"))
}
write_markdown(
  paste0(output_prefix, "_summary.md"),
  all_targets = all_targets,
  key_targets = key_targets,
  top_targets = top_targets,
  pathway_table = kegg$pathways
)

message("Wrote: ", paste0(output_prefix, "_all_targets.csv"))
message("Wrote: ", paste0(output_prefix, "_key_targets.csv"))
message("Wrote: ", paste0(output_prefix, "_summary.md"))
