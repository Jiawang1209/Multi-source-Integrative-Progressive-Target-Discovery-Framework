#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(writexl)
  library(ggvenn)
  library(jsonlite)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggnewscale)
  library(legendry)
  library(circlize)
  library(ComplexHeatmap)
})

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  parsed <- list(
    case_dir = NULL,
    swiss = NULL,
    sea = NULL,
    chembl = NULL,
    ppb2 = NULL,
    idmapping = NULL,
    ontology_n = 5L,
    kegg_n = 5L,
    liver_regex = paste(
      c(
        "liver", "hepatic", "hepatocyte", "xenobiotic", "drug",
        "lipid", "fatty", "cholesterol", "sterol", "bile", "biliary",
        "oxidative", "mitochond", "endoplasmic reticulum", "peroxisom",
        "apopt", "inflamm", "transport", "metabolic", "insulin",
        "ampk", "ppar", "tnf", "nf-kappa", "foxo", "chemokine",
        "mapk", "pi3k", "akt", "abc"
      ),
      collapse = "|"
    )
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    value <- if (i + 1 <= length(args)) args[[i + 1]] else NULL
    if (key == "--case-dir") parsed$case_dir <- value
    if (key == "--swiss") parsed$swiss <- value
    if (key == "--sea") parsed$sea <- value
    if (key == "--chembl") parsed$chembl <- value
    if (key == "--ppb2") parsed$ppb2 <- value
    if (key == "--idmapping") parsed$idmapping <- value
    if (key == "--ontology-n") parsed$ontology_n <- as.integer(value)
    if (key == "--kegg-n") parsed$kegg_n <- as.integer(value)
    if (key == "--liver-regex") parsed$liver_regex <- value
    i <- i + 2
  }

  required <- c("case_dir", "swiss", "sea", "chembl", "idmapping")
  missing <- required[vapply(required, function(x) is.null(parsed[[x]]), logical(1))]
  if (length(missing) > 0) {
    stop("Missing required arguments: ", paste(missing, collapse = ", "))
  }

  parsed
}

resolve_input_path <- function(case_dir, path_value, must_work = TRUE) {
  if (is.null(path_value)) {
    return(NULL)
  }
  candidate <- if (grepl("^/", path_value)) path_value else file.path(case_dir, path_value)
  normalizePath(candidate, winslash = "/", mustWork = must_work)
}

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

read_swiss_targets <- function(path) {
  ext <- tolower(tools::file_ext(path))
  swiss_skip <- if (ext %in% c("xlsx", "xls")) 1 else 0
  read_table_auto(path, skip = swiss_skip) %>%
    dplyr::transmute(target = str_to_upper(`Common name`)) %>%
    dplyr::filter(!is.na(target), target != "") %>%
    dplyr::distinct(target) %>%
    pull(target)
}

read_sea_targets <- function(path) {
  read_delim(path, delim = ",", show_col_types = FALSE) %>%
    dplyr::filter(str_detect(`Target ID`, "_HUMAN$")) %>%
    dplyr::transmute(target = str_to_upper(Name)) %>%
    dplyr::filter(!is.na(target), target != "") %>%
    dplyr::distinct(target) %>%
    pull(target)
}

read_ppb2_targets <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    return(character())
  }

  read_table_auto(path) %>%
    dplyr::transmute(target = str_to_upper(Symbol)) %>%
    dplyr::filter(!is.na(target), target != "") %>%
    dplyr::distinct(target) %>%
    pull(target)
}

read_idmapping <- function(path) {
  read_tsv(path, show_col_types = FALSE) %>%
    dplyr::distinct(From, To)
}

read_chembl_targets <- function(path, idmapping_tbl) {
  read_table_auto(path) %>%
    dplyr::filter(Organism == "Homo sapiens", !is.na(Accessions), Accessions != "") %>%
    dplyr::select(Accessions) %>%
    separate_rows(Accessions, sep = "\\|") %>%
    left_join(idmapping_tbl, by = c("Accessions" = "From"), relationship = "many-to-many") %>%
    dplyr::transmute(target = str_to_upper(To)) %>%
    dplyr::filter(!is.na(target), target != "") %>%
    dplyr::distinct(target) %>%
    pull(target)
}

build_go_kegg_plot <- function(go_df, kegg_df, ontology_n, kegg_n, liver_regex) {
  go_selected <- go_df %>%
    dplyr::filter(str_detect(str_to_lower(Description), str_to_lower(liver_regex))) %>%
    dplyr::arrange(ONTOLOGY, p.adjust, desc(Count)) %>%
    dplyr::group_by(ONTOLOGY) %>%
    slice_head(n = ontology_n) %>%
    dplyr::ungroup()

  fallback_go <- go_df %>%
    dplyr::arrange(ONTOLOGY, p.adjust, desc(Count)) %>%
    dplyr::group_by(ONTOLOGY) %>%
    slice_head(n = ontology_n) %>%
    dplyr::ungroup()

  if (nrow(go_selected) == 0) {
    go_selected <- fallback_go
  } else {
    missing_ontologies <- setdiff(unique(fallback_go$ONTOLOGY), unique(go_selected$ONTOLOGY))
    if (length(missing_ontologies) > 0) {
      go_selected <- bind_rows(
        go_selected,
        fallback_go %>% dplyr::filter(ONTOLOGY %in% missing_ontologies)
      )
    }
  }

  kegg_selected <- kegg_df %>%
    dplyr::filter(str_detect(str_to_lower(Description), str_to_lower(liver_regex))) %>%
    dplyr::mutate(ONTOLOGY = "KEGG") %>%
    dplyr::arrange(p.adjust, desc(Count)) %>%
    slice_head(n = kegg_n)

  if (nrow(kegg_selected) == 0) {
    kegg_selected <- kegg_df %>%
      dplyr::mutate(ONTOLOGY = "KEGG") %>%
      dplyr::arrange(p.adjust, desc(Count)) %>%
      slice_head(n = kegg_n)
  }

  bind_rows(
    go_selected %>% dplyr::select(ONTOLOGY, ID, Description, GeneRatio, BgRatio, RichFactor, FoldEnrichment, zScore, pvalue, p.adjust, qvalue, geneID, Count),
    kegg_selected %>% dplyr::select(ONTOLOGY, ID, Description, GeneRatio, BgRatio, RichFactor, FoldEnrichment, zScore, pvalue, p.adjust, qvalue, geneID, Count)
  )
}

plot_panel_c <- function(plot_df, out_path) {
  if (nrow(plot_df) == 0) {
    stop("GO_KEGG_plot dataframe is empty.")
  }

  plot_df <- plot_df %>%
    tidyr::separate(col = GeneRatio, sep = "/", into = c("n1", "n2")) %>%
    dplyr::mutate(GeneRatio = as.numeric(n1) / as.numeric(n2)) %>%
    dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "CC", "MF", "KEGG")), ordered = TRUE)) %>%
    dplyr::arrange(ONTOLOGY, Count) %>%
    dplyr::mutate(Description = str_to_sentence(Description)) %>%
    dplyr::mutate(Description = factor(Description, levels = Description, ordered = TRUE))

  p <- plot_df %>%
    ggplot() +
    geom_point(
      data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
      aes(x = Count, y = interaction(Description, ONTOLOGY, sep = "&"), fill = -log10(p.adjust), size = Count),
      shape = 21
    ) +
    scale_fill_gradient(low = "#a1d99b", high = "#238b45", name = "KEGG", guide = guide_colorbar(order = 4)) +
    ggnewscale::new_scale_fill() +
    geom_point(
      data = plot_df %>% dplyr::filter(ONTOLOGY == "MF"),
      aes(x = Count, y = interaction(Description, ONTOLOGY, sep = "&"), fill = -log10(p.adjust), size = Count),
      shape = 21
    ) +
    scale_fill_gradient(low = "#a6bddb", high = "#0570b0", name = "MF", guide = guide_colorbar(order = 3)) +
    ggnewscale::new_scale_fill() +
    geom_point(
      data = plot_df %>% dplyr::filter(ONTOLOGY == "CC"),
      aes(x = Count, y = interaction(Description, ONTOLOGY, sep = "&"), fill = -log10(p.adjust), size = Count),
      shape = 21
    ) +
    scale_fill_gradient(low = "#fcc5c0", high = "#dd3497", name = "CC", guide = guide_colorbar(order = 2)) +
    ggnewscale::new_scale_fill() +
    geom_point(
      data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
      aes(x = Count, y = interaction(Description, ONTOLOGY, sep = "&"), fill = -log10(p.adjust), size = Count),
      shape = 21
    ) +
    scale_fill_gradient(low = "#8c96c6", high = "#8c6bb1", name = "BP", guide = guide_colorbar(order = 1)) +
    guides(
      y = legendry::guide_axis_nested(
        key = "&",
        type = "bracket",
        levels_text = list(
          element_text(color = c(rep("#41ae76", sum(plot_df$ONTOLOGY == "KEGG")), rep("#225ea8", sum(plot_df$ONTOLOGY == "MF")), rep("#dd3497", sum(plot_df$ONTOLOGY == "CC")), rep("#88419d", sum(plot_df$ONTOLOGY == "BP")))),
          element_text(color = c("#41ae76", "#225ea8", "#dd3497", "#88419d"))
        )
      )
    ) +
    ggtitle(label = "GO and KEGG annotation") +
    labs(x = "Count", y = "Description") +
    scale_size(range = c(3, 7), guide = guide_legend(override.aes = list(fill = "#000000"))) +
    theme_bw() +
    theme_guide(bracket = element_line(c("#74c476", "#41b6c4", "#dd3497", "#9e9ac8"), linewidth = 2)) +
    theme(
      panel.border = element_rect(linewidth = 0.5),
      plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
      axis.text = element_text(color = "#000000", size = 11),
      axis.title = element_text(color = "#000000", size = 15),
      plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
    ) +
    coord_cartesian(clip = "off")

  ggsave(filename = out_path, plot = p, height = 10, width = 10)
}

plot_panel_d <- function(plot_df, out_path) {
  circlize_plot_df <- plot_df %>% dplyr::filter(ONTOLOGY == "KEGG")
  if (nrow(circlize_plot_df) == 0) {
    stop("No KEGG rows available for circos plot.")
  }

  select_id <- circlize_plot_df %>%
    dplyr::select(ID, geneID, FoldEnrichment) %>%
    separate_rows(geneID, sep = "/") %>%
    set_names(c("from", "to", "value")) %>%
    pull(to) %>%
    table() %>%
    as.data.frame() %>%
    set_names(c("ID", "Count")) %>%
    dplyr::filter(Count != 1)

  if (nrow(select_id) == 0) {
    select_id <- circlize_plot_df %>%
      dplyr::select(ID, geneID, FoldEnrichment) %>%
      separate_rows(geneID, sep = "/") %>%
      set_names(c("from", "to", "value")) %>%
      dplyr::distinct(to) %>%
      dplyr::transmute(ID = to, Count = 1L)
  }

  circlize_plot_df2 <- circlize_plot_df %>%
    dplyr::select(ID, geneID, FoldEnrichment) %>%
    separate_rows(geneID, sep = "/") %>%
    set_names(c("from", "to", "value")) %>%
    dplyr::arrange(to, desc(value)) %>%
    dplyr::filter(to %in% select_id$ID) %>%
    dplyr::distinct(to, .keep_all = TRUE)

  sectors <- unique(c(circlize_plot_df2$from, circlize_plot_df2$to))
  colors <- stats::setNames(grDevices::hcl.colors(length(sectors), palette = "Dynamic"), sectors)

  pdf(file = out_path, height = 12, width = 12)
  circos.par(circle.margin = c(0.5, 0.5, 0.5, 0.5))
  chordDiagram(
    circlize_plot_df2,
    grid.col = colors,
    col = "#74c476",
    annotationTrack = "grid"
  )
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[1] + 2,
      CELL_META$sector.index,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.85
    )
  }, bg.border = NA)
  circos.clear()
  dev.off()
}

opt <- parse_args(args)
case_dir <- normalizePath(opt$case_dir, winslash = "/", mustWork = TRUE)
swiss_path <- resolve_input_path(case_dir, opt$swiss, must_work = TRUE)
sea_path <- resolve_input_path(case_dir, opt$sea, must_work = TRUE)
chembl_path <- resolve_input_path(case_dir, opt$chembl, must_work = TRUE)
idmapping_path <- resolve_input_path(case_dir, opt$idmapping, must_work = TRUE)
ppb2_path <- resolve_input_path(case_dir, opt$ppb2, must_work = FALSE)

idmapping_tbl <- read_idmapping(idmapping_path)

venn <- list(
  ChEMBL = read_chembl_targets(chembl_path, idmapping_tbl),
  PPB2 = read_ppb2_targets(ppb2_path),
  Swiss = read_swiss_targets(swiss_path),
  SEA = read_sea_targets(sea_path)
)

saveRDS(venn, file = file.path(case_dir, "venn.rds"))
writeLines(jsonlite::toJSON(venn, auto_unbox = TRUE, pretty = TRUE), con = file.path(case_dir, "venn_inputs.json"))

p_venn <- ggvenn(
  venn,
  fill_color = c("#fc8d62", "#a6d854", "#8da0cb", "#bc80bd"),
  text_size = 7,
  set_name_size = 15
)
ggsave(filename = file.path(case_dir, "1.p_venn.pdf"), plot = p_venn, height = 9, width = 10)

union_id <- unique(unlist(venn))
geneid_trans <- bitr(
  geneID = union_id,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db"
)

go_result <- enrichGO(
  gene = geneid_trans$ENTREZID,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "ALL"
)
go_result_2 <- setReadable(go_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
go_df <- as.data.frame(go_result_2)
write_xlsx(go_df, path = file.path(case_dir, "2.GO_result_2.xlsx"))

kegg_result <- enrichKEGG(
  gene = geneid_trans$ENTREZID,
  organism = "hsa"
)
kegg_result_2 <- setReadable(kegg_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_df <- as.data.frame(kegg_result_2) %>%
  dplyr::mutate(ONTOLOGY = "KEGG")
write_xlsx(kegg_df %>% dplyr::select(-ONTOLOGY), path = file.path(case_dir, "2.KEGG_result_2.xlsx"))

plot_df <- build_go_kegg_plot(go_df, kegg_df %>% dplyr::select(-ONTOLOGY), opt$ontology_n, opt$kegg_n, opt$liver_regex)
write_xlsx(plot_df, path = file.path(case_dir, "GO_KEGG_plot.xlsx"))
plot_panel_c(plot_df, file.path(case_dir, "2.GO_KEGG.pdf"))
plot_panel_d(plot_df, file.path(case_dir, "3.KEGG_circos.pdf"))

key_targets <- plot_df %>%
  dplyr::filter(ONTOLOGY == "KEGG") %>%
  dplyr::select(ID, Description, geneID) %>%
  separate_rows(geneID, sep = "/") %>%
  dplyr::distinct(ID, Description, geneID)
write_csv(key_targets, file.path(case_dir, "key_targets_from_kegg.csv"))

message("Wrote: ", file.path(case_dir, "1.p_venn.pdf"))
message("Wrote: ", file.path(case_dir, "2.GO_result_2.xlsx"))
message("Wrote: ", file.path(case_dir, "2.KEGG_result_2.xlsx"))
message("Wrote: ", file.path(case_dir, "GO_KEGG_plot.xlsx"))
message("Wrote: ", file.path(case_dir, "2.GO_KEGG.pdf"))
message("Wrote: ", file.path(case_dir, "3.KEGG_circos.pdf"))
