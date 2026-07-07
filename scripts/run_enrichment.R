suppressPackageStartupMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

get_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 0) {
    return(getwd())
  }
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
}

read_config_value <- function(key, config_file) {
  if (!file.exists(config_file)) {
    return(NA_character_)
  }
  lines <- readLines(config_file, warn = FALSE)
  pattern <- paste0("^", key, "\\s*=\\s*r?[\"'](.+)[\"']")
  hit <- grep(pattern, lines, value = TRUE)
  if (length(hit) == 0) {
    return(NA_character_)
  }
  sub(pattern, "\\1", hit[[1]])
}

prepare_gene_sets <- function(species = "Homo sapiens", sources = c("H", "GO_BP")) {
  gene_sets <- list()

  if ("H" %in% sources || "Hallmark" %in% sources) {
    gene_sets[["Hallmark"]] <- msigdbr(species = species, category = "H") %>%
      dplyr::select(gs_name, gene_symbol)
  }

  if ("GO_BP" %in% sources || "GOBP" %in% sources) {
    gene_sets[["GO_BP"]] <- msigdbr(species = species, category = "C5", subcategory = "GO:BP") %>%
      dplyr::select(gs_name, gene_symbol)
  }

  if ("Reactome" %in% sources || "REACTOME" %in% sources) {
    gene_sets[["Reactome"]] <- msigdbr(species = species, category = "C2", subcategory = "CP:REACTOME") %>%
      dplyr::select(gs_name, gene_symbol)
  }

  if (length(gene_sets) == 0) {
    stop("No valid gene set source selected. Use H, GO_BP, Reactome.")
  }

  dplyr::bind_rows(gene_sets) %>%
    dplyr::distinct(gs_name, gene_symbol)
}

run_gsea <- function(all_df, term2gene, out_file, gsea_padj_cutoff = 0.05, nes_cutoff = 1) {
  rank_df <- all_df %>%
    dplyr::filter(
      !is.na(gene),
      !is.na(log2FoldChange),
      !is.na(pvalue),
      pvalue > 0
    ) %>%
    dplyr::mutate(rank_score = sign(log2FoldChange) * -log10(pvalue)) %>%
    dplyr::filter(is.finite(rank_score)) %>%
    dplyr::arrange(dplyr::desc(abs(rank_score))) %>%
    dplyr::distinct(gene, .keep_all = TRUE)

  gene_list <- rank_df$rank_score
  names(gene_list) <- rank_df$gene
  gene_list <- sort(gene_list, decreasing = TRUE)

  if (length(gene_list) < 50) {
    message("Too few ranked genes for GSEA, skip.")
    return(NULL)
  }

  res <- GSEA(
    geneList = gene_list,
    TERM2GENE = term2gene,
    pvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500,
    eps = 0,
    verbose = FALSE
  )

  res_df <- as.data.frame(res) %>%
    dplyr::filter(!is.na(p.adjust), p.adjust < gsea_padj_cutoff, abs(NES) > nes_cutoff)
  res_df$show <- TRUE
  write.csv(res_df, out_file, row.names = FALSE)
  res_df
}

run_ora <- function(genes, term2gene, out_file, min_genes = 10, ora_padj_cutoff = 0.05) {
  genes <- unique(genes)
  genes <- genes[!is.na(genes) & genes != ""]

  if (length(genes) < min_genes) {
    message("Too few genes for ORA: ", length(genes), " < ", min_genes, ", skip ", basename(out_file))
    return(NULL)
  }

  res <- enricher(
    gene = genes,
    TERM2GENE = term2gene,
    pvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500
  )

  res_df <- as.data.frame(res) %>%
    dplyr::filter(!is.na(p.adjust), p.adjust < ora_padj_cutoff)
  res_df$show <- TRUE
  write.csv(res_df, out_file, row.names = FALSE)
  res_df
}

script_dir <- get_script_dir()
project_dir <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
config_base_dir <- read_config_value("R_BASE_DIR", file.path(project_dir, "config.py"))
if (is.na(config_base_dir)) {
  config_base_dir <- "D:/lxk/project/lishan-20260613"
}

base_dir <- ifelse(length(args) >= 1, args[[1]], config_base_dir)
deg_name <- ifelse(length(args) >= 2, args[[2]], "EVT")
min_ora_genes <- ifelse(length(args) >= 3, as.integer(args[[3]]), 10)
gene_set_sources <- ifelse(length(args) >= 4, args[[4]], "H,GO_BP")
gene_set_sources <- trimws(unlist(strsplit(gene_set_sources, ",")))
gsea_padj_cutoff <- ifelse(length(args) >= 5, as.numeric(args[[5]]), 0.05)
nes_cutoff <- ifelse(length(args) >= 6, as.numeric(args[[6]]), 1)
ora_padj_cutoff <- ifelse(length(args) >= 7, as.numeric(args[[7]]), 0.05)

deg_dir <- file.path(base_dir, "DEG", deg_name)
all_file <- file.path(deg_dir, "all.csv")
sig_file <- file.path(deg_dir, "sig.csv")

all_df <- read.csv(all_file, stringsAsFactors = FALSE)
sig_df <- read.csv(sig_file, stringsAsFactors = FALSE)

term2gene <- prepare_gene_sets(sources = gene_set_sources)

message("Running enrichment for: ", deg_name)
message("Gene sets: ", paste(gene_set_sources, collapse = ", "))
message("GSEA filter: p.adjust < ", gsea_padj_cutoff, ", abs(NES) > ", nes_cutoff)
message("ORA filter: p.adjust < ", ora_padj_cutoff)

run_gsea(
  all_df = all_df,
  term2gene = term2gene,
  out_file = file.path(deg_dir, "GSEA.csv"),
  gsea_padj_cutoff = gsea_padj_cutoff,
  nes_cutoff = nes_cutoff
)

up_genes <- sig_df %>%
  dplyr::filter(change == "Up") %>%
  dplyr::pull(gene)

down_genes <- sig_df %>%
  dplyr::filter(change == "Down") %>%
  dplyr::pull(gene)

run_ora(
  genes = up_genes,
  term2gene = term2gene,
  out_file = file.path(deg_dir, "upORA.csv"),
  min_genes = min_ora_genes,
  ora_padj_cutoff = ora_padj_cutoff
)

run_ora(
  genes = down_genes,
  term2gene = term2gene,
  out_file = file.path(deg_dir, "downORA.csv"),
  min_genes = min_ora_genes,
  ora_padj_cutoff = ora_padj_cutoff
)

message("Done.")
