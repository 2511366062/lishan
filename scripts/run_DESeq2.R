suppressPackageStartupMessages({
  library(DESeq2)
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

script_dir <- get_script_dir()
project_dir <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)

read_config_value <- function(key, config_file = "config.py") {
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

config_base_dir <- read_config_value("R_BASE_DIR", file.path(project_dir, "config.py"))
if (is.na(config_base_dir)) {
  config_base_dir <- "D:/lxk/project/lishan-20260613"
}

base_dir <- ifelse(length(args) >= 1, args[[1]], config_base_dir)
deg_name <- ifelse(length(args) >= 2, args[[2]], "EVT")
padj_cutoff <- ifelse(length(args) >= 3, as.numeric(args[[3]]), 0.05)
logfc_cutoff <- ifelse(length(args) >= 4, as.numeric(args[[4]]), 0.25)

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  gsub("^_+|_+$", "", x)
}

prefix <- safe_name(deg_name)
setwd(base_dir)

deg_dir <- file.path("DEG", prefix)

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) {
    stop(paste("File not found. Tried:", paste(paths, collapse = "; ")))
  }
  hit[[1]]
}

counts_file <- first_existing(c(
  file.path(deg_dir, "counts.csv"),
  file.path(deg_dir, paste0(prefix, "_counts.csv")),
  file.path(deg_dir, paste0(prefix, "_counts_filtered.csv")),
  file.path(deg_dir, paste0(prefix, "_pseudobulk_counts_filtered.csv"))
))

metadata_file <- first_existing(c(
  file.path(deg_dir, "metadata.csv"),
  file.path(deg_dir, paste0(prefix, "_metadata.csv")),
  file.path(deg_dir, paste0(prefix, "_metadata_filtered.csv")),
  file.path(deg_dir, paste0(prefix, "_pseudobulk_metadata_filtered.csv"))
))

out_all <- file.path(deg_dir, "all.csv")
out_sig <- file.path(deg_dir, "sig.csv")


counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

rownames(metadata) <- metadata$pb_sample
metadata <- metadata[colnames(counts), , drop = FALSE]
metadata$group <- factor(metadata$group, levels = c("NP", "PE"))

print(metadata[, c("pb_sample", "sample", "group", "celltype", "n_cells", "total_counts")])
print(table(metadata$group))

if (any(table(metadata$group) < 2)) {
  stop("Each group should contain at least 2 pseudo-bulk samples after filtering.")
}


# Keep genes with at least 10 counts in at least 2 pseudo-bulk samples.
keep_genes <- rowSums(counts >= 10) >= 2
counts <- counts[keep_genes, ]
print(paste("Genes kept:", nrow(counts)))


dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData = metadata,
  design = ~ group
)

dds <- DESeq(dds)

res_raw <- results(dds, contrast = c("group", "PE", "NP"))

if (requireNamespace("ashr", quietly = TRUE)) {
  res <- lfcShrink(dds, contrast = c("group", "PE", "NP"), res = res_raw, type = "ashr")
} else {
  message("Package 'ashr' is not installed. Saving unshrunk log2FoldChange.")
  res <- res_raw
}

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
for (col in c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")) {
  if (!col %in% colnames(res_df)) {
    res_df[[col]] <- NA
  }
}

res_df <- res_df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

res_df$change <- "NS"
res_df$change[!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange > logfc_cutoff] <- "Up"
res_df$change[!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange < -logfc_cutoff] <- "Down"
res_df$show <- TRUE

res_df <- res_df[order(res_df$padj, res_df$pvalue), ]

write.csv(res_df, out_all, row.names = FALSE)

sig_df <- subset(res_df, !is.na(padj) & padj < padj_cutoff & abs(log2FoldChange) > logfc_cutoff)
write.csv(sig_df, out_sig, row.names = FALSE)

print(paste("padj cutoff:", padj_cutoff))
print(paste("abs log2FC cutoff:", logfc_cutoff))
print(paste("Saved:", out_all))
print(paste("Saved:", out_sig))
print(paste("FDR significant genes:", nrow(sig_df)))
