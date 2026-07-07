options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(limma)
  library(pheatmap)
  library(ggrepel)
  library(patchwork)
  library(cowplot)
  library(grid)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "fig", "RNA")
tab_dir <- file.path(out_dir, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "GSE190971"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "GSE75010"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "integrated"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "scRNA_enrichment"), recursive = TRUE, showWarnings = FALSE)

cols <- list(
  NP = "#8DA6BF",
  PE = "#E7A6A1",
  grey = "#6B7280",
  dark = "#1F2933",
  orange = "#C98B2E",
  teal = "#3B8C88",
  purple = "#8064A2",
  ink = "#243447"
)

theme_nature <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = cols$dark),
      axis.text = element_text(color = cols$dark, size = base_size - 1),
      axis.title = element_text(color = cols$dark, size = base_size),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 2),
      plot.subtitle = element_text(color = "#4B5563", size = base_size - 1),
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 1),
      legend.position = "right",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(color = "#222222", fill = NA, linewidth = 0.35),
      axis.line = element_blank(),
      plot.margin = margin(6, 8, 6, 8)
    )
}

save_pdf <- function(plot, filename, width = 6.2, height = 4.8) {
  path <- file.path(out_dir, filename)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(path, plot, width = width, height = height,
         device = cairo_pdf, units = "in", bg = "white")
}

cap_z <- function(x, lim = 2.5) pmax(pmin(x, lim), -lim)

dedup_matrix <- function(df, gene_col) {
  df <- df %>% filter(!is.na(.data[[gene_col]]), .data[[gene_col]] != "")
  mat <- df %>%
    mutate(.gene = toupper(.data[[gene_col]])) %>%
    select(.gene, where(is.numeric)) %>%
    group_by(.gene) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  m <- as.matrix(mat[, -1, drop = FALSE])
  rownames(m) <- mat$.gene
  m
}

limma_de <- function(expr, group, pe_label = "PE", np_label = "NP") {
  group <- factor(group, levels = c(np_label, pe_label))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- c(np_label, pe_label)
  fit <- lmFit(expr, design)
  cont <- makeContrasts(PE_vs_NP = PE - NP, levels = setNames(colnames(design), colnames(design)))
  fit2 <- eBayes(contrasts.fit(fit, cont))
  tt <- topTable(fit2, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  tt <- tt %>% rename(logFC = logFC, pvalue = P.Value, padj = adj.P.Val)
  tt
}

make_pca_plot <- function(expr, meta, title, subtitle) {
  vars <- apply(expr, 1, var, na.rm = TRUE)
  keep <- names(sort(vars, decreasing = TRUE))[seq_len(min(3000, length(vars)))]
  pc <- prcomp(t(expr[keep, , drop = FALSE]), scale. = TRUE)
  pdat <- as.data.frame(pc$x[, 1:2])
  pdat$sample <- rownames(pdat)
  pdat <- left_join(pdat, meta, by = "sample")
  ve <- (pc$sdev^2 / sum(pc$sdev^2))[1:2] * 100
  ggplot(pdat, aes(PC1, PC2, color = group, fill = group)) +
    stat_ellipse(type = "norm", geom = "polygon", alpha = 0.08, color = NA) +
    geom_point(size = 2.8, shape = 21, stroke = 0.35, color = "white") +
    scale_color_manual(values = c(NP = cols$NP, PE = cols$PE)) +
    scale_fill_manual(values = c(NP = cols$NP, PE = cols$PE)) +
    labs(title = title, subtitle = subtitle,
         x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2])) +
    theme_nature(9)
}

make_volcano <- function(de, genes, title, subtitle) {
  d <- de %>%
    mutate(sig = case_when(
      padj < 0.05 & logFC > 0 ~ "PE higher",
      padj < 0.05 & logFC < 0 ~ "NP higher",
      TRUE ~ "NS"
    ),
    label = ifelse(gene %in% genes | (padj < 0.01 & abs(logFC) > quantile(abs(logFC), 0.995, na.rm = TRUE)), gene, NA))
  ggplot(d, aes(logFC, -log10(pvalue))) +
    geom_point(aes(color = sig), size = 0.65, alpha = 0.65) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.25, color = "#9CA3AF") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.25, color = "#9CA3AF") +
    geom_text_repel(aes(label = label), size = 2.2, max.overlaps = 30, min.segment.length = 0.05,
                    box.padding = 0.25, segment.size = 0.2, na.rm = TRUE) +
    scale_color_manual(values = c("PE higher" = cols$PE, "NP higher" = cols$NP, "NS" = "#CBD5E1")) +
    labs(title = title, subtitle = subtitle, x = "log2 fold-change (PE / NP)", y = "-log10(P)") +
    theme_nature(9)
}

make_design_bar <- function(meta, title, subtitle) {
  meta %>%
    count(group, name = "n") %>%
    ggplot(aes(x = group, y = n, fill = group)) +
    geom_col(width = 0.6, color = "white", linewidth = 0.3) +
    geom_text(aes(label = n), vjust = -0.35, size = 3) +
    scale_fill_manual(values = c(NP = cols$NP, PE = cols$PE)) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "Samples") +
    coord_cartesian(ylim = c(0, max(table(meta$group)) * 1.18)) +
    theme_nature(9) +
    theme(legend.position = "none")
}

make_effect_dot <- function(de, gene_set, title, subtitle) {
  d <- de %>%
    filter(gene %in% gene_set) %>%
    mutate(direction = ifelse(logFC >= 0, "PE higher", "NP higher"),
           gene = factor(gene, levels = gene[order(logFC)]))
  ggplot(d, aes(logFC, gene, color = direction, size = -log10(pvalue))) +
    geom_vline(xintercept = 0, color = "#9CA3AF", linewidth = 0.35) +
    geom_segment(aes(x = 0, xend = logFC, yend = gene), color = "#CBD5E1", linewidth = 0.55) +
    geom_point(alpha = 0.95) +
    scale_color_manual(values = c("PE higher" = cols$PE, "NP higher" = cols$NP)) +
    scale_size_continuous(range = c(1.7, 5.5)) +
    labs(title = title, subtitle = subtitle, x = "log2 fold-change (PE / NP)", y = NULL, size = "-log10(P)") +
    theme_nature(9)
}

make_gene_box <- function(expr, meta, genes, title, subtitle, ncol = 4) {
  genes <- intersect(toupper(genes), rownames(expr))
  dat <- as.data.frame(t(expr[genes, , drop = FALSE]))
  dat$sample <- rownames(dat)
  dat <- dat %>%
    pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression") %>%
    left_join(meta, by = "sample")
  pdat <- dat %>%
    group_by(gene) %>%
    summarise(
      pvalue = tryCatch(wilcox.test(expression[group == "PE"], expression[group == "NP"])$p.value, error = function(e) NA_real_),
      y = max(expression, na.rm = TRUE) + 0.08 * diff(range(expression, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(label = case_when(
      is.na(pvalue) ~ "p = NA",
      pvalue < 0.001 ~ "p < 0.001",
      TRUE ~ sprintf("p = %.3f", pvalue)
    ))
  ggplot(dat, aes(group, expression, fill = group)) +
    geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.35, alpha = 0.72) +
    geom_jitter(aes(color = group), width = 0.14, size = 0.9, alpha = 0.72, show.legend = FALSE) +
    geom_text(data = pdat, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, size = 2.4) +
    facet_wrap(~ gene, scales = "free_y", ncol = ncol) +
    scale_fill_manual(values = c(NP = cols$NP, PE = cols$PE)) +
    scale_color_manual(values = c(NP = cols$NP, PE = cols$PE)) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "Expression") +
    coord_cartesian(clip = "off") +
    theme_nature(8) +
    theme(legend.position = "none")
}

make_signature_scores <- function(expr, meta, gene_sets) {
  scores <- lapply(names(gene_sets), function(nm) {
    genes <- intersect(toupper(gene_sets[[nm]]), rownames(expr))
    if (length(genes) == 0) return(NULL)
    z <- t(scale(t(expr[genes, , drop = FALSE])))
    data.frame(sample = colnames(expr), signature = nm, score = colMeans(z, na.rm = TRUE))
  }) %>% bind_rows()
  left_join(scores, meta, by = "sample")
}

make_signature_violin <- function(scores, title, subtitle) {
  ggplot(scores, aes(signature, score, fill = group, color = group)) +
    geom_hline(yintercept = 0, color = "#9CA3AF", linewidth = 0.25) +
    geom_violin(width = 0.75, alpha = 0.25, linewidth = 0.25, position = position_dodge(width = 0.74)) +
    geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.72, linewidth = 0.25, position = position_dodge(width = 0.74)) +
    scale_fill_manual(values = c(NP = cols$NP, PE = cols$PE)) +
    scale_color_manual(values = c(NP = cols$NP, PE = cols$PE)) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "Mean z-score") +
    theme_nature(8) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

write_heatmap <- function(expr, meta, genes, filename, title, width = 7.4, height = 5.2, cluster_rows = FALSE) {
  genes <- intersect(toupper(genes), rownames(expr))
  if (length(genes) < 2) return(invisible(NULL))
  ann <- meta %>% select(sample, group) %>% distinct() %>% arrange(group, sample) %>% as.data.frame()
  rownames(ann) <- ann$sample
  ann$sample <- NULL
  z <- t(scale(t(expr[genes, rownames(ann), drop = FALSE])))
  z[is.na(z)] <- 0
  z <- cap_z(z)
  row_effect <- rowMeans(expr[genes, rownames(ann)[ann$group == "PE"], drop = FALSE], na.rm = TRUE) -
    rowMeans(expr[genes, rownames(ann)[ann$group == "NP"], drop = FALSE], na.rm = TRUE)
  z <- z[names(sort(row_effect, decreasing = TRUE)), , drop = FALSE]
  ann_cols <- list(group = c(NP = cols$NP, PE = cols$PE))
  path <- file.path(out_dir, filename)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  pheatmap(z, filename = path, width = width, height = height,
           color = colorRampPalette(c(cols$NP, "white", cols$PE))(101),
           annotation_col = ann, annotation_colors = ann_cols,
           show_colnames = FALSE, border_color = NA, main = title,
           fontsize = 8, fontsize_row = 7, cluster_rows = cluster_rows, cluster_cols = FALSE)
}

write_corr_heatmap <- function(expr, meta, filename, title, width = 6.5, height = 5.6) {
  vars <- apply(expr, 1, var, na.rm = TRUE)
  keep <- names(sort(vars, decreasing = TRUE))[seq_len(min(2000, length(vars)))]
  cm <- cor(expr[keep, , drop = FALSE], method = "pearson", use = "pairwise.complete.obs")
  ann <- meta %>% select(sample, group) %>% distinct() %>% arrange(group, sample) %>% as.data.frame()
  rownames(ann) <- ann$sample
  ann$sample <- NULL
  ann_cols <- list(group = c(NP = cols$NP, PE = cols$PE))
  cm <- cm[rownames(ann), rownames(ann), drop = FALSE]
  path <- file.path(out_dir, filename)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  pheatmap(cm, filename = path, width = width, height = height,
           color = colorRampPalette(c("#F8FBFD", "#C9D8E6", "#657A91"))(101),
           annotation_col = ann, annotation_row = ann, annotation_colors = ann_cols,
           show_colnames = FALSE, show_rownames = FALSE, border_color = NA,
           main = title, fontsize = 8, cluster_rows = FALSE, cluster_cols = FALSE)
}

parse_sample_id_map <- function(filelist) {
  lines <- readLines(filelist, warn = FALSE)
  rows <- lapply(lines[grepl("^File\\t", lines)], function(x) {
    name <- strsplit(x, "\t", fixed = TRUE)[[1]][2]
    m <- regexec("^(GSM[0-9]+)_([0-9]+)_", name)
    hit <- regmatches(name, m)[[1]]
    if (length(hit) == 0) return(NULL)
    data.frame(geo_accession = hit[2], sample = paste0("SAMPLE_", hit[3]))
  })
  bind_rows(rows)
}

# Candidate genes from the ATAC/scRNA integration, plus PE placenta markers.
medium_path <- file.path(root, "fig", "atac_fig", "ATAC_linked_gene_sets", "medium_absATAC0.5", "medium_absATAC0.5_linked_gene_table.tsv")
strict_path <- file.path(root, "fig", "atac_fig", "ATAC_linked_gene_sets", "strict", "strict_linked_gene_table.tsv")
medium <- read_tsv(medium_path, show_col_types = FALSE)
strict <- read_tsv(strict_path, show_col_types = FALSE)
linked_genes <- unique(toupper(medium$gene))
strict_genes <- unique(toupper(strict$gene))
manual_genes <- toupper(c("FLT1", "ENG", "PGF", "VEGFA", "LEP", "HTRA4", "ADAM12", "PAPPA2", "INHBA", "SERPINE1",
                          "KRT8", "LYVE1", "PLA2G2A", "CXCL13", "TNFSF10", "LPL", "HPGD", "IL18BP", "DAPP1"))
plot_genes <- unique(c(strict_genes, linked_genes, manual_genes))

gene_sets <- list(
  `Angiogenic imbalance` = c("FLT1", "ENG", "PGF", "VEGFA", "KDR", "TEK", "ANGPT2"),
  `Inflammatory chemokine` = c("CXCL8", "CXCL10", "CXCL13", "CCL2", "TNF", "IL6", "IL1B", "TNFSF10", "IL18BP"),
  `Trophoblast stress` = c("HTRA4", "ADAM12", "PAPPA2", "INHBA", "LEP", "KRT8", "GCM1", "CGB5"),
  `Vascular lymphatic` = c("LYVE1", "PECAM1", "VWF", "KDR", "FLT1", "TEK"),
  `ATAC-linked strict` = strict_genes,
  `ATAC-linked medium` = linked_genes
)

# GSE190971 placenta RNA-seq.
g190 <- read_tsv(file.path(root, "external_data", "GSE190971", "GSE190971_Normalised_gene_counts_matrix_PLAC.txt.gz"),
                 show_col_types = FALSE)
expr190 <- dedup_matrix(g190, "Gene_Symbol")
expr190 <- log2(expr190 + 1)
meta190 <- data.frame(sample = colnames(expr190)) %>%
  mutate(group = ifelse(grepl("_PE$", sample), "PE", "NP"),
         subtype = group,
         dataset = "GSE190971")
de190 <- limma_de(expr190, meta190$group, pe_label = "PE", np_label = "NP")
write_tsv(de190, file.path(tab_dir, "GSE190971_limma_PE_vs_NP.tsv"))
write_tsv(meta190, file.path(tab_dir, "GSE190971_samples_used.tsv"))

# GSE75010 large placenta cohort.
g750 <- read_csv(file.path(root, "external_data", "GSE75010", "GSE75010_complete_dataset.csv.gz"),
                 show_col_types = FALSE)
names(g750)[1] <- "Gene_Symbol"
expr750_all <- dedup_matrix(g750, "Gene_Symbol")
map750 <- parse_sample_id_map(file.path(root, "external_data", "GSE75010", "filelist.txt"))
meta750_raw <- read_tsv(file.path(root, "external_data", "GSE75010", "GSE75010_sample_metadata.tsv"), show_col_types = FALSE)
meta750 <- map750 %>%
  left_join(meta750_raw, by = "geo_accession") %>%
  filter(sample %in% colnames(expr750_all), char_diagnosis %in% c("PE", "non-PE")) %>%
  transmute(sample, geo_accession, group = ifelse(char_diagnosis == "PE", "PE", "NP"),
            subtype = title, gestational_week = suppressWarnings(as.numeric(char_ga_week)),
            diagnosis = char_diagnosis, dataset = "GSE75010")
expr750 <- expr750_all[, meta750$sample, drop = FALSE]
de750 <- limma_de(expr750, meta750$group, pe_label = "PE", np_label = "NP")
write_tsv(de750, file.path(tab_dir, "GSE75010_limma_PE_vs_NP.tsv"))
write_tsv(meta750, file.path(tab_dir, "GSE75010_samples_used.tsv"))

# Cross-dataset summaries.
combined_effect <- full_join(
  de190 %>% select(gene, logFC_190 = logFC, p_190 = pvalue, padj_190 = padj),
  de750 %>% select(gene, logFC_750 = logFC, p_750 = pvalue, padj_750 = padj),
  by = "gene"
) %>%
  mutate(in_linked = gene %in% linked_genes, in_strict = gene %in% strict_genes)
write_tsv(combined_effect, file.path(tab_dir, "public_RNA_cross_dataset_effects.tsv"))

linked_summary <- medium %>%
  mutate(gene = toupper(gene)) %>%
  group_by(gene) %>%
  summarise(atac_log2fc = mean(atac_log2fc, na.rm = TRUE),
            scrna_log2fc = mean(scrna_log2fc, na.rm = TRUE),
            candidate = paste(unique(candidate), collapse = ";"),
            celltype = paste(unique(celltype), collapse = ";"),
            .groups = "drop") %>%
  left_join(combined_effect, by = "gene")
write_tsv(linked_summary, file.path(tab_dir, "ATAC_scRNA_public_RNA_linked_gene_validation.tsv"))

# Figures.
save_pdf(make_design_bar(meta190, "GSE190971 placenta RNA-seq design", "Placenta samples used for external PE/NP validation"),
         "GSE190971/01_design_sample_counts.pdf", 3.2, 3.2)
save_pdf(make_pca_plot(expr190, meta190, "GSE190971 placenta transcriptome PCA", "External RNA-seq: PE separates from NP across the leading expression axes"),
         "GSE190971/02_PCA_PE_NP.pdf", 5.0, 4.2)
write_corr_heatmap(expr190, meta190, "GSE190971/03_sample_correlation_heatmap_no_tree.pdf",
                   "GSE190971 sample correlation")
save_pdf(make_volcano(de190, plot_genes, "GSE190971 PE-vs-NP volcano", "Labels highlight ATAC/scRNA-linked and canonical PE genes"),
         "GSE190971/04_volcano_linked_genes.pdf", 6.2, 4.8)
top190 <- de190 %>% arrange(padj, desc(abs(logFC))) %>% slice_head(n = 55) %>% pull(gene)
write_heatmap(expr190, meta190, top190, "GSE190971/05_top_DE_gene_heatmap_no_tree.pdf",
              "GSE190971 top PE-associated genes", width = 6.8, height = 6.8)
write_heatmap(expr190, meta190, plot_genes, "GSE190971/06_ATAC_scRNA_linked_gene_heatmap_no_tree.pdf",
              "GSE190971 ATAC/scRNA-linked gene program", width = 7.2, height = 5.8)
save_pdf(make_gene_box(expr190, meta190, strict_genes, "Strict ATAC/scRNA genes in GSE190971", "Direct expression validation of the strict linked genes", ncol = 3),
         "GSE190971/07_strict_linked_gene_boxplots_with_pvalue.pdf", 6.4, 4.6)
save_pdf(make_gene_box(expr190, meta190, c("FLT1", "ENG", "PGF", "HTRA4", "ADAM12", "LEP", "PAPPA2", "INHBA"),
                       "Canonical PE and trophoblast-stress genes", "GSE190971 placenta RNA-seq", ncol = 4),
         "GSE190971/08_PE_marker_gene_boxplots_with_pvalue.pdf", 7.2, 4.8)
scores190 <- make_signature_scores(expr190, meta190, gene_sets)
write_tsv(scores190, file.path(tab_dir, "GSE190971_signature_scores.tsv"))
save_pdf(make_signature_violin(scores190, "GSE190971 PE pathway signatures", "Module scores calculated from linked and canonical PE gene sets"),
         "GSE190971/09_signature_violin.pdf", 6.8, 4.5)
save_pdf(make_effect_dot(de190, intersect(plot_genes, rownames(expr190)),
                         "GSE190971 linked-gene effect sizes", "Expression fold-changes for ATAC/scRNA-linked genes"),
         "GSE190971/10_linked_gene_effect_dotplot.pdf", 5.6, 5.8)

save_pdf(make_design_bar(meta750, "GSE75010 placenta cohort design", "Mapped GEO samples: 80 PE and 77 non-PE placentas"),
         "GSE75010/01_design_sample_counts.pdf", 3.2, 3.2)
save_pdf(make_pca_plot(expr750, meta750, "GSE75010 placenta transcriptome PCA", "Large external cohort using mapped processed expression values"),
         "GSE75010/02_PCA_PE_NP.pdf", 5.0, 4.2)
write_corr_heatmap(expr750, meta750, "GSE75010/03_sample_correlation_heatmap_no_tree.pdf",
                   "GSE75010 sample correlation", width = 6.8, height = 6.0)
save_pdf(make_volcano(de750, plot_genes, "GSE75010 PE-vs-non-PE volcano", "Large-cohort validation of linked and canonical PE genes"),
         "GSE75010/04_volcano_linked_genes.pdf", 6.2, 4.8)
top750 <- de750 %>% arrange(padj, desc(abs(logFC))) %>% slice_head(n = 60) %>% pull(gene)
write_heatmap(expr750, meta750, top750, "GSE75010/05_top_DE_gene_heatmap_no_tree.pdf",
              "GSE75010 top PE-associated genes", width = 7.0, height = 7.2)
write_heatmap(expr750, meta750, plot_genes, "GSE75010/06_ATAC_scRNA_linked_gene_heatmap_no_tree.pdf",
              "GSE75010 ATAC/scRNA-linked gene program", width = 7.2, height = 5.8)
save_pdf(make_gene_box(expr750, meta750, strict_genes, "Strict ATAC/scRNA genes in GSE75010", "Large-cohort expression validation", ncol = 3),
         "GSE75010/07_strict_linked_gene_boxplots_with_pvalue.pdf", 6.4, 4.6)
save_pdf(make_gene_box(expr750, meta750, c("FLT1", "ENG", "PGF", "HTRA4", "ADAM12", "LEP", "PAPPA2", "INHBA"),
                       "Canonical PE and trophoblast-stress genes", "GSE75010 placenta cohort", ncol = 4),
         "GSE75010/08_PE_marker_gene_boxplots_with_pvalue.pdf", 7.2, 4.8)
scores750 <- make_signature_scores(expr750, meta750, gene_sets)
write_tsv(scores750, file.path(tab_dir, "GSE75010_signature_scores.tsv"))
save_pdf(make_signature_violin(scores750, "GSE75010 PE pathway signatures", "Module scores in a large placenta PE cohort"),
         "GSE75010/09_signature_violin.pdf", 6.8, 4.5)
save_pdf(make_effect_dot(de750, intersect(plot_genes, rownames(expr750)),
                         "GSE75010 linked-gene effect sizes", "Expression fold-changes for ATAC/scRNA-linked genes"),
         "GSE75010/10_linked_gene_effect_dotplot.pdf", 5.6, 5.8)

rep <- linked_summary %>%
  filter(!is.na(logFC_190), !is.na(logFC_750)) %>%
  mutate(strict = ifelse(gene %in% strict_genes, "Strict", "Medium"),
         replicated = sign(logFC_190) == sign(logFC_750),
         label = ifelse(gene %in% strict_genes | abs(logFC_190) + abs(logFC_750) > 1.2, gene, NA))
save_pdf(
  ggplot(rep, aes(logFC_190, logFC_750, color = strict, size = -log10(pmin(p_190, p_750, na.rm = TRUE)))) +
    geom_hline(yintercept = 0, color = "#9CA3AF", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "#9CA3AF", linewidth = 0.3) +
    geom_point(alpha = 0.86) +
    geom_text_repel(aes(label = label), size = 2.4, max.overlaps = 30, na.rm = TRUE) +
    scale_color_manual(values = c(Strict = cols$purple, Medium = cols$teal)) +
    scale_size_continuous(range = c(2, 6)) +
    labs(title = "Cross-cohort replication of ATAC/scRNA-linked genes",
         subtitle = "External RNA effect sizes in GSE190971 and GSE75010",
         x = "GSE190971 log2FC (PE / NP)", y = "GSE75010 log2FC (PE / non-PE)", size = "Best -log10(P)") +
    theme_nature(9),
  "integrated/01_cross_dataset_linked_gene_replication_scatter.pdf", 5.6, 5.2
)

multi <- linked_summary %>%
  select(gene, atac_log2fc, scrna_log2fc, logFC_190, logFC_750) %>%
  filter(gene %in% unique(c(strict_genes, head(linked_genes, 18)))) %>%
  pivot_longer(-gene, names_to = "evidence", values_to = "effect") %>%
  mutate(evidence = recode(evidence, atac_log2fc = "ATAC", scrna_log2fc = "scRNA", logFC_190 = "GSE190971 RNA", logFC_750 = "GSE75010 RNA"),
         gene = factor(gene, levels = rev(unique(gene))))
save_pdf(
  ggplot(multi, aes(evidence, gene, fill = effect)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(is.na(effect), "", sprintf("%.2f", effect))), size = 2.2) +
    scale_fill_gradient2(low = cols$NP, mid = "white", high = cols$PE, midpoint = 0, na.value = "#F3F4F6") +
    labs(title = "Integrated regulatory-to-transcript validation matrix",
         subtitle = "ATAC, scRNA and two public placenta RNA cohorts shown on the same genes",
         x = NULL, y = NULL, fill = "Effect") +
    theme_nature(8) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1)),
  "integrated/02_ATAC_scRNA_public_RNA_effect_matrix.pdf", 6.7, 5.6
)

sub_scores <- scores750 %>%
  filter(signature %in% c("ATAC-linked strict", "Trophoblast stress", "Angiogenic imbalance")) %>%
  mutate(subtype_clean = gsub(", rep[0-9]+", "", subtype),
         subtype_clean = gsub("^Cont", "NP", subtype_clean),
         subtype_clean = gsub("^PE", "PE", subtype_clean))
save_pdf(
  ggplot(sub_scores, aes(reorder(subtype_clean, score, median, na.rm = TRUE), score, fill = group)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.25, alpha = 0.75) +
    geom_jitter(width = 0.14, size = 0.55, alpha = 0.5) +
    facet_wrap(~ signature, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c(NP = cols$NP, PE = cols$PE)) +
    labs(title = "GSE75010 clinical subtype structure of PE signatures",
         subtitle = "Linked and canonical programs across term/preterm and AGA/SGA categories",
         x = NULL, y = "Signature score") +
    theme_nature(7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top"),
  "GSE75010/11_subtype_signature_boxplots.pdf", 7.6, 7.2
)

ga_sig <- scores750 %>%
  filter(signature == "ATAC-linked medium") %>%
  filter(!is.na(gestational_week))
save_pdf(
  ggplot(ga_sig, aes(gestational_week, score, color = group, fill = group)) +
    geom_point(size = 1.7, alpha = 0.78) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.6, alpha = 0.15) +
    scale_color_manual(values = c(NP = cols$NP, PE = cols$PE)) +
    scale_fill_manual(values = c(NP = cols$NP, PE = cols$PE)) +
    labs(title = "ATAC-linked RNA signature versus gestational age",
         subtitle = "GSE75010 shows whether linked PE programs persist across gestational age",
         x = "Gestational age at delivery (weeks)", y = "ATAC-linked medium score") +
    theme_nature(9),
  "GSE75010/12_ATAC_linked_signature_by_gestational_age.pdf", 5.7, 4.5
)

# A compact multi-panel summary for article assembly.
pA <- make_pca_plot(expr190, meta190, "GSE190971", "Placenta RNA-seq")
pB <- make_volcano(de190, strict_genes, "GSE190971", "Strict linked genes labeled")
pC <- make_pca_plot(expr750, meta750, "GSE75010", "Placenta cohort")
pD <- make_volcano(de750, strict_genes, "GSE75010", "Strict linked genes labeled")
summary_panel <- (pA + pB) / (pC + pD) + plot_annotation(title = "Public placenta RNA validation of PE-associated regulatory genes")
save_pdf(summary_panel, "integrated/03_public_RNA_validation_summary_panel.pdf", 10.5, 8.0)

# scRNA-derived enrichment views and external RNA validation of scRNA pathways.
scrna_gsea_path <- file.path(root, "fig", "GSEA", "GSEA_dotplot_data.csv")
if (file.exists(scrna_gsea_path)) {
  scrna_gsea <- read_csv(scrna_gsea_path, show_col_types = FALSE) %>%
    mutate(
      pathway_label = ifelse(is.na(pathway_label) | pathway_label == "", Description, pathway_label),
      direction = ifelse(NES >= 0, "PE enriched", "NP enriched"),
      neg_log10_padj = -log10(p.adjust + 1e-300),
      pathway_short = gsub("Hallmark ", "", pathway_label)
    )
  write_tsv(scrna_gsea, file.path(tab_dir, "scRNA_GSEA_terms_used_for_public_RNA_validation.tsv"))

  top_gsea <- scrna_gsea %>%
    arrange(p.adjust) %>%
    group_by(celltype) %>%
    slice_head(n = 3) %>%
    ungroup() %>%
    arrange(p.adjust) %>%
    mutate(pathway_short = factor(pathway_short, levels = rev(unique(pathway_short))),
           celltype = factor(celltype, levels = unique(celltype)))
  save_pdf(
    ggplot(top_gsea, aes(celltype, pathway_short, size = neg_log10_padj, fill = NES)) +
      geom_point(shape = 21, color = "#2F3A45", linewidth = 0.22, alpha = 0.95) +
      scale_fill_gradient2(low = cols$NP, mid = "white", high = cols$PE, midpoint = 0) +
      scale_size_continuous(range = c(1.8, 5.2)) +
      labs(title = "scRNA-derived GSEA programs linked to PE",
           subtitle = "Top enriched pathways per cell type from the scRNA analysis",
           x = NULL, y = NULL, size = "-log10(adj. P)", fill = "NES") +
      theme_nature(8) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1),
            panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.25),
            panel.border = element_rect(color = "#222222", fill = NA, linewidth = 0.35)),
    "scRNA_enrichment/01_scRNA_GSEA_top_pathway_dotplot.pdf", 9.5, 6.8
  )

  gsea_mat_dat <- scrna_gsea %>%
    filter(show %in% c(TRUE, "TRUE", "True") | p.adjust < 0.05) %>%
    arrange(p.adjust) %>%
    distinct(celltype, pathway_short, .keep_all = TRUE) %>%
    group_by(pathway_short) %>%
    mutate(best = min(p.adjust, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(best) %>%
    filter(pathway_short %in% unique(pathway_short)[1:min(22, n_distinct(pathway_short))])
  if (nrow(gsea_mat_dat) > 0) {
    nes_mat <- gsea_mat_dat %>%
      select(pathway_short, celltype, NES) %>%
      pivot_wider(names_from = celltype, values_from = NES, values_fill = 0) %>%
      as.data.frame()
    rownames(nes_mat) <- nes_mat$pathway_short
    nes_mat$pathway_short <- NULL
    pheatmap(as.matrix(nes_mat),
             filename = file.path(out_dir, "scRNA_enrichment/02_scRNA_GSEA_NES_heatmap_no_tree.pdf"),
             width = 8.0, height = 6.4,
             color = colorRampPalette(c(cols$NP, "white", cols$PE))(101),
             border_color = NA, cluster_rows = FALSE, cluster_cols = FALSE,
             fontsize = 7, fontsize_row = 6.5, main = "scRNA GSEA NES matrix")
  }

  linked_upper <- toupper(linked_genes)
  overlap_dat <- scrna_gsea %>%
    mutate(core_list = strsplit(core_enrichment, "/", fixed = TRUE),
           overlap_genes = lapply(core_list, function(x) intersect(toupper(x), linked_upper)),
           overlap_n = lengths(overlap_genes),
           overlap_label = vapply(overlap_genes, function(x) paste(head(x, 6), collapse = ", "), character(1))) %>%
    filter(overlap_n > 0) %>%
    arrange(p.adjust) %>%
    slice_head(n = 30) %>%
    mutate(pathway_short = factor(pathway_short, levels = rev(unique(pathway_short))))
  if (nrow(overlap_dat) > 0) {
    write_tsv(overlap_dat %>% select(celltype, Description, NES, p.adjust, overlap_n, overlap_label),
              file.path(tab_dir, "scRNA_GSEA_overlap_with_ATAC_linked_genes.tsv"))
    save_pdf(
      ggplot(overlap_dat, aes(overlap_n, pathway_short, color = direction, size = -log10(p.adjust + 1e-300))) +
        geom_point(alpha = 0.9) +
        facet_wrap(~ celltype, scales = "free_y", ncol = 3) +
        scale_color_manual(values = c("PE enriched" = cols$PE, "NP enriched" = cols$NP)) +
        scale_size_continuous(range = c(2, 6)) +
        labs(title = "scRNA pathways carrying ATAC-linked genes",
             subtitle = "Overlap between scRNA GSEA core genes and medium-stringency ATAC-linked genes",
             x = "ATAC-linked genes in core enrichment", y = NULL, size = "-log10(adj. P)") +
        theme_nature(7),
      "scRNA_enrichment/03_scRNA_GSEA_ATAC_linked_overlap_dotplot.pdf", 8.4, 6.5
    )
  }

  selected_pathways <- scrna_gsea %>%
    filter(grepl("PLACENTA|EPITHELIAL_MESENCHYMAL|TNFA|INTERFERON|ANGIO|HYPOXIA|INFLAMMATORY|T_CELL|OXIDATIVE", Description, ignore.case = TRUE)) %>%
    arrange(p.adjust) %>%
    distinct(Description, .keep_all = TRUE) %>%
    slice_head(n = 8)
  if (nrow(selected_pathways) < 5) {
    selected_pathways <- scrna_gsea %>% arrange(p.adjust) %>% distinct(Description, .keep_all = TRUE) %>% slice_head(n = 8)
  }
  scrna_sets <- selected_pathways$core_enrichment
  names(scrna_sets) <- selected_pathways$pathway_short
  scrna_sets <- lapply(scrna_sets, function(x) unique(toupper(strsplit(x, "/", fixed = TRUE)[[1]])))
  scores190_sc <- make_signature_scores(expr190, meta190, scrna_sets)
  scores750_sc <- make_signature_scores(expr750, meta750, scrna_sets)
  write_tsv(scores190_sc, file.path(tab_dir, "GSE190971_scRNA_GSEA_core_signature_scores.tsv"))
  write_tsv(scores750_sc, file.path(tab_dir, "GSE75010_scRNA_GSEA_core_signature_scores.tsv"))
  save_pdf(
    make_signature_violin(scores190_sc, "GSE190971 validation of scRNA GSEA core programs",
                          "Scores use core-enrichment genes from scRNA PE/NP pathways"),
    "scRNA_enrichment/04_GSE190971_scRNA_GSEA_core_signature_violin.pdf", 7.4, 4.8
  )
  save_pdf(
    make_signature_violin(scores750_sc, "GSE75010 validation of scRNA GSEA core programs",
                          "Large-cohort expression scores for scRNA-derived PE/NP pathways"),
    "scRNA_enrichment/05_GSE75010_scRNA_GSEA_core_signature_violin.pdf", 7.4, 4.8
  )
  combined_sc_scores <- bind_rows(
    scores190_sc %>% mutate(dataset = "GSE190971"),
    scores750_sc %>% mutate(dataset = "GSE75010")
  )
  save_pdf(
    ggplot(combined_sc_scores, aes(group, score, fill = group, color = group)) +
      geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.3, alpha = 0.75) +
      geom_jitter(width = 0.14, size = 0.55, alpha = 0.45) +
      facet_grid(dataset ~ signature, scales = "free_y") +
      scale_fill_manual(values = c(NP = cols$NP, PE = cols$PE)) +
      scale_color_manual(values = c(NP = cols$NP, PE = cols$PE)) +
      labs(title = "Public RNA validation of scRNA-enriched pathway programs",
           subtitle = "Core GSEA genes from scRNA are scored in each public placenta RNA dataset",
           x = NULL, y = "Mean z-score") +
      theme_nature(7) +
      theme(axis.text.x = element_text(angle = 0), legend.position = "top"),
    "scRNA_enrichment/06_public_RNA_scRNA_pathway_score_grid.pdf", 12.0, 5.4
  )
}

endo_ora_path <- file.path(root, "DEG", "Endo", "downORA.csv")
if (file.exists(endo_ora_path)) {
  endo_ora <- read_csv(endo_ora_path, show_col_types = FALSE) %>%
    filter(!is.na(p.adjust), p.adjust < 0.1) %>%
    arrange(p.adjust) %>%
    slice_head(n = 18) %>%
    mutate(label = gsub("^GOBP_|^HALLMARK_", "", Description),
           label = gsub("_", " ", tools::toTitleCase(tolower(label))),
           label = factor(label, levels = rev(label)))
  if (nrow(endo_ora) > 0) {
    save_pdf(
      ggplot(endo_ora, aes(FoldEnrichment, label, size = Count, color = -log10(p.adjust + 1e-300))) +
        geom_point(alpha = 0.92) +
        scale_color_gradient(low = cols$NP, high = cols$PE) +
        scale_size_continuous(range = c(2, 7)) +
        labs(title = "scRNA Endo downregulated-gene ORA",
             subtitle = "Placenta-development and EMT-related terms connect external RNA validation to scRNA biology",
             x = "Fold enrichment", y = NULL, color = "-log10(adj. P)") +
        theme_nature(8),
      "scRNA_enrichment/07_Endo_downORA_placenta_EMT_dotplot.pdf", 6.8, 5.6
    )
  }
}

if (exists("rep") && nrow(rep) > 0 && exists("scrna_gsea")) {
  cross_path <- scrna_gsea %>%
    mutate(core_list = strsplit(core_enrichment, "/", fixed = TRUE)) %>%
    rowwise() %>%
    mutate(public_mean_effect = mean(combined_effect$logFC_750[match(intersect(toupper(core_list), combined_effect$gene), combined_effect$gene)], na.rm = TRUE),
           public_overlap = length(intersect(toupper(core_list), combined_effect$gene))) %>%
    ungroup() %>%
    filter(public_overlap >= 5, is.finite(public_mean_effect)) %>%
    arrange(p.adjust) %>%
    slice_head(n = 40)
  if (nrow(cross_path) > 0) {
    save_pdf(
      ggplot(cross_path, aes(NES, public_mean_effect, color = direction, size = public_overlap)) +
        geom_hline(yintercept = 0, color = "#9CA3AF", linewidth = 0.3) +
        geom_vline(xintercept = 0, color = "#9CA3AF", linewidth = 0.3) +
        geom_point(alpha = 0.88) +
        geom_text_repel(aes(label = pathway_short), size = 2.1, max.overlaps = 18) +
        scale_color_manual(values = c("PE enriched" = cols$PE, "NP enriched" = cols$NP)) +
        scale_size_continuous(range = c(2, 6)) +
        labs(title = "scRNA pathway NES versus public RNA effect",
             subtitle = "Average GSE75010 RNA effect among scRNA GSEA core genes",
             x = "scRNA GSEA NES", y = "Mean public RNA log2FC", size = "Core genes") +
        theme_nature(8),
      "scRNA_enrichment/08_scRNA_pathway_NES_public_RNA_effect_scatter.pdf", 6.8, 5.2
    )
  }
}

manifest <- tibble(
  file = sort(gsub(paste0("^", gsub("/", "\\\\/", out_dir), "/?"), "", list.files(out_dir, pattern = "\\.pdf$", recursive = TRUE, full.names = TRUE))),
  purpose = case_when(
    grepl("GSE190971", file) ~ "GSE190971 public placenta RNA validation",
    grepl("GSE75010", file) ~ "GSE75010 public placenta RNA validation",
    grepl("scRNA_enrichment", file) ~ "scRNA pathway enrichment and public RNA signature validation",
    grepl("integrated", file) ~ "Integrated ATAC/scRNA/public RNA validation",
    TRUE ~ "Legacy top-level RNA validation figure"
  )
)
write_tsv(manifest, file.path(out_dir, "RNA_public_validation_figure_manifest.tsv"))

message("Generated ", nrow(manifest), " PDF figures in ", out_dir)
