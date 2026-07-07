suppressPackageStartupMessages({
  library(ggplot2)
  library(ChIPseeker)
  library(GenomicFeatures)
  library(txdbmaker)
  library(GenomicRanges)
  library(rtracklayer)
  library(clusterProfiler)
  library(enrichplot)
  library(ReactomePA)
  library(KEGGREST)
})

root <- normalizePath("D:/lxk/project/lishan-20260613", winslash = "/", mustWork = TRUE)
atac <- file.path(root, "analysis/ATAC")
out <- file.path(root, "fig/ATAC_windowsR")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cache <- file.path(atac, "annotation/windowsR_cache")
dir.create(cache, recursive = TRUE, showWarnings = FALSE)

gtf <- normalizePath("C:/Users/admin/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu24.04LTS_79rhkp1fndgsc/LocalState/rootfs/home/vapor/database/source/human/genes.filtered.gtf",
                     winslash = "/", mustWork = FALSE)
if (!file.exists(gtf)) {
  gtf <- "//wsl.localhost/Ubuntu-24.04/home/vapor/database/source/human/genes.filtered.gtf"
}
if (!file.exists(gtf)) stop("Cannot find genes.filtered.gtf from Windows R.")

txdb_file <- file.path(cache, "genes.filtered.sqlite")
if (!file.exists(txdb_file)) {
  message("Building TxDb from GTF: ", gtf)
  txdb <- txdbmaker::makeTxDbFromGFF(gtf, format = "gtf")
  saveDb(txdb, txdb_file)
} else {
  txdb <- loadDb(txdb_file)
}

message("Importing GTF gene symbols")
gtf_gr <- import(gtf)
gene_rows <- gtf_gr[gtf_gr$type == "gene"]
gene_map <- unique(data.frame(
  gene_id = sub("\\..*$", "", mcols(gene_rows)$gene_id),
  gene_name = ifelse(is.na(mcols(gene_rows)$gene_name), sub("\\..*$", "", mcols(gene_rows)$gene_id), mcols(gene_rows)$gene_name),
  stringsAsFactors = FALSE
))
gene_map <- gene_map[!duplicated(gene_map$gene_id), ]
symbol_to_ensembl <- gene_map$gene_id
names(symbol_to_ensembl) <- gene_map$gene_name

make_peak_bed <- function(df, path) {
  write.table(df[, c("chrom", "start", "end")], path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}

desc <- read.delim(file.path(atac, "counts/consensus_peak_descriptive_PE_vs_NP.tsv"), check.names = FALSE)
sets <- list(
  consensus = desc,
  PE_open = desc[desc$candidate == "PE_open", ],
  PE_closed = desc[desc$candidate == "PE_closed", ]
)

anno_list <- list()
for (nm in names(sets)) {
  bed <- file.path(cache, paste0(nm, ".bed"))
  make_peak_bed(sets[[nm]], bed)
  peak <- readPeakFile(bed)
  anno <- annotatePeak(
    peak,
    TxDb = txdb,
    tssRegion = c(-3000, 3000),
    annoDb = NULL,
    verbose = FALSE
  )
  anno_df <- as.data.frame(anno)
  anno_df$set <- nm
  anno_list[[nm]] <- anno_df
  write.csv(anno_df, file.path(cache, paste0(nm, "_ChIPseeker_annotation.csv")), row.names = FALSE)

  pdf(file.path(out, paste0("ChIPseeker_", nm, "_annotation_pie.pdf")), width = 4.2, height = 4.2)
  print(plotAnnoPie(anno))
  dev.off()

  pdf(file.path(out, paste0("ChIPseeker_", nm, "_annotation_bar.pdf")), width = 5.2, height = 3.6)
  print(plotAnnoBar(anno))
  dev.off()

  pdf(file.path(out, paste0("ChIPseeker_", nm, "_distance_to_TSS.pdf")), width = 4.8, height = 3.6)
  print(plotDistToTSS(anno, title = paste0(nm, " peak distance to TSS")))
  dev.off()
}

anno_all <- do.call(rbind, anno_list)
write.csv(anno_all, file.path(out, "ChIPseeker_all_peak_annotations.csv"), row.names = FALSE)

pdf(file.path(out, "ChIPseeker_compare_annotation_bar.pdf"), width = 6.8, height = 4.0)
anno_tab <- as.data.frame(table(anno_all$set, anno_all$annotation))
colnames(anno_tab) <- c("set", "annotation", "count")
anno_tab <- anno_tab[anno_tab$count > 0, ]
print(
  ggplot(anno_tab, aes(set, count, fill = annotation)) +
    geom_col(position = "fill", width = 0.72, color = "white", linewidth = 0.1) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = NULL, y = "Fraction of peaks", fill = "Annotation") +
    theme_classic(base_size = 8) +
    theme(legend.position = "right")
)
dev.off()

linked <- read.delim(file.path(atac, "annotation/atac_scrna_linked_gene_summary.tsv"), check.names = FALSE)
scrna <- read.delim(file.path(atac, "annotation/scrna_sig_genes_all.tsv"), check.names = FALSE)
candidate_genes <- unique(c(linked$gene, scrna$gene[abs(scrna$log2FoldChange) > 1]))
candidate_genes <- candidate_genes[candidate_genes %in% names(symbol_to_ensembl)]

symbol_to_entrez_kegg <- function(symbols) {
  out_rows <- list()
  for (sym in unique(symbols)) {
    res <- tryCatch(keggFind("genes", paste0("hsa ", sym)), error = function(e) character())
    if (length(res) == 0) next
    ids <- sub("^hsa:", "", names(res))
    out_rows[[sym]] <- data.frame(SYMBOL = sym, ENTREZID = ids, stringsAsFactors = FALSE)
  }
  unique(do.call(rbind, out_rows))
}

kegg_ids <- symbol_to_entrez_kegg(candidate_genes)
write.csv(kegg_ids, file.path(out, "candidate_gene_kegg_entrez_mapping.csv"), row.names = FALSE)

if (!is.null(kegg_ids) && nrow(kegg_ids) >= 3) {
  kk <- enrichKEGG(gene = unique(kegg_ids$ENTREZID), organism = "hsa", pvalueCutoff = 0.2, qvalueCutoff = 0.5)
  if (!is.null(kk) && nrow(as.data.frame(kk)) > 0) {
    write.csv(as.data.frame(kk), file.path(out, "clusterProfiler_KEGG_enrichment.csv"), row.names = FALSE)
    pdf(file.path(out, "clusterProfiler_KEGG_dotplot.pdf"), width = 6.2, height = 4.2)
    print(dotplot(kk, showCategory = min(15, nrow(as.data.frame(kk)))) + ggtitle("KEGG enrichment of ATAC-scRNA linked genes"))
    dev.off()
    pdf(file.path(out, "clusterProfiler_KEGG_barplot.pdf"), width = 6.2, height = 4.2)
    print(barplot(kk, showCategory = min(15, nrow(as.data.frame(kk)))) + ggtitle("KEGG enrichment of ATAC-scRNA linked genes"))
    dev.off()
  }
}

reactome_ids <- kegg_ids
if (requireNamespace("org.Hs.eg.db", quietly = TRUE) && !is.null(reactome_ids) && nrow(reactome_ids) >= 3) {
  rr <- enrichPathway(gene = unique(reactome_ids$ENTREZID), organism = "human", pvalueCutoff = 0.2, readable = FALSE)
  if (!is.null(rr) && nrow(as.data.frame(rr)) > 0) {
    write.csv(as.data.frame(rr), file.path(out, "ReactomePA_pathway_enrichment.csv"), row.names = FALSE)
    pdf(file.path(out, "ReactomePA_dotplot.pdf"), width = 6.2, height = 4.2)
    print(dotplot(rr, showCategory = min(15, nrow(as.data.frame(rr)))) + ggtitle("Reactome enrichment of ATAC-scRNA linked genes"))
    dev.off()
  }
} else {
  message("Skip ReactomePA: org.Hs.eg.db is not installed in this Windows R library.")
}

message("Wrote Windows R ATAC figures to: ", out)
