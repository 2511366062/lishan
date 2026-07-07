suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
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

as_bool <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

clean_pathway <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("^GOBP_", "", x)
  x <- gsub("^GO_", "", x)
  x <- gsub("^REACTOME_", "", x)
  x <- gsub("_", " ", x)
  str_to_title(tolower(x))
}

script_dir <- get_script_dir()
project_dir <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
config_base_dir <- read_config_value("R_BASE_DIR", file.path(project_dir, "config.py"))
if (is.na(config_base_dir)) {
  config_base_dir <- "D:/lxk/project/lishan-20260613"
}

base_dir <- ifelse(length(args) >= 1, args[[1]], config_base_dir)
top_n <- ifelse(length(args) >= 2, as.integer(args[[2]]), 60)

deg_root <- file.path(base_dir, "DEG")
fig_dir <- file.path(base_dir, "fig", "GSEA")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

deg_dirs <- list.dirs(deg_root, recursive = FALSE, full.names = TRUE)
deg_dirs <- deg_dirs[basename(deg_dirs) != "scripts"]

gsea_list <- lapply(deg_dirs, function(dir) {
  file <- file.path(dir, "GSEA.csv")
  if (!file.exists(file)) {
    return(NULL)
  }
  df <- read.csv(file, stringsAsFactors = FALSE)
  if (nrow(df) == 0) {
    return(NULL)
  }
  if (!"show" %in% colnames(df)) {
    df$show <- TRUE
  }
  df <- df[as_bool(df$show), , drop = FALSE]
  if (nrow(df) == 0) {
    return(NULL)
  }
  df$celltype <- basename(dir)
  df
})

gsea_df <- bind_rows(gsea_list)
if (nrow(gsea_df) == 0) {
  stop("No shown GSEA pathways found.")
}

gsea_df <- gsea_df %>%
  filter(!is.na(NES), !is.na(p.adjust)) %>%
  mutate(
    pathway = ifelse(!is.na(Description) & Description != "", Description, ID),
    pathway_label = clean_pathway(pathway),
    neg_log10_padj = -log10(p.adjust),
    direction = ifelse(NES >= 0, "PE enriched", "NP enriched")
  )

pathway_order <- gsea_df %>%
  group_by(pathway_label) %>%
  summarise(
    max_abs_nes = max(abs(NES), na.rm = TRUE),
    min_padj = min(p.adjust, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(max_abs_nes), min_padj) %>%
  slice_head(n = top_n) %>%
  pull(pathway_label)

plot_df <- gsea_df %>%
  filter(pathway_label %in% pathway_order)

cell_order <- c(
  "DSC", "Decidualized_DSC", "Contractile_DSC", "Fibroblast",
  "EVT", "Trophoblast",
  "Endo", "Vascular_Endo", "ACKR1_Endo", "Lymphatic_Endo",
  "Myeloid", "HLAII_Mye", "FOLR2_Macro", "Infla_Mono", "APOE_Macro",
  "T", "CD8_cytoT", "TCF7_T", "Activated_T", "S100A4_T"
)
cell_order <- c(cell_order[cell_order %in% unique(plot_df$celltype)],
                setdiff(unique(plot_df$celltype), cell_order))

plot_df$celltype <- factor(plot_df$celltype, levels = cell_order)
plot_df$pathway_label <- factor(plot_df$pathway_label, levels = rev(pathway_order))

max_abs_nes <- max(abs(plot_df$NES), na.rm = TRUE)
size_max <- max(plot_df$neg_log10_padj, na.rm = TRUE)
size_breaks <- pretty(c(0, size_max), n = 4)
size_breaks <- size_breaks[size_breaks > 0]

height <- max(5.5, 0.23 * length(unique(plot_df$pathway_label)) + 1.6)
width <- max(8.5, 0.36 * length(unique(plot_df$celltype)) + 3.2)

p <- ggplot(plot_df, aes(x = celltype, y = pathway_label)) +
  geom_point(
    aes(size = neg_log10_padj, fill = NES),
    shape = 21,
    color = "grey18",
    stroke = 0.25,
    alpha = 0.95
  ) +
  scale_fill_gradient2(
    low = "#7FB89A",
    mid = "white",
    high = "#D98AA0",
    midpoint = 0,
    limits = c(-max_abs_nes, max_abs_nes),
    name = "NES"
  ) +
  scale_size_continuous(
    range = c(1.8, 7.5),
    breaks = size_breaks,
    name = expression(-log[10]("FDR"))
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 8.5),
    axis.text.y = element_text(color = "black", size = 8.3),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(8, 12, 8, 8)
  ) +
  guides(
    fill = guide_colorbar(order = 1, barheight = unit(38, "mm"), barwidth = unit(4, "mm")),
    size = guide_legend(order = 2, override.aes = list(fill = "grey70"))
  )

pdf_file <- file.path(fig_dir, "GSEA_dotplot.pdf")
png_file <- file.path(fig_dir, "GSEA_dotplot.png")

ggsave(pdf_file, p, width = width, height = height, device = cairo_pdf)
ggsave(png_file, p, width = width, height = height, dpi = 400)

write.csv(plot_df, file.path(fig_dir, "GSEA_dotplot_data.csv"), row.names = FALSE)

message("Saved: ", pdf_file)
message("Saved: ", png_file)
message("Saved: ", file.path(fig_dir, "GSEA_dotplot_data.csv"))

single_dir <- file.path(fig_dir, "by_cell")
dir.create(single_dir, recursive = TRUE, showWarnings = FALSE)

for (ct in unique(plot_df$celltype)) {
  single_df <- plot_df %>%
    dplyr::filter(celltype == ct) %>%
    dplyr::arrange(p.adjust, dplyr::desc(abs(NES))) %>%
    dplyr::distinct(pathway_label, .keep_all = TRUE) %>%
    dplyr::arrange(dplyr::desc(NES))

  if (nrow(single_df) == 0) {
    next
  }

  single_df$pathway_label <- factor(single_df$pathway_label, levels = rev(single_df$pathway_label))
  single_max_abs_nes <- max(abs(single_df$NES), na.rm = TRUE)
  single_size_max <- max(single_df$neg_log10_padj, na.rm = TRUE)
  single_size_breaks <- pretty(c(0, single_size_max), n = 4)
  single_size_breaks <- single_size_breaks[single_size_breaks > 0]

  p_single <- ggplot(single_df, aes(x = NES, y = pathway_label)) +
    geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35) +
    geom_segment(
      aes(x = 0, xend = NES, y = pathway_label, yend = pathway_label),
      color = "grey82",
      linewidth = 0.35
    ) +
    geom_point(
      aes(size = neg_log10_padj, fill = NES),
      shape = 21,
      color = "grey18",
      stroke = 0.25,
      alpha = 0.95
    ) +
    scale_fill_gradient2(
      low = "#7FB89A",
      mid = "white",
      high = "#D98AA0",
      midpoint = 0,
      limits = c(-single_max_abs_nes, single_max_abs_nes),
      name = "NES"
    ) +
    scale_size_continuous(
      range = c(2.2, 7.8),
      breaks = single_size_breaks,
      name = expression(-log[10]("FDR"))
    ) +
    labs(
      title = as.character(ct),
      x = "Normalized enrichment score",
      y = NULL
    ) +
    theme_classic(base_size = 9) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(color = "black", size = 8.5),
      axis.text.y = element_text(color = "black", size = 8.3),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      plot.margin = margin(8, 12, 8, 8)
    ) +
    guides(
      fill = guide_colorbar(order = 1, barheight = unit(34, "mm"), barwidth = unit(4, "mm")),
      size = guide_legend(order = 2, override.aes = list(fill = "grey70"))
    )

  single_height <- max(3.2, 0.28 * nrow(single_df) + 1.4)
  single_width <- 6.2
  ct_safe <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(ct))

  ggsave(
    file.path(single_dir, paste0(ct_safe, "_GSEA.pdf")),
    p_single,
    width = single_width,
    height = single_height,
    device = cairo_pdf
  )
  ggsave(
    file.path(single_dir, paste0(ct_safe, "_GSEA.png")),
    p_single,
    width = single_width,
    height = single_height,
    dpi = 400
  )
}

message("Saved per-cell plots to: ", single_dir)
