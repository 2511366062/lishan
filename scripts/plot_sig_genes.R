suppressPackageStartupMessages({
  library(ggplot2)
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

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  gsub("^_+|_+$", "", x)
}

as_bool <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

p_to_stars <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}

make_half_violin_data <- function(plot_df, group_order, width = 0.28) {
  out <- list()

  for (i in seq_along(group_order)) {
    g <- group_order[[i]]
    values <- plot_df$expression[as.character(plot_df$group) == g]
    values <- values[is.finite(values)]

    if (length(unique(values)) < 2) {
      next
    }

    d <- density(values, na.rm = TRUE)
    scaled <- d$y / max(d$y) * width

    out[[g]] <- data.frame(
      group = g,
      x = c(rep(i, length(d$x)), i + rev(scaled)),
      y = c(d$x, rev(d$x))
    )
  }

  if (length(out) == 0) {
    return(data.frame(group = character(), x = numeric(), y = numeric()))
  }

  do.call(rbind, out)
}

script_dir <- get_script_dir()
project_dir <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
config_base_dir <- read_config_value("R_BASE_DIR", file.path(project_dir, "config.py"))
if (is.na(config_base_dir)) {
  config_base_dir <- "D:/lxk/project/lishan-20260613"
}

base_dir <- ifelse(length(args) >= 1, args[[1]], config_base_dir)
deg_name <- ifelse(length(args) >= 2, args[[2]], "EVT")

deg_dir <- file.path(base_dir, "DEG", deg_name)
fig_dir <- file.path(base_dir, "fig", "DEG", deg_name)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

counts <- read.csv(file.path(deg_dir, "counts.csv"), row.names = 1, check.names = FALSE)
metadata <- read.csv(file.path(deg_dir, "metadata.csv"), stringsAsFactors = FALSE)
sig <- read.csv(file.path(deg_dir, "sig.csv"), stringsAsFactors = FALSE)

metadata <- metadata[match(colnames(counts), metadata$pb_sample), , drop = FALSE]

lib_size <- colSums(counts)
log_cpm <- log2(t(t(counts) / lib_size * 1e6) + 1)

if (!"show" %in% colnames(sig)) {
  sig$show <- TRUE
}

genes <- sig$gene[as_bool(sig$show)]
genes <- genes[genes %in% rownames(log_cpm)]

group_order <- c("NP", "PE")
group_order <- c(group_order[group_order %in% unique(metadata$group)],
                 setdiff(unique(metadata$group), group_order))
metadata$group <- factor(metadata$group, levels = group_order)

group_colors <- NULL
if ("group_color" %in% colnames(metadata)) {
  color_df <- unique(metadata[, c("group", "group_color")])
  color_df <- color_df[!is.na(color_df$group_color) & grepl("^#", color_df$group_color), ]
  if (nrow(color_df) > 0) {
    group_colors <- setNames(color_df$group_color, as.character(color_df$group))
  }
}

for (gene in genes) {
  row <- sig[match(gene, sig$gene), , drop = FALSE]
  direction <- tolower(ifelse("change" %in% colnames(row), row$change, "NS"))
  if (!direction %in% c("up", "down")) {
    direction <- "ns"
  }

  plot_df <- metadata
  plot_df$expression <- as.numeric(log_cpm[gene, plot_df$pb_sample])
  plot_df$x_center <- as.numeric(plot_df$group)
  set.seed(1)
  plot_df$x_point <- plot_df$x_center - 0.27 + runif(nrow(plot_df), -0.035, 0.035)
  violin_df <- make_half_violin_data(plot_df, group_order)
  violin_df$group <- factor(violin_df$group, levels = group_order)

  padj <- ifelse("padj" %in% colnames(row), row$padj, NA)
  logfc <- ifelse("log2FoldChange" %in% colnames(row), row$log2FoldChange, NA)
  y_max <- max(plot_df$expression, na.rm = TRUE)
  y_bar <- y_max + max(abs(y_max) * 0.08, 0.2)
  y_text <- y_bar + max(abs(y_max) * 0.04, 0.1)

  p <- ggplot() +
    geom_polygon(
      data = violin_df,
      aes(x = x, y = y, fill = group, group = group),
      alpha = 0.28,
      color = NA
    ) +
    geom_boxplot(
      data = plot_df,
      aes(x = x_center, y = expression, color = group, group = group),
      width = 0.25,
      outlier.shape = NA,
      fill = "white",
      linewidth = 0.45
    ) +
    geom_point(
      data = plot_df,
      aes(x = x_point, y = expression, color = group, fill = group),
      shape = 21,
      size = 2.3,
      stroke = 0.7,
      alpha = 0.85
    ) +
    annotate("segment", x = 1, xend = 1, y = y_bar, yend = y_text, linewidth = 0.35) +
    annotate("segment", x = 1, xend = 2, y = y_text, yend = y_text, linewidth = 0.35) +
    annotate("segment", x = 2, xend = 2, y = y_text, yend = y_bar, linewidth = 0.35) +
    annotate("text", x = 1.5, y = y_text, label = p_to_stars(padj), vjust = -0.25) +
    labs(
      title = gene,
      subtitle = sprintf("log2FC=%.2f, padj=%.2g", logfc, padj),
      x = NULL,
      y = "log2(CPM + 1)"
    ) +
    scale_x_continuous(
      breaks = seq_along(group_order),
      labels = group_order,
      limits = c(0.55, length(group_order) + 0.45),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.position = "none",
      axis.text.x = element_text(size = 11)
    )

  if (!is.null(group_colors)) {
    p <- p +
      scale_color_manual(values = group_colors, drop = FALSE) +
      scale_fill_manual(values = group_colors, drop = FALSE)
  }

  out_file <- file.path(fig_dir, paste0(safe_name(gene), "_", direction, ".pdf"))
  ggsave(out_file, p, width = 3.2, height = 4.2, device = cairo_pdf)
}

message("DEG folder: ", deg_dir)
message("Output folder: ", fig_dir)
message("Genes plotted: ", length(genes))
