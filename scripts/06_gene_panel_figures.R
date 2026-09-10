#!/usr/bin/env Rscript
# Raw methylation (top) / DML dots (bottom) per gene
# Three genes side by side
# Figure 2: hsp70.1, hsp70.2, hsp70.3
# Figure 3: b2m, cyp1a, actb2

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "config.R"))
CpG_FILE    <- CPG_FILTERED
RESULT_FILE <- DSS_RESULTS
PLOT_DIR    <- FIGURE_DIR
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

COMP_ORDER <- c('Fipronil_Low_7d','Fipronil_High_7d',
                'Cyfluthrin_Low_7d','Cyfluthrin_High_7d')
COMP_LABELS <- c(
  'Fipronil_Low_7d'   = 'Fipronil Low 7d',
  'Fipronil_High_7d'  = 'Fipronil High 7d',
  'Cyfluthrin_Low_7d' = 'Cyfluthrin Low 7d',
  'Cyfluthrin_High_7d'= 'Cyfluthrin High 7d'
)
GROUP_ORDER <- c('Control 7d','Fipronil Low 7d','Fipronil High 7d',
                 'Cyfluthrin Low 7d','Cyfluthrin High 7d')
GROUP_COLOURS <- c(
  'Control 7d'        = "#555555",
  'Fipronil Low 7d'   = "#a6cee3",
  'Fipronil High 7d'  = "#1f78b4",
  'Cyfluthrin Low 7d' = "#fdbf6f",
  'Cyfluthrin High 7d'= "#ff7f00"
)

parse_group <- function(s) {
  if (grepl("-K-",s) & (grepl("-7-",s)|grepl("7$",s)) &
      !grepl("-96-",s) & !grepl("96$",s))                                return("Control 7d")
  if (grepl("-170-",s) & (grepl("-7-",s)|grepl("7$",s)) &
      !grepl("-96-",s) & !grepl("96$",s))                                return("Fipronil Low 7d")
  if (grepl("-1700-",s) & (grepl("-7-",s)|grepl("7$",s)) &
      !grepl("-96-",s) & !grepl("96$",s))                                return("Fipronil High 7d")
  if (grepl("-18-",s) & !grepl("-180-|-180$",s) &
      (grepl("-7-",s)|grepl("7$",s)) &
      !grepl("-96-",s) & !grepl("96$",s))                                return("Cyfluthrin Low 7d")
  if (grepl("-180-|-180$",s) & (grepl("-7-",s)|grepl("7$",s)) &
      !grepl("-96-",s) & !grepl("96$",s))                                return("Cyfluthrin High 7d")
  return(NA_character_)
}

theme_pub <- function(base_size = 15) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor  = element_blank(),
      strip.background  = element_blank(),
      strip.text        = element_text(face = "bold", size = base_size),
      axis.title        = element_text(size = base_size),
      axis.text         = element_text(size = base_size - 1),
      plot.title        = element_text(size = base_size + 1, face = "bold"),
      legend.text       = element_text(size = base_size - 1),
      legend.title      = element_text(size = base_size, face = "bold"),
      legend.position   = "bottom"
    )
}

# Load
message(">>> Loading data...")
raw <- fread(CpG_FILE)[coverage >= 5]
raw[, meth_percent := 100 * meth / coverage]
meta <- data.table(sample = unique(raw$sample))
meta[, group := sapply(sample, parse_group)]
meta <- meta[!is.na(group)]
raw  <- merge(raw, meta, by = "sample")
raw[, group := factor(group, levels = GROUP_ORDER)]

dml <- fread(RESULT_FILE)[comparison %in% COMP_ORDER]
dml[, comparison := factor(comparison, levels = COMP_ORDER)]

# Raw panel — y from 60 to 100 to reduce height
make_raw_panel <- function(gene, title_label) {
  g_raw <- raw[amplicon == gene]
  if (nrow(g_raw) == 0) return(NULL)
  pos_map <- data.table(pos = sort(unique(g_raw$pos)))[, cpg_index := .I]
  g_raw   <- merge(g_raw, pos_map, by = "pos")

  ggplot(g_raw, aes(x = cpg_index, y = meth_percent,
                    group = sample, color = group)) +
    geom_line(alpha = 0.7, linewidth = 0.6) +
    geom_point(size = 1.5, alpha = 0.8) +
    scale_color_manual(values = GROUP_COLOURS, name = "Group") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(limits = c(0, 100),
                       breaks = seq(0, 100, 25)) +
    facet_wrap(~group, nrow = 2, ncol = 3) +
    theme_pub() +
    theme(legend.position = "none") +
    labs(title = title_label, x = "CpG index", y = "Methylation (%)")
}

# DML panel
make_dml_panel <- function(gene) {
  g_dml <- dml[amplicon == gene]
  if (nrow(g_dml) == 0) return(NULL)
  pos_map <- data.table(pos = sort(unique(g_dml$pos)))[, cpg_index := .I]
  g_dml   <- merge(g_dml, pos_map, by = "pos")
  g_dml[, sig       := fdr < 0.05]
  g_dml[, comp_label := factor(COMP_LABELS[as.character(comparison)],
                                levels = COMP_LABELS)]

  ggplot(g_dml, aes(x = cpg_index, y = diff * 100, color = sig)) +
    geom_hline(yintercept = 0, linetype = "dotted",
               color = "grey50", linewidth = 0.5) +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_color_manual(values = c("FALSE"="black","TRUE"="red"),
                       labels = c("FALSE"="FDR >= 0.05","TRUE"="FDR < 0.05"),
                       name = "") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(limits = c(-10, 10), breaks = seq(-10, 10, 5)) +
    facet_wrap(~comp_label, nrow = 4, ncol = 1) +
    theme_pub() +
    theme(legend.position = "bottom") +
    labs(x = "CpG Index", y = "Delta Methylation %")
}

build_figure <- function(genes, out_file, fig_width = 28, fig_height = 16) {
  message(">>> Building ", out_file, "...")
  panel_list <- list()
  for (i in seq_along(genes)) {
    gene        <- genes[i]
    letter_raw  <- letters[(i-1)*2 + 1]
    letter_dml  <- letters[(i-1)*2 + 2]
    p_raw <- make_raw_panel(gene, paste0(letter_raw, "   ", gene,
                                          " - methylation per sample"))
    p_dml <- make_dml_panel(gene)
    if (is.null(p_raw) | is.null(p_dml)) next
    p_dml <- p_dml + labs(title = paste0(letter_dml, "   ", gene, " - DML"))
    panel_list[[gene]] <- p_raw / p_dml +
      plot_layout(heights = c(1, 2))
  }
  if (length(panel_list) == 0) return(invisible(NULL))
  combined <- wrap_plots(panel_list, nrow = 1) +
    plot_annotation(theme = theme(plot.margin = margin(5,5,5,5)))
  pdf(file.path(PLOT_DIR, out_file), width = fig_width, height = fig_height)
  print(combined)
  dev.off()
  message("  Saved: ", out_file)
}

build_figure(c('hsp70.1','hsp70.2','hsp70.3'), "Figure2_hsp70.pdf",
             fig_width = 28, fig_height = 16)
build_figure(c('b2m','cyp1a','actb2'), "Figure3_b2m_cyp1a_actb2.pdf",
             fig_width = 28, fig_height = 16)

message("Done.")
