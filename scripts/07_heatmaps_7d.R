#!/usr/bin/env Rscript
# Supplementary Figure: CpG methylation heatmaps — 7d only
# Shows ONLY positions with >=5x coverage in all 5 groups (no gaps, no "-")
# Treatment labels placed once on the right of the combined figure
# Wide cells, compact height

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "config.R"))
RAW_FILE    <- CPG_FILTERED
RESULT_FILE <- DSS_RESULTS
PLOT_DIR    <- FIGURE_DIR
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

GENE_ORDER <- c('hsp70.1','hsp70.2','hsp70.3','actb2','b2m','cyp1a')

GROUP_ORDER_7D <- c(
  'Control 7d',
  'Fipronil Low 7d', 'Fipronil High 7d',
  'Cyfluthrin Low 7d', 'Cyfluthrin High 7d'
)
# Short labels for y axis
GROUP_SHORT <- c(
  'Control 7d'         = 'Control',
  'Fipronil Low 7d'    = 'Fip Low',
  'Fipronil High 7d'   = 'Fip High',
  'Cyfluthrin Low 7d'  = 'Cyt Low',
  'Cyfluthrin High 7d' = 'Cyt High'
)

parse_group_7d <- function(s) {
  if (grepl("-K-", s)) {
    if (grepl("-24-",s)|grepl("-96-",s)|grepl("96$",s)) return(NA_character_)
    if (grepl("-7-", s)|grepl("7$",  s))                return("Control 7d")
    return(NA_character_)
  }
  if (grepl("-170-",  s)&(grepl("-7-",s)|grepl("7$",s))&!grepl("-96-",s)&!grepl("96$",s)) return("Fipronil Low 7d")
  if (grepl("-1700-", s)&(grepl("-7-",s)|grepl("7$",s))&!grepl("-96-",s)&!grepl("96$",s)) return("Fipronil High 7d")
  if (grepl("-18-",s)&!grepl("-180-|-180$",s)&(grepl("-7-",s)|grepl("7$",s))&!grepl("-96-",s)&!grepl("96$",s)) return("Cyfluthrin Low 7d")
  if (grepl("-180-|-180$",s)&(grepl("-7-",s)|grepl("7$",s))&!grepl("-96-",s)&!grepl("96$",s)) return("Cyfluthrin High 7d")
  return(NA_character_)
}

# ======================================
# LOAD
# ======================================
message(">>> Loading data...")
raw <- fread(RAW_FILE)[coverage >= 5]
raw[, meth_percent := 100 * meth / coverage]
meta <- data.table(sample = unique(raw$sample))
meta[, group := sapply(sample, parse_group_7d)]
raw  <- merge(raw, meta[!is.na(group)], by = "sample")
raw[, group := factor(group, levels = GROUP_ORDER_7D)]

dml      <- fread(RESULT_FILE)
comps_7d <- c('Fipronil_Low_7d','Fipronil_High_7d',
               'Cyfluthrin_Low_7d','Cyfluthrin_High_7d')
dss_pos  <- unique(dml[comparison %in% comps_7d, .(amplicon, pos)])

# ======================================
# HEATMAP FUNCTION
# show_y: whether to show treatment labels on y axis
# ======================================
make_heatmap <- function(gene, panel_letter, show_y = FALSE) {
  message("  ", gene)
  gene_pos <- dss_pos[amplicon == gene]
  g_raw    <- merge(raw[amplicon == gene],
                    gene_pos[, .(pos)], by = "pos")

  # Group mean per position (uses available samples)
  g_mean <- g_raw[, .(meth = mean(meth_percent, na.rm = TRUE)),
                  by = .(group, pos)]

  # Re-index ALL DSS positions 1..N (keep all 42 for cyp1a etc.)
  pos_map <- data.table(pos = sort(unique(gene_pos$pos)))[, cpg_index := .I]

  # Complete grid: every DSS position x every group
  full_grid <- CJ(group = factor(GROUP_ORDER_7D, levels = rev(GROUP_ORDER_7D)),
                  pos   = pos_map$pos)
  g_mean <- merge(full_grid, g_mean, by = c("group","pos"), all.x = TRUE)
  g_mean <- merge(g_mean, pos_map, by = "pos")

  g_mean[, group := factor(group, levels = rev(GROUP_ORDER_7D))]
  n_pos <- nrow(pos_map)

  p <- ggplot(g_mean, aes(x = cpg_index, y = group)) +
    geom_tile(aes(fill = meth), color = "white", linewidth = 0.4) +
    scale_fill_gradientn(
      colors = c("#2166AC","#6BAED6","#F7FBFF","#FDAE6B","#D94801"),
      limits = c(0, 100),
      na.value = "grey85",
      name   = "Mean methylation (%)",
      breaks = c(0, 25, 50, 75, 100),
      guide  = guide_colorbar(barwidth = 6, barheight = 0.5,
                               title.position = "top", title.hjust = 0.5)) +
    scale_x_continuous(
      breaks = if (n_pos <= 15) seq(1, n_pos)
               else if (n_pos <= 25) seq(2, n_pos, 2)
               else seq(4, n_pos, 4),
      expand = expansion(add = c(0.5, 0.5))) +
    scale_y_discrete(labels = GROUP_SHORT[rev(GROUP_ORDER_7D)]) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid    = element_blank(),
      panel.border  = element_rect(color = "grey70", linewidth = 0.4),
      axis.text.x   = element_text(size = 7.5),
      axis.title.x  = element_text(size = 8),
      axis.title.y  = element_blank(),
      axis.ticks    = element_line(color = "grey70", linewidth = 0.3),
      plot.title    = element_text(size = 9, face = "bold",
                                   margin = margin(b = 1)),
      plot.subtitle = element_text(size = 7.5, color = "grey50",
                                   margin = margin(b = 3))
    ) +
    labs(title    = paste0(panel_letter, "   ", gene),
         subtitle = paste0(n_pos, " CpG sites (DSS-tested)"),
         x        = "CpG index")

  # Only show y-axis treatment labels on left-column panels
  if (show_y) {
    p <- p + theme(axis.text.y = element_text(size = 8, hjust = 1,
                                              face = "bold"))
  } else {
    p <- p + theme(axis.text.y = element_blank(),
                   axis.ticks.y = element_blank())
  }
  p
}

# ======================================
# BUILD — show y labels only on left column (genes 1 and 4)
# ======================================
message(">>> Building panels...")
panels <- list()
for (i in seq_along(GENE_ORDER)) {
  show_y <- i %in% c(1, 4)   # left column in 2x3 grid
  panels[[GENE_ORDER[i]]] <- make_heatmap(GENE_ORDER[i], letters[i], show_y)
}

message(">>> Combining...")
combined <- wrap_plots(panels, nrow = 2, ncol = 3,
                       guides = "collect") &
  theme(legend.position = "bottom")

combined <- combined +
  plot_annotation(
    title    = "Supplementary Figure S2 — CpG methylation heatmaps (7-day exposures)",
    subtitle = paste0("Group mean CpG methylation (%) at all DSS-tested positions ",
                      "(matches Figures 2–3). Grey = not covered at ≥5× in that group."),
    theme = theme(
      plot.title    = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 8, color = "grey40"))
  )

# Wide, compact height — fits A4 portrait width
pdf(file.path(PLOT_DIR, "FigS_MethylationHeatmaps_7d.pdf"),
    width = 8.27, height = 6)
print(combined)
dev.off()

png(file.path(PLOT_DIR, "FigS_MethylationHeatmaps_7d.png"),
    width = 8.27, height = 6, units = "in", res = 350)
print(combined)
dev.off()

message("Saved: FigS_MethylationHeatmaps_7d.pdf + .png")
message("Done.")
