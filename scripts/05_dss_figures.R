#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
})

# ======================================
# 1. PATHS
# ======================================
source(file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "config.R"))
RESULT_FILE <- DSS_RESULTS
RAW_FILE    <- CPG_FILTERED
PLOT_DIR    <- FIGURE_DIR

dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

message(">>> Loading Data...")
dml <- fread(RESULT_FILE)
raw <- fread(RAW_FILE)

# ======================================
# 2. DATA HARMONIZATION
# dml$amplicon is now gene name (e.g. "hsp70.1") from master_bed_filtered
# raw$gene_name is also gene name -> join directly on amplicon + pos
# ======================================
message(">>> Mapping DML results to gene names...")

pos_to_gene   <- unique(raw[, .(amplicon, pos, gene_name)])
dml_annotated <- merge(dml, pos_to_gene, by = c("amplicon", "pos"), all.x = TRUE)

# ======================================
# 3. METADATA & CLEAN LABELS
# Same sample-prefix-aware parser as DSS_v7.R
# ======================================
parse_sample_metadata <- function(sample_names) {
  dt <- data.table(sample = sample_names)

  dt[, compound := fcase(
    grepl("^F", sample) & grepl("-K-", sample),                                               "Control",
    grepl("^F", sample) & grepl("-170-", sample),                                             "Fipronil",
    grepl("^F", sample) & grepl("-1700-", sample),                                            "Fipronil",
    grepl("^S", sample) & grepl("-K-", sample),                                               "Control",
    grepl("^S", sample) & grepl("-18-|18-", sample) & !grepl("-180-|-180$", sample),          "Cyfluthrin",
    grepl("^S", sample) & grepl("-180-|-180$", sample),                                       "Cyfluthrin",
    default = "Unknown"
  )]

  dt[, dose := fcase(
    grepl("-K-", sample),                                                                      "Control",
    grepl("-170-", sample) & !grepl("-1700-", sample),                                        "Low",
    grepl("-1700-", sample),                                                                   "High",
    grepl("^S", sample) & grepl("-18-|18-", sample) & !grepl("-180-|-180$|180-", sample),     "Low",
    grepl("^S", sample) & grepl("-180-|-180$|180-", sample),                                  "High",
    default = "Unknown"
  )]

  dt[, timepoint := fcase(
    grepl("^F.*-K-24",    sample),                                    "24h",
    grepl("^F.*-K-96",    sample),                                    "96h",
    grepl("^F.*-K-7",     sample),                                    "7d",
    grepl("^F.*-170-7|^F.*-1700-7", sample),                         "7d",
    grepl("^F.*-1700-24", sample),                                    "24h",
    grepl("^S.*-K-96",    sample),                                    "96h",
    grepl("^S.*-K-7",     sample),                                    "7d",
    grepl("^S.*-96$|-96-", sample),                                   "96h",
    grepl("^S.*-7$", sample) | grepl("-7-[0-9]|[0-9]-7$", sample),  "7d",
    default = "Unknown"
  )]

  # Fallback second pass for any S* still Unknown
  s_mask <- grepl("^S", dt$sample) & dt$timepoint == "Unknown"
  if (any(s_mask)) {
    dt[s_mask, timepoint := fcase(
      grepl("96",      sample), "96h",
      grepl("-7-|-7$", sample), "7d",
      default = "Unknown"
    )]
  }

  dt[, pub_label := fcase(
    dose == "Control",                         paste0("Control ", timepoint),
    compound == "Fipronil"   & dose == "Low",  paste0("Fipronil Low ", timepoint),
    compound == "Fipronil"   & dose == "High", paste0("Fipronil High ", timepoint),
    compound == "Cyfluthrin" & dose == "Low",  paste0("Cyfluthrin Low ", timepoint),
    compound == "Cyfluthrin" & dose == "High", paste0("Cyfluthrin High ", timepoint),
    default = "Unknown"
  )]

  return(dt)
}

raw[, meth_percent := 100 * (meth / coverage)]
raw_plot <- raw[coverage >= 5]
meta     <- parse_sample_metadata(unique(raw_plot$sample))

unknown_labels <- meta[pub_label == "Unknown", sample]
if (length(unknown_labels) > 0)
  warning(length(unknown_labels), " sample(s) got Unknown label: ",
          paste(unknown_labels, collapse = ", "))

raw_plot <- merge(raw_plot, meta, by = "sample")

# ======================================
# 4. GLOBAL PLOT CONSTANTS
# Fixed across all genes so figures are directly comparable
# ======================================

# All 6 comparisons in a fixed display order
ALL_COMPARISONS <- c(
  "Fipronil_Low_7d",
  "Fipronil_High_7d",
  "Cyfluthrin_Low_96h",
  "Cyfluthrin_High_96h",
  "Cyfluthrin_Low_7d",
  "Cyfluthrin_High_7d"
)

Y_MIN     <- -10  # covers b2m Fipronil_Low_7d which reaches -8.3%
Y_MAX     <-  10
PDF_W     <- 12   # fixed PDF width (inches)
PDF_H_P1  <- 14   # fixed PDF height for line plot page (6 panels * ~2in each + margins)
PDF_H_P2  <-  7   # fixed PDF height for heatmap page

# ======================================
# 5. PLOTTING LOOP
# ======================================
genes_to_plot <- unique(dml_annotated[!is.na(gene_name)]$gene_name)
message(">>> Genes to plot: ", paste(sort(genes_to_plot), collapse = ", "))

for (gn in genes_to_plot) {

  message(">>> Creating Publication Analysis for: ", gn)

  g_dml <- dml_annotated[gene_name == gn][order(pos)]
  g_raw <- raw_plot[gene_name == gn]

  if (nrow(g_raw) == 0) {
    message("  No raw data for ", gn, " — skipping")
    next
  }

  pos_map <- data.table(pos = sort(unique(g_dml$pos)))[, cpg_index := .I]
  g_dml   <- merge(g_dml, pos_map, by = "pos")
  g_raw   <- merge(g_raw, pos_map, by = "pos")  # inner join: only DSS-tested positions

  # Force comparison to ordered factor BEFORE passing to ggplot
  # This ensures facet panels appear in ALL_COMPARISONS order for every gene
  g_dml[, comparison := factor(comparison, levels = ALL_COMPARISONS)]

  # --- Page 1: Delta Methylation ---
  pdf(file.path(PLOT_DIR, paste0(gn, "_DML.pdf")), width = PDF_W, height = PDF_H_P1)

  p1 <- ggplot(g_dml, aes(x = cpg_index, y = diff * 100, color = fdr < 0.05)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "black", alpha = 0.6) +
    geom_line(color = "grey", alpha = 0.4) +
    geom_point(size = 2.5) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"),
                       labels = c("FALSE" = "FDR >= 0.05", "TRUE" = "FDR < 0.05"),
                       name   = "") +
    scale_y_continuous(limits = c(Y_MIN, Y_MAX),
                       breaks = seq(Y_MIN, Y_MAX, by = 5)) +
    facet_wrap(~comparison, ncol = 1, drop = FALSE) +  # drop=FALSE keeps empty panels
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text        = element_text(size = 9, face = "bold"),
          legend.position   = "bottom",
          panel.grid.minor  = element_blank()) +
    labs(title    = "Differential Methylation Analysis",
         subtitle = paste("Gene:", gn),
         y        = "Delta Methylation %",
         x        = "CpG Index")
  print(p1)
  dev.off()

  # --- Page 2: Raw methylation heatmap ---
  pdf(file.path(PLOT_DIR, paste0(gn, "_Heatmap.pdf")), width = PDF_W, height = PDF_H_P2)

  # Samples within same pub_label averaged into one row
  mat_wide <- dcast(g_raw, pub_label ~ cpg_index,
                    value.var     = "meth_percent",
                    fun.aggregate = mean)
  mat_plot           <- as.matrix(mat_wide[, -1])
  rownames(mat_plot) <- mat_wide$pub_label

  pheatmap(mat_plot,
           cluster_cols    = FALSE,
           cluster_rows    = FALSE,
           main            = paste("Raw Methylation Heatmap:", gn),
           color           = colorRampPalette(c("blue", "white", "red"))(100),
           na_col          = "grey90",
           border_color    = "white",
           fontsize_row    = 10,
           display_numbers = FALSE)

  dev.off()
  message("  Saved: ", gn, "_DML.pdf + ", gn, "_Heatmap.pdf")
}

message("\n>>> SUCCESS. Plots saved to: ", PLOT_DIR)
message(">>> Each gene produces two files: *_DML.pdf and *_Heatmap.pdf")
