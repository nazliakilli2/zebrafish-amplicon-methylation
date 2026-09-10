#!/usr/bin/env Rscript

# ============================================================
# Region-level (per-amplicon) differential methylation summary
#
# Standalone companion to DSS_v7.R. It does NOT modify any existing
# result. Rationale: the six amplicons are PRE-DEFINED target regions,
# so DMR *discovery* (callDMR) is not meaningful here — each short
# amplicon collapses to at most one region. What matters per region is
# (1) the effect size and (2) whether it changed.
#
# This script therefore treats each amplicon as the unit and tests it
# at the level of BIOLOGICAL REPLICATES, which is the honest inference:
#   - collapse every sample to ONE coverage-weighted methylation value
#     per amplicon (sum meth / sum coverage across the amplicon's CpGs),
#   - compare treatment vs timepoint-matched control samples with a
#     Welch t-test on those per-sample values.
# This deliberately uses the sample (n = 2-3), NOT the read count, as the
# unit — so it avoids the over-powering / correlated-CpG inflation that
# a read-count test on summed coverage would produce.
#
# Outputs:
#   DSS_results_v7/DSS_region_summary.tsv        - per gene x comparison table
#   plots_v10_publication/DSS_region_effect_tile.{pdf,png}  - main summary tile
#   plots_v10_publication/DSS_region_dotplot.{pdf,png}      - per-sample data
#
# Input: the same filtered per-sample CpG count table used by DSS_v7.R.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ======================================
# 1. PATHS & SETTINGS
# ======================================
source(file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "config.R"))
INPUT_FILE <- CPG_FILTERED
OUT_DIR    <- RESULTS_DIR
PLOT_DIR   <- FIGURE_DIR

MIN_COV  <- 10   # per-CpG coverage floor (matches DSS_v7.R)
MIN_NCPG <- 3    # ignore an amplicon in a sample if fewer CpGs pass the floor

dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

# ======================================
# 2. DISPLAY CONSTANTS (match Figures 1-3)
# ======================================
GENE_ORDER <- c('hsp70.1','hsp70.2','hsp70.3','actb2','b2m','cyp1a')

COMP_ORDER <- c('Fipronil_Low_7d','Fipronil_High_7d',
                'Cyfluthrin_Low_96h','Cyfluthrin_High_96h',
                'Cyfluthrin_Low_7d','Cyfluthrin_High_7d')
COMP_LABELS <- c(
  'Fipronil_Low_7d'    = 'Fip Low 7d',
  'Fipronil_High_7d'   = 'Fip High 7d',
  'Cyfluthrin_Low_96h' = 'Cyt Low 96h',
  'Cyfluthrin_High_96h'= 'Cyt High 96h',
  'Cyfluthrin_Low_7d'  = 'Cyt Low 7d',
  'Cyfluthrin_High_7d' = 'Cyt High 7d'
)

# Groups (for the per-sample dotplot), controls first
GROUP_ORDER <- c('Control_24h','Control_96h','Control_7d',
                 'Fipronil_Low_7d','Fipronil_High_24h','Fipronil_High_7d',
                 'Cyfluthrin_Low_96h','Cyfluthrin_High_96h',
                 'Cyfluthrin_Low_7d','Cyfluthrin_High_7d')
GROUP_COLOURS <- c(
  'Control_24h'='#999999','Control_96h'='#777777','Control_7d'='#555555',
  'Fipronil_Low_7d'='#a6cee3','Fipronil_High_24h'='#7fbfff','Fipronil_High_7d'='#1f78b4',
  'Cyfluthrin_Low_96h'='#fb9a99','Cyfluthrin_High_96h'='#e31a1c',
  'Cyfluthrin_Low_7d'='#fdbf6f','Cyfluthrin_High_7d'='#ff7f00'
)

theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", size = base_size - 1),
      axis.text        = element_text(size = base_size - 2),
      axis.title       = element_text(size = base_size - 1),
      plot.title       = element_text(face = "bold", size = base_size),
      plot.subtitle    = element_text(size = base_size - 2, color = "grey40"),
      legend.text      = element_text(size = base_size - 2),
      legend.title     = element_text(size = base_size - 1, face = "bold")
    )
}

stars <- function(p) fifelse(is.na(p), "",
                    fifelse(p < 0.001, "***",
                    fifelse(p < 0.01,  "**",
                    fifelse(p < 0.05,  "*", ""))))

# ======================================
# 3. METADATA PARSER (identical to DSS_v7.R)
# ======================================
parse_sample_metadata <- function(sample_names) {
  dt <- data.table(sample = sample_names)

  dt[, compound := fcase(
    grepl("^F", sample) & grepl("-K-", sample),                                          "Control",
    grepl("^F", sample) & grepl("-170-", sample),                                        "Fipronil",
    grepl("^F", sample) & grepl("-1700-", sample),                                       "Fipronil",
    grepl("^S", sample) & grepl("-K-", sample),                                          "Control",
    grepl("^S", sample) & grepl("-18-|18-", sample) & !grepl("-180-|-180$", sample),     "Cyfluthrin",
    grepl("^S", sample) & grepl("-180-|-180$", sample),                                  "Cyfluthrin",
    default = "Unknown"
  )]

  dt[, dose := fcase(
    grepl("-K-", sample),                                                                "Control",
    grepl("-170-", sample) & !grepl("-1700-", sample),                                   "Low",
    grepl("-1700-", sample),                                                             "High",
    grepl("^S", sample) & grepl("-18-|18-", sample) & !grepl("-180-|-180$|180-", sample),"Low",
    grepl("^S", sample) & grepl("-180-|-180$|180-", sample),                             "High",
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

  s_mask <- grepl("^S", dt$sample) & dt$timepoint == "Unknown"
  if (any(s_mask)) {
    dt[s_mask, timepoint := fcase(
      grepl("96",       sample), "96h",
      grepl("-7-|-7$",  sample), "7d",
      default = "Unknown"
    )]
  }

  dt[, group := fcase(
    dose == "Control",                          paste0("Control_",         timepoint),
    compound == "Fipronil"   & dose == "Low",  paste0("Fipronil_Low_",    timepoint),
    compound == "Fipronil"   & dose == "High", paste0("Fipronil_High_",   timepoint),
    compound == "Cyfluthrin" & dose == "Low",  paste0("Cyfluthrin_Low_",  timepoint),
    compound == "Cyfluthrin" & dose == "High", paste0("Cyfluthrin_High_", timepoint),
    default = "Unknown"
  )]

  dt[]
}

# ======================================
# 4. LOAD & COLLAPSE TO PER-SAMPLE, PER-AMPLICON METHYLATION
# ======================================
if (!file.exists(INPUT_FILE))
  stop("Filtered CpG counts not found: ", INPUT_FILE)

message(">>> Loading filtered CpG counts: ", INPUT_FILE)
raw <- fread(INPUT_FILE)

need <- c("chr","pos","sample","meth","coverage","amplicon")
miss <- setdiff(need, names(raw))
if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))

raw <- raw[coverage >= MIN_COV]

meta <- parse_sample_metadata(unique(raw$sample))
unknown <- meta[group %like% "Unknown", sample]
if (length(unknown) > 0)
  warning(length(unknown), " sample(s) unparsed and dropped:\n  ",
          paste(unknown, collapse = "\n  "))
meta <- meta[!(group %like% "Unknown")]
raw  <- merge(raw, meta, by = "sample")

# Coverage-weighted amplicon methylation per sample (the region value)
reg <- raw[, .(meth      = sum(meth),
               coverage  = sum(coverage),
               n_cpg     = uniqueN(pos)),
           by = .(sample, amplicon, group, compound, dose, timepoint)]
reg <- reg[n_cpg >= MIN_NCPG]
reg[, region_meth := 100 * meth / coverage]

message(">>> Per-sample x amplicon region values: ", nrow(reg), " rows")

# ======================================
# 5. REGION-LEVEL TEST: treatment vs timepoint-matched control
# ======================================
targets <- COMP_ORDER
res_list <- list()

for (tr in targets) {
  tp   <- sub(".*_", "", tr)          # trailing 24h / 96h / 7d
  ctrl <- paste0("Control_", tp)

  for (amp in GENE_ORDER) {
    t_vals <- reg[amplicon == amp & group == tr,   region_meth]
    c_vals <- reg[amplicon == amp & group == ctrl, region_meth]

    if (length(t_vals) >= 2 && length(c_vals) >= 2) {
      p <- tryCatch(t.test(t_vals, c_vals)$p.value, error = function(e) NA_real_)
      res_list[[paste(amp, tr, sep = "__")]] <- data.table(
        gene       = amp,
        comparison = tr,
        n_treat    = length(t_vals),
        n_ctrl     = length(c_vals),
        mean_treat = mean(t_vals),
        mean_ctrl  = mean(c_vals),
        diff_pp    = mean(t_vals) - mean(c_vals),
        pval       = p
      )
    }
  }
}

if (length(res_list) == 0)
  stop("No comparisons had >=2 samples in both groups. Check metadata parsing.")

res <- rbindlist(res_list)
res[, fdr := p.adjust(pval, method = "BH"), by = comparison]   # per comparison, like DSS_v7
res[, sig := stars(fdr)]
res[, gene       := factor(gene, levels = GENE_ORDER)]
res[, comparison := factor(comparison, levels = COMP_ORDER)]

# ======================================
# 6. WRITE TABLE
# ======================================
out <- res[order(comparison, gene),
           .(gene, comparison,
             n_treat, n_ctrl,
             mean_treat = round(mean_treat, 2),
             mean_ctrl  = round(mean_ctrl, 2),
             diff_pp    = round(diff_pp, 2),
             pval       = signif(pval, 3),
             fdr        = signif(fdr, 3),
             sig)]
fwrite(out, file.path(OUT_DIR, "DSS_region_summary.tsv"), sep = "\t")
message(">>> Wrote: ", file.path(OUT_DIR, "DSS_region_summary.tsv"),
        "  (", nrow(out), " gene x comparison tests)")

# ======================================
# 7. FIGURE 1 (main): region effect-size tile
# fill = mean methylation change (pp); label = effect + significance stars
# ======================================
message(">>> Building region effect tile...")

grid <- CJ(gene       = factor(GENE_ORDER, levels = GENE_ORDER),
           comparison = factor(COMP_ORDER, levels = COMP_ORDER))
tile <- merge(grid, res, by = c("gene","comparison"), all.x = TRUE)
tile[, tested := !is.na(diff_pp)]
tile[, label  := fifelse(tested,
                         paste0(sprintf("%+.1f", diff_pp), sig),
                         "n/a")]

lim <- max(abs(tile$diff_pp), na.rm = TRUE)
lim <- ifelse(is.finite(lim) && lim > 0, lim, 1)

p_tile <- ggplot(tile, aes(x = comparison, y = gene, fill = diff_pp)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label,
                color = tested & abs(diff_pp) > lim * 0.6),
            size = 3.1, fontface = "bold", show.legend = FALSE) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#D6604D",
                       midpoint = 0, limits = c(-lim, lim),
                       na.value = "grey92",
                       name = "Δ methylation\n(treat − ctrl, pp)") +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "grey25")) +
  scale_x_discrete(labels = COMP_LABELS[COMP_ORDER]) +
  scale_y_discrete(limits = rev(GENE_ORDER)) +
  theme_pub() +
  theme(axis.text.x  = element_text(angle = 35, hjust = 1),
        panel.grid   = element_blank(),
        panel.border = element_blank()) +
  labs(title    = "Region-level methylation change per gene and treatment",
       subtitle = "Coverage-weighted amplicon methylation, treatment vs matched control (Welch t-test on replicates).\nCell = Δ in percentage points; * FDR<0.05  ** <0.01  *** <0.001; n/a = <2 samples.",
       x = NULL, y = NULL)

ggsave(file.path(PLOT_DIR, "DSS_region_effect_tile.pdf"), p_tile, width = 8.5, height = 5)
ggsave(file.path(PLOT_DIR, "DSS_region_effect_tile.png"), p_tile, width = 8.5, height = 5,
       dpi = 350)
message("  Saved: DSS_region_effect_tile.{pdf,png}")

# ======================================
# 8. FIGURE 2 (companion): per-sample amplicon methylation by group
# Shows the actual replicates behind every cell (honest about n)
# ======================================
message(">>> Building per-sample dotplot...")

reg[, group := factor(group, levels = GROUP_ORDER)]
dot <- reg[!is.na(group)]
dot[, gene := factor(amplicon, levels = GENE_ORDER)]

p_dot <- ggplot(dot, aes(x = group, y = region_meth, fill = group)) +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = "crossbar", width = 0.6,
               color = "grey30", linewidth = 0.3, fatten = 0) +
  geom_jitter(width = 0.12, size = 2, shape = 21, color = "grey30",
              stroke = 0.3, alpha = 0.9) +
  scale_fill_manual(values = GROUP_COLOURS, drop = FALSE, guide = "none") +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  theme_pub(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
  labs(title    = "Amplicon methylation per sample (coverage-weighted)",
       subtitle = "Each point = one biological replicate; bar = group mean",
       x = NULL, y = "Region methylation (%)")

ggsave(file.path(PLOT_DIR, "DSS_region_dotplot.pdf"), p_dot, width = 12, height = 7)
ggsave(file.path(PLOT_DIR, "DSS_region_dotplot.png"), p_dot, width = 12, height = 7,
       dpi = 350)
message("  Saved: DSS_region_dotplot.{pdf,png}")

# ======================================
# 9. CONSOLE SUMMARY (self-check on a local run)
# ======================================
message("\n>>> Region effect table (Δ pp, FDR stars):")
print(dcast(res, gene ~ comparison,
            value.var = "diff_pp", fun.aggregate = function(x) round(x[1], 1)))
message("\n>>> Significant regions (FDR<0.05):")
print(res[fdr < 0.05, .(gene, comparison, diff_pp = round(diff_pp,1),
                        n_treat, n_ctrl, fdr = signif(fdr,3))][order(comparison, gene)])

message("\n>>> Done.")
message("    Table:   ", file.path(OUT_DIR, "DSS_region_summary.tsv"))
message("    Figures: ", file.path(PLOT_DIR, "DSS_region_effect_tile.pdf"))
message("             ", file.path(PLOT_DIR, "DSS_region_dotplot.pdf"))
