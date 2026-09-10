#!/usr/bin/env Rscript

# ============================================================
# SUPPLEMENTARY robustness figure: the hypomethylation signal does not
# depend on DSS smoothing or on treating correlated CpGs as independent.
#
# Standalone companion to DSS_v7.R / DSS_region_summary.R. It does NOT
# modify any existing result and re-runs no model.
#
# Motivation. In the main analysis DSS is run with smoothing = TRUE, which
# borrows information across neighbouring CpGs. Because CpGs a few bp apart
# are biologically co-methylated, the per-CpG DML *counts* are not
# independent findings. Two objections follow, and this figure answers both
# WITHOUT any modelling assumption:
#
#   Panel A - RAW SIGNAL (no model, no smoothing).
#       At each CpG, delta = mean(treatment methylation %) -
#       mean(control methylation %), computed directly from the pooled
#       count table. If hypomethylation is already visible here, it was not
#       manufactured by smoothing.  (This is a cleaner "without smoothing"
#       demonstration than a smoothing = FALSE DSS run, which still models.)
#
#   Panel B - INDEPENDENT UNIT (amplicon level).
#       Each pool is collapsed to ONE coverage-weighted methylation value
#       per amplicon; treatment and control pools are shown per locus. Here
#       every amplicon is a single observation, so the pattern cannot depend
#       on treating within-amplicon CpGs as independent.
#
# NOTE on replication: each sequencing sample is a POOL of ~10 fish, so a
# point in Panel B is a group-level average over many animals; the
# experimental unit for inference remains the pool (n per group is small).
# This figure is descriptive (effect sizes / direction), not a second
# significance test.
#
# Outputs (plots_v10_publication/):
#   DSS_supp_robustness.{pdf,png}        - combined Panel A + Panel B
#   DSS_supp_rawdelta_perCpG.{pdf,png}   - Panel A alone
#   DSS_supp_amplicon_pools.{pdf,png}    - Panel B alone
# Outputs (DSS_results_v7/):
#   DSS_supp_rawdelta_summary.tsv        - per gene x comparison raw-delta table
#   DSS_supp_panelA_by_amplicon.tsv      - Panel A per-amplicon table (summed
#                                          over the 4 comparisons): CpGs studied
#                                          vs raw hypo/hyper/no-change + percents
#   DSS_supp_panelA_counts.tsv           - Panel A count table: CpGs studied vs
#                                          below/above the ±DELTA_THRESH cutoff,
#                                          per comparison, with totals + percents
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

MIN_COV     <- 10  # per-CpG coverage floor (matches DSS_v7.R so Panel A uses the
                   # same sites the DSS analysis saw)
MIN_NCPG    <- 3   # ignore an amplicon in a pool if fewer CpGs pass the floor
MIN_SAMP_A  <- 1   # Panel A: min pools per group required at a CpG to draw its
                   # raw delta (1 = show maximum breadth; set 2 to require both).
DELTA_THRESH <- 2.5  # |Δ| (pp) needed to colour a CpG hypo/hyper in Panel A;
                   # smaller changes are drawn grey ("no change").

dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

# ======================================
# 2. DISPLAY CONSTANTS (match Figures 1-3 and DSS_region_summary.R)
# ======================================
GENE_ORDER <- c('hsp70.1','hsp70.2','hsp70.3','actb2','b2m','cyp1a')

# The four 7-day treatment comparisons reported in the paper, vs 7d control.
CTRL_7D <- 'Control_7d'
COMP_7D <- c('Fipronil_Low_7d','Fipronil_High_7d',
             'Cyfluthrin_Low_7d','Cyfluthrin_High_7d')
COMP_LABELS <- c(
  'Control_7d'         = 'Control 7d',
  'Fipronil_Low_7d'    = 'Fip Low 7d',
  'Fipronil_High_7d'   = 'Fip High 7d',
  'Cyfluthrin_Low_7d'  = 'Cyt Low 7d',
  'Cyfluthrin_High_7d' = 'Cyt High 7d'
)
GROUP_ORDER_B  <- c(CTRL_7D, COMP_7D)
GROUP_COLOURS  <- c(
  'Control_7d'         = '#555555',
  'Fipronil_Low_7d'    = '#a6cee3',
  'Fipronil_High_7d'   = '#1f78b4',
  'Cyfluthrin_Low_7d'  = '#fdbf6f',
  'Cyfluthrin_High_7d' = '#ff7f00'
)

# Direction colours: hypo / no-change / hyper, consistent with the diverging
# fill used elsewhere (blue = hypo, red = hyper, grey = below threshold).
DIR_COLOURS <- c(hypo = '#2166AC', none = '#BDBDBD', hyper = '#D6604D')

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

# ======================================
# 3. METADATA PARSER (identical to DSS_v7.R / DSS_region_summary.R)
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
# 4. LOAD & ATTACH METADATA
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

# keep only the 7d groups relevant to this figure
raw <- raw[group %in% GROUP_ORDER_B]
raw[, meth_pct := 100 * meth / coverage]

if (nrow(raw) == 0)
  stop("No 7d samples found after parsing. Check metadata parsing / group names.")

# ======================================
# 5. PANEL A - RAW per-CpG delta (no model, no smoothing)
#    group methylation at a CpG = mean of per-pool methylation % (pools
#    weighted equally); delta = treatment - control at the SAME CpG.
# ======================================
message(">>> Panel A: raw per-CpG deltas...")

cpg_grp <- raw[, .(grp_meth = mean(meth_pct), n_samp = .N),
               by = .(chr, pos, amplicon, group)]

ctrl_cpg <- cpg_grp[group == CTRL_7D,
                    .(chr, pos, amplicon, c_meth = grp_meth, n_c = n_samp)]

dA_list <- list()
for (tr in COMP_7D) {
  t_cpg <- cpg_grp[group == tr,
                   .(chr, pos, amplicon, t_meth = grp_meth, n_t = n_samp)]
  m <- merge(t_cpg, ctrl_cpg, by = c("chr","pos","amplicon"))
  m <- m[n_t >= MIN_SAMP_A & n_c >= MIN_SAMP_A]
  if (nrow(m) == 0) next
  m[, `:=`(comparison = tr, delta = t_meth - c_meth)]
  dA_list[[tr]] <- m
}
if (length(dA_list) == 0)
  stop("Panel A: no CpGs with matched treatment/control coverage. Check groups.")
dA <- rbindlist(dA_list)

# rank each CpG 5'->3' within its amplicon, then spread the ranks evenly across
# the full panel width (0-1) so every amplicon FILLS its rectangle regardless of
# CpG count -- low-CpG amplicons are not "chopped off" on the left.
idx <- unique(dA[, .(amplicon, pos)])[order(amplicon, pos)]
idx[, cpg_index := seq_len(.N), by = amplicon]
idx[, n_amp     := .N,          by = amplicon]
idx[, xpos      := (cpg_index - 0.5) / n_amp]    # evenly spaced within (0, 1)
dA <- merge(dA, idx, by = c("amplicon","pos"))

dA[, gene       := factor(amplicon,   levels = GENE_ORDER)]
dA[, comparison := factor(comparison, levels = COMP_7D)]

# effect-size-thresholded direction: only |Δ| > DELTA_THRESH counts as a change
dA[, dir := fcase(delta < -DELTA_THRESH, "hypo",
                  delta >  DELTA_THRESH, "hyper",
                  default = "none")]
dA[, dir := factor(dir, levels = c("hypo","none","hyper"))]

# per gene x comparison raw-delta summary (self-check + supporting table)
sumA <- dA[, .(n_cpg        = .N,
               n_hypo       = sum(delta < -DELTA_THRESH),
               n_hyper      = sum(delta >  DELTA_THRESH),
               n_nochange   = sum(abs(delta) <= DELTA_THRESH),
               mean_delta   = round(mean(delta), 2),
               median_delta = round(median(delta), 2)),
           by = .(gene, comparison)][order(comparison, gene)]
fwrite(sumA, file.path(OUT_DIR, "DSS_supp_rawdelta_summary.tsv"), sep = "\t")
message("    Wrote: ", file.path(OUT_DIR, "DSS_supp_rawdelta_summary.tsv"))

# ---- Panel A COUNT TABLE (for the supplement) --------------------------------
# "Of all CpGs studied, how many fall below / above the +/-DELTA_THRESH cutoff?"
# Reported per comparison, with an all-genes total per comparison and a grand
# total across the four 7d comparisons. Percentages are of the CpGs tested in
# that row (n_cpg). Counts are CpG x comparison observations, i.e. the points
# actually drawn in Panel A.
count_block <- function(d, group_cols) {
  d[, .(n_cpg      = .N,
        n_hypo     = sum(delta < -DELTA_THRESH),
        n_hyper    = sum(delta >  DELTA_THRESH),
        n_nochange = sum(abs(delta) <= DELTA_THRESH)),
    by = group_cols]
}

add_pct <- function(d) {
  d[, `:=`(
    pct_hypo     = round(100 * n_hypo     / n_cpg, 1),
    pct_hyper    = round(100 * n_hyper    / n_cpg, 1),
    pct_nochange = round(100 * n_nochange / n_cpg, 1)
  )][]
}

# ---- (a) PER-AMPLICON table (summed over the four 7d comparisons) ------------
# This is the supplementary-table view: for each amplicon, how many CpGs were
# studied and how many are raw-hypo / raw-hyper / no-change at |Δ| > DELTA_THRESH.
tbl_amp <- count_block(dA, "gene")[order(gene)]
tbl_amp[, gene := as.character(gene)]
tbl_tot <- count_block(dA, character(0))
tbl_tot[, gene := "TOTAL"]
panelA_by_amplicon <- add_pct(rbindlist(list(tbl_amp, tbl_tot), use.names = TRUE))
setcolorder(panelA_by_amplicon,
            c("gene","n_cpg","n_hypo","n_hyper","n_nochange",
              "pct_hypo","pct_hyper","pct_nochange"))
fwrite(panelA_by_amplicon, file.path(OUT_DIR, "DSS_supp_panelA_by_amplicon.tsv"), sep = "\t")
message("    Wrote: ", file.path(OUT_DIR, "DSS_supp_panelA_by_amplicon.tsv"))

# ---- (b) full breakdown: per comparison x gene, per-comparison, grand total --
tbl_gene <- count_block(dA, c("comparison","gene"))[order(comparison, gene)]
tbl_gene[, gene := as.character(gene)]

tbl_comp <- count_block(dA, "comparison")[order(comparison)]
tbl_comp[, gene := "ALL genes"]

tbl_all  <- count_block(dA, character(0))
tbl_all[, `:=`(comparison = "ALL 7d comparisons", gene = "ALL genes")]

panelA_counts <- rbindlist(list(tbl_gene, tbl_comp, tbl_all), use.names = TRUE)
panelA_counts[, comparison := as.character(comparison)]
setcolorder(panelA_counts, c("comparison","gene","n_cpg",
                             "n_hypo","n_hyper","n_nochange"))
add_pct(panelA_counts)
fwrite(panelA_counts, file.path(OUT_DIR, "DSS_supp_panelA_counts.tsv"), sep = "\t")
message("    Wrote: ", file.path(OUT_DIR, "DSS_supp_panelA_counts.tsv"))

pA <- ggplot(dA, aes(xpos, delta, color = dir)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.3) +
  geom_hline(yintercept = c(-DELTA_THRESH, DELTA_THRESH),
             linetype = "dotted", color = "grey70", linewidth = 0.3) +
  geom_point(size = 1.05, alpha = 0.85) +
  facet_grid(gene ~ comparison,
             labeller = labeller(comparison = COMP_LABELS)) +
  scale_color_manual(values = DIR_COLOURS,
                     breaks = c("hypo","none","hyper"),
                     labels = c(sprintf("hypomethylated (Δ < −%g)", DELTA_THRESH),
                                sprintf("no change (|Δ| ≤ %g)", DELTA_THRESH),
                                sprintf("hypermethylated (Δ > +%g)", DELTA_THRESH)),
                     name   = "Raw effect (pp)",
                     drop   = FALSE) +
  scale_x_continuous(breaks = NULL, expand = expansion(mult = 0.04)) +
  theme_pub(base_size = 10) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "top") +
  labs(title    = "A  Raw per-CpG methylation change (no model, no smoothing)",
       subtitle = sprintf("Δ = mean treatment − mean control methylation at each CpG, from pooled counts. Coloured only when |Δ| > %g pp; hypomethylation dominates.", DELTA_THRESH),
       x = "CpG position within amplicon (5′ → 3′)",
       y = "Δ methylation (percentage points)")

ggsave(file.path(PLOT_DIR, "DSS_supp_rawdelta_perCpG.pdf"), pA, width = 9.5, height = 8)
ggsave(file.path(PLOT_DIR, "DSS_supp_rawdelta_perCpG.png"), pA, width = 9.5, height = 8, dpi = 350)
message("    Saved: DSS_supp_rawdelta_perCpG.{pdf,png}")

# ======================================
# 6. PANEL B - amplicon-level methylation, one dot per pool
#    coverage-weighted methylation per pool (identical collapse to
#    DSS_region_summary.R), treatment and control pools shown per locus.
# ======================================
message(">>> Panel B: amplicon-level pools...")

reg <- raw[, .(meth = sum(meth), coverage = sum(coverage), n_cpg = uniqueN(pos)),
           by = .(sample, amplicon, group)]
reg <- reg[n_cpg >= MIN_NCPG]
reg[, region_meth := 100 * meth / coverage]
reg[, group := factor(group, levels = GROUP_ORDER_B)]
reg[, gene  := factor(amplicon, levels = GENE_ORDER)]

pB <- ggplot(reg, aes(group, region_meth, fill = group)) +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = "crossbar", width = 0.6,
               color = "grey30", linewidth = 0.3, fatten = 0) +
  geom_point(shape = 21, size = 2.2, color = "grey20", stroke = 0.3,
             position = position_jitter(width = 0.1, height = 0)) +
  facet_wrap(~ gene, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = GROUP_COLOURS, drop = FALSE, guide = "none") +
  scale_x_discrete(labels = COMP_LABELS[GROUP_ORDER_B]) +
  theme_pub(base_size = 10) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8)) +
  labs(title    = "B  Amplicon-level methylation per pool (the independent unit)",
       subtitle = "Each point = one pooled sample (~10 fish); bar = group mean. The locus-level shift does not depend on within-amplicon CpG independence.",
       x = NULL, y = "Coverage-weighted amplicon methylation (%)")

ggsave(file.path(PLOT_DIR, "DSS_supp_amplicon_pools.pdf"), pB, width = 9.5, height = 6)
ggsave(file.path(PLOT_DIR, "DSS_supp_amplicon_pools.png"), pB, width = 9.5, height = 6, dpi = 350)
message("    Saved: DSS_supp_amplicon_pools.{pdf,png}")

# ======================================
# 7. COMBINED SUPPLEMENTARY FIGURE (Panel A over Panel B)
# ======================================
if (requireNamespace("patchwork", quietly = TRUE)) {
  suppressPackageStartupMessages(library(patchwork))
  combined <- pA / pB + plot_layout(heights = c(2.1, 1.2))
  ggsave(file.path(PLOT_DIR, "DSS_supp_robustness.pdf"), combined, width = 10, height = 13)
  ggsave(file.path(PLOT_DIR, "DSS_supp_robustness.png"), combined, width = 10, height = 13, dpi = 300)
  message("    Saved combined: DSS_supp_robustness.{pdf,png}")
} else {
  message("    patchwork not installed - combined figure skipped; ",
          "use the two individual panels (install.packages('patchwork') to combine).")
}

# ======================================
# 8. CONSOLE SUMMARY (self-check on a local run)
# ======================================
message("\n>>> Raw per-CpG delta summary (Panel A; threshold ±", DELTA_THRESH, " pp):")
print(sumA)
message("\n>>> Panel A PER-AMPLICON table (raw, no model; |Δ| > ", DELTA_THRESH, " pp):")
print(panelA_by_amplicon)
message("\n>>> Panel A full breakdown (per comparison x gene):")
print(panelA_counts)

message("\n>>> Done.")
message("    Tables:  ", file.path(OUT_DIR, "DSS_supp_panelA_by_amplicon.tsv"))
message("             ", file.path(OUT_DIR, "DSS_supp_panelA_counts.tsv"))
message("             ", file.path(OUT_DIR, "DSS_supp_rawdelta_summary.tsv"))
message("    Figures: ", file.path(PLOT_DIR, "DSS_supp_robustness.pdf"))
message("             ", file.path(PLOT_DIR, "DSS_supp_rawdelta_perCpG.pdf"))
message("             ", file.path(PLOT_DIR, "DSS_supp_amplicon_pools.pdf"))
