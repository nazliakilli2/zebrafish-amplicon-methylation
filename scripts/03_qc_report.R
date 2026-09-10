#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

# ======================================
# 1. PATHS
# ======================================
source(file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "config.R"))
CpG_FILE <- CPG_FILTERED
CHH_FILE <- CHH_FILTERED
CHG_FILE <- CHG_FILTERED
# QC_DIR is provided by config.R
dir.create(QC_DIR, showWarnings = FALSE, recursive = TRUE)

# ======================================
# 2. SHARED THEME
# Clean publication theme - no ggpubr dependency
# ======================================
theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor  = element_blank(),
      strip.background  = element_blank(),
      strip.text        = element_text(face = "bold", size = base_size - 1),
      legend.text       = element_text(size = base_size - 2),
      legend.title      = element_text(size = base_size - 1, face = "bold"),
      axis.text         = element_text(size = base_size - 2),
      axis.title        = element_text(size = base_size - 1),
      plot.title        = element_text(face = "bold", size = base_size),
      plot.subtitle     = element_text(size = base_size - 2, color = "grey40")
    )
}

# ======================================
# 3. COLOUR PALETTE
# Fixed group order: Control first, then Fipronil, then Cyfluthrin
# ======================================
GROUP_COLOURS <- c(
  "Control_Control_24h"  = "#999999",
  "Control_Control_96h"  = "#777777",
  "Control_Control_7d"   = "#555555",
  "Fipronil_Low_7d"      = "#a6cee3",
  "Fipronil_High_24h"    = "#7fbfff",
  "Fipronil_High_7d"     = "#1f78b4",
  "Cyfluthrin_Low_96h"   = "#fb9a99",
  "Cyfluthrin_High_96h"  = "#e31a1c",
  "Cyfluthrin_Low_7d"    = "#fdbf6f",
  "Cyfluthrin_High_7d"   = "#ff7f00"
)

GROUP_ORDER <- names(GROUP_COLOURS)

# ======================================
# 4. METADATA PARSER
# Same sample-prefix-aware parser as DSS_v7.R
# ======================================
parse_meta <- function(sample_names) {
  res <- data.table(sample = sample_names)

  res[, Treatment := fcase(
    grepl("^F", sample) & grepl("-K-", sample),                                               "Control",
    grepl("^F", sample) & grepl("-170-", sample),                                             "Fipronil",
    grepl("^F", sample) & grepl("-1700-", sample),                                            "Fipronil",
    grepl("^S", sample) & grepl("-K-", sample),                                               "Control",
    grepl("^S", sample) & grepl("-18-|18-", sample) & !grepl("-180-|-180$", sample),          "Cyfluthrin",
    grepl("^S", sample) & grepl("-180-|-180$", sample),                                       "Cyfluthrin",
    default = "Unknown"
  )]

  res[, Dose := fcase(
    grepl("-K-", sample),                                                                      "Control",
    grepl("-170-", sample) & !grepl("-1700-", sample),                                        "Low",
    grepl("-1700-", sample),                                                                   "High",
    grepl("^S", sample) & grepl("-18-|18-", sample) & !grepl("-180-|-180$|180-", sample),     "Low",
    grepl("^S", sample) & grepl("-180-|-180$|180-", sample),                                  "High",
    default = "Unknown"
  )]

  res[, Duration := fcase(
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

  # Second-pass fallback for any S* still Unknown
  s_mask <- grepl("^S", res$sample) & res$Duration == "Unknown"
  if (any(s_mask)) {
    res[s_mask, Duration := fcase(
      grepl("96",      sample), "96h",
      grepl("-7-|-7$", sample), "7d",
      default = "Unknown"
    )]
  }

  res[, Full_Group := paste(Treatment, Dose, Duration, sep = "_")]

  bad <- res[Full_Group %like% "Unknown", sample]
  if (length(bad) > 0)
    warning(length(bad), " sample(s) could not be parsed:\n  ", paste(bad, collapse = "\n  "))

  res[]
}

# ======================================
# 5. LOAD DATA
# ======================================
message(">>> Loading CpG data...")
cpg <- fread(CpG_FILE)
cpg[, context := "CpG"]
cpg[, meth_percent := 100 * meth / coverage]
cpg <- cpg[coverage >= 5]  # consistent threshold throughout

meta <- parse_meta(unique(cpg$sample))
message(">>> Sample groups detected:")
print(meta[, .N, by = Full_Group])

cpg <- merge(cpg, meta, by = "sample")
cpg[, Full_Group := factor(Full_Group, levels = GROUP_ORDER)]

load_context <- function(path, ctx) {
  if (!file.exists(path)) {
    message("  [SKIP] ", ctx, " not found: ", path)
    return(NULL)
  }
  dt <- fread(path)
  dt[, context := ctx]
  dt[, meth_percent := 100 * meth / coverage]
  dt <- dt[coverage >= 5]
  dt <- merge(dt, meta, by = "sample")
  dt[, Full_Group := factor(Full_Group, levels = GROUP_ORDER)]
  dt
}

message(">>> Loading CHH and CHG data...")
chh <- load_context(CHH_FILE, "CHH")
chg <- load_context(CHG_FILE, "CHG")

ctx_list <- list(cpg)
if (!is.null(chh)) ctx_list <- c(ctx_list, list(chh))
if (!is.null(chg)) ctx_list <- c(ctx_list, list(chg))
all_ctx <- rbindlist(ctx_list, use.names = TRUE, fill = TRUE)
all_ctx[, Full_Group := factor(Full_Group, levels = GROUP_ORDER)]

# ======================================
# FIGURE 1: PCA on CpG methylation
# ======================================
message(">>> Fig 1: PCA...")

cpg[, site_id := paste(chr, pos, sep = "_")]
pca_wide <- dcast(cpg, sample ~ site_id, value.var = "meth_percent", fun.aggregate = mean)
pca_mat  <- as.matrix(pca_wide[, -1])
rownames(pca_mat) <- pca_wide$sample          # critical: name rows by sample
pca_mat  <- pca_mat[, colSums(is.na(pca_mat)) == 0]  # keep sites in ALL samples

if (ncol(pca_mat) < 2) {
  message("  [WARN] Too few complete sites for PCA - skipping.")
} else {
  pca_res <- prcomp(pca_mat, scale. = TRUE)
  var_exp <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

  pca_df <- data.frame(pca_res$x[, 1:min(4, ncol(pca_res$x))],
                       sample = pca_wide$sample)
  pca_df <- merge(pca_df, meta, by = "sample")
  pca_df$Full_Group <- factor(pca_df$Full_Group, levels = GROUP_ORDER)

  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2,
                               color = Full_Group, shape = Treatment,
                               label = sample)) +
    geom_point(size = 4, alpha = 0.9) +
    geom_text(size = 2.5, vjust = -0.8, hjust = 0.5, show.legend = FALSE) +
    scale_color_manual(values = GROUP_COLOURS, drop = FALSE) +
    scale_shape_manual(values = c(Control = 15, Fipronil = 17, Cyfluthrin = 19)) +
    labs(
      title = "PCA - CpG methylation (all amplicons)",
      x     = paste0("PC1 (", var_exp[1], "%)"),
      y     = paste0("PC2 (", var_exp[2], "%)"),
      color = "Group", shape = "Compound"
    ) +
    theme_pub()

  pdf(file.path(QC_DIR, "Fig1_PCA.pdf"), width = 10, height = 8)
  print(p_pca)
  dev.off()
  message("  Saved: Fig1_PCA.pdf")
}

# ======================================
# FIGURE 2: Sample-sample correlation heatmap
# ======================================
message(">>> Fig 2: Correlation heatmap...")

cor_mat <- cor(t(pca_mat), use = "pairwise.complete.obs")

ann_df <- merge(
  data.frame(sample = rownames(cor_mat), stringsAsFactors = FALSE),
  as.data.frame(meta[, .(sample, Full_Group, Treatment)]),
  by = "sample", all.x = TRUE
)
ann_df <- ann_df[match(rownames(cor_mat), ann_df$sample), ]

row_ann <- data.frame(
  Group    = ann_df$Full_Group,
  Compound = ann_df$Treatment,
  row.names = rownames(cor_mat),
  stringsAsFactors = TRUE
)
row_ann <- row_ann[complete.cases(row_ann), , drop = FALSE]
cor_mat <- cor_mat[rownames(row_ann), rownames(row_ann)]

groups_present    <- levels(droplevels(row_ann$Group))
compounds_present <- levels(droplevels(row_ann$Compound))
all_compound_cols <- c(Control = "#999999", Fipronil = "#1f78b4", Cyfluthrin = "#e31a1c")
extra_groups      <- setdiff(groups_present, names(GROUP_COLOURS))
extra_cols        <- setNames(
  colorRampPalette(c("#b2df8a", "#33a02c", "#cab2d6", "#6a3d9a"))(max(length(extra_groups), 1))[seq_along(extra_groups)],
  extra_groups
)
ann_colours <- list(
  Group    = c(GROUP_COLOURS, extra_cols)[groups_present],
  Compound = all_compound_cols[compounds_present]
)

pdf(file.path(QC_DIR, "Fig2_SampleCorrelation_Heatmap.pdf"), width = 10, height = 9)
pheatmap(
  cor_mat,
  annotation_row    = row_ann,
  annotation_col    = row_ann,
  annotation_colors = ann_colours,
  color             = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
  border_color      = NA,
  main              = "Sample-sample Pearson correlation (CpG methylation)",
  fontsize_row      = 8,
  fontsize_col      = 8
)
dev.off()
message("  Saved: Fig2_SampleCorrelation_Heatmap.pdf")

# ======================================
# FIGURE 3: Coverage distribution per sample per context
# ======================================
message(">>> Fig 3: Coverage distributions...")

pdf(file.path(QC_DIR, "Fig3_Coverage_by_Context.pdf"), width = 12, height = 6)
for (ctx in c("CpG", "CHG", "CHH")) {
  sub <- all_ctx[context == ctx & coverage >= 5]
  if (nrow(sub) == 0) next

  p_cov <- ggplot(sub, aes(x = reorder(sample, coverage, FUN = median),
                             y = coverage, fill = Full_Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85) +
    scale_fill_manual(values = GROUP_COLOURS, drop = FALSE) +
    coord_flip() +
    scale_y_continuous(limits = c(0, quantile(sub$coverage, 0.99))) +
    labs(
      title = paste0(ctx, " - read coverage per sample"),
      x     = "",
      y     = "Read depth",
      fill  = "Group"
    ) +
    theme_pub()
  print(p_cov)
}
dev.off()
message("  Saved: Fig3_Coverage_by_Context.pdf")

# ======================================
# FIGURE 4: Mean methylation per context per group
# ======================================
message(">>> Fig 4: Mean methylation by context...")

ctx_summary <- all_ctx[, .(
  mean_meth   = mean(meth_percent, na.rm = TRUE),
  median_meth = median(meth_percent, na.rm = TRUE),
  n_sites     = .N
), by = .(sample, Full_Group, Treatment, Dose, Duration, context)]

p_ctx <- ggplot(ctx_summary, aes(x = Full_Group, y = mean_meth, fill = Full_Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, aes(color = Full_Group)) +
  scale_fill_manual(values  = GROUP_COLOURS, drop = FALSE) +
  scale_color_manual(values = GROUP_COLOURS, drop = FALSE) +
  facet_wrap(~context, scales = "free_y", ncol = 1) +
  theme_pub() +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(
    title = "Mean methylation % by methylation context and treatment group",
    x     = "",
    y     = "Mean methylation (%)"
  )

pdf(file.path(QC_DIR, "Fig4_MeanMethylation_by_Context.pdf"), width = 10, height = 10)
print(p_ctx)
dev.off()
message("  Saved: Fig4_MeanMethylation_by_Context.pdf")

# ======================================
# FIGURE 5: Covered sites per context per sample
# ======================================
message(">>> Fig 5: Site counts per context...")

site_counts <- all_ctx[, .(n_sites = uniqueN(paste(chr, pos))),
                        by = .(sample, Full_Group, context)]

p_sites <- ggplot(site_counts, aes(x = reorder(sample, n_sites),
                                    y = n_sites, fill = Full_Group)) +
  geom_col(alpha = 0.85) +
  scale_fill_manual(values = GROUP_COLOURS, drop = FALSE) +
  coord_flip() +
  facet_wrap(~context, scales = "free_x") +
  labs(
    title = "Number of covered CpG sites (>=5x) per sample",
    x     = "",
    y     = "N sites",
    fill  = "Group"
  ) +
  theme_pub()

pdf(file.path(QC_DIR, "Fig5_SiteCount_by_Context.pdf"), width = 12, height = 7)
print(p_sites)
dev.off()
message("  Saved: Fig5_SiteCount_by_Context.pdf")

# ======================================
# FIGURE 6: Per-gene methylation profiles - one PDF per gene
# ======================================
message(">>> Fig 6: Per-gene methylation profiles...")

genes <- sort(unique(cpg$gene_name[!is.na(cpg$gene_name) & cpg$gene_name != ""]))

for (gn in genes) {
  g_dt <- cpg[gene_name == gn]
  if (nrow(g_dt) == 0) next

  pos_map <- data.table(pos = sort(unique(g_dt$pos)))[, cpg_index := .I]
  g_dt    <- merge(g_dt, pos_map, by = "pos")

  p_line <- ggplot(g_dt, aes(x = cpg_index, y = meth_percent,
                               group = sample, color = Full_Group)) +
    geom_line(alpha = 0.6, linewidth = 0.7) +
    geom_point(size = 1.8, alpha = 0.8) +
    scale_color_manual(values = GROUP_COLOURS, drop = FALSE) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    facet_wrap(~Full_Group, ncol = 3, drop = FALSE) +
    ylim(0, 100) +
    theme_pub() +
    theme(legend.position = "none") +
    labs(
      title = paste0(gn, " - CpG methylation per sample"),
      x     = "CpG index",
      y     = "Methylation (%)"
    )

  pdf(file.path(QC_DIR, paste0("Fig6_", gn, "_MethylationProfile.pdf")),
      width = 14, height = 8)
  print(p_line)
  dev.off()
}
message("  Saved: Fig6_*_MethylationProfile.pdf (one per gene)")

# ======================================
# FIGURE 7: Non-CpG methylation (CHH) per sample
# Note: CHH methylation in vertebrate somatic tissue is biologically real
# and NOT a reliable bisulfite conversion proxy without spike-in controls.
# This plot shows non-CpG methylation levels for characterisation purposes.
# ======================================
if (!is.null(chh)) {
  message(">>> Fig 7: Non-CpG (CHH) methylation...")

  conv_proxy <- chh[, .(
    mean_CHH   = mean(meth_percent, na.rm = TRUE),
    median_CHH = median(meth_percent, na.rm = TRUE)
  ), by = .(sample, Full_Group)]

  p_conv <- ggplot(conv_proxy, aes(x = reorder(sample, mean_CHH),
                                    y = mean_CHH, fill = Full_Group)) +
    geom_col(alpha = 0.85) +
    scale_fill_manual(values = GROUP_COLOURS, drop = FALSE) +
    coord_flip() +
    labs(
      title = "Non-CpG methylation (CHH context) per sample",
      x     = "",
      y     = "Mean CHH methylation (%)",
      fill  = "Group"
    ) +
    theme_pub()

  pdf(file.path(QC_DIR, "Fig7_CHH_Methylation.pdf"), width = 10, height = 7)
  print(p_conv)
  dev.off()
  message("  Saved: Fig7_CHH_Methylation.pdf")
}

# ======================================
# FIGURE 8: Per-amplicon coverage uniformity
# ======================================
message(">>> Fig 8: Amplicon coverage uniformity...")

amp_cov <- cpg[, .(
  mean_cov   = mean(coverage),
  median_cov = as.double(median(coverage)),
  pct_ge10   = mean(coverage >= 10) * 100
), by = .(sample, amplicon, Full_Group)]

p_amp <- ggplot(amp_cov, aes(x = amplicon, y = median_cov, fill = Full_Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, aes(color = Full_Group), alpha = 0.6) +
  scale_fill_manual(values  = GROUP_COLOURS, drop = FALSE) +
  scale_color_manual(values = GROUP_COLOURS, drop = FALSE) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_y_log10(labels = scales::comma) +
  coord_flip() +
  labs(
    title    = "Read coverage per amplicon (log scale)",
    subtitle = "Dashed line = 10x minimum threshold for DSS",
    x        = "",
    y        = "Median read depth (log10)",
    fill     = "Group", color = "Group"
  ) +
  theme_pub()

pdf(file.path(QC_DIR, "Fig8_AmpliconCoverage.pdf"), width = 12, height = 7)
print(p_amp)
dev.off()
message("  Saved: Fig8_AmpliconCoverage.pdf")

# ======================================
# FIGURE 9: Methylation distribution per sample (violin + density)
# Shows bimodality expected in CpG methylation (most sites ~0% or ~100%)
# ======================================
message(">>> Fig 9: Methylation distribution per sample...")

p_violin <- ggplot(cpg, aes(x = reorder(sample, meth_percent, FUN = median),
                              y = meth_percent, fill = Full_Group)) +
  geom_violin(alpha = 0.7, scale = "width", linewidth = 0.3) +
  geom_boxplot(width = 0.08, outlier.shape = NA, fill = "white", alpha = 0.8) +
  scale_fill_manual(values = GROUP_COLOURS, drop = FALSE) +
  coord_flip() +
  labs(
    title = "CpG methylation distribution per sample",
    x     = "",
    y     = "Methylation (%)",
    fill  = "Group"
  ) +
  theme_pub()

pdf(file.path(QC_DIR, "Fig9_MethylationDistribution.pdf"), width = 10, height = 8)
print(p_violin)
dev.off()
message("  Saved: Fig9_MethylationDistribution.pdf")

# ======================================
# FIGURE 10: Within-group reproducibility
# Coefficient of variation (CV) per CpG site across samples in same group.
# Low CV = consistent methylation across biological replicates.
# ======================================
message(">>> Fig 10: Within-group reproducibility (CV per site)...")

repro <- cpg[, .(
  mean_meth = mean(meth_percent, na.rm = TRUE),
  sd_meth   = sd(meth_percent,   na.rm = TRUE),
  n         = .N
), by = .(amplicon, pos, Full_Group)]

repro[, cv := ifelse(mean_meth > 0, (sd_meth / mean_meth) * 100, NA_real_)]
repro <- repro[n >= 2]  # need at least 2 samples to compute SD
repro[, Full_Group := factor(Full_Group, levels = GROUP_ORDER)]

# Shortened x-axis labels to prevent overlap in facets
repro[, Group_short := gsub("Control_Control_", "Ctrl_", as.character(Full_Group))]
repro[, Group_short := gsub("Fipronil_",   "Fip_",  Group_short)]
repro[, Group_short := gsub("Cyfluthrin_", "Cyt_",  Group_short)]
repro[, Group_short := factor(Group_short,
        levels = gsub("Control_Control_", "Ctrl_",
                 gsub("Fipronil_",   "Fip_",
                 gsub("Cyfluthrin_", "Cyt_", GROUP_ORDER))))]

p_cv <- ggplot(repro, aes(x = Group_short, y = cv, fill = Full_Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.4) +
  scale_fill_manual(values = GROUP_COLOURS, drop = FALSE, guide = "none") +
  facet_wrap(~amplicon, ncol = 3, scales = "free_y") +
  coord_cartesian(ylim = NULL) +
  theme_pub() +
  theme(axis.text.x  = element_text(angle = 55, hjust = 1, size = 7),
        legend.position = "none") +
  labs(
    title    = "Within-group reproducibility per CpG site",
    subtitle = "Coefficient of variation (CV%) across biological replicates",
    x        = "",
    y        = "CV (%)"
  )

pdf(file.path(QC_DIR, "Fig10_WithinGroup_CV.pdf"), width = 14, height = 8)
print(p_cv)
dev.off()
message("  Saved: Fig10_WithinGroup_CV.pdf")

# ======================================
# FIGURE 11: Sample-level outlier detection
# Z-score of mean methylation per sample per amplicon.
# Samples with |Z| > 2 in any amplicon are flagged.
# ======================================
message(">>> Fig 11: Outlier detection...")

sample_amp_meth <- cpg[, .(
  mean_meth = mean(meth_percent, na.rm = TRUE)
), by = .(sample, amplicon, Full_Group)]

# Z-score within each amplicon
sample_amp_meth[, z_score := scale(mean_meth)[, 1], by = amplicon]
sample_amp_meth[, outlier := abs(z_score) > 2]
sample_amp_meth[, Full_Group := factor(Full_Group, levels = GROUP_ORDER)]

# Print flagged samples to console
flagged <- sample_amp_meth[outlier == TRUE]
if (nrow(flagged) > 0) {
  message("  [WARN] Potential outlier samples (|Z| > 2):")
  print(flagged[, .(sample, amplicon, Full_Group, mean_meth = round(mean_meth, 1), z_score = round(z_score, 2))])
} else {
  message("  No outliers detected (|Z| > 2 threshold)")
}

p_outlier <- ggplot(sample_amp_meth,
                     aes(x = amplicon, y = z_score,
                         color = Full_Group, shape = outlier, label = ifelse(outlier, sample, ""))) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0,        linetype = "solid",  color = "grey80") +
  geom_jitter(width = 0.2, size = 3, alpha = 0.85) +
  ggplot2::geom_text(size = 2.5, vjust = -1, hjust = 0.5, show.legend = FALSE) +
  scale_color_manual(values = GROUP_COLOURS, drop = FALSE) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 8),
                     labels = c("FALSE" = "Normal", "TRUE" = "Outlier (|Z|>2)"),
                     name   = "") +
  labs(
    title    = "Sample outlier detection per amplicon",
    subtitle = "Z-score of mean CpG methylation; dashed lines = +/-2 SD",
    x        = "Amplicon",
    y        = "Z-score",
    color    = "Group"
  ) +
  theme_pub()

pdf(file.path(QC_DIR, "Fig11_OutlierDetection.pdf"), width = 12, height = 7)
print(p_outlier)
dev.off()
message("  Saved: Fig11_OutlierDetection.pdf")

# ======================================
# SUMMARY TABLE - for Methods section
# ======================================
message(">>> Writing summary table...")

summary_tbl <- all_ctx[, .(
  n_sites    = uniqueN(paste(chr, pos)),
  mean_cov   = round(mean(coverage), 1),
  median_cov = as.double(median(coverage)),
  mean_meth  = round(mean(meth_percent, na.rm = TRUE), 2),
  sd_meth    = round(sd(meth_percent,   na.rm = TRUE), 2)
), by = .(sample, Full_Group, amplicon, context)]

fwrite(summary_tbl, file.path(QC_DIR, "QC_summary_table.tsv"), sep = "\t")
message("  Saved: QC_summary_table.tsv")

# Outlier summary
fwrite(sample_amp_meth[, .(sample, amplicon, Full_Group, mean_meth = round(mean_meth,1),
                            z_score = round(z_score,2), outlier)],
       file.path(QC_DIR, "QC_outlier_zscore.tsv"), sep = "\t")
message("  Saved: QC_outlier_zscore.tsv")

message("\n>>> QC complete. All outputs in: ", QC_DIR)
message("    Fig1_PCA.pdf")
message("    Fig2_SampleCorrelation_Heatmap.pdf")
message("    Fig3_Coverage_by_Context.pdf")
message("    Fig4_MeanMethylation_by_Context.pdf")
message("    Fig5_SiteCount_by_Context.pdf")
message("    Fig6_*_MethylationProfile.pdf")
if (!is.null(chh)) message("    Fig7_CHH_Methylation.pdf")
message("    Fig8_AmpliconCoverage.pdf")
message("    Fig9_MethylationDistribution.pdf")
message("    Fig10_WithinGroup_CV.pdf")
message("    Fig11_OutlierDetection.pdf")
message("    QC_summary_table.tsv")
message("    QC_outlier_zscore.tsv")
