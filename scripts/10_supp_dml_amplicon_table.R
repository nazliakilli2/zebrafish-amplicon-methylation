#!/usr/bin/env Rscript

# ============================================================
# SUPPLEMENTARY TABLE: per-amplicon DML counts
#   For each amplicon: how many CpG sites were TESTED, and how many are
#   significantly HYPO- or HYPER-methylated (differentially methylated loci,
#   DMLs). These are the significance-based numbers behind the main-text
#   Table 2 / Results (e.g. "582 tested CpG sites, 178 DMLs, 99% hypo").
#
# Standalone. Reads ONLY the per-CpG DSS output (DSS_DML_results.tsv) already
# produced and validated by DSS_v7.R. Re-runs no model, modifies no result,
# writes new table files only.
#
# DIRECTION CONVENTION (must match DSS_v7.R):
#   DMLtest() was called with group1 = treatment, group2 = control, so DSS
#   reports  diff = mu1 - mu2 = (treatment methylation) - (control methylation).
#     diff < 0  ->  treatment LOWER than control  ->  HYPO-methylation
#     diff > 0  ->  treatment HIGHER than control ->  HYPER-methylation
#   A CpG is a DML when fdr < FDR_THRESH (BH-adjusted within each comparison,
#   exactly as DSS_v7.R computed it).
#
# COUNTING UNIT:
#   "Tested" counts CpG x comparison tests (one CpG contributes once per
#   comparison it is tested in). Summed over the four 7-day comparisons this
#   reproduces the main-text totals (e.g. hsp70 cluster 66 CpGs x 4 = 264
#   tests; 4 comparisons x all amplicons = 582). Effect sizes (diff) are on the
#   0-1 methylation scale in DSS and are reported here in percentage points.
#
# Outputs (DSS_results_v7/):
#   DSS_DML_amplicon_table.tsv          - per amplicon (summed over comparisons)
#                                         + a genomic hsp70-cluster row + TOTAL
#   DSS_DML_amplicon_by_comparison.tsv  - per amplicon x comparison breakdown
# ============================================================

suppressPackageStartupMessages(library(data.table))

# ======================================
# 1. PATHS & SETTINGS
# ======================================
source(file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "config.R"))
OUT_DIR     <- RESULTS_DIR
INPUT_FILE  <- DSS_RESULTS

FDR_THRESH  <- 0.05   # a CpG is a DML when fdr < FDR_THRESH (matches DSS_v7.R)

# Restrict to the comparisons reported in the paper. Set to NULL to summarise
# every comparison present in the file instead.
COMPARISONS <- c('Fipronil_Low_7d','Fipronil_High_7d',
                 'Cyfluthrin_Low_7d','Cyfluthrin_High_7d')

GENE_ORDER  <- c('hsp70.1','hsp70.2','hsp70.3','actb2','b2m','cyp1a')

# ======================================
# 2. LOAD
# ======================================
if (!file.exists(INPUT_FILE)) stop("DSS DML results not found: ", INPUT_FILE)
message(">>> Loading DSS per-CpG results: ", INPUT_FILE)
dml <- fread(INPUT_FILE)

need <- c("amplicon","comparison","diff","fdr")
miss <- setdiff(need, names(dml))
if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))

# a handful of CpGs can have NA fdr (model returned NA); they were tested but
# cannot be called significant -- keep them in "tested", never in hypo/hyper.
n_na_fdr <- dml[is.na(fdr), .N]
if (n_na_fdr > 0)
  message("    Note: ", n_na_fdr, " tested CpG(s) have NA fdr; counted as tested, not as DML.")

if (!is.null(COMPARISONS)) {
  missing_comp <- setdiff(COMPARISONS, unique(dml$comparison))
  if (length(missing_comp) > 0)
    warning("Requested comparison(s) not in file: ", paste(missing_comp, collapse = ", "))
  dml <- dml[comparison %in% COMPARISONS]
}
if (nrow(dml) == 0) stop("No rows left after comparison filter. Check COMPARISONS.")

# order amplicons; keep any unexpected amplicon name at the end rather than dropping it
lev <- c(GENE_ORDER, setdiff(unique(dml$amplicon), GENE_ORDER))
dml[, amplicon := factor(amplicon, levels = lev)]

# ======================================
# 3. COUNTER
#    tested / hypo / hyper / DML counts + DML effect-size summary (pp)
# ======================================
summarise_dml <- function(d, by) {
  d[, {
    is_dml  <- !is.na(fdr) & fdr < FDR_THRESH
    is_hypo <- is_dml & diff < 0
    is_hyp  <- is_dml & diff > 0
    dml_pp  <- 100 * diff[is_dml]          # effect sizes of significant CpGs (pp)
    .(
      n_tested   = .N,
      n_DML      = sum(is_dml),
      n_hypo     = sum(is_hypo),
      n_hyper    = sum(is_hyp),
      pct_DML    = round(100 * sum(is_dml)  / .N, 1),
      pct_hypo   = round(100 * sum(is_hypo) / .N, 1),
      pct_hyper  = round(100 * sum(is_hyp)  / .N, 1),
      # effect sizes among DMLs only (percentage points); NA when no DML
      mean_delta_pp   = if (any(is_dml)) round(mean(dml_pp), 2)   else NA_real_,
      median_delta_pp = if (any(is_dml)) round(median(dml_pp), 2) else NA_real_,
      min_delta_pp    = if (any(is_dml)) round(min(dml_pp), 2)    else NA_real_,
      max_delta_pp    = if (any(is_dml)) round(max(dml_pp), 2)    else NA_real_
    )
  }, by = by]
}

# ---- (a) per amplicon x comparison ------------------------------------------
by_comp <- summarise_dml(dml, c("amplicon","comparison"))
setorder(by_comp, amplicon, comparison)
fwrite(by_comp, file.path(OUT_DIR, "DSS_DML_amplicon_by_comparison.tsv"), sep = "\t")
message("    Wrote: ", file.path(OUT_DIR, "DSS_DML_amplicon_by_comparison.tsv"))

# ---- (b) per amplicon (summed over the comparisons above) -------------------
per_amp <- summarise_dml(dml, "amplicon")
setorder(per_amp, amplicon)
per_amp[, amplicon := as.character(amplicon)]

# genomic hsp70 cluster row (the three tandem hsp70 amplicons pooled), matching
# the way the main text groups them; only added if those amplicons are present.
hsp <- c('hsp70.1','hsp70.2','hsp70.3')
if (all(hsp %in% dml$amplicon)) {
  cl <- summarise_dml(dml[amplicon %in% hsp], character(0))
  cl[, amplicon := "hsp70_cluster (1+2+3)"]
  setcolorder(cl, names(per_amp))
} else cl <- NULL

# grand total across all amplicons
tot <- summarise_dml(dml, character(0))
tot[, amplicon := "TOTAL"]
setcolorder(tot, names(per_amp))

per_amp_out <- rbindlist(c(list(per_amp), if (!is.null(cl)) list(cl), list(tot)),
                         use.names = TRUE)
fwrite(per_amp_out, file.path(OUT_DIR, "DSS_DML_amplicon_table.tsv"), sep = "\t")
message("    Wrote: ", file.path(OUT_DIR, "DSS_DML_amplicon_table.tsv"))

# ======================================
# 4. CONSOLE SUMMARY (self-check on a local run)
# ======================================
message("\n>>> Comparisons included: ",
        if (is.null(COMPARISONS)) "ALL in file" else paste(COMPARISONS, collapse = ", "))
message(">>> DML defined as fdr < ", FDR_THRESH,
        "  (diff<0 = hypomethylation, diff>0 = hypermethylation)\n")
message(">>> Per-amplicon DML table (summed over comparisons):")
print(per_amp_out)
message("\n>>> Per amplicon x comparison:")
print(by_comp)

message("\n>>> Done.")
message("    Tables: ", file.path(OUT_DIR, "DSS_DML_amplicon_table.tsv"))
message("            ", file.path(OUT_DIR, "DSS_DML_amplicon_by_comparison.tsv"))
