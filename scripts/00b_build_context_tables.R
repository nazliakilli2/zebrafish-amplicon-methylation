#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# =============================================================================
# 00b_build_context_tables.R
#
# Aggregates the per-sample Bismark methylation-extractor output into the
# combined per-context count tables that scripts/01_filter_to_amplicons.R
# consumes. This is the bridge between:
#   scripts/00_bismark_align_extract.sh   (per-sample *_context_*.txt.gz files)
#   scripts/01_filter_to_amplicons.R      (reads data/tables_fixed_final/all_sites_*.tsv)
#
# For every cytosine it sums methylated / unmethylated calls per genomic
# position and per sample, then writes one table per sequence context.
# =============================================================================

# ======================================
# 1. PATHS
# ======================================
source(file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "config.R"))

# Directory holding the per-sample Bismark context files
# (CpG_context_*.txt.gz / CHG_context_*.txt.gz / CHH_context_*.txt.gz).
# Defaults to data/bismark_extraction; override with BISMARK_DIR if the
# extraction output lives elsewhere.
IN_DIR  <- Sys.getenv("BISMARK_DIR", unset = file.path(DATA_DIR, "bismark_extraction"))

# Combined tables land here; this matches BASE_IN in 01_filter_to_amplicons.R.
OUT_DIR <- file.path(DATA_DIR, "tables_fixed_final")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ======================================
# 2. FIND CONTEXT FILES
# ======================================
# Matches the Bismark methylation-extractor naming style, e.g.
#   CpG_context_S13-7-18-7_..._bismark_bt2_pe.txt.gz
files <- list.files(
  path = IN_DIR,
  pattern = "context_.*_bismark_bt2.*\\.txt\\.gz$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(files) == 0) {
  stop(paste("No context files found in:", IN_DIR,
             "\nCheck that Bismark has finished, or set BISMARK_DIR to the extraction folder."))
}

cat("Found", length(files), "context files. Starting aggregation...\n")

# ======================================
# 3. PARSER FUNCTION
# ======================================
parse_context_file <- function(file_path) {
  fname <- basename(file_path)

  # Determine context
  ctx <- if (grepl("CpG_context", fname)) "CpG"
         else if (grepl("CHG_context", fname)) "CHG"
         else if (grepl("CHH_context", fname)) "CHH"
         else "Unknown"

  # Clean Sample Name: Removes the context prefix and the Bismark suffixes
  # Example: "CHH_context_S13-7-18-7_..._bismark_bt2_pe.txt.gz" -> "S13-7-18-7"
  sname <- gsub("^(CpG|CHG|CHH)_context_", "", fname)
  sname <- gsub("_S[0-9]+_R[0-9]+_.*", "", sname) # Strips sequencer tags
  sname <- gsub("\\.namesorted.*", "", sname)    # Strips sorting tags

  message("  Processing Sample: ", sname, " [Context: ", ctx, "]")

  # Read raw Bismark extraction (Skip the 1-line header)
  # Columns: 1:ReadID, 2:Strand, 3:Chr, 4:Pos, 5:Call
  dt <- fread(file_path, sep = "\t", skip = 1, header = FALSE,
              col.names = c("id", "strand", "chr", "pos", "call"))

  if (nrow(dt) == 0) return(NULL)

  # Aggregate counts per genomic position
  # Uppercase (Z,C,H,X) = Methylated; Lowercase (z,c,h,x) = Unmethylated
  res <- dt[, .(
    meth     = sum(grepl("[ZCHX]", call)),
    unmeth   = sum(grepl("[zchx]", call)),
    coverage = .N
  ), by = .(chr, pos)]

  res[, `:=`(
    percent = 100 * (meth / coverage),
    sample  = sname,
    context = ctx
  )]

  return(res)
}

# ======================================
# 4. EXECUTION
# ======================================
all_list <- lapply(files, parse_context_file)
all_sites <- rbindlist(all_list, use.names = TRUE, fill = TRUE)

if (nrow(all_sites) == 0) {
  stop("Aggregation produced 0 rows. Check if input files are empty.")
}

setorder(all_sites, chr, pos)

# ======================================
# 5. SAVE OUTPUTS
# ======================================
cat("\n=== Saving Results to:", OUT_DIR, "===\n")

# A. Save Context-Specific Tables (all_sites_CpG.tsv, etc.)
for (ctx in unique(all_sites$context)) {
  fwrite(all_sites[context == ctx],
         file.path(OUT_DIR, paste0("all_sites_", ctx, ".tsv")),
         sep = "\t")
}

# B. Save CpG-only count table (chr, pos, sample, meth, coverage)
fwrite(all_sites[context == "CpG", .(chr, pos, sample, meth, coverage)],
       file.path(OUT_DIR, "all_sites_CpG_counts.tsv"),
       sep = "\t")

# C. Save QC Summary
# Use this to verify CHH/CHG methylation is low (bisulfite conversion efficiency).
qc_summary <- all_sites[, .(
  mean_meth   = mean(percent),
  median_meth = median(percent),
  mean_cov    = mean(coverage),
  total_sites = .N
), by = .(sample, context)]

fwrite(qc_summary, file.path(OUT_DIR, "qc_summary.tsv"), sep = "\t")

cat("\n>>> DONE. Please check qc_summary.tsv for conversion efficiency.\n")
