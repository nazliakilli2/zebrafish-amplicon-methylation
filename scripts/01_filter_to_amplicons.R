#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ======================================
# PATHS
# ======================================
source(file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "config.R"))
BASE_IN  <- file.path(DATA_DIR, "tables_fixed_final")  # raw per-context extraction tables
OUT_DIR  <- DATA_DIR                                   # filtered tables land here
# BED_FILE is provided by config.R

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ======================================
# LOAD BED (amplicon boundaries)
# ======================================
bed <- fread(BED_FILE, col.names = c("chr", "amplicon", "start", "end"))
bed[, chr   := as.character(chr)]
bed[, start := as.integer(start)]
bed[, end   := as.integer(end)]
setkey(bed, chr, start, end)
cat("Loaded BED:", nrow(bed), "amplicon regions\n")

# ======================================
# PROCESS EACH CONTEXT
# ======================================
contexts <- c("CpG", "CHH", "CHG")

for (ctx in contexts) {

  in_file  <- file.path(BASE_IN,  paste0("all_sites_", ctx, ".tsv"))
  out_file <- file.path(OUT_DIR,  paste0("all_sites_", ctx, "_counts.master_bed.tsv"))
  out_filt <- file.path(OUT_DIR,  paste0("all_sites_", ctx, "_counts.master_bed.filtered.tsv"))

  if (!file.exists(in_file)) {
    cat("  [SKIP]", ctx, "- file not found:", in_file, "\n")
    next
  }

  cat("\n>>> Processing context:", ctx, "\n")
  dt <- fread(in_file)

  # Rename 'percent' to 'meth_percent' for consistency if present
  if ("percent" %in% names(dt) && !"meth_percent" %in% names(dt)) {
    setnames(dt, "percent", "meth_percent")
  }

  # Compute meth_percent if not present
  if (!"meth_percent" %in% names(dt)) {
    dt[, meth_percent := fifelse(coverage > 0, 100 * meth / coverage, NA_real_)]
  }

  # BED-style coordinates for foverlaps
  dt[, chr  := as.character(chr)]
  dt[, pos0 := as.integer(pos) - 1L]
  dt[, pos1 := pos0]
  setkey(dt, chr, pos0, pos1)

  # Overlap with amplicon BED
  filtered <- foverlaps(
    dt, bed,
    by.x    = c("chr", "pos0", "pos1"),
    by.y    = c("chr", "start", "end"),
    nomatch = 0L
  )

  # amplicon column doubles as gene_name for downstream scripts
  filtered[, gene_name := amplicon]

  # Clean up helper columns
  filtered[, c("pos0", "pos1") := NULL]

  cat("  Sites after BED overlap:", nrow(filtered), "\n")
  fwrite(filtered, out_file, sep = "\t")
  cat("  Wrote:", out_file, "\n")

  # ---- TABSAT-style adaptive coverage filter ----
  filtered <- filtered[coverage >= 5]
  filtered <- filtered[, {
    m      <- mean(coverage, na.rm = TRUE)
    s      <- sd(coverage,   na.rm = TRUE)
    cutoff <- max(5, floor(m - s))
    .SD[coverage >= cutoff]
  }, by = .(sample, amplicon)]

  cat("  Sites after coverage filter:", nrow(filtered), "\n")
  fwrite(filtered, out_filt, sep = "\t")
  cat("  Wrote:", out_filt, "\n")
}

cat("\n>>> DONE. All context files written to:", OUT_DIR, "\n")
cat(">>> Update QC script paths to:\n")
cat("    CHH_FILE <- \"", file.path(OUT_DIR, "all_sites_CHH_counts.master_bed.filtered.tsv"), "\"\n", sep="")
cat("    CHG_FILE <- \"", file.path(OUT_DIR, "all_sites_CHG_counts.master_bed.filtered.tsv"), "\"\n", sep="")
