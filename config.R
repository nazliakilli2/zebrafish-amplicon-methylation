# =============================================================================
# config.R  --  central paths and settings for the amplicon methylation pipeline
#
# Every script in scripts/ sources this file, so all paths live in one place
# and no absolute, machine-specific path appears anywhere in the code.
#
# Run scripts from the repository root, e.g.:
#
#     Rscript scripts/02_dss_dmltest.R
#
# By default inputs are read from ./data and ./reference, and outputs are
# written to ./results and ./figures (all relative to the repository root).
# To keep data and outputs somewhere else, point PROJECT_ROOT at that location
# before running:
#
#     PROJECT_ROOT=/path/to/project Rscript scripts/02_dss_dmltest.R
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = ".")

# ---- directories ----
DATA_DIR    <- file.path(PROJECT_ROOT, "data")       # input count tables
REF_DIR     <- file.path(PROJECT_ROOT, "reference")  # amplicon coordinates
RESULTS_DIR <- file.path(PROJECT_ROOT, "results")    # DSS + summary tables
FIGURE_DIR  <- file.path(PROJECT_ROOT, "figures")    # publication figures
QC_DIR      <- file.path(RESULTS_DIR, "qc")          # QC report outputs

# ---- reference ----
# amplicon coordinates, columns: chr  amplicon  start  end
BED_FILE <- file.path(REF_DIR, "master_amplicons.bed")

# ---- inputs: per-context, amplicon- and coverage-filtered count tables ----
# produced by scripts/01_filter_to_amplicons.R
CPG_FILTERED <- file.path(DATA_DIR, "all_sites_CpG_counts.master_bed.filtered.tsv")
CHH_FILTERED <- file.path(DATA_DIR, "all_sites_CHH_counts.master_bed.filtered.tsv")
CHG_FILTERED <- file.path(DATA_DIR, "all_sites_CHG_counts.master_bed.filtered.tsv")

# ---- key result: per-CpG DSS DMLtest output ----
# produced by scripts/02_dss_dmltest.R
DSS_RESULTS <- file.path(RESULTS_DIR, "DSS_DML_results.tsv")

# ---- ensure output directories exist ----
for (.d in c(RESULTS_DIR, FIGURE_DIR, QC_DIR)) {
  dir.create(.d, showWarnings = FALSE, recursive = TRUE)
}
rm(.d)
