#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(DSS)
  library(bsseq)
})

# ======================================
# 1. PATHS & SETTINGS
# ======================================
source(file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "config.R"))
INPUT_FILE <- CPG_FILTERED
OUT_DIR    <- RESULTS_DIR

# FROM v6: 10x is more conservative and appropriate for DSS smoothing
MIN_COV <- 10

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

message(">>> Loading master BED filtered counts...")
raw <- fread(INPUT_FILE)

# FROM v5: hard stop if expected columns are missing — catches wrong input file immediately
required_cols <- c("chr", "pos", "sample", "meth", "coverage", "amplicon")
missing_cols  <- setdiff(required_cols, names(raw))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

raw <- raw[coverage >= MIN_COV]

# ======================================
# 2. METADATA PARSER
# ======================================
# FROM v5: sample-prefix-aware patterns handle both F* and S* naming conventions
parse_sample_metadata <- function(sample_names) {
  dt <- data.table(sample = sample_names)

  dt[, compound := fcase(
    grepl("^F", sample) & grepl("-K-", sample),                                                "Control",
    grepl("^F", sample) & grepl("-170-", sample),                                              "Fipronil",
    grepl("^F", sample) & grepl("-1700-", sample),                                             "Fipronil",
    grepl("^S", sample) & grepl("-K-", sample),                                                "Control",
    grepl("^S", sample) & grepl("-18-|18-", sample) & !grepl("-180-|-180$", sample),           "Cyfluthrin",
    grepl("^S", sample) & grepl("-180-|-180$", sample),                                        "Cyfluthrin",
    default = "Unknown"
  )]

  dt[, dose := fcase(
    grepl("-K-", sample),                                                                       "Control",
    grepl("-170-", sample) & !grepl("-1700-", sample),                                         "Low",
    grepl("-1700-", sample),                                                                    "High",
    grepl("^S", sample) & grepl("-18-|18-", sample) & !grepl("-180-|-180$|180-", sample),      "Low",
    grepl("^S", sample) & grepl("-180-|-180$|180-", sample),                                   "High",
    default = "Unknown"
  )]

  dt[, timepoint := fcase(
    grepl("^F.*-K-24",    sample),                                     "24h",
    grepl("^F.*-K-96",    sample),                                     "96h",
    grepl("^F.*-K-7",     sample),                                     "7d",
    grepl("^F.*-170-7|^F.*-1700-7", sample),                          "7d",
    grepl("^F.*-1700-24", sample),                                     "24h",
    grepl("^S.*-K-96",    sample),                                     "96h",
    grepl("^S.*-K-7",     sample),                                     "7d",
    grepl("^S.*-96$|-96-", sample),                                    "96h",
    grepl("^S.*-7$", sample) | grepl("-7-[0-9]|[0-9]-7$", sample),   "7d",
    default = "Unknown"
  )]

  # FROM v5: second-pass fallback for any S* samples still Unknown after main parser
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

meta <- parse_sample_metadata(unique(raw$sample))

# FROM v5: warn and exclude any samples that could not be parsed
unknown_samples <- meta[group == "Unknown", sample]
if (length(unknown_samples) > 0) {
  warning(
    length(unknown_samples), " sample(s) could not be parsed and will be excluded:\n  ",
    paste(unknown_samples, collapse = "\n  ")
  )
}
meta <- meta[group != "Unknown"]

raw <- merge(raw, meta, by = "sample")

message(">>> Sample groups detected:")
print(meta[, .N, by = group])

# ======================================
# 3. COMPARISON TARGETS
# ======================================
targets <- c(
  "Fipronil_Low_7d",
  "Fipronil_High_7d",
  "Cyfluthrin_Low_96h",
  "Cyfluthrin_High_96h",
  "Cyfluthrin_Low_7d",
  "Cyfluthrin_High_7d"
)

all_results <- list()

# ======================================
# 4. DSS LOOP
# FDR applied after, globally per comparison (across all amplicons)
# ======================================
for (amp in unique(raw$amplicon)) {
  message(">>> Processing amplicon: ", amp)

  dg           <- raw[amplicon == amp]
  sample_names <- unique(dg$sample)

  # FROM v6: skip amplicons with too few CpG sites for meaningful smoothing
  if (uniqueN(dg$pos) < 5) {
    message("  Skipping ", amp, ": too few sites (", uniqueN(dg$pos), ")")
    next
  }

  # FROM v5: build BSseq once per amplicon using all samples present
  sample_list <- lapply(sample_names, function(s) {
    as.data.frame(dg[sample == s, .(chr, pos, N = coverage, X = meth)])
  })
  names(sample_list) <- sample_names
  BSobj <- makeBSseqData(sample_list, sampleNames = sample_names)

  for (tr in targets) {
    tp_match <- sub(".*_", "", tr)
    ctrl_grp <- paste0("Control_", tp_match)

    s_treat <- meta[group == tr          & sample %in% sample_names, sample]
    s_ctrl  <- meta[group == ctrl_grp   & sample %in% sample_names, sample]

    if (length(s_treat) >= 2 && length(s_ctrl) >= 2) {
      message("  Testing: ", tr, " vs ", ctrl_grp,
              " (treat n=", length(s_treat), ", ctrl n=", length(s_ctrl), ")")

      # FROM v6: smoothing=TRUE with span=200 bp, appropriate for amplicon size
      dml_test <- tryCatch(
        DMLtest(
          BSobj,
          group1         = s_treat,
          group2         = s_ctrl,
          smoothing      = TRUE,
          smoothing.span = 200
        ),
        error = function(e) {
          message("    Skipping due to error: ", e$message)
          NULL
        }
      )

      if (!is.null(dml_test)) {
        dml_test             <- as.data.table(dml_test)
        dml_test[, amplicon   := amp]
        dml_test[, comparison := tr]
        all_results[[paste(amp, tr, sep = "__")]] <- dml_test
      }

    } else {
      message("  Skipping: ", tr,
              " (treat n=", length(s_treat), ", ctrl n=", length(s_ctrl), ")")
    }
  }
}

# ======================================
# 5. FDR CORRECTION & SAVE
# BH correction within each comparison, across all amplicons.
# Each comparison is an independent biological contrast — this is
# more honest than correcting 5 CpGs at a time per amplicon.
# ======================================
if (length(all_results) == 0) {
  stop("No successful DSS comparisons were produced.")
}

final_dml <- rbindlist(all_results, use.names = TRUE, fill = TRUE)
final_dml[, fdr := p.adjust(pval, method = "BH"), by = comparison]

fwrite(final_dml, file.path(OUT_DIR, "DSS_DML_results.tsv"), sep = "\t")

message("\n>>> DONE")
message(">>> Results saved to: ", file.path(OUT_DIR, "DSS_DML_results.tsv"))
message(">>> Total CpG tests across all comparisons: ", nrow(final_dml))
message(">>> Significant (FDR < 0.05) per comparison:")
print(final_dml[, .(n_sig = sum(fdr < 0.05, na.rm = TRUE), n_total = .N), by = comparison])
