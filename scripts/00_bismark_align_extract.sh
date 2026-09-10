#!/usr/bin/env bash
set -euo pipefail

# =====================================================
# Bismark WGBS Post-alignment Pipeline (Corrected)
# - paired-end BAMs assumed
# - extracts methylation with Bismark
# - builds global and context-specific summaries
# - creates BigWig tracks from coverage / cytosine report
# - generates Bismark HTML reports if available
# =====================================================

# -----------------------------
# 0) USER SETTINGS
# -----------------------------
# Paths default to the repository layout but can be overridden from the
# environment, e.g.  BAM_FOLDER=/data/bam GENOME_FOLDER=/ref/genome ./00_bismark_align_extract.sh
PROJECT_ROOT="${PROJECT_ROOT:-.}"
BAM_FOLDER="${BAM_FOLDER:-$PROJECT_ROOT/data/bam}"
GENOME_FOLDER="${GENOME_FOLDER:-$PROJECT_ROOT/reference/genome}"
ANALYSIS_FOLDER="${ANALYSIS_FOLDER:-$PROJECT_ROOT/results/bismark}"
TRACKS_FOLDER="${TRACKS_FOLDER:-$PROJECT_ROOT/results/tracks}"
LOG_FOLDER="$ANALYSIS_FOLDER/logs"

# reference FASTA inside genome folder
GENOME_FASTA="$GENOME_FOLDER/genome.fa"

# resources
THREADS=8
BISMARK_MULTICORE=4

# coverage threshold for optional filtered tracks/summaries
MIN_COVERAGE=5

# set to "true" if BAMs are paired-end and need name sorting
PAIRED_END="true"

mkdir -p "$ANALYSIS_FOLDER" "$TRACKS_FOLDER" "$LOG_FOLDER"

# -----------------------------
# 1) CHECK DEPENDENCIES
# -----------------------------
need_cmd() {
    command -v "$1" >/dev/null 2>&1 || {
        echo "[ERROR] Required program not found in PATH: $1" >&2
        exit 1
    }
}

for cmd in samtools awk sort gzip zcat bedGraphToBigWig bismark_methylation_extractor; do
    need_cmd "$cmd"
done

# optional but recommended
HAS_BISMARK2REPORT="false"
HAS_BISMARK2SUMMARY="false"
if command -v bismark2report >/dev/null 2>&1; then HAS_BISMARK2REPORT="true"; fi
if command -v bismark2summary >/dev/null 2>&1; then HAS_BISMARK2SUMMARY="true"; fi

# -----------------------------
# 2) CHROM SIZES
# -----------------------------
CHRSIZES="$GENOME_FOLDER/chrom.sizes"

if [[ ! -f "$CHRSIZES" ]]; then
    echo "[INFO] Generating chromosome sizes from $GENOME_FASTA"
    samtools faidx "$GENOME_FASTA"
    cut -f1,2 "${GENOME_FASTA}.fai" > "$CHRSIZES"
fi

# -----------------------------
# 3) METHYLATION EXTRACTION
# -----------------------------
echo "[INFO] Starting methylation extraction"

shopt -s nullglob
BAMS=("$BAM_FOLDER"/*.bam)

if [[ ${#BAMS[@]} -eq 0 ]]; then
    echo "[ERROR] No BAM files found in $BAM_FOLDER" >&2
    exit 1
fi

for BAM in "${BAMS[@]}"; do
    SAMPLE=$(basename "$BAM" .bam)
    echo "[INFO] Processing sample: $SAMPLE"

    INPUT_BAM="$BAM"

    if [[ "$PAIRED_END" == "true" ]]; then
        NAME_SORTED="$ANALYSIS_FOLDER/${SAMPLE}.name_sorted.bam"
        if [[ ! -f "$NAME_SORTED" ]]; then
            echo "[INFO] Name-sorting BAM for paired-end extractor input: $SAMPLE"
            samtools sort -n -@ "$THREADS" -o "$NAME_SORTED" "$BAM"
        fi
        INPUT_BAM="$NAME_SORTED"
    fi

    # Bismark extraction
    # --comprehensive pools strands into context-dependent files
    # --bedGraph writes bedGraph + coverage output
    # --cytosine_report + --CX creates genome-wide cytosine report incl. non-CG contexts
    LOGFILE="$LOG_FOLDER/${SAMPLE}.methylation_extractor.log"

    echo "[INFO] Running bismark_methylation_extractor for $SAMPLE"
    bismark_methylation_extractor \
        --gzip \
        --bedGraph \
        --comprehensive \
        --cytosine_report \
        --CX \
        --multicore "$BISMARK_MULTICORE" \
        --genome_folder "$GENOME_FOLDER" \
        --output "$ANALYSIS_FOLDER" \
        "$INPUT_BAM" \
        >"$LOGFILE" 2>&1
done

# -----------------------------
# 4) NORMALIZE OUTPUT FILENAMES
# -----------------------------
# This block tries to map outputs back to a simple sample name.
# Bismark names files from the input BAM basename; if input was name-sorted,
# outputs will contain ".name_sorted". We create symlink-like copies with cleaner names.
echo "[INFO] Normalizing output names"

for BAM in "${BAMS[@]}"; do
    SAMPLE=$(basename "$BAM" .bam)

    # candidate prefixes
    PREFIX1="${SAMPLE}"
    PREFIX2="${SAMPLE}.name_sorted"

    # helper: if file with PREFIX2 exists and PREFIX1 does not, copy it
    for pattern in \
        "*.bismark.cov.gz" \
        "*.bedGraph.gz" \
        "*.CX_report.txt.gz" \
        "*M-bias.txt" \
        "*splitting_report.txt" \
        "CpG_context_*.txt.gz" \
        "CHG_context_*.txt.gz" \
        "CHH_context_*.txt.gz"
    do
        for f in "$ANALYSIS_FOLDER"/${pattern}; do
            [[ -e "$f" ]] || continue
            base=$(basename "$f")
            if [[ "$base" == *"$PREFIX2"* ]]; then
                clean="${base//$PREFIX2/$PREFIX1}"
                if [[ ! -e "$ANALYSIS_FOLDER/$clean" ]]; then
                    cp -n "$f" "$ANALYSIS_FOLDER/$clean"
                fi
            fi
        done
    done
done

# -----------------------------
# 5) GLOBAL CONTEXT-SPECIFIC SUMMARY
# Uses CX_report columns:
# 1 chr, 2 pos, 3 strand, 4 methylated, 5 unmethylated, 6 context, 7 trinuc
# -----------------------------
echo "[INFO] Calculating global methylation summary"

QC_OUT="$ANALYSIS_FOLDER/global_methylation_summary.tsv"
echo -e "Sample\tContext\tMethylated\tUnmethylated\tTotal\tPercent_Methylation" > "$QC_OUT"

CX_REPORTS=("$ANALYSIS_FOLDER"/*CX_report.txt.gz)
if [[ ${#CX_REPORTS[@]} -gt 0 ]]; then
    for CX in "${CX_REPORTS[@]}"; do
        SAMPLE=$(basename "$CX" .CX_report.txt.gz)
        zcat "$CX" | awk -v sample="$SAMPLE" '
            BEGIN { OFS="\t" }
            {
                context = $6
                meth[context] += $4
                unmeth[context] += $5
            }
            END {
                for (c in meth) {
                    total = meth[c] + unmeth[c]
                    perc = (total > 0 ? (meth[c] / total) * 100 : 0)
                    printf "%s\t%s\t%d\t%d\t%d\t%.4f\n",
                           sample, c, meth[c], unmeth[c], total, perc
                }
            }
        ' >> "$QC_OUT"
    done
else
    echo "[WARN] No CX reports found; skipping global methylation summary"
fi

# -----------------------------
# 6) FILTERED GLOBAL SUMMARY (coverage >= MIN_COVERAGE)
# -----------------------------
echo "[INFO] Calculating filtered methylation summary (coverage >= $MIN_COVERAGE)"

QC_FILT_OUT="$ANALYSIS_FOLDER/global_methylation_summary.cov${MIN_COVERAGE}.tsv"
echo -e "Sample\tContext\tMethylated\tUnmethylated\tTotal\tPercent_Methylation" > "$QC_FILT_OUT"

if [[ ${#CX_REPORTS[@]} -gt 0 ]]; then
    for CX in "${CX_REPORTS[@]}"; do
        SAMPLE=$(basename "$CX" .CX_report.txt.gz)
        zcat "$CX" | awk -v sample="$SAMPLE" -v mincov="$MIN_COVERAGE" '
            BEGIN { OFS="\t" }
            {
                cov = $4 + $5
                if (cov >= mincov) {
                    context = $6
                    meth[context] += $4
                    unmeth[context] += $5
                }
            }
            END {
                for (c in meth) {
                    total = meth[c] + unmeth[c]
                    perc = (total > 0 ? (meth[c] / total) * 100 : 0)
                    printf "%s\t%s\t%d\t%d\t%d\t%.4f\n",
                           sample, c, meth[c], unmeth[c], total, perc
                }
            }
        ' >> "$QC_FILT_OUT"
    done
fi

# -----------------------------
# 7) PER-CONTEXT TEXT TABLES FROM CX REPORT
# -----------------------------
echo "[INFO] Writing per-context cytosine tables"

if [[ ${#CX_REPORTS[@]} -gt 0 ]]; then
    for CX in "${CX_REPORTS[@]}"; do
        SAMPLE=$(basename "$CX" .CX_report.txt.gz)

        zcat "$CX" | awk 'BEGIN{OFS="\t"} $6=="CG"  {print $0}' \
            | gzip -c > "$ANALYSIS_FOLDER/${SAMPLE}.CpG.CX_report.txt.gz"

        zcat "$CX" | awk 'BEGIN{OFS="\t"} $6=="CHG" {print $0}' \
            | gzip -c > "$ANALYSIS_FOLDER/${SAMPLE}.CHG.CX_report.txt.gz"

        zcat "$CX" | awk 'BEGIN{OFS="\t"} $6=="CHH" {print $0}' \
            | gzip -c > "$ANALYSIS_FOLDER/${SAMPLE}.CHH.CX_report.txt.gz"
    done
fi

# -----------------------------
# 8) BIGWIGS FROM COVERAGE FILE (.cov.gz)
# The coverage file has no context column.
# We therefore:
#   - make "combined" BigWig from .cov.gz methylation %
#   - make context-specific BigWigs from CX_report using weighted methylation %
# -----------------------------
echo "[INFO] Creating BigWig tracks"

COV_FILES=("$ANALYSIS_FOLDER"/*.bismark.cov.gz)
if [[ ${#COV_FILES[@]} -gt 0 ]]; then
    for COV in "${COV_FILES[@]}"; do
        SAMPLE=$(basename "$COV" .bismark.cov.gz)

        # combined methylation% track from Bismark coverage file
        BEDGRAPH="$ANALYSIS_FOLDER/${SAMPLE}.combined.bedGraph"
        SORTED_BG="$ANALYSIS_FOLDER/${SAMPLE}.combined.sorted.bedGraph"

        zcat "$COV" | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $3, $4}' > "$BEDGRAPH"
        sort -k1,1 -k2,2n "$BEDGRAPH" > "$SORTED_BG"
        bedGraphToBigWig "$SORTED_BG" "$CHRSIZES" "$TRACKS_FOLDER/${SAMPLE}.combined.bw"
        rm -f "$BEDGRAPH" "$SORTED_BG"
    done
else
    echo "[WARN] No .bismark.cov.gz files found"
fi

# -----------------------------
# 9) CONTEXT-SPECIFIC BIGWIGS FROM CX REPORT
# For each cytosine:
#   methylation% = 100 * meth / (meth + unmeth)
# Creates:
#   sample.CG.bw
#   sample.CHG.bw
#   sample.CHH.bw
# and filtered versions with coverage >= MIN_COVERAGE
# -----------------------------
echo "[INFO] Creating context-specific BigWigs from CX reports"

if [[ ${#CX_REPORTS[@]} -gt 0 ]]; then
    for CX in "${CX_REPORTS[@]}"; do
        SAMPLE=$(basename "$CX" .CX_report.txt.gz)

        for CONTEXT in CG CHG CHH; do
            BG="$ANALYSIS_FOLDER/${SAMPLE}.${CONTEXT}.bedGraph"
            BG_SORT="$ANALYSIS_FOLDER/${SAMPLE}.${CONTEXT}.sorted.bedGraph"

            zcat "$CX" | awk -v ctx="$CONTEXT" '
                BEGIN { OFS="\t" }
                $6 == ctx {
                    cov = $4 + $5
                    if (cov > 0) {
                        perc = 100 * $4 / cov
                        print $1, $2-1, $2, perc
                    }
                }
            ' > "$BG"

            if [[ -s "$BG" ]]; then
                sort -k1,1 -k2,2n "$BG" > "$BG_SORT"
                bedGraphToBigWig "$BG_SORT" "$CHRSIZES" "$TRACKS_FOLDER/${SAMPLE}.${CONTEXT}.bw"
            fi

            rm -f "$BG" "$BG_SORT"

            # filtered by minimum coverage
            BG_F="$ANALYSIS_FOLDER/${SAMPLE}.${CONTEXT}.cov${MIN_COVERAGE}.bedGraph"
            BG_F_SORT="$ANALYSIS_FOLDER/${SAMPLE}.${CONTEXT}.cov${MIN_COVERAGE}.sorted.bedGraph"

            zcat "$CX" | awk -v ctx="$CONTEXT" -v mincov="$MIN_COVERAGE" '
                BEGIN { OFS="\t" }
                $6 == ctx {
                    cov = $4 + $5
                    if (cov >= mincov) {
                        perc = 100 * $4 / cov
                        print $1, $2-1, $2, perc
                    }
                }
            ' > "$BG_F"

            if [[ -s "$BG_F" ]]; then
                sort -k1,1 -k2,2n "$BG_F" > "$BG_F_SORT"
                bedGraphToBigWig "$BG_F_SORT" "$CHRSIZES" "$TRACKS_FOLDER/${SAMPLE}.${CONTEXT}.cov${MIN_COVERAGE}.bw"
            fi

            rm -f "$BG_F" "$BG_F_SORT"
        done
    done
fi

# -----------------------------
# 10) OPTIONAL: Bismark HTML REPORTS
# bismark2report: per-sample report
# bismark2summary: directory-wide summary
# -----------------------------
echo "[INFO] Generating Bismark reports where possible"

if [[ "$HAS_BISMARK2REPORT" == "true" ]]; then
    for BAM in "${BAMS[@]}"; do
        SAMPLE=$(basename "$BAM" .bam)
        (
            cd "$ANALYSIS_FOLDER"
            # bismark2report auto-discovers matching files in many setups
            bismark2report --dir "$ANALYSIS_FOLDER" >/dev/null 2>&1 || true
        )
    done
else
    echo "[WARN] bismark2report not found; skipping per-sample HTML reports"
fi

if [[ "$HAS_BISMARK2SUMMARY" == "true" ]]; then
    (
        cd "$ANALYSIS_FOLDER"
        bismark2summary >/dev/null 2>&1 || true
    )
else
    echo "[WARN] bismark2summary not found; skipping global HTML summary"
fi

# -----------------------------
# 11) SIMPLE QC TABLE FROM COVERAGE FILES
# -----------------------------
echo "[INFO] Building simple per-sample coverage QC table"

COV_QC="$ANALYSIS_FOLDER/coverage_qc_summary.tsv"
echo -e "Sample\tSites\tMean_Coverage\tMean_Methylation_Percent" > "$COV_QC"

if [[ ${#COV_FILES[@]} -gt 0 ]]; then
    for COV in "${COV_FILES[@]}"; do
        SAMPLE=$(basename "$COV" .bismark.cov.gz)
        zcat "$COV" | awk -v sample="$SAMPLE" '
            BEGIN { OFS="\t" }
            {
                cov = $5 + $6
                n++
                sum_cov += cov
                sum_meth += $4
            }
            END {
                mean_cov  = (n > 0 ? sum_cov / n : 0)
                mean_meth = (n > 0 ? sum_meth / n : 0)
                printf "%s\t%d\t%.4f\t%.4f\n", sample, n, mean_cov, mean_meth
            }
        ' >> "$COV_QC"
    done
fi

# -----------------------------
# 12) DONE
# -----------------------------
echo "[INFO] Pipeline complete"
echo "[INFO] Analysis folder: $ANALYSIS_FOLDER"
echo "[INFO] Tracks folder:   $TRACKS_FOLDER"
echo "[INFO] Global summary:  $QC_OUT"
echo "[INFO] Filtered summary:$QC_FILT_OUT"
echo "[INFO] Coverage QC:     $COV_QC"
