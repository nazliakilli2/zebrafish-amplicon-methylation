# Amplicon bisulfite methylation analysis — zebrafish pesticide exposure

Analysis code for a targeted (amplicon) bisulfite-sequencing study of DNA
methylation in adult zebrafish (*Danio rerio*, GRCz11 / Ensembl) exposed to
sublethal fipronil and cyfluthrin. Six amplicons are assayed — three heat-shock
loci (*hsp70.1*, *hsp70.2*, *hsp70.3*), a reference locus (*actb2*), an immune
locus (*b2m*), and a xenobiotic-response locus (*cyp1a*) — and differential
methylation is called per CpG with DSS.

This repository contains the analysis scripts only. It reproduces the tables and
figures from the per-sample methylation count tables; it does not redistribute
the raw sequencing data.

## Repository layout

```
config.R              central paths and settings (sourced by every script)
reference/
  master_amplicons.bed   amplicon coordinates: chr  amplicon  start  end
scripts/
  00_bismark_align_extract.sh   alignment + methylation extraction (Bismark)
  00b_build_context_tables.R    aggregate per-sample extractor output into per-context tables
  01_filter_to_amplicons.R      restrict counts to amplicons, coverage filter
  02_dss_dmltest.R              DSS beta-binomial DMLtest (the main analysis)
  03_qc_report.R                coverage / methylation QC
  04_figure1_context_overview.R CpG/CHH/CHG methylation overview figure
  05_dss_figures.R              DSS result figures
  06_gene_panel_figures.R       per-gene panel figure
  07_heatmaps_7d.R              7-day methylation heatmaps
  08_dss_region_summary.R       amplicon-level region summary
  09_supp_nomodel_robustness.R  model-free per-CpG robustness check (supplement)
  10_supp_dml_amplicon_table.R  per-amplicon DML count table (supplement)
data/        input count tables         (contents git-ignored)
results/     DSS results + tables        (contents git-ignored)
figures/     publication figures         (contents git-ignored)
```

`data/`, `results/`, and `figures/` are tracked as empty directories; their
contents are not committed.

## Requirements

- **R** (≥ 4.0) with:
  - CRAN: `data.table`, `ggplot2`, `patchwork`, `pheatmap`, `RColorBrewer`
  - Bioconductor: `DSS`, `bsseq`
- **Bismark** and **samtools** (only for `00_bismark_align_extract.sh`, the
  upstream alignment step).

```r
install.packages(c("data.table", "ggplot2", "patchwork", "pheatmap", "RColorBrewer"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DSS", "bsseq"))
```

## Configuration

All paths are defined once in `config.R`. By default, inputs are read from
`./data` and `./reference` and outputs are written to `./results` and
`./figures`, relative to the repository root. Run scripts from the repository
root:

```sh
Rscript scripts/02_dss_dmltest.R
```

To keep data and outputs elsewhere, set `PROJECT_ROOT`:

```sh
PROJECT_ROOT=/path/to/project Rscript scripts/02_dss_dmltest.R
```

## Running the pipeline

1. **Align and extract methylation** (optional; skip if you already have the
   per-context tables). Edit the settings block or set the environment variables
   at the top of `scripts/00_bismark_align_extract.sh`, then run it. This
   produces the per-sample Bismark methylation-extractor files
   (`CpG_context_*.txt.gz`, `CHG_context_*.txt.gz`, `CHH_context_*.txt.gz`).
2. **Build per-context tables.** Point `BISMARK_DIR` at the folder holding those
   per-sample files (default `data/bismark_extraction`) and run
   `Rscript scripts/00b_build_context_tables.R`. This aggregates the per-sample
   calls into `data/tables_fixed_final/` (`all_sites_CpG.tsv`,
   `all_sites_CHH.tsv`, `all_sites_CHG.tsv`) plus a `qc_summary.tsv` for checking
   bisulfite conversion. If you already have those `all_sites_*.tsv` tables, just
   place them under `data/tables_fixed_final/` and skip this step.
3. **Filter to amplicons.** Run `Rscript scripts/01_filter_to_amplicons.R`. This
   reads `data/tables_fixed_final/all_sites_{CpG,CHH,CHG}.tsv` and writes the
   amplicon- and coverage-filtered tables into `data/`:
   `all_sites_{CpG,CHH,CHG}_counts.master_bed.filtered.tsv`.
4. **Call differential methylation.** `Rscript scripts/02_dss_dmltest.R` →
   `results/DSS_DML_results.tsv` (per-CpG DSS DMLtest, smoothing enabled,
   BH-FDR within each comparison).
5. **QC, figures, and supplementary tables** (each reads the filtered tables
   and/or the DSS results; run in any order):
   `03_qc_report.R`, `04`–`07` (figures), `08_dss_region_summary.R`,
   `09_supp_nomodel_robustness.R`, `10_supp_dml_amplicon_table.R`.

## Methods notes

- **Direction convention:** DSS `DMLtest()` is called with `group1 = treatment`,
  `group2 = control`, so `diff = mu1 - mu2 = treatment - control`. A negative
  `diff` is hypomethylation (treatment lower than control); positive is
  hypermethylation.
- **Coverage floor:** CpG sites are kept at ≥10× coverage for the DSS analysis
  (`MIN_COV` in `02_dss_dmltest.R`).
- **Model-free robustness check:** `09_supp_nomodel_robustness.R` recomputes the
  per-CpG change without any model or smoothing (mean treatment − mean control
  at each site, classified by a fixed effect-size threshold), as a supplement
  demonstrating the directional signal does not depend on DSS smoothing. These
  counts use a display threshold, not FDR, and are reported separately from the
  DSS-based counts.

## Amplicons

| gene    | chr | start    | end      |
|---------|-----|----------|----------|
| hsp70.1 | 3   | 26113353 | 26113669 |
| hsp70.2 | 3   | 26120181 | 26120519 |
| hsp70.3 | 3   | 26126142 | 26126458 |
| actb2   | 3   | 40528408 | 40528826 |
| b2m     | 4   | 12794094 | 12794494 |
| cyp1a   | 18  | 5594770  | 5595161  |

## Data and code availability

- **Code.** This repository contains all analysis code. For a citable, permanent
  archive, the repository is linked to [Zenodo](https://zenodo.org): logging in
  to Zenodo with GitHub, enabling this repository there, and publishing a GitHub
  Release causes Zenodo to archive that snapshot and mint a DOI (a version DOI
  for each release and a concept DOI covering all versions). Cite the concept
  DOI in the manuscript:

  > *Code availability:* Analysis code is available at
  > https://github.com/nazliakilli2/zebrafish-amplicon-methylation and archived at Zenodo
  > (https://doi.org/10.5281/zenodo.22696675).

- **Sequencing data.** Raw bisulfite-sequencing reads are deposited in
  <!-- SRA / GEO / ENA --> under accession <!-- accession number -->. The
  per-sample methylation count tables that this code takes as input are provided
  as supplementary data / in the same archive.

## Reproducibility

Record the exact package versions used, so the analysis can be reproduced. After
loading the analysis packages in R, save the session information:

```r
writeLines(capture.output(sessionInfo()), "SESSIONINFO.txt")
```

Commit the resulting `SESSIONINFO.txt` alongside the code.

## Citing this repository

Citation metadata is in `CITATION.cff`; GitHub shows a "Cite this repository"
button from it. Update the `doi:`, `version:`, and `date-released:` fields when
the Zenodo DOI is minted and the release is tagged.

## License

Released under the MIT License (see `LICENSE`). Copyright (c) 2026 Nazlı Akıllı.
