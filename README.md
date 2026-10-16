# OSMR breast cancer: Figure 1 and Figure 4 analysis code

Code for **Figure 1 (page 1, panels a–j)** and **Figure 4 (page 4, panels a–e)** in `20260930_Main_Figures_SPark.pdf`, prepared on 2026-09-30.

This repository contains analysis scripts, plotting code and documentation only. Input data, clinical tables, expression matrices, cached objects, gene-set databases, statistical results and figure files are **not included**. All inputs and outputs must reside outside this repository.

## Manuscript

Analysis code for the in silico components of:

> Park S, Nguyen DT, Cho M, In S, Kim J, Hurh S, Park YS, Choi J, Kim W-K, Park BJ, Hwang J-I.
> *Clinical significance of OSMR expression and divergent growth responses to oncostatin M in breast cancer cells.*
> Manuscript submitted to the Journal of Experimental & Clinical Cancer Research, 2026.
> DOI: [to be added on acceptance]

Figures 2, 3, 5 and 6 report wet-laboratory measurements and have no analysis code. Supplementary Figures S1 and S2 are quantitative RT-PCR measurements and likewise have no analysis code. Supplementary Table S1 lists primer sequences. The in silico content of the study is therefore Figure 1 and Figure 4, which this repository covers in full.

## Data availability

No data are distributed here. Obtain each input from its original source:

| Input | Source |
|---|---|
| TCGA breast invasive carcinoma (Nature 2012) expression and clinical tables | cBioPortal export, as described in [docs/INPUTS.md](docs/INPUTS.md) |
| Breast cancer single-cell atlas (Wu et al.) | NCBI GEO, accession **GSE176078** |
| Cell-line RNA sequencing generated for this study (acute and sustained-expression TPM tables) | NCBI GEO, accession **[GEO accession to be added]**; raw reads under SRA **[SRA accession to be added]** |

Place the downloaded files outside this repository and point the configuration at them, as described in [docs/INPUTS.md](docs/INPUTS.md).

## Contents

| Directory / file | Purpose |
|---|---|
| `figure1/bulk.R` | TCGA expression, survival, correlations, MCP-counter, Hallmark scores and Metascape input |
| `figure1/sc_preprocess.R` | Optional upstream count-matrix QC, normalization, Harmony and UMAP |
| `figure1/copykat.R` | Optional sample-wise CopyKAT with matched reference cells |
| `figure1/single_cell.R` | Single-cell panels, OSM potential, epithelial CopyKAT labels, NATMI-style specificity and joint regression |
| `figure4/preprocess_and_fit.R` | Supplied TPM matrices, duplicate symbols, metadata and limma-trend |
| `figure4/analyze.R` | Baseline Welch/Holm tests, moderated effects/CI and Hallmark GSEA |
| `figure4/plot.R` | Figure 4a–e and panel exports; panel e spacing/grid revisions and EGFR nominal-P footer |
| `R/` | External path configuration and gene-set loading |
| `scripts/` | Entry point, environment check and Git data guard |
| `docs/` | Inputs, panel mapping, methods, source provenance, environment and validation |

## Run

Use R with the packages listed in [docs/ENVIRONMENT.md](docs/ENVIRONMENT.md). From this directory:

```sh
Rscript scripts/check_environment.R
cp config.example.R /absolute/path/outside/repository/osmr_config.R
# Edit the copied configuration: paths and output_dir must be absolute.
Rscript scripts/run.R figure1-bulk /absolute/path/outside/repository/osmr_config.R
Rscript scripts/run.R figure1-single-cell /absolute/path/outside/repository/osmr_config.R
Rscript scripts/run.R figure4 /absolute/path/outside/repository/osmr_config.R
```

Each step requires only its own inputs. Use the archived processed Seurat object and sample-wise CopyKAT prediction files to reproduce the published single-cell comparisons. To reconstruct these inputs from the public count matrix instead, run `figure1-preprocess`, then `figure1-copykat`, then `figure1-single-cell`. Leave `scrna_rds` and `copykat_dir` as `NULL` to use the objects generated under the external output root. CopyKAT is a separate, computationally expensive step and is never launched implicitly.

`figure4-analyze` prepares TPMs, refits the eight displayed contrasts and runs enrichment without rendering. `figure4-plot` renders existing `figure4/rds/figure4_objects.rds` in the configured output directory. The complete `figure4` command runs both. Figure 4 starts from **gene-level TPM spreadsheets**; limma does not calculate TPM.

Each completed run writes its R session information and configuration under the **external** `provenance/` directory. An optional external `hallmark_rds` preserves a particular MSigDB release. When absent, `msigdbr` uses its installed release/cache and may need network access to obtain its database. The used gene-set definitions are saved outside Git.

## Figure boundaries

See [docs/PANEL_MAP.md](docs/PANEL_MAP.md). The corrected `20260721_Figures_redraw/Figure2.R` master is the authoritative source for Figure 1f–j. Panels h–j consistently use CopyKAT labels transferred only within the original epithelial compartment.

Removed Figure 4 panels f/g and unrelated exploratory analyses are excluded. CLDN4/CLDN7 remain in the baseline test family because removing them changes the Holm adjustment for the displayed genes.

Figure 1e's right-hand network was produced in **Metascape**. This repository generates the original selection rule's top-500 gene input; upload and network export are separate steps documented in [docs/INPUTS.md](docs/INPUTS.md). The supplied PDF is manually assembled artwork. These scripts reproduce the underlying analyses and panel exports, not every manual layout/font/color edit in that PDF.

## Commit code only

This directory is the repository root; do not initialize Git in the parent research workspace. A strict `.gitignore` permits only the listed code/documentation paths. The pre-commit hook checks **all staged blobs**, including files added with `git add -f`, and rejects unexpected paths, symlinks, binary content and oversized files.

For a new checkout, activate the hook once:

```sh
git config core.hooksPath .githooks
git add .
python3 scripts/check_commit.py
git diff --cached --stat
# Review, then make the desired commit.
```

The prepared local repository has the hook enabled, with no commit, remote or publication created. Input data and generated results are excluded even if someone accidentally places them here; runtime configuration also rejects input/output paths inside this repository. `config.local.R` is ignored. To add a new code file later, update both the explicit `.gitignore` allowlist and `ALLOWED` in `scripts/check_commit.py` after reviewing that file.

Review [docs/VALIDATION.md](docs/VALIDATION.md) for the executed checks and the upstream steps not rerun during code packaging.

## License

MIT License. See [LICENSE](LICENSE).
