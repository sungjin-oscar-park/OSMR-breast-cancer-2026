# Packaging validation — 2026-09-30

All validation inputs, logs, intermediate objects, rendered figures and numerical comparison reports were written outside this repository. No input datasets or derived result tables are included here.

## Executed checks

| Check | Result |
|---|---|
| R parsing | All included R scripts parse successfully |
| Figure 1 bulk route | Executed from archived cBioPortal expression/clinical files; panels a–e and top-500 Metascape input generated |
| Figure 1 single-cell route | Executed from archived Seurat object and per-sample CopyKAT predictions; panels f–j generated |
| Figure 1i/j | Epithelial counts, receptor detection and standardized regression coefficients agree with the submitted panel values |
| Figure 1h | Re-executed successfully; calculation matches the corrected master’s `run_fig2c()` after panel renaming and score export |
| Figure 4 complete route | Re-executed from the two supplied TPM spreadsheets through models, enrichment and figure exports |
| Figure 4 moderated statistics | All gene-by-contrast effects, SEs, CI bounds, nominal P and BH FDR values exactly match the archived objects for all eight displayed contrasts |
| Figure 4 baseline correction | Complete 18-test family matches the original objects |
| Figure 4 plotted statistics | DEG counts and CDKN1A/EGFR effect tables match the prior Figure 4 source-data exports |
| Figure 4 enrichment | Displayed IL6–JAK–STAT3 NES and adjusted P values match the archived output |
| Figure 4 violin widths | Original cached bandwidths independently reconstructed from TPM observations, within numerical precision |
| Render inspection | Regenerated Figure 4 and Figure 1h visually inspected; Figure 4e retains raised symbols and no x/y grid |
| Git exclusions | Representative expression, clinical, object, image, PDF, raw-read and local-config paths ignored |
| Pre-commit guard | Code-only index accepted; force-added data path, symlink, binary disguise and oversized-file probes rejected in a separate temporary repository |
| Original source integrity | Source scripts are unchanged; original SHA-256 values are recorded in PROVENANCE.md |

## Authoritative single-cell source

The corrected `20260721_Figures_redraw/Figure2.R` master is the user-confirmed analysis source. The Figure 1h function follows that master's `run_fig2c()` directly; panels h–j share the corrected epithelial-only CopyKAT transfer. The superseded all-lineage calculation and its comparison export have been removed from the package. The manuscript and supplied PDF were not modified during this code correction.

**2026-10-15.** The package was re-executed end to end against the archived inputs. Twenty of the twenty-one values reported in the manuscript reproduced exactly; the only exception was Figure 1h, whose published panel (0.047 fibroblasts, 0.041 epithelium) came from the superseded all-lineage calculation. The manuscript Results text has since been updated to the epithelial-only values produced by this code (0.060 fibroblasts, 0.027 epithelium, below 0.001 T cells). Manuscript, figure and code therefore describe the same calculation.

**2026-10-16.** The reference artwork is now `20261015_Main_Figures_SPark.pdf`, whose Figure 1h panel prints the redrawn 0.060, 0.027 and <0.001. A page-by-page text comparison with `20260930_Main_Figures_SPark.pdf`, the file the re-execution was checked against, found that panel and the Figure 2b LIFR amplicon label (484 bp to 504 bp) to be the only changes; the Figure 4 page shows no difference in extracted text. The 2026-10-15 value checks therefore still hold for the current PDF.

## Boundaries of validation

- Raw-count QC/Harmony/UMAP reconstruction and the expensive sample-wise CopyKAT runs were included and syntax-checked, but not rerun. Panel validation used their archived output files. Software versions can alter reconstructed embeddings/copy-number calls.
- Metascape upload and network generation were not rerun. The top-500 input was regenerated; the archived code does not establish the original external network settings/database version.
- Full-page manual artwork assembly was not reconstructed. The supplied PDF remains the visual reference; panel source exports preserve the available plotting code.
- The historical Cox model emits a convergence warning for its sex covariate. Its definition was preserved rather than silently altered. Tied ranks and fgsea probability-floor warnings likewise remain visible in the external logs. Numerical reproduction does not resolve these statistical assumptions.

No Git commit, remote connection or upload was performed. Only the code repository was initialized and its local pre-commit hook enabled.
