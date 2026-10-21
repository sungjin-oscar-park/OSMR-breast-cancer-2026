# Mapping to 20261021_Main_Figures.pdf

The page and letter assignments below use the final PDF, not the historical source filenames.

| Final panel | Analysis and output | Included source |
|---|---|---|
| 1a | IL-6-family bulk heatmap, OSM-sorted samples | `figure1/bulk.R`, `run_fig1a` |
| 1b | OSMR, IL6ST and IL27RA Kaplan–Meier/Cox analyses | `figure1/bulk.R`, `run_fig1b` |
| 1c | OSM–receptor Spearman correlations | `figure1/bulk.R`, `run_fig1c` |
| 1d | OSMR/IL6ST high–low MCP-counter signature differences | `figure1/bulk.R`, `run_fig1d` |
| 1e, left | OSMR-stratified Hallmark rank-score heatmap | `figure1/bulk.R`, `run_fig1e` |
| 1e, right | Ranked upregulated top-500 input; network built externally in Metascape | `figure1/bulk.R`, `run_fig1e_metascape_deg` |
| 1f | Major cell-type UMAP | `figure1/single_cell.R`, `run_fig1f` |
| 1g | IL-6-family mean expression by recipient/source compartment | `figure1/single_cell.R`, `run_fig1g` |
| 1h | Macrophage-OSM × receptor-complex potential using the corrected master’s epithelial-only CopyKAT transfer | `figure1/single_cell.R`, `run_fig1h` |
| 1i | CopyKAT-resolved epithelial UMAP, receptor detection and NATMI-style edge specificity | `figure1/single_cell.R`, `run_fig1i` |
| 1j | Joint standardized receptor associations with the OSM-response module | `figure1/single_cell.R`, `run_fig1j` |
| 4a | Untreated OSMR, IL6ST, CDKN1A and EGFR; 18-test Holm family | `figure4/preprocess_and_fit.R`, `analyze.R`, `plot.R` |
| 4b | Differential-expression breadth across eight comparisons | Same Figure 4 pipeline |
| 4c | IL6–JAK–STAT3 Hallmark enrichment, moderated-t ranking | Same Figure 4 pipeline |
| 4d | CDKN1A observations and within-line response statistics | Same Figure 4 pipeline |
| 4e | CDKN1A/EGFR effects, pointwise CI, FDR symbols and EGFR nominal P | Same Figure 4 pipeline |

Historical remapping: `20260721_Figures_redraw/Figure1.R` old a/b/d/e/f become final 1a/b/c/d/e; `Figure2.R` old a–e become final 1f–j. The old clinical-variable landscape (old Figure 1c), additional immune/Hallmark panels, Figure 4 RTK-regulator panels, interaction contrasts and the separate RNA-seq-to-TCGA exploratory bridge are not included.

The statistical baseline family intentionally retains CLDN4 and CLDN7, although they are not displayed, so the original Holm correction is preserved. Figure 4b retains the source code's absolute effect cutoff **0.585** (on the mean log2(TPM + 1) difference scale), which produces the displayed counts; it is not silently changed to 0.5.

Panel exports retain source plotting styles and the latest scripted Figure 4e formatting. The assembled PDF contains manual artwork edits beyond those scripts; this repository does not claim exact pixel-level reconstruction of the full pages.
