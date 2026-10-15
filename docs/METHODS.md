# Analysis definitions retained from the figure source

## Figure 1: bulk tumors

The 14-gene panel comprises OSM, OSMR, IL6ST, IL6, IL6R, IL11, IL11RA, IL27, IL27RA, LIF, LIFR, CLCF1, CNTFR and CTF1. High is expression at or above the cohort median; low is below it. Survival Cox models adjust for age group (<65/≥65) and sex, not stage or intrinsic subtype. The archived workflow does not test proportional hazards. Spearman correlation is used for the expression comparisons.

MCP-counter is called through immunedeconv on the archived gene-wise z-scores. Outputs are exploratory relative signatures, not cell fractions. The plotted high–low effect is an arithmetic mean-score difference; nominal Wilcoxon P values determine point sizes, and BH-adjusted P values are exported. Packaging fixes a local assignment bug in the historical plotting function and stops if MCP-counter fails instead of substituting another deconvolution method.

Hallmark scores use the custom within-sample formula `1 - mean(descending gene ranks within the set) / number of ranked genes`, requiring at least 15 overlapping genes. This is not the GSVA/ssGSEA algorithm. Two-sided Welch tests and BH correction select OSMR-associated Hallmarks for the heatmap. The Metascape input procedure is in [INPUTS.md](INPUTS.md).

## Figure 1: single cells

Optional upstream processing uses >200 and <7,500 detected genes, mitochondrial fraction <20%, Seurat log normalization, 2,000 variable features, PCA, Harmony by sample and UMAP on Harmony dimensions 1–20. The original public atlas had already undergone doublet filtering; this reconstruction adds no new DoubletFinder step.

CopyKAT runs separately per sample on epithelial query cells and matched T-cell/CAF/endothelial references. Samples with fewer than 20 query cells are skipped; fewer than 10 reference cells triggers the original automatic-reference setting. Parameters are `id.type = "S"`, `cell.line = "no"`, `ngene.chr = 5`, `win.size = 25`, `KS.cut = 0.1`, `distance = "euclidean"`. Only the original epithelial compartment receives malignant (aneuploid) or non-malignant (diploid) labels. Panels h–j all use this corrected labeling from the authoritative `Figure2.R` master.

The Figure 1h potential is:

```text
mean_macrophages(OSM) × mean_target_cells(sqrt((OSMR + LIFR) × IL6ST))
```

Here expression is the Seurat normalized `data` layer; the square root is calculated per target cell before averaging. It is an expression-based potential, not a measured signaling flux. The epithelial target pool combines CopyKAT-resolved malignant and non-malignant epithelial cells. CAFs and T cells retain their original lineage labels, including cells used as CopyKAT references.

Figure 1i NATMI-style edge specificity is computed directly from expression, without loading a cached LIANA result:

```text
[mean_source(OSM) / sum_celltypes mean(OSM)] ×
[mean_target(receptor complex) / sum_celltypes mean(receptor complex)]
```

Receptor-complex expression is the arithmetic mean of IL6ST and OSMR or LIFR. The normalization universe is macrophages, malignant epithelium, non-malignant epithelium, CAFs and T cells. The displayed edges select the OSMR complex in malignant and LIFR complex in non-malignant epithelium; their display does not establish absence of other edges. Receptor detection uses raw count >0; compartment comparisons use chi-square tests.

The study-defined module contains SOCS3, CCL2, CXCL1, CXCL8, IL6, MMP1, MMP3, MMP9, VEGFA, ICAM1, CD44, SAA1, SAA2, STAT3, MYC, BCL2L1, MCL1 and BIRC5. Seurat AddModuleScore subtracts matched control-gene expression from the mean signature expression. Figure 1j fits `z(module) ~ z(log1p(OSMR counts)) + z(log1p(IL6ST counts))` in malignant epithelial cells. Bars are joint standardized coefficients; error bars are coefficient ±1.96 × SE. The coefficient difference uses the model covariance matrix and a two-sided t test. The cell-level model does not adjust for donor clustering; these are exploratory associations, not causal receptor effects.

## Figure 4

Duplicate gene symbols are averaged on the supplied TPM scale before log2(TPM + 1). Filtering requires positive expression in at least three acute or two sustained libraries across the respective experiment. Group-based limma models include no batch covariates; `eBayes(trend = TRUE)` supplies moderated statistics for six acute and two sustained within-line contrasts.

Effects and confidence intervals describe differences in mean log2(TPM + 1), not directly calculated log2 ratios of raw TPM. CIs use moderated standard errors and total degrees of freedom and are pointwise/unadjusted. Figure 4b requires genome-wide BH FDR <0.05 and absolute effect >0.585. Figure 4e significance symbols use genome-wide FDR; its footer additionally shows nominal two-sided moderated-t EGFR P values.

Baseline Welch tests compare the three cell lines for OSMR, IL6ST, CDKN1A, EGFR, CLDN4 and CLDN7. Holm correction is applied jointly over all 18 comparisons, including the two genes not displayed. These comparisons remain exploratory given unresolved library/batch provenance.

Hallmark enrichment uses moderated t-statistics, retaining the largest absolute statistic for duplicate identifiers. fgsea uses set sizes 10–500, `nPermSimple = 10000`, seed 1 and BH adjustment over tested sets within each contrast. The packaged runner uses one enrichment worker. An NES represents a treatment contrast, not an individual control sample. MSigDB definitions are supplied externally or fetched by the installed msigdbr package; release changes can affect results.

Violin bandwidths are reconstructed from the historical formula `max(0.055, 1.06 × SD(within-group residuals) × 3^(-0.2))`, pooling the acute and sustained experiments for each gene. Densities are descriptive; raw observations and mean/SD overlays are exported separately.
