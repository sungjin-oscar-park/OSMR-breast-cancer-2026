# External inputs and output locations

No files described in the input table are distributed in this repository.

| Configuration key | Required content | Used by |
|---|---|---|
| `tcga_dir` | `data_mrna_mirna_merged_zscores.txt` and `clinical_data.tsv`, from the cBioPortal TCGA breast invasive carcinoma Nature 2012 export | Figure 1a–e |
| `scrna_raw_dir` | Public GSE176078/Wu et al. count files: `count_matrix_sparse.mtx`, `count_matrix_genes.tsv`, `count_matrix_barcodes.tsv`, `metadata.csv` | Optional single-cell reconstruction |
| `scrna_rds` | `BRCA_scRNASeq_processed.rds`: Seurat RNA counts/data, UMAP and annotations; optional archived `OSM_Response` score | Figure 1f–j |
| `copykat_dir` | Per-sample `*_ck_copykat_prediction.txt` with `cell.names` and `copykat.pred` columns | Figure 1h–j |
| `tpm_acute` | Original `20260111_gene_expression_TPM.xlsx` or identically structured gene-level TPM table | Figure 4 |
| `tpm_sustained` | Original `20260111_gene_expression_TPM2.xlsx` or identically structured gene-level TPM table | Figure 4 |
| `hallmark_rds` | Optional external data frame returned by human Hallmark `msigdbr`, including `gs_name`, `gene_symbol`, preferably `db_version` | Figure 1e and 4c |
| `output_dir` | Dedicated writable directory outside this repository | All steps |

## TCGA schema

The expression table has gene symbols in the first column, optional `Entrez_Gene_Id`, then sample columns containing gene-wise z-scores. The clinical table requires `Sample ID`, `Overall Survival (Months)`, numeric `Overall Survival Status` (0 censored, 1 death), `Diagnosis Age`, `Sex`, `ER Status`, `PR Status`, `HER2 Status`, and `Converted Stage`. Expression and clinical records are matched by sample ID. Cohort selection is defined by these archived input files; no live portal query or substituted cohort is performed.

## Single-cell schema

Raw metadata rows match cell barcodes. Required annotations are `sample` (or `orig.ident`), `celltype_major` and `celltype_minor`. The processed object additionally carries `Analysis_CellType`. This annotation identifies the original epithelial compartment; CopyKAT calls from reference T cells, CAFs or endothelial cells must not be relabeled as epithelial. CopyKAT labels are transferred only within the original epithelial compartment, following the corrected `Figure2.R` master. Figure 1h pools resolved malignant and non-malignant epithelial cells as its epithelial target group; Figure 1i/j use the same corrected labels. Undefined epithelial calls remain unresolved and are excluded from those resolved epithelial groups.

With `recompute_module = FALSE`, a stored `OSM_Response` column is retained to reproduce the archived analysis. If absent, or if explicitly requested, the study-defined module is recalculated using `AddModuleScore` (seed 1, 24 bins, 100 control genes per feature). This score is an expression summary, not an independently validated OSM-specific classifier.

The preprocessing script preserves original lineage annotations instead of applying CopyKAT predictions to every lineage. UMAP coordinates reconstructed with a different Seurat/Harmony version need not be pixel-identical to the archived coordinates.

## TPM spreadsheets

The first column contains gene symbols; numeric columns contain TPM observations. Duplicate symbols are averaged on the TPM scale before transformation. Sample column names encode metadata:

- Acute: cell-line tokens `BT549`, `MDA231`, `MDA468`; condition tokens `_NT`, `_6h`, `_24h`.
- Sustained: `MB231` or `MB468`; `_V` for vector and `_OSM` for sustained OSM expression.

All available supplied libraries are retained under the historical design. No batch effect is inferred from naming blocks. Library independence and batch provenance require the experimental record. Raw-read quantification scripts were not available among the figure source files, so this repository starts at the supplied TPM matrices.

## Metascape network (Figure 1e, right)

1. Run `figure1-bulk`.
2. Use `figure1/bulk/Figure1e_OSMR_High_vs_Low_Top500_for_Metascape.xlsx` under `output_dir`.
3. Submit only the `Metascape_Input` sheet's gene list, with human gene identifiers. Other sheets contain the statistical audit trail, not additional gene lists.
4. Export the enriched-term network in Metascape. Keep its session/settings and exports outside this repository.

Selection: OSMR-high versus OSMR-low, two-sided Welch tests, BH adjustment across successfully tested nonzero genes, FDR < 0.05 and mean z-score difference > 0.5; select the first 500 genes in ascending FDR then nominal P order. The historical list retains OSMR itself. The archived code does not specify the Metascape database release, network filtering or manual layout; those settings cannot be reconstructed from this script alone.

## Output structure

```text
output_dir/
  figure1/bulk/           # a-e, including Metascape input workbook
  figure1/preprocessed/   # optional reconstructed Seurat object
  figure1/copykat/        # optional sample-wise CopyKAT results
  figure1/single_cell/    # f-j and associated statistics
  figure4/rds/            # fitted objects and panel input object
  figure4/tables/         # statistical exports
  figure4/panels/         # a-e PDFs, SVGs and PNGs
  figure4/source_data/    # observations and plotted summaries
  figure4/qa/             # plot checks and session information
  provenance/            # per-step session/configuration and gene-set snapshot
```
