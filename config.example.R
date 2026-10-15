# Copy outside the repository (or to ignored config.local.R), then edit paths.
# Use absolute paths. NULL inputs are needed only for the corresponding step.
list(
  output_dir = "/absolute/path/to/osmr-analysis-output",
  tcga_dir = NULL,         # cBioPortal z-score matrix and clinical_data.tsv
  scrna_raw_dir = NULL,    # count_matrix_sparse.mtx, genes, barcodes, metadata
  scrna_rds = NULL,        # processed Seurat object; optional after preprocess
  copykat_dir = NULL,      # *_ck_copykat_prediction.txt; optional after copykat
  tpm_acute = NULL,        # supplied acute gene-level TPM .xlsx
  tpm_sustained = NULL,    # supplied sustained-expression gene-level TPM .xlsx
  hallmark_rds = NULL,    # optional external msigdbr human Hallmark data.frame
  # FALSE uses the archived OSM_Response column if present; TRUE recalculates.
  recompute_module = FALSE,
  workers = 1L,           # CopyKAT cores; increase for an adequately sized host
  font_family = "Arial"
)
