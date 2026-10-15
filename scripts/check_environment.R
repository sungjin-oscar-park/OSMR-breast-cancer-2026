#!/usr/bin/env Rscript
packages <- list(
  common = c("dplyr", "tidyr", "tibble", "ggplot2", "msigdbr"),
  figure1_bulk = c("readr", "reshape2", "ggpubr", "pheatmap", "RColorBrewer",
    "survival", "survminer", "patchwork", "cowplot", "gridExtra", "ragg",
    "immunedeconv", "MCPcounter", "writexl"),
  figure1_single_cell = c("Seurat", "pheatmap", "gtable", "patchwork", "ragg", "RColorBrewer"),
  optional_single_cell_upstream = c("Matrix", "harmony", "copykat"),
  figure4 = c("readxl", "limma", "fgsea", "BiocParallel", "svglite", "jsonlite"))
cat(R.version.string, "\n")
missing <- character()
for (group in names(packages)) {
  cat("\n", group, "\n", sep = "")
  for (p in packages[[group]]) {
    ok <- requireNamespace(p, quietly = TRUE)
    cat(sprintf("  %-18s %s\n", p, if (ok) as.character(packageVersion(p)) else "MISSING"))
    if (!ok) missing <- c(missing, p)
  }
}
if (length(missing)) {
  cat("\nInstall the packages required for your chosen step; see docs/ENVIRONMENT.md.\n")
  quit(status = 1L)
}
