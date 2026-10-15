# Gene-level TPMs are supplied inputs; limma models log2(TPM + 1), not TPM generation.
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tibble); library(limma)
})
read_tpm <- function(path, sustained = FALSE) {
  raw <- readxl::read_excel(path)
  gene_col <- names(raw)[1]
  expr <- raw %>% group_by(across(all_of(gene_col))) %>%
    summarise(across(where(is.numeric), mean), .groups = "drop") %>%
    column_to_rownames(gene_col) %>% as.matrix()
  if (!ncol(expr) || any(!is.finite(expr)) || any(expr < 0)) stop("Invalid TPM matrix.")
  meta <- data.frame(Sample = colnames(expr), stringsAsFactors = FALSE)
  if (!sustained) {
    meta <- meta %>% mutate(
      CellLine = case_when(grepl("BT549", Sample) ~ "BT549",
                           grepl("MDA231", Sample) ~ "MDA231",
                           grepl("MDA468", Sample) ~ "MDA468", TRUE ~ NA_character_),
      Condition = case_when(grepl("_NT", Sample) ~ "NT",
                            grepl("_6h", Sample) ~ "6h",
                            grepl("_24h", Sample) ~ "24h", TRUE ~ NA_character_),
      Treatment = ifelse(Condition == "NT", "NT", "OSM"))
  } else {
    meta <- meta %>% mutate(
      CellLine = case_when(grepl("MB231", Sample) ~ "MB231",
                           grepl("MB468", Sample) ~ "MB468", TRUE ~ NA_character_),
      Treatment = case_when(grepl("_OSM", Sample) ~ "OSM",
                            grepl("_V", Sample) ~ "Vector", TRUE ~ NA_character_),
      Condition = Treatment)
  }
  if (anyNA(meta)) stop("Unrecognized sample names; see docs/INPUTS.md.")
  rownames(meta) <- meta$Sample
  list(expr_mat = expr, log_expr = log2(expr + 1), metadata = meta)
}
fit_tpm <- function(obj, contrasts, min_libraries) {
  meta <- obj$metadata
  meta$Group <- factor(paste(meta$CellLine, meta$Condition, sep = "_"))
  design <- model.matrix(~ 0 + Group, data = meta)
  colnames(design) <- levels(meta$Group)
  if (qr(design)$rank < ncol(design) || nrow(design) <= ncol(design))
    stop("Design has no residual degrees of freedom or is rank deficient.")
  keep <- rowSums(obj$log_expr > 0) >= min_libraries
  fit <- lmFit(obj$log_expr[keep, , drop = FALSE], design)
  cm <- makeContrasts(contrasts = contrasts, levels = design)
  colnames(cm) <- names(contrasts)
  eBayes(contrasts.fit(fit, cm), trend = TRUE)
}
acute <- read_tpm(required_input(cfg, "tpm_acute"))
oe <- read_tpm(required_input(cfg, "tpm_sustained"), sustained = TRUE)
# Only the eight within-line comparisons used in Figure 4 are retained.
acute_contrasts <- c(
  BT549_6h_vs_NT = "BT549_6h - BT549_NT", BT549_24h_vs_NT = "BT549_24h - BT549_NT",
  MDA231_6h_vs_NT = "MDA231_6h - MDA231_NT", MDA231_24h_vs_NT = "MDA231_24h - MDA231_NT",
  MDA468_6h_vs_NT = "MDA468_6h - MDA468_NT", MDA468_24h_vs_NT = "MDA468_24h - MDA468_NT")
sustained_contrasts <- c(MB231_Response = "MB231_OSM - MB231_Vector",
                         MB468_Response = "MB468_OSM - MB468_Vector")
de <- list(fit2 = fit_tpm(acute, acute_contrasts, 3L))
deoe <- list(fit2 = fit_tpm(oe, sustained_contrasts, 2L))
OUT <- output_path(cfg, "figure4")
dir.create(file.path(OUT, "rds"), showWarnings = FALSE)
saveRDS(list(acute = acute, sustained = oe, fit_acute = de$fit2, fit_sustained = deoe$fit2),
        file.path(OUT, "rds/tpm_and_models.rds"))
