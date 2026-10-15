# Pathway definitions are external inputs; no measured data are bundled.
load_hallmarks <- function(cfg) {
  if (!is.null(cfg$hallmark_rds)) {
    hs <- readRDS(required_input(cfg, "hallmark_rds"))
  } else {
    hs <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
  }
  stopifnot(is.data.frame(hs), all(c("gs_name", "gene_symbol") %in% names(hs)))
  hs <- hs[grepl("^HALLMARK_", hs$gs_name), , drop = FALSE]
  stopifnot(length(unique(hs$gs_name)) == 50L)
  # Preserve the precise definitions/version used for each run outside Git.
  provenance <- output_path(cfg, "provenance")
  saveRDS(hs, file.path(provenance, "hallmark_definitions.rds"))
  if ("db_version" %in% names(hs))
    writeLines(unique(hs$db_version), file.path(provenance, "msigdb_version.txt"))
  hs
}
