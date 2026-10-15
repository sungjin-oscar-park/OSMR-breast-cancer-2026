# Optional, computationally expensive upstream CopyKAT run, separately per sample.
suppressPackageStartupMessages({ library(Seurat); library(copykat) })
set.seed(42)
seu <- readRDS(required_input(cfg, "scrna_rds", file.path(cfg$output_dir, "figure1/preprocessed/BRCA_scRNASeq_processed.rds")))
res_dir <- output_path(cfg, "figure1/copykat")
# CopyKAT writes files in its working directory.
setwd(res_dir)
# --- REFINED COPYKAT: Sample-wise with References ---
samples <- unique(seu$sample)
all_ck_results <- list()

for (s in samples) {
  print(paste("Processing CopyKAT for sample:", s))
  
  # 1. Define Query (Epithelial candidates)
  # We use the paper's annotation as a base, which already filtered immune/stromal contamination
  query_cells <- WhichCells(seu, expression = (sample == s) & (Analysis_CellType == "Epithelial"))
  
  # 2. Define Reference (Normal cells from the same sample)
  ref_cells <- WhichCells(seu, expression = (sample == s) & (Analysis_CellType %in% c("T-cell", "CAF", "Endothelial")))
  
  if (length(query_cells) < 20) {
    print(paste("Skipping sample", s, "due to low epithelial count."))
    next
  }
  
  if (length(ref_cells) < 10) {
    print(paste("Sample", s, "has low reference cells. Attempting auto-ref discovery by copykat..."))
    # If no local reference, run with id.type='S' and let copykat find its own or we can use a pooled reference if necessary.
    # For now, let's try to run it without explicit ref list and let it discover if ref is too small.
    target_cells <- query_cells
    ref_param <- NULL
  } else {
    target_cells <- c(query_cells, ref_cells)
    ref_param <- ref_cells
  }
  
  # Subset counts
  counts_sub <- GetAssayData(seu, layer = "counts")[, target_cells]
  
  # Run CopyKAT
  ck_res <- tryCatch({
    copykat(rawmat = counts_sub, 
            id.type = "S", 
            cell.line = "no", 
            ngene.chr = 5, 
            win.size = 25, 
            KS.cut = 0.1, 
            sam.name = paste0(s, "_ck"), 
            norm.cell.names = ref_param,
            distance = "euclidean",
            n.cores = cfg$workers)
  }, error = function(e) {
    print(paste("CopyKAT failed for sample", s, ":", e$message))
    stop("CopyKAT failed; inspect the external output directory before resuming: ", e$message)
  })
  
  if (!is.null(ck_res)) {
    all_ck_results[[s]] <- ck_res$prediction
  }
}


if (!length(all_ck_results)) stop("No CopyKAT prediction was generated.")
