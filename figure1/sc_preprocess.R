# Optional upstream reconstruction from the published count matrix and metadata.
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(harmony)
})
set.seed(42)
data_dir <- required_input(cfg, "scrna_raw_dir")
res_dir <- output_path(cfg, "figure1/preprocessed")
print("Loading data...")
barcodes <- read.table(file.path(data_dir, "count_matrix_barcodes.tsv"), stringsAsFactors = FALSE)
genes <- read.table(file.path(data_dir, "count_matrix_genes.tsv"), stringsAsFactors = FALSE)
matrix_data <- readMM(file.path(data_dir, "count_matrix_sparse.mtx"))

rownames(matrix_data) <- genes$V1
colnames(matrix_data) <- barcodes$V1

metadata <- read.csv(file.path(data_dir, "metadata.csv"), row.names = 1)

# Ensure sample column exists
if(!"sample" %in% colnames(metadata)) {
  metadata$sample <- metadata$orig.ident
}

print("Creating Seurat object...")
seu <- CreateSeuratObject(counts = matrix_data, meta.data = metadata, project = "BRCA_OSM", min.cells = 3, min.features = 200)

cat("Initial cells:\n", ncol(seu), "\n")

# Free memory
rm(matrix_data)
gc()

print("Subsetting to desired subtypes for processing if necessary, or processing all...")
# 2. QC
seu[["percent.mito"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.rb"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")

# filter low-quality
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mito < 20)
cat("Cells after standard QC:\n", ncol(seu), "\n")

# Skipping DoubletFinder because original paper already filtered them, and Seurat v5 merge causes reduction conflicts.
cat("Cells ready for Integration:\n", ncol(seu), "\n")


# 4. Normalization, Integration (Harmony)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu)

print("Running Harmony...")
seu <- RunHarmony(seu, group.by.vars = "sample")

seu <- RunUMAP(seu, reduction = "harmony", dims = 1:20)

# Minor Cell Type Annotation adjustment to match the requested targets
seu$Analysis_CellType <- "Other"
seu$Analysis_CellType[seu$celltype_major == "Myeloid" & grepl("Macrophage", seu$celltype_minor, ignore.case=T)] <- "Macrophage"
seu$Analysis_CellType[seu$celltype_major == "CAFs"] <- "CAF"
seu$Analysis_CellType[seu$celltype_major %in% c("Cancer Epithelial", "Normal Epithelial", "Epithelial")] <- "Epithelial"
seu$Analysis_CellType[seu$celltype_major == "T-cells"] <- "T-cell"
seu$Analysis_CellType[seu$celltype_major == "Endothelial"] <- "Endothelial"


# Keep the original epithelial annotation; CopyKAT labels are applied only
# within this compartment in single_cell.R.
saveRDS(seu, file.path(res_dir, "BRCA_scRNASeq_processed.rds"))
