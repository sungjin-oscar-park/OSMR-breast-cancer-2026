# Final Figure 1f-j. Processed Seurat object and CopyKAT predictions are external.
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(tidyr)
  library(pheatmap); library(grid); library(gtable)
  library(patchwork); library(ragg); library(RColorBrewer)
})

RDS_INPUT <- required_input(cfg, "scrna_rds", file.path(cfg$output_dir, "figure1/preprocessed/BRCA_scRNASeq_processed.rds"))
DATA_DIR <- required_input(cfg, "copykat_dir", file.path(cfg$output_dir, "figure1/copykat"))
OUT_DIR <- output_path(cfg, "figure1/single_cell")

# Nature-style, colour-blind-safe visual system. Cell-type colours are fixed by
# biological identity and reused wherever those identities recur.
FONT_FAMILY <- cfg$font_family
BASE_FONT_SIZE <- 6
GEOM_TEXT_SIZE <- BASE_FONT_SIZE / ggplot2::.pt
CELL_TYPE_COLORS <- c(
  "CAF" = "#CC79A7",
  "Epithelial" = "#E69F00",
  "Macrophage" = "#009E73",
  "Other" = "#999999",
  "T cell" = "#56B4E9"
)
EPITHELIAL_STATE_COLORS <- c(
  "Malignant Epithelial" = "#D55E00",
  "Non-malignant Epithelial" = "#0072B2",
  "Malignant\nEpithelial" = "#D55E00",
  "Non-malignant\nEpithelial" = "#0072B2"
)
RECEPTOR_COLORS <- c("OSMR" = "#D55E00", "IL6ST" = "#009E73")

theme_nature <- function(base_size = BASE_FONT_SIZE) {
  theme_classic(base_size = base_size, base_family = FONT_FAMILY) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks.length = unit(1.4, "mm"),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = "black"),
      panel.spacing = unit(3.2, "mm"),
      plot.margin = margin(4, 5, 4, 4)
    )
}

save_pub_plot <- function(plot, stem, width, height, dpi = 600) {
  ggsave(file.path(OUT_DIR, paste0(stem, ".pdf")), plot = plot,
         device = cairo_pdf, width = width, height = height, units = "in",
         bg = "white")
  ragg::agg_tiff(file.path(OUT_DIR, paste0(stem, ".tiff")),
                 width = width, height = height, units = "in", res = dpi,
                 compression = "lzw", background = "white")
  print(plot)
  dev.off()
}

save_pub_grid <- function(draw_fun, stem, width, height, dpi = 600) {
  cairo_pdf(file.path(OUT_DIR, paste0(stem, ".pdf")),
            width = width, height = height, family = FONT_FAMILY,
            bg = "white")
  draw_fun()
  dev.off()
  ragg::agg_tiff(file.path(OUT_DIR, paste0(stem, ".tiff")),
                 width = width, height = height, units = "in", res = dpi,
                 compression = "lzw", background = "white")
  draw_fun()
  dev.off()
}

# =============================================================================
# 0.  Shared data loading + CopyKAT-refined annotation (run once)
# =============================================================================
cat("Loading processed Seurat object (once)...\n")
seu <- readRDS(RDS_INPUT)

# ---- CopyKAT refinement of EPITHELIAL CELLS ONLY ------------------------------
# CopyKAT prediction files contain other lineages as well. Restricting the label
# transfer here prevents diploid T cells/CAFs/endothelial cells from being
# mislabelled as "Normal_Epithelial" in Fig. 1i/j.
ck_files <- list.files(DATA_DIR, pattern = "_ck_copykat_prediction.txt$", full.names = TRUE)
if (length(ck_files) > 0) {
  ck <- do.call(rbind, unname(lapply(ck_files, function(f){ d <- read.table(f, header = TRUE, sep = "\t"); rownames(d) <- d$cell.names; d })))
  seu$CT_Refined <- as.character(seu$Analysis_CellType)
  epithelial_cells <- Cells(seu)[as.character(seu$Analysis_CellType) == "Epithelial"]
  mc_epi <- intersect(epithelial_cells, rownames(ck))
  ck_pred <- as.character(ck[mc_epi, "copykat.pred"])
  defined <- ck_pred %in% c("aneuploid", "diploid")
  seu@meta.data[mc_epi[defined], "CT_Refined"] <- ifelse(
    ck_pred[defined] == "aneuploid",
    "Malignant_Epithelial",
    "Normal_Epithelial"
  )
  cat(
    "CopyKAT epithelial-only labels:",
    sum(ck_pred == "aneuploid"), "malignant;",
    sum(ck_pred == "diploid"), "non-malignant;",
    sum(!defined), "undefined epithelial cells retained as Epithelial\n"
  )
} else {
  stop("No CopyKAT prediction files found in copykat_dir.")
}

# ---- OSM response score (Fig 2e): use stored column, else recompute ----
if (isTRUE(cfg$recompute_module) || !("OSM_Response" %in% colnames(seu@meta.data))) {
  osm_sig <- c("SOCS3","CCL2","CXCL1","CXCL8","IL6","MMP1","MMP3","MMP9","VEGFA","ICAM1",
               "CD44","SAA1","SAA2","STAT3","MYC","BCL2L1","MCL1","BIRC5")
  osm_sig <- intersect(osm_sig, rownames(seu))
  seu <- AddModuleScore(seu, features = list(osm_sig), name = "OSM_Response_Score", seed = 1, nbin = 24, ctrl = 100)
  seu$OSM_Response <- seu$OSM_Response_Score1
}

get_stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

# =============================================================================
# 1f.  Major cell-type UMAP (5 categories; Endothelial merged into Other)
# =============================================================================
run_fig1f <- function() {
  ct <- rep("Other", ncol(seu))
  ct[seu$celltype_major == "Myeloid" & grepl("Macrophage", seu$celltype_minor, ignore.case = TRUE)] <- "Macrophage"
  ct[seu$celltype_major == "CAFs"] <- "CAF"
  ct[seu$celltype_major %in% c("Cancer Epithelial","Normal Epithelial","Epithelial")] <- "Epithelial"
  ct[seu$celltype_major == "T-cells"] <- "T cell"
  seu$Fig1f_CellType <- factor(ct, levels = c("CAF","Epithelial","Macrophage","Other","T cell"))
  set.seed(2026)
  p <- DimPlot(
    seu, reduction = "umap", group.by = "Fig1f_CellType",
    pt.size = 0.25, raster = TRUE, raster.dpi = c(600, 600), shuffle = TRUE
  ) +
    scale_color_manual(values = CELL_TYPE_COLORS, name = "Cell type") +
    theme_nature() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.position = "right"
    ) +
    guides(color = guide_legend(override.aes = list(size = 2.8, alpha = 1))) +
    labs(title = NULL, x = "UMAP 1", y = "UMAP 2")
  save_pub_plot(p, "Figure1f_UMAP_MajorTypes", width = 4.6, height = 3.5)
  cat("Fig1f saved\n")
}

# =============================================================================
# 1g.  IL6-family ligand/receptor mean-expression heatmap (10 genes x 4 cell types)
# =============================================================================
run_fig1g <- function() {
  ct <- rep("Other", ncol(seu))
  ct[seu$celltype_major == "Myeloid" & grepl("Macrophage", seu$celltype_minor, ignore.case = TRUE)] <- "Macrophage"
  ct[seu$celltype_major == "CAFs"] <- "CAF"
  ct[seu$celltype_major %in% c("Cancer Epithelial","Normal Epithelial","Epithelial")] <- "Epithelial"
  ct[seu$celltype_major == "T-cells"] <- "T cell"
  seu$Fig1g_CT <- ct
  genes <- c("LIF","OSMR","LIFR","IL11","CNTFR","IL6ST","CLCF1","CNTF","OSM","IL6R")
  celltypes <- c("CAF","Macrophage","Epithelial","T cell")
  sub <- subset(seu, cells = colnames(seu)[seu$Fig1g_CT %in% celltypes])
  sub$Fig1g_CT <- factor(sub$Fig1g_CT, levels = celltypes)
  avg <- as.matrix(AverageExpression(sub, features = genes, group.by = "Fig1g_CT", layer = "data")$RNA)[genes, celltypes]
  zmat <- t(scale(t(avg))); zmat[zmat > 2] <- 2; zmat[zmat < -2] <- -2
  row_order <- c("OSM","LIF","IL11","CLCF1","CNTF","OSMR","LIFR","IL6ST","CNTFR","IL6R")
  zmat <- zmat[row_order, ]; avg <- avg[row_order, ]
  num_labels <- matrix(ifelse(avg < 0.005, "<0.01", sprintf("%.2f", avg)),
                       nrow = nrow(avg), dimnames = dimnames(avg))
  row_anno <- data.frame(Class = factor(c(rep("Ligand", 5), rep("Receptor", 5)),
                                        levels = c("Ligand", "Receptor")))
  rownames(row_anno) <- row_order
  ph <- pheatmap(
    zmat,
    color = colorRampPalette(c("#3B4CC0", "#F7F7F7", "#B40426"))(100),
    breaks = seq(-2, 2, length.out = 101),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    gaps_row = 5,
    annotation_row = row_anno,
    annotation_colors = list(Class = c("Ligand" = "#999999", "Receptor" = "#332288")),
    annotation_names_row = FALSE,
    display_numbers = num_labels,
    number_color = "black",
    fontsize_number = BASE_FONT_SIZE,
    fontsize = BASE_FONT_SIZE,
    angle_col = 45,
    border_color = "white",
    main = "",
    legend_breaks = c(-2, 0, 2),
    legend_labels = c("\u22122", "0", "2"),
    silent = TRUE
  )
  row_names_index <- which(ph$gtable$layout$name == "row_names")
  if (length(row_names_index) == 1) {
    ph$gtable$grobs[[row_names_index]]$gp <- gpar(
      fontfamily = FONT_FAMILY,
      fontsize = BASE_FONT_SIZE,
      fontface = "italic"
    )
  }
  save_pub_grid(
    function() { grid.newpage(); grid.draw(ph$gtable) },
    "Figure1g_IL6family_heatmap", width = 4.7, height = 4.2
  )
  cat("Fig1g saved\n")
}

# =============================================================================
# 1h.  Total heterodimeric OSM signaling-potential barchart (Epithelial/CAF/T cell)
#       Score = Mean(OSM in Macrophages) * sqrt((OSMR + LIFR) * IL6ST) per cell, averaged
# =============================================================================
run_fig1h <- function() {
  tg <- case_when(seu$CT_Refined == "CAF" ~ "CAFs",
                  seu$CT_Refined %in% c("Malignant_Epithelial","Normal_Epithelial") ~ "Epithelial Cells",
                  seu$CT_Refined == "T-cell" ~ "T-cells", TRUE ~ "Others")
  seu$Target_Group <- tg
  mac_cells <- which(seu$Analysis_CellType == "Macrophage")
  osm_prod <- mean(GetAssayData(seu, layer = "data")["OSM", mac_cells])
  st <- subset(seu, subset = Target_Group %in% c("CAFs","Epithelial Cells","T-cells"))
  ed <- FetchData(st, vars = c("LIFR","OSMR","IL6ST","Target_Group"), layer = "data")
  ed$Total_Complex <- sqrt((ed$OSMR + ed$LIFR) * ed$IL6ST)
  ed$Total_Potential <- osm_prod * ed$Total_Complex
  tp <- ed %>% group_by(Target_Group) %>% summarise(Mean = mean(Total_Potential), .groups = "drop")
  disp <- c("Epithelial Cells"="Epithelial","CAFs"="CAF","T-cells"="T cell")
  tp$lab <- disp[tp$Target_Group]
  tp$lab <- factor(tp$lab, levels = c("T cell", "Epithelial", "CAF"))
  p <- ggplot(tp, aes(Mean, lab, colour = lab)) +
    geom_segment(aes(x = 0, xend = Mean, yend = lab),
                 colour = "grey78", linewidth = 0.7) +
    geom_point(size = 3.1) +
    geom_text(aes(label = sprintf("%.3f", Mean)), hjust = -0.45,
              colour = "black", size = GEOM_TEXT_SIZE, family = FONT_FAMILY) +
    scale_colour_manual(values = CELL_TYPE_COLORS, guide = "none") +
    scale_x_continuous(
      limits = c(0, max(tp$Mean) * 1.28),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      x = expression("Total " * italic(OSM) * " signalling potential (a.u.)"),
      y = NULL
    ) +
    theme_nature()
  save_pub_plot(p, "Figure1h_Total_OSM_Signaling_Potential", width = 3.5, height = 2.5)
  write.csv(tp, file.path(OUT_DIR, "Figure1h_potential_scores.csv"), row.names = FALSE)
  cat("Fig1h saved. values:\n"); print(as.data.frame(tp[, c("lab","Mean")]))
}

# =============================================================================
# 1i.  Malignant/Non-malignant UMAP + NATMI-style edge specificity + status barcharts
# =============================================================================
run_fig1i <- function(file_prefix = "") {
  sub <- subset(seu, subset = CT_Refined %in% c("Malignant_Epithelial","Normal_Epithelial"))
  stopifnot(all(as.character(sub$Analysis_CellType) == "Epithelial"))
  sub$Predicted <- factor(ifelse(sub$CT_Refined == "Malignant_Epithelial", "Malignant Epithelial", "Non-malignant Epithelial"),
                          levels = c("Malignant Epithelial","Non-malignant Epithelial"))
  set.seed(2026)
  p_umap <- DimPlot(
    sub, group.by = "Predicted", reduction = "umap", pt.size = 0.55,
    alpha = 1,  # 100% opacity; ggplot/Seurat alpha is expressed on a 0-1 scale.
    raster = TRUE, raster.dpi = c(600, 600), shuffle = TRUE
  ) +
    scale_color_manual(values = EPITHELIAL_STATE_COLORS, name = NULL) +
    theme_nature() +
    theme(
      axis.text = element_blank(), axis.ticks = element_blank(),
      axis.line = element_blank(), legend.position = "bottom",
      legend.direction = "vertical"
    ) +
    guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    labs(title = NULL, x = "UMAP 1", y = "UMAP 2")

  # NATMI-style expression specificity dotplot — NATMI edge specificity computed DIRECTLY from
  # expression (no CSV read). Complex receptor expr = mean of subunits; specificity
  # over the same 5 cell types LIANA subset to. CopyKAT labels are transferred
  # only within the original epithelial compartment.
  liana_cts <- c("Macrophage","Malignant_Epithelial","CAF","Normal_Epithelial","T-cell")
  liana_cts <- intersect(liana_cts, unique(seu$CT_Refined))
  dmat <- GetAssayData(seu, layer = "data")
  mean_ct <- function(gene) sapply(liana_cts, function(c) mean(dmat[gene, seu$CT_Refined == c]))
  osm_m <- mean_ct("OSM"); il6st_m <- mean_ct("IL6ST"); lifr_m <- mean_ct("LIFR"); osmr_m <- mean_ct("OSMR")
  comp_LIFR <- (il6st_m + lifr_m) / 2; comp_OSMR <- (il6st_m + osmr_m) / 2
  edge_spec <- function(lig, comp, tgt) as.numeric((lig["Macrophage"]/sum(lig)) * (comp[tgt]/sum(comp)))
  spec_LIFR <- edge_spec(osm_m, comp_LIFR, "Normal_Epithelial")
  spec_OSMR <- edge_spec(osm_m, comp_OSMR, "Malignant_Epithelial")
  dp <- data.frame(
    interaction = c("OSM \u2192\nIL6ST-LIFR","OSM \u2192\nIL6ST-OSMR"),
    target = c("Non-malignant\nEpithelial","Malignant\nEpithelial"),
    spec = c(spec_LIFR, spec_OSMR))
  dp$interaction <- factor(dp$interaction, levels = c("OSM \u2192\nIL6ST-OSMR","OSM \u2192\nIL6ST-LIFR"))
  dp$target <- factor(dp$target, levels = c("Malignant\nEpithelial","Non-malignant\nEpithelial"))
  write.csv(dp, file.path(OUT_DIR, "Figure1i_edge_specificity.csv"), row.names = FALSE)
  dp$label <- sprintf("%.3f", dp$spec)
  p_dot <- ggplot(dp, aes(target, interaction, colour = target)) +
    geom_point(aes(size = spec), alpha = 0.95) +
    geom_text(aes(label = label), nudge_y = -0.20, colour = "black",
              size = GEOM_TEXT_SIZE, family = FONT_FAMILY) +
    scale_colour_manual(values = c(
      "Malignant\nEpithelial" = EPITHELIAL_STATE_COLORS[["Malignant\nEpithelial"]],
      "Non-malignant\nEpithelial" = EPITHELIAL_STATE_COLORS[["Non-malignant\nEpithelial"]]
    ), guide = "none") +
    scale_size_continuous(range = c(4.0, 5.2), guide = "none") +
    scale_x_discrete(drop = FALSE, expand = expansion(mult = c(0.25, 0.25))) +
    scale_y_discrete(drop = FALSE, expand = expansion(mult = c(0.34, 0.34))) +
    labs(
      x = "Target cell",
      y = expression(italic(OSM) ~ "\u2192 receptor complex"),
      subtitle = "Macrophage source; NATMI edge specificity"
    ) +
    theme_nature() +
    theme(
      plot.subtitle = element_text(size = BASE_FONT_SIZE, colour = "grey30"),
      axis.text.x = element_text(size = BASE_FONT_SIZE),
      axis.text.y = element_text(size = BASE_FONT_SIZE, face = "italic"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.minor = element_blank()
    )

  # status distributions (compute directly from counts for full independence)
  epi <- subset(seu, subset = CT_Refined %in% c("Malignant_Epithelial","Normal_Epithelial"))
  stopifnot(all(as.character(epi$Analysis_CellType) == "Epithelial"))
  cnt <- FetchData(epi, vars = c("OSMR","IL6ST"), layer = "counts")
  cnt$CT <- factor(ifelse(epi$CT_Refined == "Malignant_Epithelial", "Malignant\nEpithelial", "Non-malignant\nEpithelial"),
                   levels = c("Malignant\nEpithelial","Non-malignant\nEpithelial"))
  status_stats <- function(gene) {
    st <- data.frame(CT = cnt$CT, S = ifelse(cnt[[gene]] > 0, paste0(gene,"+"), paste0(gene,"-")))
    tab <- as.data.frame(table(st$CT, st$S)); colnames(tab) <- c("CT","S","Count")
    tab <- tab %>% group_by(CT) %>% mutate(Pct = Count/sum(Count)*100) %>% ungroup()
    tab$S <- factor(tab$S, levels = c(paste0(gene,"+"), paste0(gene,"-")))
    pval <- chisq.test(table(st$CT, st$S))$p.value
    list(tab = tab, pval = pval)
  }
  status_bar <- function(gene, poscol, stat) {
    tab <- stat$tab
    pval <- stat$pval
    tab$full_lab <- ifelse(tab$Pct >= 7, sprintf("%.1f%%", tab$Pct), "")
    p_lab <- if (pval < 2.2e-16) expression(italic(P) < 2.2 %*% 10^-16) else
      paste0("P = ", format.pval(pval, digits = 2))
    g <- ggplot(tab, aes(CT, Pct, fill = S)) +
      geom_col(width = 0.68, colour = "white", linewidth = 0.25) +
      geom_text(aes(label = full_lab), position = position_stack(vjust = 0.5),
                colour = "white", size = GEOM_TEXT_SIZE, family = FONT_FAMILY) +
      scale_fill_manual(values = setNames(c(poscol, "#999999"), c(paste0(gene,"+"), paste0(gene,"-"))),
                        name = NULL) +
      scale_y_continuous(
        limits = c(0, 100), breaks = c(0, 50, 100),
        labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0))
      ) +
      labs(title = gene, subtitle = p_lab, x = NULL, y = "Cells") +
      theme_nature() +
      theme(
        plot.title = element_text(face = "bold.italic", size = BASE_FONT_SIZE),
        plot.subtitle = element_text(size = BASE_FONT_SIZE),
        axis.text.x = element_text(size = BASE_FONT_SIZE),
        legend.text = element_text(size = BASE_FONT_SIZE, face = "italic"),
        legend.position = "bottom"
      )
    cat(gene, "chi-sq p =", format.pval(pval), "\n"); g
  }
  osmr_stat <- status_stats("OSMR")
  il6st_stat <- status_stats("IL6ST")
  p_osmr  <- status_bar("OSMR", RECEPTOR_COLORS[["OSMR"]], osmr_stat)
  p_il6st <- status_bar("IL6ST", RECEPTOR_COLORS[["IL6ST"]], il6st_stat)

  status_out <- bind_rows(
    osmr_stat$tab %>% transmute(
      Gene = "OSMR",
      Predicted_Type = gsub("\n", " ", as.character(CT), fixed = TRUE),
      Status = as.character(S),
      Count = Count,
      Percent = round(Pct, 4),
      Chi_square_P = osmr_stat$pval
    ),
    il6st_stat$tab %>% transmute(
      Gene = "IL6ST",
      Predicted_Type = gsub("\n", " ", as.character(CT), fixed = TRUE),
      Status = as.character(S),
      Count = Count,
      Percent = round(Pct, 4),
      Chi_square_P = il6st_stat$pval
    )
  )
  status_file <- file.path(
    OUT_DIR,
    paste0(file_prefix, "Figure1i_OSMR_IL6ST_Status_Distribution.txt")
  )
  write.table(
    status_out, status_file, sep = "\t", row.names = FALSE,
    col.names = TRUE, quote = FALSE
  )

  fig <- p_umap + p_osmr + p_il6st + p_dot +
    plot_layout(nrow = 1, widths = c(0.72, 1.00, 1.00, 1.28))
  save_pub_plot(fig, paste0(file_prefix, "Figure1i_full"), width = 7.2, height = 2.65)
  cat("Fig1i and status table saved with prefix: ", file_prefix, "\n", sep = "")
}

# =============================================================================
# 1j.  OSMR vs IL6ST standardized effect-size bar chart (Malignant Epithelial)
#       Standardized multiple regression: OSM_Response_Z ~ OSMR_Z + IL6ST_Z
# =============================================================================
run_fig1j <- function(file_prefix = "") {
  mal <- subset(seu, subset = CT_Refined == "Malignant_Epithelial")
  stopifnot(all(as.character(mal$Analysis_CellType) == "Epithelial"))
  ed <- FetchData(mal, vars = c("OSMR","IL6ST"), layer = "counts")
  df <- data.frame(R = as.numeric(scale(mal$OSM_Response)),
                   OSMR_Z = as.numeric(scale(log1p(ed$OSMR))),
                   IL6ST_Z = as.numeric(scale(log1p(ed$IL6ST))))
  m <- lm(R ~ OSMR_Z + IL6ST_Z, data = df); sm <- summary(m)
  b <- coef(m); se <- sm$coefficients[, "Std. Error"]; vc <- vcov(m)
  se_diff <- sqrt(vc["OSMR_Z","OSMR_Z"] + vc["IL6ST_Z","IL6ST_Z"] - 2*vc["OSMR_Z","IL6ST_Z"])
  tstat <- (b["OSMR_Z"] - b["IL6ST_Z"]) / se_diff
  pval <- 2 * pt(-abs(tstat), df = m$df.residual)
  cd <- data.frame(Receptor = factor(c("OSMR","IL6ST"), levels = c("OSMR","IL6ST")),
                   Beta = c(b["OSMR_Z"], b["IL6ST_Z"]), SE = c(se["OSMR_Z"], se["IL6ST_Z"]))
  cd$CI_Low <- cd$Beta - 1.96 * cd$SE
  cd$CI_High <- cd$Beta + 1.96 * cd$SE
  bracket_y <- max(cd$CI_High) + 0.010
  bracket_tick <- 0.004
  y_max <- bracket_y + 0.014
  p <- ggplot(cd, aes(Receptor, Beta, fill = Receptor)) +
    geom_col(width = 0.58, colour = "black", linewidth = 0.25) +
    geom_errorbar(
      aes(ymin = CI_Low, ymax = CI_High),
      width = 0.16, linewidth = 0.50, colour = "black"
    ) +
    annotate(
      "segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y,
      linewidth = 0.35, colour = "grey35"
    ) +
    annotate(
      "segment", x = c(1, 2), xend = c(1, 2),
      y = bracket_y, yend = bracket_y - bracket_tick,
      linewidth = 0.35, colour = "grey35"
    ) +
    annotate(
      "text", x = 1.5, y = bracket_y + 0.002,
      label = get_stars(pval), vjust = 0, size = GEOM_TEXT_SIZE,
      family = FONT_FAMILY, fontface = "bold"
    ) +
    scale_fill_manual(values = RECEPTOR_COLORS, guide = "none") +
    scale_y_continuous(
      limits = c(0, y_max),
      breaks = seq(0, 0.20, by = 0.05),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      x = "Receptor unit",
      y = "Standardized effect size"
    ) +
    theme_nature() +
    theme(
      axis.text.x = element_text(face = "italic", size = BASE_FONT_SIZE),
      plot.margin = margin(7, 5, 4, 5)
    )
  save_pub_plot(p, paste0(file_prefix, "Figure1j_Receptor_EffectSize"), width = 3.5, height = 2.5)
  result_out <- cd %>%
    mutate(
      N_malignant_epithelial = nrow(df),
      OSMR_minus_IL6ST_P = pval
    ) %>%
    select(
      Receptor, N_malignant_epithelial, Beta, SE, CI_Low, CI_High,
      OSMR_minus_IL6ST_P
    )
  write.table(
    result_out,
    file.path(OUT_DIR, paste0(file_prefix, "Figure1j_Receptor_EffectSize_results.txt")),
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  cat(sprintf("Fig1j saved. OSMR beta=%.4f  IL6ST beta=%.4f  diff p=%.2e  %s\n",
              b["OSMR_Z"], b["IL6ST_Z"], pval, get_stars(pval)))
}


run_fig1f()
run_fig1g()
run_fig1h()
run_fig1i()
run_fig1j()
