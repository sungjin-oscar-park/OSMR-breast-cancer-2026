# Final Figure 1a-e: TCGA BRCA 2012; see docs/PANEL_MAP.md.
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(tidyr)
  library(reshape2); library(ggpubr); library(ggplot2)
  library(pheatmap); library(RColorBrewer)
  library(survival); library(survminer)
  library(patchwork); library(msigdbr); library(grid)
  library(cowplot); library(gridExtra); library(ragg)
})

DATA_DIR <- required_input(cfg, "tcga_dir")
OUT_DIR <- output_path(cfg, "figure1/bulk")
set.seed(1)

# =============================================================================
# 0.  Shared data loading
# =============================================================================
target_genes <- c("OSMR","OSM","IL6ST","IL6","CLCF1","CNTFR","CTF1",
                  "IL6R","IL11","IL11RA","IL27","IL27RA","LIF","LIFR")

gene_expr_raw <- read_tsv(file.path(DATA_DIR, "data_mrna_mirna_merged_zscores.txt"), show_col_types = FALSE)
if ("Entrez_Gene_Id" %in% colnames(gene_expr_raw)) gene_expr_raw <- gene_expr_raw %>% select(-Entrez_Gene_Id)
gnc <- colnames(gene_expr_raw)[1]
gene_symbols_raw <- as.character(gene_expr_raw[[gnc]])
gnu <- make.names(gene_expr_raw[[gnc]], unique = TRUE)
gene_expr_data <- gene_expr_raw %>% select(-all_of(gnc)) %>% as.data.frame()
rownames(gene_expr_data) <- gnu

clinical_data <- read_tsv(file.path(DATA_DIR, "clinical_data.tsv"), show_col_types = FALSE) %>%
  rename(OS = `Overall Survival (Months)`, OS.IDX = `Overall Survival Status`)

matched_samples   <- intersect(colnames(gene_expr_data), clinical_data$`Sample ID`)
available_genes   <- intersect(target_genes, rownames(gene_expr_data))
gene_expr_matched <- gene_expr_data[available_genes, matched_samples, drop = FALSE]
clinical_matched  <- clinical_data[clinical_data$`Sample ID` %in% matched_samples, ]
full_mat          <- as.matrix(gene_expr_data[, matched_samples])   # all genes x matched samples
cat("Loaded:", length(available_genes), "target genes,", length(matched_samples), "samples\n")

get_stars <- function(p) ifelse(is.na(p), "ns",
  ifelse(p < 0.0001, "****", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))))


write_excel_workbook <- function(sheets, path) {
  writexl::write_xlsx(sheets, path = path)
}

run_fig1e_metascape_deg <- function(
    target_gene = "OSMR",
    fdr_cutoff = 0.05,
    expression_difference_cutoff = 0.5,
    top_n = 500) {

  mat <- full_mat
  storage.mode(mat) <- "numeric"
  symbols <- gene_symbols_raw

  keep_nonzero <- rowSums(abs(mat), na.rm = TRUE) > 0
  mat <- mat[keep_nonzero, , drop = FALSE]
  symbols <- symbols[keep_nonzero]
  mat[!is.finite(mat)] <- 0

  target_idx <- which(symbols == target_gene)
  if (length(target_idx) == 0) stop(target_gene, " was not found in the expression matrix.")
  if (length(target_idx) > 1) {
    warning("Multiple rows matched ", target_gene, "; the first row is being used for grouping.")
    target_idx <- target_idx[1]
  }

  target_expression <- as.numeric(mat[target_idx, ])
  target_median <- median(target_expression, na.rm = TRUE)
  group <- ifelse(target_expression >= target_median, "High", "Low")
  high_idx <- which(group == "High")
  low_idx <- which(group == "Low")

  mean_high <- rowMeans(mat[, high_idx, drop = FALSE], na.rm = TRUE)
  mean_low <- rowMeans(mat[, low_idx, drop = FALSE], na.rm = TRUE)
  p_value <- vapply(seq_len(nrow(mat)), function(i) {
    high_values <- as.numeric(mat[i, high_idx])
    low_values <- as.numeric(mat[i, low_idx])
    high_values <- high_values[!is.na(high_values)]
    low_values <- low_values[!is.na(low_values)]
    if (length(high_values) < 3 || length(low_values) < 3) return(NA_real_)
    tryCatch(
      t.test(high_values, low_values, alternative = "two.sided", var.equal = FALSE)$p.value,
      error = function(e) NA_real_
    )
  }, numeric(1))

  tested <- !is.na(p_value)
  deg_all <- data.frame(
    Gene_symbol = symbols[tested],
    Mean_High = mean_high[tested],
    Mean_Low = mean_low[tested],
    Expression_Difference = mean_high[tested] - mean_low[tested],
    P_value = p_value[tested],
    stringsAsFactors = FALSE
  )
  deg_all$FDR_BH <- p.adjust(deg_all$P_value, method = "BH")

  significant_deg <- deg_all %>%
    filter(FDR_BH < fdr_cutoff) %>%
    arrange(FDR_BH, P_value)

  metascape_candidates <- significant_deg %>%
    filter(Expression_Difference > expression_difference_cutoff)

  if (nrow(metascape_candidates) < top_n) {
    stop(
      "Only ", nrow(metascape_candidates), " genes passed FDR < ", fdr_cutoff,
      " and Expression_Difference > ", expression_difference_cutoff,
      "; cannot create a Top ", top_n, " list."
    )
  }

  top_deg <- metascape_candidates %>%
    slice_head(n = top_n) %>%
    mutate(Rank = row_number(), .before = 1)

  # Keep the first worksheet to a single gene-list column so that the workbook
  # can be uploaded directly to Metascape. Remaining sheets provide the audit
  # trail and should not be selected as additional gene lists.
  metascape_input <- setNames(
    data.frame(top_deg$Gene_symbol, check.names = FALSE),
    paste0(target_gene, "_high_Top", top_n)
  )

  parameters <- data.frame(
    Parameter = c(
      "Cohort", "Target_gene", "Group_definition", "High_n", "Low_n",
      "Target_median_z_score", "Expression_unit", "Statistical_test",
      "Multiple_testing_method", "FDR_correction_family",
      "Significant_DEG_rule", "Metascape_candidate_rule", "Ranking_rule",
      "Top_n", "Target_gene_retained", "Raw_gene_rows",
      "Nonzero_genes_tested", "Significant_DEGs", "Metascape_candidates",
      "Generated_at"
    ),
    Value = c(
      "TCGA BRCA Nature 2012",
      target_gene,
      paste0(target_gene, " z-score >= median: High; < median: Low"),
      length(high_idx),
      length(low_idx),
      format(target_median, digits = 12),
      "cBioPortal mRNA/miRNA merged expression z-score",
      "Independent two-sided Welch t-test",
      "Benjamini-Hochberg",
      paste0(nrow(deg_all), " successfully tested non-zero genes"),
      paste0("FDR_BH < ", fdr_cutoff),
      paste0("FDR_BH < ", fdr_cutoff, " and Expression_Difference > ",
             expression_difference_cutoff),
      "Ascending unrounded FDR_BH, then P_value",
      top_n,
      "Yes (historical-list reproduction)",
      nrow(full_mat),
      nrow(deg_all),
      nrow(significant_deg),
      nrow(metascape_candidates),
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    ),
    stringsAsFactors = FALSE
  )

  output_file <- file.path(
    OUT_DIR,
    paste0("Figure1e_", target_gene, "_High_vs_Low_Top", top_n, "_for_Metascape.xlsx")
  )
  write_excel_workbook(
    list(
      Metascape_Input = metascape_input,
      Top500_DEG_with_statistics = top_deg,
      All_significant_DEGs = significant_deg,
      Parameters = parameters
    ),
    output_file
  )

  cat(
    "Fig1e Metascape input saved:", output_file, "\n",
    "  Tested genes:", nrow(deg_all),
    "| FDR-significant:", nrow(significant_deg),
    "| candidates:", nrow(metascape_candidates),
    "| exported:", nrow(top_deg), "\n"
  )
  invisible(list(
    file = output_file,
    metascape_input = metascape_input,
    top_deg = top_deg,
    significant_deg = significant_deg,
    parameters = parameters
  ))
}

# =============================================================================
# 1a.  IL6-family expression heatmap, samples sorted by OSM expression
run_fig1a <- function() {
  osmr_expr <- as.numeric(gene_expr_matched["OSMR", ]); names(osmr_expr) <- colnames(gene_expr_matched)
  osmr_groups <- ifelse(osmr_expr >= median(osmr_expr, na.rm = TRUE), "High", "Low")
  osmr_factor <- factor(osmr_groups, levels = c("Low","High"))

  ann <- clinical_matched %>%
    select(`Sample ID`, `ER Status`, `PR Status`, `HER2 Status`, `Converted Stage`, `Diagnosis Age`) %>%
    column_to_rownames("Sample ID")
  ann$OSMR_Expression <- osmr_factor[rownames(ann)]
  ann$Stage2 <- case_when(
    ann$`Converted Stage` %in% c("Stage I","Stage IIA","Stage IIB") ~ "Early",
    ann$`Converted Stage` %in% c("Stage IIIA","Stage IIIB","Stage IV") ~ "Advanced", TRUE ~ NA_character_)
  ann$Age2 <- case_when(ann$`Diagnosis Age` >= 65 ~ "Late", ann$`Diagnosis Age` < 65 ~ "Early", TRUE ~ NA_character_)
  ann_final <- ann %>% select(OSMR_Expression, `ER Status`, `PR Status`, `HER2 Status`, Stage2, Age2)

  osm_expr <- as.numeric(gene_expr_matched["OSM", ]); names(osm_expr) <- colnames(gene_expr_matched)
  ann_final$OSM_Expression <- factor(ifelse(osm_expr >= median(osm_expr), "High", "Low"), levels = c("Low","High"))

  cmap <- function(v, pos="#9e0142", neg="#5e4fa2", oth="#b2abd2") {
    u <- unique(v[!is.na(v)]); m <- setNames(rep(oth, length(u)), u)
    if ("Positive" %in% u) m["Positive"] <- pos; if ("Negative" %in% u) m["Negative"] <- neg; m }
  ann_colors <- list(
    OSMR_Expression = c("Low"="#487cb7","High"="#d23625"),
    OSM_Expression  = c("Low"="#abdda4","High"="#fdae61"),
    `ER Status` = cmap(ann_final$`ER Status`), `PR Status` = cmap(ann_final$`PR Status`),
    `HER2 Status` = cmap(ann_final$`HER2 Status`),
    Stage2 = c("Early"="#2166ac","Advanced"="#d6604d"), Age2 = c("Early"="#74add1","Late"="#f46d43"))

  heatmap_colors <- colorRampPalette(c("#276419","#ffffff","#8e0152"))(100)
  top_priority <- c("OSM","OSMR","IL6ST","LIFR","IL6","IL6R","LIF","CNTFR","CTF1","IL11","IL11RA","IL27","IL27RA","CLCF1")
  gene_order <- c(intersect(top_priority, rownames(gene_expr_matched)),
                  setdiff(rownames(gene_expr_matched), top_priority))
  hd <- t(scale(t(gene_expr_matched[gene_order, ]))); hd[hd > 3] <- 3; hd[hd < -3] <- -3

  sample_order_osm <- names(sort(osm_expr))    # <-- FINAL: sort by OSM expression
  pheatmap(hd[, sample_order_osm], color = heatmap_colors,
           annotation_col = ann_final[sample_order_osm, ], annotation_colors = ann_colors,
           show_rownames = TRUE, show_colnames = FALSE, cluster_rows = FALSE, cluster_cols = FALSE,
           main = "IL6 Family Heatmap (Sorted by OSM Expression)",
           filename = file.path(OUT_DIR, "Figure1a_Sample_Heatmap_IL6Family_by_OSM.pdf"), width = 12, height = 8)
  cat("Fig1a saved (OSM-sorted heatmap)\n")
}

# =============================================================================
# 1b.  Per-gene Kaplan-Meier survival (OSMR / IL6ST / IL27RA), adj. Age+Sex Cox
# =============================================================================
run_fig1b <- function() {
  osmr_expr <- as.numeric(gene_expr_matched["OSMR", ]); names(osmr_expr) <- colnames(gene_expr_matched)
  osmr_groups <- ifelse(osmr_expr >= median(osmr_expr, na.rm = TRUE), "High", "Low")
  clusters_info <- data.frame(Sample_ID = names(osmr_groups),
                              Cluster = ifelse(osmr_groups == "High", 2, 1))
  cwc <- clinical_matched %>%
    left_join(clusters_info, by = c("Sample ID" = "Sample_ID")) %>%
    mutate(Age2 = case_when(`Diagnosis Age` >= 65 ~ "Late", `Diagnosis Age` < 65 ~ "Early", TRUE ~ NA_character_))

  key_genes <- c("OSMR","IL6ST","IL27RA"); plots <- list(); stats <- list()
  for (g in key_genes) {
    te <- as.numeric(gene_expr_matched[g, ]); names(te) <- colnames(gene_expr_matched)
    grp <- ifelse(te >= median(te, na.rm = TRUE), "High", "Low")
    sd <- cwc %>%
      mutate(Target = factor(grp[`Sample ID`], levels = c("Low","High"))) %>%
      filter(!is.na(OS) & !is.na(OS.IDX) & !is.na(Target)) %>%
      mutate(OS_years = as.numeric(OS)/12, OS_event = as.numeric(OS.IDX))
    fit <- survfit(Surv(OS_years, OS_event) ~ Target, data = sd)
    cd <- sd %>% filter(!is.na(Age2) & !is.na(Sex)) %>%
      mutate(Target = relevel(factor(Target), ref="Low"),
             Age2 = relevel(factor(Age2), ref="Early"), Sex = relevel(factor(Sex), ref="Female"))
    cox <- coxph(Surv(OS_years, OS_event) ~ Target + Age2 + Sex, data = cd)
    cs <- summary(cox)
    hr <- round(cs$conf.int["TargetHigh","exp(coef)"], 2); p <- signif(cs$coefficients["TargetHigh","Pr(>|z|)"], 3)
    lr <- survdiff(Surv(OS_years, OS_event) ~ Target, data = sd)
    lrp <- signif(1 - pchisq(lr$chisq, length(lr$n) - 1), 3)
    stats[[g]] <- c(HR = hr, cox_p = p, logrank_p = lrp,
                    High = sum(sd$Target=="High"), Low = sum(sd$Target=="Low"))
    kp <- ggsurvplot(fit, data = sd, palette = c("#487cb7","#d23625"),
                     risk.table = FALSE, pval = FALSE, conf.int = TRUE, surv.median.line = "hv",
                     break.time.by = 2, xlim = c(0,20), legend.labs = c("Low","High"),
                     legend.title = g, xlab = "Time (years)", ylab = "Survival probability",
                     ggtheme = theme_classic(base_size = 12))
    kp$plot <- kp$plot + annotate("text", x = 1, y = 0.08,
                 label = paste0("adj. HR = ", hr, "\np = ", p), hjust = 0, size = 4, fontface = "italic")
    plots[[g]] <- kp$plot
  }
  fn <- file.path(OUT_DIR, "Figure1b.jpg")
  ragg::agg_jpeg(fn, width = 15, height = 4, units = "in", res = 600)
  gridExtra::grid.arrange(grobs = plots, nrow = 1)
  dev.off()
  write.csv(do.call(rbind, stats), file.path(OUT_DIR, "Figure1b_survival_statistics.csv"))
  cat("Fig1b saved. stats:\n"); print(do.call(rbind, stats))
}

run_fig1c <- function() {
  pairs <- c("OSMR","IL6ST","IL27RA","LIFR")
  osm <- as.numeric(gene_expr_matched["OSM", ])
  pl <- lapply(pairs, function(g) {
    y <- as.numeric(gene_expr_matched[g, ]); df <- data.frame(OSM = osm, Y = y)
    ct <- cor.test(osm, y, method = "spearman")
    ggplot(df, aes(OSM, Y)) +
      geom_point(color = "grey25", alpha = 0.45, size = 0.9) +
      geom_smooth(method = "lm", se = TRUE, color = "#e8502a", fill = "#fde0c8", linewidth = 0.7) +
      annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
               label = sprintf("rho = %.2f\np = %.2g", ct$estimate, ct$p.value), size = 3.5) +
      labs(x = "OSM expression (z)", y = paste0(g, " expression (z)")) +
      theme_classic(base_size = 11) + theme(axis.ticks = element_blank())
  })
  fig <- (pl[[1]] | pl[[2]]) / (pl[[3]] | pl[[4]])
  ggsave(file.path(OUT_DIR, "Figure1c.jpg"), fig, width = 6.5, height = 6, dpi = 600, device = ragg::agg_jpeg)
  cat("Fig1c saved\n")
}

# =============================================================================
# 1d.  Immune-cell deconvolution (MCP-counter via immunedeconv)
#       Exploratory relative signatures from input z-scores; nominal Wilcoxon P.
#       Requires immunedeconv with MCP-counter; errors if unavailable.
# =============================================================================
run_fig1d <- function() {
  if (!requireNamespace("immunedeconv", quietly = TRUE))
    stop("Figure 1d requires immunedeconv and MCPcounter; no alternative method is substituted.")
  gene_expr_matrix <- full_mat
  gs <- rowSums(abs(gene_expr_matrix), na.rm = TRUE)
  gene_expr_matrix <- gene_expr_matrix[gs > 0, ]; gene_expr_matrix[!is.finite(gene_expr_matrix)] <- 0
  primary_method <- "mcp_counter"
  immune_result <- immunedeconv::deconvolute(gene_expr_matrix, "mcp_counter")
  immune_df <- as.data.frame(immune_result); cell_types <- immune_df$cell_type
  immune_mat <- t(as.matrix(immune_df[, -1])); colnames(immune_mat) <- cell_types
  immune_mat[!is.finite(immune_mat)] <- 0

  summ_file <- file.path(OUT_DIR, "Immune_Cell_Deconvolution_Summary.txt")
  cat("Gene","Cell_Type","P_Value","Mean_High","Mean_Low","Fold_Change","P_Adj",
      file = summ_file, sep = "\t", append = FALSE); cat("\n", file = summ_file, append = TRUE)

  target_deconv_genes <- c("OSMR","IL6ST"); common_ct_order <- NULL; common_ct_colors <- NULL
  for (target_gene in target_deconv_genes) {
    if (!(target_gene %in% rownames(gene_expr_matched))) next
    g_expr <- as.numeric(gene_expr_matched[target_gene, ]); names(g_expr) <- colnames(gene_expr_matched)
    g_groups <- ifelse(g_expr >= median(g_expr, na.rm = TRUE), "High", "Low")
    common_samples <- intersect(rownames(immune_mat), names(g_groups)); if (length(common_samples) < 10) next
    gg <- factor(g_groups[common_samples], levels = c("Low","High")); im <- immune_mat[common_samples, ]

    st <- data.frame(cell_type=character(), p_value=numeric(), mean_high=numeric(), mean_low=numeric(), fold_change=numeric())
    for (cell in colnames(im)) {
      hv <- as.numeric(im[gg=="High", cell]); lv <- as.numeric(im[gg=="Low", cell])
      if (all(hv==hv[1]) && all(lv==lv[1])) next
      if (length(hv)>=3 & length(lv)>=3) {
        tr <- tryCatch(wilcox.test(hv, lv), error = function(e) NULL)
        if (!is.null(tr)) { mh<-mean(hv,na.rm=TRUE); ml<-mean(lv,na.rm=TRUE)
          st <- rbind(st, data.frame(cell_type=cell, p_value=tr$p.value, mean_high=mh, mean_low=ml,
                                     fold_change=ifelse(ml==0,NA,mh/ml))) }
      }
    }
    if (nrow(st)==0) next
    st$p_adj <- p.adjust(st$p_value, "fdr"); st <- st[order(st$p_value), ]
    write.table(cbind(Gene=target_gene, st), sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE, file=summ_file, append=TRUE)

    if (is.null(common_ct_order) && target_gene == "OSMR") {
      os <- st[order(st$mean_high, decreasing = TRUE), ]
      common_ct_order  <- as.character(head(os$cell_type, 10))
      pc <- brewer.pal(min(12, length(common_ct_order)+2), "Paired")
      common_ct_colors <- setNames(pc[seq_along(common_ct_order)], common_ct_order)
    }
    pc <- st[st$cell_type %in% common_ct_order, ]; pc <- pc[match(common_ct_order, pc$cell_type), ]; pc <- na.omit(pc)
    pc$log2FC <- if (primary_method=="mcp_counter") pc$mean_high - pc$mean_low
                 else log2((pc$mean_high+1e-6)/(pc$mean_low+1e-6))
    pc$log10P <- -log10(pc$p_value + 1e-10)
    pc$cell_type <- factor(pc$cell_type, levels = common_ct_order)

    pdd <- pc; pdd$cell_type <- factor(pdd$cell_type, levels=rev(common_ct_order))
    p_right <- ggplot(pdd, aes(log2FC, cell_type)) +
      geom_vline(xintercept=0, color="grey30") + geom_point(aes(size=log10P, color=cell_type)) +
      scale_color_manual(values=common_ct_colors, guide="none") + scale_size_continuous(range=c(2,8)) +
      labs(title=paste0(target_gene, " High vs Low"), x="Mean signature difference", y="") +
      theme_classic() + theme(axis.text.y=element_text(size=9))
    ggsave(file.path(OUT_DIR, paste0("Figure1d_Immune_Cell_Comparison_", target_gene, ".pdf")), p_right, width=5, height=6)
    cat("Fig1d", target_gene, "saved\n")
  }
}

# =============================================================================
# 1e.  Hallmark rank_score heatmap + OSMR-high Top 500 Metascape network input
# =============================================================================
hallmark_rank_score <- function(expr_mat, gene_sets){
  ns <- ncol(expr_mat); scores <- matrix(0, length(gene_sets), ns,
    dimnames=list(names(gene_sets), colnames(expr_mat)))
  rn_all <- rownames(expr_mat)
  for (i in seq_along(gene_sets)) {
    sg <- intersect(gene_sets[[i]], rn_all)
    if (length(sg) >= 15) { idx <- rn_all %in% sg
      for (j in 1:ns) { rk <- rank(-expr_mat[,j], ties.method="average"); scores[i,j] <- 1 - mean(rk[idx])/length(rk) } }
  }
  scores
}
run_fig1e <- function() {
  ger <- read_tsv(file.path(DATA_DIR, "data_mrna_mirna_merged_zscores.txt"), show_col_types = FALSE)
  gcol <- colnames(ger)[1]; if ("Entrez_Gene_Id" %in% colnames(ger)) ger <- ger %>% select(-Entrez_Gene_Id)
  gnu2 <- make.names(ger[[gcol]], unique = TRUE)
  rn <- make.unique(gsub("\\.\\d+$", "", gnu2))
  mat <- ger %>% select(-all_of(gcol)) %>% as.data.frame(); rownames(mat) <- rn; mat <- as.matrix(mat)
  ms <- intersect(colnames(mat), clinical_matched$`Sample ID`); mat <- mat[, ms]
  gs <- rowSums(abs(mat), na.rm = TRUE); mat <- mat[gs > 0, ]; mat[!is.finite(mat)] <- 0

  hg <- load_hallmarks(cfg)
  gene_col <- if ("gene_symbol" %in% colnames(hg)) "gene_symbol" else "human_gene_symbol"
  hallmark_list <- split(hg[[gene_col]], hg$gs_name)

  hallmark_mapping <- data.frame(
    Geneset = c(
      "HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_WNT_BETA_CATENIN_SIGNALING","HALLMARK_HEDGEHOG_SIGNALING",
      "HALLMARK_NOTCH_SIGNALING","HALLMARK_PHOSPHATIDYLINOSITOL_3_KINASE_AKT_MTOR_SIGNALING","HALLMARK_MTORC1_SIGNALING",
      "HALLMARK_KRAS_SIGNALING_UP","HALLMARK_KRAS_SIGNALING_DN","HALLMARK_TGF_BETA_SIGNALING",
      "HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_IL6_JAK_STAT3_SIGNALING",
      "HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_ALLOGRAFT_REJECTION","HALLMARK_INFLAMMATORY_RESPONSE",
      "HALLMARK_COMPLEMENT","HALLMARK_COAGULATION","HALLMARK_ADIPOGENESIS","HALLMARK_BILE_ACID_METABOLISM",
      "HALLMARK_CHOLESTEROL_HOMEOSTASIS","HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_GLYCOLYSIS",
      "HALLMARK_HEME_METABOLISM","HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_PEROXISOME",
      "HALLMARK_XENOBIOTIC_METABOLISM","HALLMARK_PROTEIN_SECRETION","HALLMARK_ANDROGEN_RESPONSE",
      "HALLMARK_ESTROGEN_RESPONSE_EARLY","HALLMARK_ESTROGEN_RESPONSE_LATE","HALLMARK_P53_PATHWAY",
      "HALLMARK_UV_RESPONSE_UP","HALLMARK_UV_RESPONSE_DN","HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
      "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY","HALLMARK_DNA_REPAIR","HALLMARK_G2M_CHECKPOINT",
      "HALLMARK_E2F_TARGETS","HALLMARK_MYC_TARGETS_V1","HALLMARK_MYC_TARGETS_V2",
      "HALLMARK_MITOTIC_SPINDLE","HALLMARK_SPERMATOGENESIS","HALLMARK_APICAL_JUNCTION",
      "HALLMARK_APICAL_SURFACE","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_ANGIOGENESIS",
      "HALLMARK_MYOGENESIS","HALLMARK_APOPTOSIS","HALLMARK_HYPOXIA","HALLMARK_PANCREAS_BETA_CELLS"),
    Category = c(rep("Signaling",9), rep("Immune",8), rep("Metabolic",10), rep("Hormone",3),
                 rep("Stress/DNA Damage",5), rep("Proliferation",6), rep("Development",6),
                 "Cell Death","Stress/DNA Damage","Development"),
    stringsAsFactors = FALSE)
  category_colors <- c("Signaling"="#E41A1C","Immune"="#377EB8","Metabolic"="#4DAF4A","Hormone"="#984EA3",
                       "Stress/DNA Damage"="#FF7F00","Proliferation"="#A65628","Development"="#999999","Cell Death"="#F781BF")
  cluster_colors <- c("1"="#E08214","2"="#4393C3","3"="#BABABA")

  cat("Fig1e: running rank-based rank_score (50 Hallmark sets x", ncol(mat), "samples)...\n")
  rank_score_scores <- hallmark_rank_score(mat, hallmark_list)
  rank_score_scores <- rank_score_scores[rowSums(rank_score_scores != 0) > 0, ]

  build_rank_score_heatmap <- function(gene_name) {
    ge <- as.numeric(mat[gene_name, ]); names(ge) <- colnames(mat)
    grp <- ifelse(ge >= median(ge, na.rm = TRUE), "High", "Low"); ss <- rank_score_scores[, names(grp)]
    res <- data.frame(Geneset = rownames(ss), stringsAsFactors = FALSE)
    res$P <- sapply(seq_len(nrow(ss)), function(i) tryCatch(t.test(ss[i, grp=="High"], ss[i, grp=="Low"])$p.value, error=function(e) NA))
    res$FDR <- p.adjust(res$P, "fdr"); sig <- res$Geneset[which(res$FDR < 0.05)]
    hd <- ss[sig, ]
    ra <- hallmark_mapping[match(sig, hallmark_mapping$Geneset), ]
    row_ann <- data.frame(Category = ra$Category, row.names = gsub("_", " ", gsub("HALLMARK_", "", sig)))
    rownames(hd) <- gsub("_", " ", gsub("HALLMARK_", "", sig))
    cc <- hclust(as.dist(1 - cor(hd)), method = "complete"); clus <- cutree(cc, k = 3)
    sizes <- sort(table(clus)); lab <- rep("3", length(clus)); names(lab) <- names(clus)
    lab[clus == as.integer(names(sizes)[1])] <- "1"; lab[clus == as.integer(names(sizes)[2])] <- "2"
    col_ann <- data.frame(Cluster = factor(lab[colnames(hd)], levels = c("1","2","3")),
                          `expression` = factor(grp[colnames(hd)], levels = c("High","Low")),
                          row.names = colnames(hd), check.names = FALSE)
    colnames(col_ann)[2] <- paste0(gene_name, " expression")
    ann_cols <- list(Category = category_colors, Cluster = cluster_colors)
    ann_cols[[paste0(gene_name, " expression")]] <- c("High"="#d73027","Low"="#4575b4")
    fn <- file.path(OUT_DIR, paste0("Figure1e_", gene_name, "_rank_score_heatmap.jpg"))
    pheatmap(hd, annotation_col = col_ann, annotation_row = row_ann, annotation_colors = ann_cols,
             cluster_rows = TRUE, cluster_cols = cc, scale = "row",
             color = colorRampPalette(c("#4575b4","white","#d73027"))(100), breaks = seq(-3, 3, length.out = 101),
             show_colnames = FALSE, fontsize_row = 7, main = paste0(gene_name, " expression"),
             filename = fn, width = 10, height = max(6, nrow(hd) * 0.28))
    cat(sprintf("Fig1e %s: %d significant sets -> %s\n", gene_name, length(sig), basename(fn)))
  }
  build_rank_score_heatmap("OSMR")
  run_fig1e_metascape_deg(target_gene = "OSMR")
}


run_fig1a()
run_fig1b()
run_fig1c()
run_fig1d()
run_fig1e()
