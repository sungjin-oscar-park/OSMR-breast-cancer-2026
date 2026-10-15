# Figure 4a-e statistics only: baseline tests, moderated effects/CI and Hallmark GSEA.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(tibble); library(fgsea) })
OUT <- output_path(cfg, "figure4")
dir.create(file.path(OUT, "tables"), showWarnings = FALSE)
set.seed(1)
LINE  <- c(BT549 = "BT-549", MDA231 = "MDA-MB-231", MDA468 = "MDA-MB-468",
           MB231 = "MDA-MB-231", MB468 = "MDA-MB-468")
SHORT <- c(`BT-549` = "BT-549", `MDA-MB-231` = "MB-231", `MDA-MB-468` = "MB-468")

BASE <- c("OSMR", "IL6ST", "CDKN1A", "EGFR", "CLDN4", "CLDN7")
nt   <- acute$metadata %>% filter(Condition == "NT")
bl   <- as.data.frame(t(acute$log_expr[BASE, nt$Sample, drop = FALSE])) %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "Gene", values_to = "expr") %>%
  left_join(nt, by = "Sample") %>%
  mutate(Line = LINE[CellLine])

bl_stat <- bl %>% group_by(Gene) %>%
  summarise(pairs = list({
    g <- split(expr, Line)
    cb <- combn(names(g), 2, simplify = FALSE)
    bind_rows(lapply(cb, function(k) {
      tt <- t.test(g[[k[1]]], g[[k[2]]])
      data.frame(A = k[1], B = k[2], diff = unname(diff(rev(tt$estimate))),
                 p = tt$p.value)
    }))
  }), .groups = "drop") %>% unnest(pairs) %>%
  mutate(p_holm = p.adjust(p, method = "holm"))

bl_mean <- bl %>% group_by(Gene, Line) %>%
  summarise(mean = mean(expr), sd = sd(expr), n = n(), .groups = "drop")

write.csv(bl,      file.path(OUT, "tables", "B1_baseline_values.csv"),  row.names = FALSE)
write.csv(bl_mean, file.path(OUT, "tables", "B2_baseline_means.csv"),   row.names = FALSE)
write.csv(bl_stat, file.path(OUT, "tables", "B3_baseline_pairwise.csv"), row.names = FALSE)
cat("== baseline means (log2 TPM+1) ==\n")
print(bl_mean %>% select(-sd, -n) %>% pivot_wider(names_from = Line, values_from = mean) %>%
        as.data.frame(), digits = 3)

##===========================================================================
## A2. Effect estimates with 95% CI   (Figure 5f pairing)
##===========================================================================
ci_from_fit <- function(fit, coefs, labels) {
  bind_rows(lapply(seq_along(coefs), function(i) {
    cf <- coefs[i]
    se <- fit$stdev.unscaled[, cf] * sqrt(fit$s2.post)
    tcrit <- qt(0.975, df = fit$df.total)
    data.frame(Gene = rownames(fit), Contrast = labels[i],
               logFC = fit$coefficients[, cf],
               se = se, lo = fit$coefficients[, cf] - tcrit * se,
               hi = fit$coefficients[, cf] + tcrit * se,
               P = fit$p.value[, cf],
               FDR = p.adjust(fit$p.value[, cf], "BH"),
               row.names = NULL)
  }))
}

acute_coefs <- c("BT549_6h_vs_NT", "BT549_24h_vs_NT", "MDA231_6h_vs_NT",
                 "MDA231_24h_vs_NT", "MDA468_6h_vs_NT", "MDA468_24h_vs_NT")
ci_acute <- ci_from_fit(de$fit2, acute_coefs, acute_coefs)
oe_coefs <- colnames(deoe$fit2)
ci_oe    <- ci_from_fit(deoe$fit2, oe_coefs, oe_coefs)
ci_all   <- bind_rows(ci_acute, ci_oe)
saveRDS(ci_all, file.path(OUT, "rds", "ci_all.rds"))

egfr_ci <- ci_all %>% filter(Gene %in% c("EGFR", "CDKN1A"))
write.csv(egfr_ci, file.path(OUT, "tables", "B4_EGFR_CDKN1A_effect_CI.csv"), row.names = FALSE)
cat("\n== EGFR effect estimates with 95% CI ==\n")
print(egfr_ci %>% filter(Gene == "EGFR") %>%
        select(Contrast, logFC, lo, hi, FDR) %>% as.data.frame(), digits = 3)

hs <- load_hallmarks(cfg)
pw <- split(hs$gene_symbol, hs$gs_name)

run_gsea <- function(d, label) {
  d <- d %>% filter(!is.na(Gene)) %>% group_by(Gene) %>%
    slice_max(abs(logFC / se), n = 1, with_ties = FALSE) %>% ungroup()
  stats <- setNames(d$logFC / d$se, d$Gene)
  fgsea(pw, stats, minSize = 10, maxSize = 500, nPermSimple = 10000,
        BPPARAM = BiocParallel::SerialParam()) %>%
    as.data.frame() %>% mutate(Contrast = label,
                               leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
}
gsea <- bind_rows(lapply(unique(ci_all$Contrast), function(k)
  run_gsea(ci_all %>% filter(Contrast == k), k))) %>%
  mutate(Hallmark = sub("^HALLMARK_", "", pathway),
         CellLine = sub("_.*", "", Contrast),
         Line = LINE[CellLine],
         Time = ifelse(grepl("Vector|Response|OE", Contrast), "autocrine",
                       sub(".*_(\\d+h)_vs_NT", "\\1", Contrast)))
write.csv(gsea, file.path(OUT, "tables", "B10_hallmark_gsea.csv"), row.names = FALSE)

cat("\n== IL6_JAK_STAT3 across all contrasts ==\n")
print(gsea %>% filter(Hallmark == "IL6_JAK_STAT3_SIGNALING") %>%
        select(Line, Time, NES, padj) %>% as.data.frame(), digits = 3)

saveRDS(list(bl = bl, bl_mean = bl_mean, bl_stat = bl_stat, ci_all = ci_all,
             gsea = gsea, log_acute = acute$log_expr, meta_acute = acute$metadata,
             log_oe = oe$log_expr, meta_oe = oe$metadata, LINE = LINE, SHORT = SHORT),
        file.path(OUT, "rds/figure4_objects.rds"))
