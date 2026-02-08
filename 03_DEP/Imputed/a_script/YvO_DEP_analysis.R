################################################################################
#
#   YvO Differential Expression of Proteins — Full Analysis Pipeline
#
#   Sections:
#     0. Setup & packages
#     1. Data loading & metadata construction
#     2. proteoDA: DAList → limma model → contrast extraction
#     3. Xiao et al. π-value scoring
#     4. Hallmark pathway enrichment (fgsea)
#     5. Volcano + GSEA flanking plots
#     6. GO:BP enrichment
#     7. Summary stacked bar plots (proteins & pathways per contrast)
#     8. Euler diagrams (protein-level & pathway-level overlap)
#     9. Custom UpSet — protein intersections + GO:BP set sizes
#    10. Composite summary figure
#    11. Export & session info
#
#   Reference: Xiao et al. (2014) Bioinformatics 30(6):801-807
#              doi:10.1093/bioinformatics/btt598
#
#   Run from: 03_DEP/Imputed/a_script/
#   Usage:    Rscript YvO_DEP_analysis.R
#
################################################################################

cat("=== YvO DEP Analysis Pipeline ===\n\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 0.  SETUP & PACKAGES
# ═══════════════════════════════════════════════════════════════════════════════

cat(">> 0 -- Setup\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(proteoDA)
  library(fgsea)
  library(msigdbr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(eulerr)
  library(patchwork)
  library(ggrepel)
  library(scales)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
base_dir    <- normalizePath(file.path(dirname(getwd()), "..", ".."), mustWork = TRUE)
DATA_FILE   <- file.path(base_dir, "02_Imputation", "c_data", "01_imputed.csv")
REPORT_DIR  <- file.path(base_dir, "03_DEP", "Imputed", "b_reports")
DATA_DIR    <- file.path(base_dir, "03_DEP", "Imputed", "c_data")
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR,   recursive = TRUE, showWarnings = FALSE)

cat("   Base directory:", base_dir, "\n")

# ── Aesthetics ────────────────────────────────────────────────────────────────
theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  strip.text = element_text(face = "bold")))

CONTRAST_META <- tribble(
  ~contrast,          ~label,                    ~color,
  "Training_Young",   "Training (Young)",        "#2166AC",
  "Training_Old",     "Training (Old)",          "#B2182B",
  "Aging",            "Aging (Old - Young)",     "#4DAF4A",
  "Interaction",      "Age x Training",          "#984EA3"
)

UP_COL <- "#D73027"
DN_COL <- "#4575B4"
NS_COL <- "grey75"

# ═══════════════════════════════════════════════════════════════════════════════
# 1.  DATA LOADING & METADATA CONSTRUCTION
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n>> 1 -- Load & Validate Data\n")

df <- read_csv(DATA_FILE, show_col_types = FALSE)
cat(sprintf("   Loaded: %d proteins x %d columns\n", nrow(df), ncol(df)))

ann_cols   <- c("uniprot_id", "protein", "gene", "description")
ann        <- df[, ann_cols]
samp_names <- setdiff(names(df), ann_cols)
mat        <- as.matrix(df[, samp_names])
rownames(mat) <- ann$uniprot_id

cat(sprintf("   Matrix: %d proteins x %d samples\n", nrow(mat), ncol(mat)))

meta <- tibble(sample_id = samp_names) |>
  mutate(
    prefix   = str_extract(sample_id, "^[A-Z]+"),
    subj_num = str_extract(sample_id, "S\\d+"),
    time     = str_extract(sample_id, "(Pre|Post)$"),
    age      = if_else(str_detect(prefix, "^O"), "Old", "Young"),
    subject  = paste0(prefix, "_", subj_num),
    group    = paste(age, time, sep = "_")
  )

meta$age   <- factor(meta$age,  levels = c("Young", "Old"))
meta$time  <- factor(meta$time, levels = c("Pre", "Post"))
meta$group <- factor(meta$group,
                     levels = c("Young_Pre", "Young_Post",
                                "Old_Pre",   "Old_Post"))

cat("\n   Sample distribution:\n")
print(table(meta$age, meta$time))
stopifnot(identical(colnames(mat), meta$sample_id))

# ═══════════════════════════════════════════════════════════════════════════════
# 2.  proteoDA: DAList -> LIMMA MODEL -> CONTRASTS
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n>> 2 -- Limma Model & Contrasts\n")

meta_df <- as.data.frame(meta)
rownames(meta_df) <- meta$sample_id
dal <- DAList(
  data       = mat,
  annotation = as.data.frame(ann),
  metadata   = meta_df
)

design <- model.matrix(~ 0 + group, data = meta)
colnames(design) <- levels(meta$group)

dupcor <- duplicateCorrelation(mat, design, block = meta$subject)
cat(sprintf("   Within-subject correlation: %.3f\n", dupcor$consensus.correlation))

fit <- lmFit(mat, design,
             block = meta$subject,
             correlation = dupcor$consensus.correlation)

contrast_matrix <- makeContrasts(
  Training_Young = Young_Post - Young_Pre,
  Training_Old   = Old_Post   - Old_Pre,
  Aging          = Old_Pre    - Young_Pre,
  Interaction    = (Old_Post - Old_Pre) - (Young_Post - Young_Pre),
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

contrast_names <- colnames(contrast_matrix)

results_list <- map(contrast_names, function(cname) {
  tt <- topTable(fit2, coef = cname, number = Inf, sort.by = "none")
  tt$uniprot_id <- rownames(tt)
  tt$contrast   <- cname
  left_join(tt, ann, by = "uniprot_id")
}) |> set_names(contrast_names)

all_results <- bind_rows(results_list)
write_csv(all_results, file.path(DATA_DIR, "limma_all_contrasts.csv"))

cat("\n   Contrast summary (FDR < 0.05):\n")
all_results |>
  group_by(contrast) |>
  summarise(n_sig = sum(adj.P.Val < 0.05),
            n_up  = sum(adj.P.Val < 0.05 & logFC > 0),
            n_dn  = sum(adj.P.Val < 0.05 & logFC < 0),
            .groups = "drop") |>
  print()

# ═══════════════════════════════════════════════════════════════════════════════
# 3.  XIAO ET AL. pi-VALUE SCORING
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n>> 3 -- Pi-Value Scoring\n")

compute_pi <- function(logFC, pval) {
  phi       <- abs(logFC)
  neglog10p <- -log10(pval)
  pi_score  <- phi * neglog10p
  Pi_score  <- 10^(-pi_score)
  tibble(phi = phi, neglog10p = neglog10p,
         pi_score = pi_score, Pi_score = Pi_score)
}

results_list <- map(results_list, function(tt) {
  pi_out <- compute_pi(tt$logFC, tt$P.Value)
  bind_cols(tt, pi_out) |>
    mutate(
      sig_pi  = Pi_score < 0.05,
      sig_fdr = adj.P.Val < 0.05,
      direction = case_when(
        sig_pi & logFC > 0 ~ "Up",
        sig_pi & logFC < 0 ~ "Down",
        TRUE               ~ "NS"
      ),
      shape_class = if_else(sig_fdr, "FDR < 0.05", "FDR >= 0.05")
    )
})

pi_results <- bind_rows(results_list)
write_csv(pi_results, file.path(DATA_DIR, "limma_pi_scored.csv"))

for (cn in contrast_names) {
  tt <- results_list[[cn]]
  cat(sprintf("   %s: %d up, %d down (Pi < 0.05)\n",
              cn, sum(tt$direction == "Up"), sum(tt$direction == "Down")))
}

# ═══════════════════════════════════════════════════════════════════════════════
# 4.  HALLMARK ENRICHMENT (fgsea)
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n>> 4 -- Hallmark Enrichment (fgsea)\n")

hallmark_df <- msigdbr(species = "Homo sapiens", collection = "H") |>
  dplyr::select(gs_name, gene_symbol)
hallmark_sets <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)
cat(sprintf("   Hallmark gene sets: %d\n", length(hallmark_sets)))

fgsea_results <- map(contrast_names, function(cname) {
  tt <- results_list[[cname]]
  ranks <- setNames(sign(tt$logFC) * tt$pi_score, tt$gene)
  ranks <- ranks[!is.na(names(ranks)) & names(ranks) != ""]
  ranks <- ranks[!duplicated(names(ranks))]
  ranks <- sort(ranks, decreasing = TRUE)
  res <- fgsea(pathways = hallmark_sets, stats = ranks,
               minSize = 10, maxSize = 500, nPermSimple = 10000)
  res$contrast <- cname
  res
}) |> set_names(contrast_names)

all_fgsea <- bind_rows(fgsea_results) |>
  mutate(
    pathway_short = str_remove(pathway, "^HALLMARK_") |>
      str_replace_all("_", " ") |> str_to_title(),
    direction = if_else(NES > 0, "Up", "Down"),
    sig       = padj < 0.05
  )

write_csv(all_fgsea |> dplyr::select(-leadingEdge),
          file.path(DATA_DIR, "hallmark_fgsea_results.csv"))

cat(sprintf("   Significant pathways: %d (across all contrasts)\n", sum(all_fgsea$sig)))

# ═══════════════════════════════════════════════════════════════════════════════
# 5.  VOLCANO + GSEA FLANKING PLOTS
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n>> 5 -- Volcano + GSEA Flanking Plots\n")

# Volcano function ─────────────────────────────────────────────────────────────
plot_pi_volcano <- function(tt, contrast_label) {

  top_genes <- tt |>
    filter(sig_pi) |>
    group_by(direction) |>
    slice_max(pi_score, n = 10) |>
    ungroup()

  n_up <- sum(tt$direction == "Up")
  n_dn <- sum(tt$direction == "Down")

  xmax <- max(abs(tt$logFC), na.rm = TRUE) + 0.5
  ymax <- max(tt$neglog10p, na.rm = TRUE) * 1.05
  curve_x <- seq(0.01, xmax, length.out = 500)
  curve_y <- (-log10(0.05)) / curve_x
  curve_df <- bind_rows(
    tibble(x =  curve_x, y = curve_y),
    tibble(x = -curve_x, y = curve_y)
  ) |> filter(y <= ymax, y >= 0)

  ggplot(tt, aes(x = logFC, y = neglog10p)) +
    geom_point(data = filter(tt, direction == "NS"),
               aes(shape = shape_class),
               color = NS_COL, size = 1.1, alpha = 0.35) +
    geom_point(data = filter(tt, direction == "Up"),
               aes(shape = shape_class),
               color = UP_COL, size = 1.6, alpha = 0.7) +
    geom_point(data = filter(tt, direction == "Down"),
               aes(shape = shape_class),
               color = DN_COL, size = 1.6, alpha = 0.7) +
    geom_line(data = curve_df, aes(x = x, y = y),
              color = "grey30", linetype = "dashed", linewidth = 0.45,
              inherit.aes = FALSE) +
    geom_text_repel(data = top_genes, aes(label = gene),
                    size = 2.6, max.overlaps = 25, segment.size = 0.25,
                    min.segment.length = 0.1, box.padding = 0.3) +
    scale_shape_manual(
      values = c("FDR < 0.05" = 15, "FDR >= 0.05" = 16),
      name   = NULL
    ) +
    annotate("text", x = xmax * 0.7, y = ymax * 0.95,
             label = paste0("Up: ", n_up), color = UP_COL,
             fontface = "bold", size = 4) +
    annotate("text", x = -xmax * 0.7, y = ymax * 0.95,
             label = paste0("Dn: ", n_dn), color = DN_COL,
             fontface = "bold", size = 4) +
    labs(title    = contrast_label,
         subtitle = expression("Dashed curve:" ~~ Pi == 0.05 ~~ "(Xiao et al. 2014)"),
         x = expression(log[2]~"Fold Change"),
         y = expression(-log[10]~italic(p))) +
    theme(legend.position  = "bottom",
          plot.title    = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"))
}

# GSEA flanking bar function ───────────────────────────────────────────────────
plot_gsea_flank <- function(fgsea_contrast_df, dir, bar_col, n_top = 10) {
  sig_df <- fgsea_contrast_df |>
    filter(sig, direction == dir) |>
    arrange(desc(abs(NES))) |>
    head(n_top)

  if (nrow(sig_df) == 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = paste0("No significant\n", tolower(dir), " pathways"),
                      size = 3, color = "grey50") +
             theme_void())
  }

  sig_df$pathway_short <- factor(sig_df$pathway_short,
                                 levels = rev(sig_df$pathway_short))

  ggplot(sig_df, aes(x = NES, y = pathway_short)) +
    geom_col(fill = bar_col, width = 0.7) +
    geom_text(aes(label = sprintf("%.2f", NES),
                  hjust = if_else(NES > 0, -0.1, 1.1)),
              size = 2.5, fontface = "bold") +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    labs(x = "NES", y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.text.y  = element_text(size = 7),
          plot.margin  = margin(5, 5, 5, 5))
}

# Build one row per contrast: [Down GSEA | Volcano | Up GSEA] ─────────────────
volcano_gsea_rows <- map(contrast_names, function(cname) {
  lbl <- CONTRAST_META$label[CONTRAST_META$contrast == cname]
  fgsea_cn <- all_fgsea |> filter(contrast == cname)

  p_dn   <- plot_gsea_flank(fgsea_cn, "Down", DN_COL)
  p_volc <- plot_pi_volcano(results_list[[cname]], lbl)
  p_up   <- plot_gsea_flank(fgsea_cn, "Up",   UP_COL)

  (p_dn | p_volc | p_up) + plot_layout(widths = c(2, 3, 2))
})

p_volcano_gsea <- wrap_plots(volcano_gsea_rows, ncol = 1) +
  plot_annotation(
    title = "Volcano Plots with Hallmark GSEA Flanking (FDR < 0.05)",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15))
  )

ggsave(file.path(REPORT_DIR, "01_volcano_gsea.pdf"),
       p_volcano_gsea, width = 20, height = 28)
ggsave(file.path(REPORT_DIR, "01_volcano_gsea.png"),
       p_volcano_gsea, width = 20, height = 28, dpi = 300)
cat("   Saved: 01_volcano_gsea.pdf\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 6.  GO:BP ENRICHMENT
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n>> 6 -- GO:BP Enrichment\n")

go_bp_results <- map(contrast_names, function(cname) {
  tt <- results_list[[cname]]
  bg <- tt$gene[!is.na(tt$gene) & tt$gene != ""]

  map(c("Up", "Down"), function(dir) {
    genes <- tt |> filter(direction == dir) |> pull(gene)
    genes <- genes[!is.na(genes) & genes != ""]
    if (length(genes) < 3) return(NULL)

    res <- tryCatch(
      enrichGO(gene = genes, universe = bg, OrgDb = org.Hs.eg.db,
               keyType = "SYMBOL", ont = "BP",
               pAdjustMethod = "BH", pvalueCutoff = 0.05,
               qvalueCutoff = 0.1, minGSSize = 5, maxGSSize = 500),
      error = function(e) NULL
    )
    if (is.null(res) || nrow(res@result |> filter(p.adjust < 0.05)) == 0) return(NULL)

    res@result |> filter(p.adjust < 0.05) |>
      mutate(contrast = cname, direction = dir)
  }) |> compact() |> bind_rows()
}) |> compact() |> bind_rows()

write_csv(go_bp_results, file.path(DATA_DIR, "go_bp_ora_results.csv"))
cat(sprintf("   GO:BP enrichment: %d significant terms\n", nrow(go_bp_results)))

# ═══════════════════════════════════════════════════════════════════════════════
# 7.  SUMMARY STACKED BAR PLOTS
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n>> 7 -- Summary Bar Plots\n")

protein_counts <- map_dfr(results_list, function(tt) {
  tibble(contrast = tt$contrast[1],
         Up   = sum(tt$direction == "Up"),
         Down = sum(tt$direction == "Down"))
}) |>
  pivot_longer(c(Up, Down), names_to = "direction", values_to = "count") |>
  mutate(
    contrast  = factor(contrast, levels = contrast_names),
    direction = factor(direction, levels = c("Up", "Down")),
    y_val     = if_else(direction == "Down", -count, count)
  )

p_prot_bar <- ggplot(protein_counts, aes(x = contrast, y = y_val, fill = direction)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_text(aes(label = count,
                vjust = if_else(direction == "Up", -0.3, 1.3)),
            size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("Up" = UP_COL, "Down" = DN_COL)) +
  scale_x_discrete(labels = CONTRAST_META$label) +
  labs(title = expression("Significant Proteins per Contrast" ~~ (Pi < 0.05)),
       x = NULL, y = "Number of proteins", fill = NULL) +
  coord_flip() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"))

pathway_counts <- all_fgsea |>
  filter(sig) |>
  group_by(contrast, direction) |>
  summarise(count = n(), .groups = "drop") |>
  complete(contrast = contrast_names, direction = c("Up", "Down"),
           fill = list(count = 0)) |>
  mutate(
    contrast  = factor(contrast, levels = contrast_names),
    direction = factor(direction, levels = c("Up", "Down")),
    y_val     = if_else(direction == "Down", -count, count)
  )

p_path_bar <- ggplot(pathway_counts, aes(x = contrast, y = y_val, fill = direction)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_text(aes(label = count,
                vjust = if_else(direction == "Up", -0.3, 1.3)),
            size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("Up" = UP_COL, "Down" = DN_COL)) +
  scale_x_discrete(labels = CONTRAST_META$label) +
  labs(title = "Significant Hallmark Pathways per Contrast (FDR < 0.05)",
       x = NULL, y = "Number of pathways", fill = NULL) +
  coord_flip() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(REPORT_DIR, "04_protein_counts_bar.pdf"),  p_prot_bar, width = 8, height = 5)
ggsave(file.path(REPORT_DIR, "04_pathway_counts_bar.pdf"),  p_path_bar, width = 8, height = 5)
ggsave(file.path(REPORT_DIR, "04_protein_counts_bar.png"),  p_prot_bar, width = 8, height = 5, dpi = 300)
ggsave(file.path(REPORT_DIR, "04_pathway_counts_bar.png"),  p_path_bar, width = 8, height = 5, dpi = 300)
cat("   Saved: 04_protein_counts_bar.pdf, 04_pathway_counts_bar.pdf\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 8.  EULER DIAGRAMS
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n>> 8 -- Euler Diagrams\n")

sig_genes <- map(results_list, function(tt) {
  list(Up   = tt$gene[tt$direction == "Up"],
       Down = tt$gene[tt$direction == "Down"])
})

up_gene_lists <- map(sig_genes, "Up") |> set_names(contrast_names)
dn_gene_lists <- map(sig_genes, "Down") |> set_names(contrast_names)
up_gene_lists <- up_gene_lists[lengths(up_gene_lists) > 0]
dn_gene_lists <- dn_gene_lists[lengths(dn_gene_lists) > 0]

up_path_lists <- list(); dn_path_lists <- list()
for (cname in contrast_names) {
  up_p <- all_fgsea |> filter(contrast == cname, sig, direction == "Up") |> pull(pathway)
  dn_p <- all_fgsea |> filter(contrast == cname, sig, direction == "Down") |> pull(pathway)
  if (length(up_p) > 0) up_path_lists[[cname]] <- up_p
  if (length(dn_p) > 0) dn_path_lists[[cname]] <- dn_p
}

make_euler_plot <- function(item_lists, title_text) {
  if (length(item_lists) < 2) {
    return(ggplot() + annotate("text", x = .5, y = .5, label = "< 2 sets", size = 3) +
             theme_void() + ggtitle(title_text))
  }
  cols <- CONTRAST_META$color[match(names(item_lists), CONTRAST_META$contrast)]
  fit  <- euler(item_lists)
  plot(fit,
       quantities = list(type = "counts", cex = 0.9),
       fills      = list(fill = cols, alpha = 0.45),
       labels     = list(cex = 0.65),
       main       = list(label = title_text, cex = 0.9))
}

pdf(file.path(REPORT_DIR, "03_euler_proteins_up.pdf"), width = 8, height = 6)
if (length(up_gene_lists) >= 2) plot(make_euler_plot(up_gene_lists, "Up-Regulated Proteins")) else
  { plot.new(); text(0.5, 0.5, "< 2 contrasts with up-regulated proteins") }
dev.off()

pdf(file.path(REPORT_DIR, "03_euler_proteins_down.pdf"), width = 8, height = 6)
if (length(dn_gene_lists) >= 2) plot(make_euler_plot(dn_gene_lists, "Down-Regulated Proteins")) else
  { plot.new(); text(0.5, 0.5, "< 2 contrasts with down-regulated proteins") }
dev.off()

pdf(file.path(REPORT_DIR, "03_euler_pathways_up.pdf"), width = 8, height = 6)
if (length(up_path_lists) >= 2) plot(make_euler_plot(up_path_lists, "Up Hallmark Pathways")) else
  { plot.new(); text(0.5, 0.5, "< 2 contrasts with up pathways") }
dev.off()

pdf(file.path(REPORT_DIR, "03_euler_pathways_down.pdf"), width = 8, height = 6)
if (length(dn_path_lists) >= 2) plot(make_euler_plot(dn_path_lists, "Down Hallmark Pathways")) else
  { plot.new(); text(0.5, 0.5, "< 2 contrasts with down pathways") }
dev.off()

cat("   Saved: 03_euler_*.pdf\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 9.  CUSTOM UPSET — PROTEIN INTERSECTIONS + GO:BP SET SIZES
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n>> 9 -- UpSet Plot\n")

# 9a. Build protein membership matrix ─────────────────────────────────────────
upset_df <- tibble(gene = results_list[[1]]$gene)
for (cn in contrast_names) {
  upset_df[[cn]] <- results_list[[cn]]$sig_pi
}

# Keep only proteins significant in at least one contrast
upset_df <- upset_df |> filter(if_any(all_of(contrast_names), ~.))

# Assign direction: use contrast with largest |logFC| among significant
upset_df$direction <- map_chr(seq_len(nrow(upset_df)), function(i) {
  best_fc <- 0
  for (cn in contrast_names) {
    if (upset_df[[cn]][i]) {
      tt <- results_list[[cn]]
      idx <- which(tt$gene == upset_df$gene[i])
      if (length(idx) > 0 && abs(tt$logFC[idx[1]]) > abs(best_fc)) {
        best_fc <- tt$logFC[idx[1]]
      }
    }
  }
  if (best_fc > 0) "Up" else "Down"
})

# Create pattern string
upset_df$pattern <- apply(upset_df[contrast_names], 1, function(x) {
  paste(as.integer(x), collapse = "-")
})

# 9b. Intersection counts ─────────────────────────────────────────────────────
isect <- upset_df |>
  group_by(pattern, direction) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(pattern) |>
  mutate(total = sum(n)) |>
  ungroup()

# Sort by total size, keep top 15
pattern_order <- isect |>
  distinct(pattern, total) |>
  arrange(desc(total)) |>
  head(15) |>
  pull(pattern)

isect <- isect |> filter(pattern %in% pattern_order)
isect$pattern <- factor(isect$pattern, levels = pattern_order)
isect$direction <- factor(isect$direction, levels = c("Up", "Down"))

# Totals for labels
isect_totals <- isect |>
  group_by(pattern) |>
  summarise(total = sum(n), .groups = "drop")
isect_totals$pattern <- factor(isect_totals$pattern, levels = pattern_order)

# 9c. Dot matrix data ─────────────────────────────────────────────────────────
dot_data <- map_dfr(pattern_order, function(pat) {
  bits <- as.integer(strsplit(pat, "-")[[1]])
  tibble(pattern = pat, contrast = contrast_names, active = as.logical(bits))
})
dot_data$pattern  <- factor(dot_data$pattern, levels = pattern_order)
dot_data$contrast <- factor(dot_data$contrast, levels = rev(contrast_names))

# Connecting lines between active dots
line_data <- dot_data |>
  filter(active) |>
  group_by(pattern) |>
  summarise(ymin = min(as.numeric(contrast)),
            ymax = max(as.numeric(contrast)),
            .groups = "drop") |>
  filter(ymin != ymax)
line_data$pattern <- factor(line_data$pattern, levels = pattern_order)

# 9d. GO:BP set sizes ─────────────────────────────────────────────────────────
go_bp_set_counts <- go_bp_results |>
  group_by(contrast, direction) |>
  summarise(n = n(), .groups = "drop") |>
  complete(contrast = contrast_names, direction = c("Up", "Down"),
           fill = list(n = 0))
go_bp_set_counts$contrast <- factor(go_bp_set_counts$contrast,
                                    levels = rev(contrast_names))
go_bp_set_counts$direction <- factor(go_bp_set_counts$direction,
                                     levels = c("Down", "Up"))

# 9e. Build the 3 panels ──────────────────────────────────────────────────────

# Panel 1: Intersection bars (stacked Up/Down)
p_isect <- ggplot(isect, aes(x = pattern, y = n, fill = direction)) +
  geom_col(width = 0.7) +
  geom_text(data = isect_totals, aes(x = pattern, y = total, label = total),
            inherit.aes = FALSE, vjust = -0.3, size = 3, fontface = "bold") +
  scale_fill_manual(values = c("Up" = UP_COL, "Down" = DN_COL), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(y = "Protein count") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x  = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "top",
        legend.key.size = unit(0.35, "cm"))

# Panel 2: Dot matrix
p_dots <- ggplot() +
  geom_segment(data = line_data,
               aes(x = pattern, xend = pattern, y = ymin, yend = ymax),
               linewidth = 0.8, color = "grey30") +
  geom_point(data = dot_data |> filter(!active),
             aes(x = pattern, y = contrast), color = "grey80", size = 2.5) +
  geom_point(data = dot_data |> filter(active),
             aes(x = pattern, y = contrast), color = "grey20", size = 3.5) +
  scale_y_discrete(labels = setNames(CONTRAST_META$label, CONTRAST_META$contrast)) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x  = element_blank(),
        axis.title   = element_blank(),
        axis.ticks   = element_blank(),
        panel.grid   = element_blank())

# Panel 3: GO:BP pathway counts per contrast (stacked Up/Down)
p_gobp_sets <- ggplot(go_bp_set_counts, aes(x = n, y = contrast, fill = direction)) +
  geom_col(width = 0.6) +
  geom_text(data = go_bp_set_counts |>
              group_by(contrast) |>
              summarise(total = sum(n), .groups = "drop"),
            aes(x = total, y = contrast, label = total),
            inherit.aes = FALSE, hjust = -0.2, size = 3, fontface = "bold") +
  scale_fill_manual(values = c("Up" = UP_COL, "Down" = DN_COL), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = "GO:BP pathways", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_blank(),
        panel.grid.minor = element_blank())

# 9f. Combine with patchwork ──────────────────────────────────────────────────
# Layout:  [  empty  | intersection bars ]
#          [ GO:BP   | dot matrix        ]
design_layout <- "
#BBBBBBB
ACCCCCCC
"

p_upset_custom <- p_gobp_sets + p_isect + p_dots +
  plot_layout(design = design_layout, heights = c(3, 2)) +
  plot_annotation(
    title = "Protein Intersection UpSet with GO:BP Pathway Counts",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )

ggsave(file.path(REPORT_DIR, "05_upset_stacked.pdf"),
       p_upset_custom, width = 14, height = 9)
ggsave(file.path(REPORT_DIR, "05_upset_stacked.png"),
       p_upset_custom, width = 14, height = 9, dpi = 300)
cat("   Saved: 05_upset_stacked.pdf\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 10. COMPOSITE FIGURE: BARS + EULER INSETS
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n>> 10 -- Composite Figure\n")

wrap_euler_gg <- function(item_lists, title_text) {
  if (length(item_lists) < 2) {
    return(ggplot() + annotate("text", x = .5, y = .5, label = "< 2 sets", size = 3) +
             theme_void() + ggtitle(title_text))
  }
  cols <- CONTRAST_META$color[match(names(item_lists), CONTRAST_META$contrast)]
  fit  <- euler(item_lists)
  p    <- plot(fit,
               quantities = list(type = "counts", cex = 0.8),
               fills = list(fill = cols, alpha = 0.45),
               labels = list(cex = 0.55),
               main = list(label = title_text, cex = 0.8))
  wrap_elements(p)
}

e_up_prot <- wrap_euler_gg(up_gene_lists,  "Up Proteins")
e_dn_prot <- wrap_euler_gg(dn_gene_lists,  "Down Proteins")
e_up_path <- wrap_euler_gg(up_path_lists,  "Up Pathways")
e_dn_path <- wrap_euler_gg(dn_path_lists,  "Down Pathways")

composite <- ((p_prot_bar | e_up_prot | e_dn_prot) /
  (p_path_bar | e_up_path | e_dn_path)) +
  plot_layout(widths = c(3, 1.5, 1.5)) +
  plot_annotation(
    title = "Differential Protein & Pathway Summary",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  )

ggsave(file.path(REPORT_DIR, "07_composite_summary.pdf"),
       composite, width = 18, height = 10)
ggsave(file.path(REPORT_DIR, "07_composite_summary.png"),
       composite, width = 18, height = 10, dpi = 300)

cat("   Saved: 07_composite_summary.pdf\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 11. SAVE & SESSION INFO
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n>> 11 -- Export & Session Info\n")

save(results_list, all_fgsea, go_bp_results, fit2, meta, ann, mat,
     contrast_matrix, design,
     file = file.path(DATA_DIR, "DEP_workspace.RData"))

cat(sprintf("   Proteins tested: %d\n", nrow(mat)))
for (cname in contrast_names) {
  tt <- results_list[[cname]]
  cat(sprintf("   %s: %d sig (Pi < 0.05) -- %d up, %d down\n",
              cname, sum(tt$sig_pi),
              sum(tt$direction == "Up"), sum(tt$direction == "Down")))
}
cat(sprintf("   Hallmark pathways significant: %d (across all contrasts)\n",
            sum(all_fgsea$sig)))
cat(sprintf("   GO:BP terms significant: %d\n", nrow(go_bp_results)))

cat("\n>> Session Info\n")
sessionInfo()

cat("\n=== YvO DEP Analysis Complete ===\n")
