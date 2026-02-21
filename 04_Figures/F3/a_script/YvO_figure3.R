# ============================================================================
# Figure 3: Rejuvenation — Does Training Reverse Aging?
# ============================================================================
# Study: YvO — DIA-MS Proteomics of Skeletal Muscle
# Design: 2x2 mixed (Age x Time) with repeated measures
#
# Central question: Does resistance training in Old adults reverse
# age-related proteomic changes?
#
# Panels:
#   A — Side-by-side volcanos (Aging | Training Old)
#   B — Reversal scatter (logFC Aging vs logFC Training Old)
#   C — RRHO2 threshold-free reversal map (Cahill et al. 2018)
#   D — mitch 2D pathway enrichment (Kaspi & Ziemann 2020)
#   E — Reversal classification (diverging lollipop)
#   F — Pathway enrichment of reversal categories (Yu et al. 2012)
#
# References:
#   Melov et al. 2007 PLOS ONE (rejuvenation of skeletal muscle)
#   Lee-Odegard et al. 2025 npj Aging
#   Cahill et al. 2018, Scientific Reports 8:9588 (RRHO2)
#   Kaspi & Ziemann 2020, BMC Genomics 21:447 (mitch)
#   Yu et al. 2012, OMICS 16:284-287 (clusterProfiler)
#   Xiao et al. 2014 (pi-score methodology)
# ============================================================================

# === 1. PACKAGES ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(grid)
  library(RRHO2)
  library(mitch)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(fgsea)
  library(ggExtra)
  library(png)
})

# === 2. SEED ================================================================

set.seed(42)

# === 3. PATH RESOLUTION =====================================================

.script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)),
                         error = function(e) {
                           args <- commandArgs(trailingOnly = FALSE)
                           f <- grep("^--file=", args, value = TRUE)
                           if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1])))
                           else normalizePath(".")
                         })
BASE_DIR <- normalizePath(file.path(.script_dir, "..", "..", ".."))
FIG_DIR  <- normalizePath(file.path(.script_dir, ".."))
RPT_DIR  <- file.path(FIG_DIR, "b_reports")
DAT_DIR  <- file.path(FIG_DIR, "c_data")
SUP_DIR  <- file.path(RPT_DIR, "supplementary")

# === 4. DIRECTORY CREATION ==================================================

dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUP_DIR, recursive = TRUE, showWarnings = FALSE)

# === 5. CANONICAL CONSTANTS =================================================

CONTRAST_COLORS <- c(Aging = "#4CAF50", Training_Young = "#E05A4E",
                     Training_Old = "#5DA5DA", Interaction = "#9B7FBF")
AGE_COLORS <- c(Young = "#4393C3", Old = "#D6604D")
DIR_COLORS <- c(Up = "#D6604D", Down = "#4393C3")
DB_COLORS  <- c(Hallmark = "#AA336A", "GO:BP" = "#00796B",
                "GO:CC" = "#26A69A", "GO:MF" = "#CD5C5C")
KEY_TEXT  <- 2.2
KEY_TITLE <- 2.3

THEME_PUB <- theme_bw(base_size = 8) +
  theme(plot.title       = element_text(face = "bold", size = 9),
        plot.subtitle    = element_text(size = 6.5, color = "grey30", face = "italic"),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 6.5),
        legend.key.size  = unit(3, "mm"))

# === 6. HELPER FUNCTIONS ====================================================

clean_pathway_name <- function(name) {
  name %>%
    str_remove("^HALLMARK_") %>%
    str_remove("^GOBP_") %>%
    str_remove("^GOCC_") %>%
    str_remove("^GOMF_") %>%
    str_replace_all("_", " ") %>%
    str_to_title() %>%
    str_replace("Mtorc1", "mTORC1") %>%
    str_replace("Tnfa", "TNFa") %>%
    str_replace("Myc ", "MYC ") %>%
    str_replace("E2f ", "E2F ") %>%
    str_replace("P53", "p53") %>%
    str_replace("Kras", "KRAS") %>%
    str_replace("Dna ", "DNA ") %>%
    str_replace("Uv ", "UV ") %>%
    str_replace("Il2", "IL-2") %>%
    str_replace("Il6", "IL-6") %>%
    str_replace("Ifn", "IFN") %>%
    str_replace("Nfkb", "NF-kB") %>%
    str_replace("Tgf", "TGF") %>%
    str_replace("Pi3k", "PI3K") %>%
    str_replace("Akt", "AKT") %>%
    str_replace("Mtor", "mTOR") %>%
    str_replace("Oxidative Phosphorylation", "OXPHOS") %>%
    str_trunc(45, ellipsis = "...")
}

sig_stars <- function(padj) {
  case_when(padj < 0.001 ~ "***",
            padj < 0.01  ~ "**",
            padj < 0.05  ~ "*",
            TRUE         ~ "")
}

# === 7. DATA LOADING ========================================================

message("Loading data...")
dep_df <- read_csv(file.path(BASE_DIR, "03_DEP", "c_data", "combined_results.csv"),
                   show_col_types = FALSE)
fgsea_all <- read_csv(file.path(BASE_DIR, "04_Figures", "F2", "c_data", "fgsea_tstat_all_v2.csv"),
                      show_col_types = FALSE)

# Validation
stopifnot(nrow(dep_df) > 2000)
stopifnot(all(c("gene", "logFC_Aging", "t_Aging",
                "logFC_Training_Old", "t_Training_Old",
                "sig_pi_Aging", "sig_pi_Training_Old") %in% names(dep_df)))
message(sprintf("Loaded %d proteins, %d fGSEA results", nrow(dep_df), nrow(fgsea_all)))

# === 8. PANEL A — Side-by-Side Volcanos (Aging | Training Old) ==============

message("Building Panel A: side-by-side volcanos (Aging | Training Old)...")

make_volcano <- function(ctr, fill_color) {
  # --- 1. Extract & rename columns for this contrast ---
  col_logFC  <- paste0("logFC_", ctr)
  col_pval   <- paste0("P.Value_", ctr)
  col_pi     <- paste0("pi_score_", ctr)
  col_sig    <- paste0("sig_pi_", ctr)

  vdf <- dep_df %>%
    dplyr::select(gene,
           logFC  = all_of(col_logFC),
           P.Value = all_of(col_pval),
           pi_score = all_of(col_pi),
           sig_pi   = all_of(col_sig)) %>%
    filter(!is.na(logFC), !is.na(P.Value)) %>%
    mutate(
      neg_log10p = -log10(P.Value),
      direction  = case_when(
        sig_pi == 1 & logFC > 0 ~ "Up",
        sig_pi == 1 & logFC < 0 ~ "Down",
        TRUE                     ~ "NS"
      )
    )

  # --- 2. Summary counts for corner annotations ---
  n_up   <- sum(vdf$direction == "Up",   na.rm = TRUE)
  n_down <- sum(vdf$direction == "Down", na.rm = TRUE)
  med_lfc_up   <- median(abs(vdf$logFC[vdf$direction == "Up"]),   na.rm = TRUE)
  med_lfc_down <- median(abs(vdf$logFC[vdf$direction == "Down"]), na.rm = TRUE)

  # --- 3. Top 6 DEPs by |logFC| among significant ---
  top_genes <- vdf %>%
    filter(sig_pi == 1) %>%
    arrange(desc(abs(logFC))) %>%
    slice_head(n = 6)

  # --- 4. Pathway inset: Hallmark, padj < 0.05 ---
  pw_df <- fgsea_all %>%
    filter(contrast == ctr, database == "Hallmark", padj < 0.05)

  pw_up <- pw_df %>%
    filter(NES > 0) %>%
    arrange(desc(NES)) %>%
    slice_head(n = 3) %>%
    mutate(label = clean_pathway_name(pathway))

  pw_down <- pw_df %>%
    filter(NES < 0) %>%
    arrange(NES) %>%
    slice_head(n = 3) %>%
    mutate(label = clean_pathway_name(pathway))

  # Axis ranges (needed for positioning pathway labels)
  y_max_est <- max(vdf$neg_log10p, na.rm = TRUE)
  x_max_est <- max(abs(vdf$logFC), na.rm = TRUE)

  # --- 5. Build ggplot ---
  # Strip title: display "Aging" or "Training (Old)"
  strip_title <- if (ctr == "Aging") "Aging" else paste0("Training (", str_extract(ctr, "Young|Old"), ")")

  p <- ggplot(vdf, aes(x = logFC, y = neg_log10p)) +
    # Reference lines
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey40", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey40", linewidth = 0.3) +
    # Points
    geom_point(aes(color = direction), size = 0.5, alpha = 0.4) +
    scale_color_manual(values = c(Up = "#D6604D", Down = "#4393C3", NS = "grey70")) +
    # Gene labels
    geom_text_repel(
      data = top_genes,
      aes(label = gene),
      size       = KEY_TEXT,
      max.overlaps  = 15,
      segment.size  = 0.2,
      fontface      = "italic",
      min.segment.length = 0
    ) +
    # Corner stat annotations — top-right for Up
    annotate("text",
             x     = x_max_est * 0.95,
             y     = y_max_est * 0.97,
             label = if (n_up > 0)
               paste0("n(Up) = ", n_up, "  med|logFC| = ",
                      sprintf("%.2f", med_lfc_up))
             else "n(Up) = 0",
             hjust = 1, vjust = 1, size = KEY_TITLE, color = "#D6604D") +
    # Corner stat annotations — top-left for Down
    annotate("text",
             x     = -x_max_est * 0.95,
             y     = y_max_est * 0.97,
             label = if (n_down > 0)
               paste0("n(Down) = ", n_down, "  med|logFC| = ",
                      sprintf("%.2f", med_lfc_down))
             else "n(Down) = 0",
             hjust = 0, vjust = 1, size = KEY_TITLE, color = "#4393C3") +
    # Title & axes
    labs(
      title = strip_title,
      x     = expression(log[2]~fold~change),
      y     = expression(-log[10]~italic(P))
    ) +
    THEME_PUB +
    theme(legend.position = "none")

  # --- 6. Add pathway inset labels ---
  if (nrow(pw_up) > 0) {
    pw_up_label <- paste(pw_up$label, collapse = "\n")
    p <- p + annotate("label",
                       x     = x_max_est * 0.60,
                       y     = y_max_est * 0.15,
                       label = pw_up_label,
                       hjust = 1, vjust = 0,
                       size  = 1.8,
                       fill  = alpha("white", 0.85),
                       label.padding = unit(1.5, "pt"),
                       color = "#D6604D")
  }

  if (nrow(pw_down) > 0) {
    pw_down_label <- paste(pw_down$label, collapse = "\n")
    p <- p + annotate("label",
                       x     = -x_max_est * 0.60,
                       y     = y_max_est * 0.15,
                       label = pw_down_label,
                       hjust = 0, vjust = 0,
                       size  = 1.8,
                       fill  = alpha("white", 0.85),
                       label.padding = unit(1.5, "pt"),
                       color = "#4393C3")
  }

  return(p)
}

# --- Panel A assembly: shared axis limits across Aging and Training_Old ---
volc_xlim <- max(abs(c(dep_df$logFC_Aging, dep_df$logFC_Training_Old)),
                 na.rm = TRUE) * 1.05
volc_ylim <- max(-log10(c(dep_df$P.Value_Aging, dep_df$P.Value_Training_Old)),
                 na.rm = TRUE)
volc_ylim <- min(volc_ylim, 15)   # cap at 15

pA_left  <- make_volcano("Aging", CONTRAST_COLORS["Aging"]) +
  coord_cartesian(xlim = c(-volc_xlim, volc_xlim), ylim = c(0, volc_ylim))
pA_right <- make_volcano("Training_Old", CONTRAST_COLORS["Training_Old"]) +
  coord_cartesian(xlim = c(-volc_xlim, volc_xlim), ylim = c(0, volc_ylim)) +
  theme(axis.title.y = element_blank())

pA <- (pA_left | pA_right) + plot_layout(widths = c(1, 1))

# --- Panel A test save ---
ggsave(file.path(RPT_DIR, "test_panelA.pdf"), pA,
       width = 250, height = 120, units = "mm")
message("Panel A test saved")

# === 9. PANEL B — Reversal Scatter (logFC Aging vs logFC Training Old) ======

message("Building Panel B: reversal scatter...")

# --- 1. Prepare data ---
scatter_df <- dep_df %>%
  transmute(gene,
            logFC_A  = logFC_Aging,
            logFC_TO = logFC_Training_Old,
            sig_A    = sig_pi_Aging,
            sig_TO   = sig_pi_Training_Old) %>%
  filter(!is.na(logFC_A), !is.na(logFC_TO)) %>%
  mutate(
    reversed = (logFC_A > 0 & logFC_TO < 0) | (logFC_A < 0 & logFC_TO > 0),
    quadrant = case_when(
      logFC_A > 0 & logFC_TO < 0 ~ "Reversed",
      logFC_A < 0 & logFC_TO > 0 ~ "Reversed",
      TRUE ~ "Exacerbated"
    )
  )

# --- 2. Compute stats ---
cor_r   <- cor(scatter_df$logFC_A, scatter_df$logFC_TO, use = "complete.obs")
cor_rho <- cor(scatter_df$logFC_A, scatter_df$logFC_TO, use = "complete.obs",
               method = "spearman")

# --- 3. Quadrant counts ---
q_detail <- scatter_df %>%
  mutate(q = case_when(
    logFC_A > 0 & logFC_TO < 0 ~ "Q2",
    logFC_A < 0 & logFC_TO > 0 ~ "Q4",
    logFC_A > 0 & logFC_TO > 0 ~ "Q1",
    TRUE ~ "Q3"
  )) %>%
  count(q) %>%
  deframe()

# --- 4. Axis range for quadrant shading ---
axis_lim <- max(abs(c(scatter_df$logFC_A, scatter_df$logFC_TO)), na.rm = TRUE) * 1.15

# --- 5. Build scatter plot ---
# Quadrant shading: Q2 (x>0, y<0) and Q4 (x<0, y>0) = REVERSED = teal
#                   Q1 (x>0, y>0) and Q3 (x<0, y<0) = EXACERBATED = orange
pB_base <- ggplot(scatter_df, aes(x = logFC_A, y = logFC_TO)) +
  # Reversed quadrants (Q2, Q4) in teal
  annotate("rect", xmin = 0, xmax =  axis_lim, ymin = -axis_lim, ymax = 0,
           fill = "#00897B", alpha = 0.06) +
  annotate("rect", xmin = -axis_lim, xmax = 0, ymin = 0, ymax =  axis_lim,
           fill = "#00897B", alpha = 0.06) +
  # Exacerbated quadrants (Q1, Q3) in orange
  annotate("rect", xmin = 0, xmax =  axis_lim, ymin = 0, ymax =  axis_lim,
           fill = "#FF8F00", alpha = 0.06) +
  annotate("rect", xmin = -axis_lim, xmax = 0, ymin = -axis_lim, ymax = 0,
           fill = "#FF8F00", alpha = 0.06) +
  # Reference lines
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  # Anti-diagonal reference line (perfect reversal: y = -x)
  geom_abline(slope = -1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  # NS points first
  geom_point(data = filter(scatter_df, sig_A == 0 & sig_TO == 0),
             size = 0.5, alpha = 0.3, color = "grey60") +
  # Significant proteins overlaid as diamonds
  geom_point(data = filter(scatter_df, sig_A == 1 | sig_TO == 1),
             aes(color = reversed),
             shape = 18, size = 2.0, alpha = 0.85) +
  scale_color_manual(values = c(`TRUE` = "#00897B", `FALSE` = "#FF8F00"),
                     guide = "none") +
  # Stats annotation — upper-left
  annotate("text",
           x = -axis_lim * 0.95, y = axis_lim * 0.95,
           label = sprintf("r = %.2f\nrho = %.2f", cor_r, cor_rho),
           hjust = 0, vjust = 1,
           size = KEY_TITLE, fontface = "bold") +
  # Quadrant count annotations
  annotate("text", x =  axis_lim, y =  axis_lim,
           label = paste("Exac. n =", q_detail["Q1"]),
           hjust = 1.1, vjust = 1.5, size = 2.0, color = "#FF8F00") +
  annotate("text", x =  axis_lim, y = -axis_lim,
           label = paste("Rev. n =", q_detail["Q2"]),
           hjust = 1.1, vjust = -0.5, size = 2.0, color = "#00897B") +
  annotate("text", x = -axis_lim, y = -axis_lim,
           label = paste("Exac. n =", q_detail["Q3"]),
           hjust = -0.1, vjust = -0.5, size = 2.0, color = "#FF8F00") +
  annotate("text", x = -axis_lim, y =  axis_lim,
           label = paste("Rev. n =", q_detail["Q4"]),
           hjust = -0.1, vjust = 1.5, size = 2.0, color = "#00897B") +
  # Axis labels and title
  labs(
    title = "Protein Reversal Map",
    subtitle = sprintf("Dashed = perfect reversal | r = %.2f", cor_r),
    x = expression(log[2]*FC ~ "(Aging)"),
    y = expression(log[2]*FC ~ "(Training Old)")
  ) +
  coord_cartesian(xlim = c(-axis_lim, axis_lim),
                  ylim = c(-axis_lim, axis_lim)) +
  THEME_PUB +
  theme(legend.position = "none")

# --- 6. Label top reversed + exacerbated proteins ---
top_reversed <- scatter_df %>%
  filter(sig_A == 1 | sig_TO == 1, reversed) %>%
  mutate(dev = abs(logFC_A) + abs(logFC_TO)) %>%
  arrange(desc(dev)) %>%
  slice_head(n = 8)

top_exacerbated <- scatter_df %>%
  filter(sig_A == 1 | sig_TO == 1, !reversed) %>%
  mutate(dev = abs(logFC_A) + abs(logFC_TO)) %>%
  arrange(desc(dev)) %>%
  slice_head(n = 4)

top_labels_B <- bind_rows(top_reversed, top_exacerbated)

pB_base <- pB_base +
  geom_text_repel(
    data = top_labels_B,
    aes(label = gene),
    size             = KEY_TEXT,
    max.overlaps     = 15,
    segment.size     = 0.2,
    fontface         = "italic",
    min.segment.length = 0
  )

# --- 7. Export reversal scatter data ---
write_csv(scatter_df, file.path(DAT_DIR, "fig3_reversal_scatter.csv"))
message(sprintf("Panel B: r = %.3f, rho = %.3f, %d reversed, %d exacerbated",
                cor_r, cor_rho,
                sum(scatter_df$reversed), sum(!scatter_df$reversed)))

# --- Test save ---
ggsave(file.path(RPT_DIR, "test_panelB.pdf"), pB_base,
       width = 150, height = 140, units = "mm")
message("Panel B test saved")

# === 10. PANEL C — RRHO2 (Aging vs Training Old) ===========================

message("Building Panel C: RRHO2 reversal map...")

# --- 1. Prepare ranked lists ---
rrho_list1 <- dep_df %>%
  dplyr::select(gene, t = t_Aging) %>%
  dplyr::filter(!is.na(t)) %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  as.data.frame()

rrho_list2 <- dep_df %>%
  dplyr::select(gene, t = t_Training_Old) %>%
  dplyr::filter(!is.na(t)) %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  as.data.frame()

# Intersect to shared genes
shared_genes <- intersect(rrho_list1$gene, rrho_list2$gene)
rrho_list1 <- rrho_list1 %>% dplyr::filter(gene %in% shared_genes)
rrho_list2 <- rrho_list2 %>% dplyr::filter(gene %in% shared_genes)

# --- 2. Run RRHO2 ---
rrho_obj <- RRHO2_initialize(
  list1 = rrho_list1,
  list2 = rrho_list2,
  labels = c("Aging", "Training (Old)"),
  log10.ind = TRUE,
  method = "hyper",
  boundary = 0.04
)

# --- 3. Export RRHO2 matrix ---
write.csv(rrho_obj$hypermat, file.path(DAT_DIR, "fig3_rrho2_matrix.csv"))
message(sprintf("RRHO2 matrix: %d x %d", nrow(rrho_obj$hypermat), ncol(rrho_obj$hypermat)))

# --- 4. Test save ---
pdf(file.path(RPT_DIR, "test_panelC.pdf"), width = 7, height = 5)
par(mar = c(2, 2, 2, 1))
RRHO2_heatmap(rrho_obj)
dev.off()
message("Panel C test saved")

# --- 5. Capture RRHO2 heatmap for patchwork ---
pC <- wrap_elements(full = ~ {
  RRHO2_heatmap(rrho_obj)
})

# === 11. PANEL D — mitch 2D Pathway Enrichment =============================

message("Building Panel D: mitch 2D pathway enrichment...")

# --- 1. Build t-stat matrix ---
tstat_mat <- dep_df %>%
  dplyr::select(gene, Aging = t_Aging, Training_Old = t_Training_Old) %>%
  dplyr::filter(!is.na(Aging), !is.na(Training_Old)) %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# --- 2. Fetch gene sets from msigdbr (as named list for mitch) ---
hallmark_gs <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol)
gobp_gs <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP") %>%
  dplyr::select(gs_name, gene_symbol)

gs_df <- bind_rows(hallmark_gs, gobp_gs)
genesets <- split(gs_df$gene_symbol, gs_df$gs_name)

# --- 3. Run mitch ---
mitch_res <- mitch_calc(x = tstat_mat, genesets = genesets,
                        priority = "effect", cores = 1, resrows = Inf)

# --- 4. Process results for plotting ---
mitch_df <- mitch_res$enrichment_result %>%
  mutate(
    padj_A   = p.adjust(p.Aging, method = "BH"),
    padj_TO  = p.adjust(p.Training_Old, method = "BH"),
    sig_A    = padj_A < 0.05,
    sig_TO   = padj_TO < 0.05,
    sig_cat  = case_when(
      sig_A & sig_TO ~ "Both",
      sig_A | sig_TO ~ "One",
      TRUE ~ "Neither"
    ),
    pathway_clean = clean_pathway_name(set)
  )

# --- 5. Identify pathways to label ---
label_keywords <- c("ribosom", "oxphos", "oxidative phosph", "mtorc1",
                     "extracellular matrix", "myogenes", "glycoly",
                     "proteasome", "translation", "unfolded protein",
                     "mitotic spindle", "muscle", "respiratory chain")
label_df <- mitch_df %>%
  filter(sig_cat != "Neither") %>%
  filter(str_detect(tolower(set), paste(label_keywords, collapse = "|"))) %>%
  bind_rows(
    mitch_df %>% filter(sig_cat == "Both") %>%
      slice_max(abs(s.Aging) + abs(s.Training_Old), n = 6)
  ) %>%
  distinct(set, .keep_all = TRUE) %>%
  slice_head(n = 12)

# --- 6. Build scatter plot ---
pw_cor <- cor(mitch_df$s.Aging, mitch_df$s.Training_Old)

pD <- ggplot(mitch_df, aes(x = s.Aging, y = s.Training_Old)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  # Identity line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40", linewidth = 0.3) +
  # Anti-diagonal reference line (perfect reversal)
  geom_abline(slope = -1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  # Points by significance category
  geom_point(data = filter(mitch_df, sig_cat == "Neither"),
             aes(size = setSize), color = "grey75", alpha = 0.3) +
  geom_point(data = filter(mitch_df, sig_cat == "One"),
             aes(size = setSize), color = "#FF8F00", alpha = 0.5) +
  geom_point(data = filter(mitch_df, sig_cat == "Both"),
             aes(size = setSize), color = "#2E7D32", alpha = 0.8) +
  scale_size_continuous(range = c(0.3, 2.5), name = "Gene set size",
                        guide = guide_legend(override.aes = list(alpha = 0.8))) +
  # Labels for key pathways
  geom_text_repel(data = label_df,
                  aes(label = pathway_clean),
                  size = KEY_TEXT, max.overlaps = 20,
                  segment.size = 0.2, fontface = "italic",
                  min.segment.length = 0) +
  # Correlation annotation
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = sprintf("r = %.2f", pw_cor),
           size = KEY_TITLE, fontface = "bold") +
  labs(x = "Enrichment (Aging)",
       y = "Enrichment (Training Old)",
       title = "Pathway Reversal (mitch 2D)") +
  THEME_PUB +
  theme(legend.position = "bottom",
        legend.key.size = unit(2, "mm"),
        legend.text = element_text(size = 5))

# --- 7. Export mitch results ---
write_csv(mitch_df, file.path(DAT_DIR, "fig3_mitch_2d_results.csv"))
message(sprintf("mitch: %d pathways, %d sig both, %d sig one, r = %.3f",
                nrow(mitch_df),
                sum(mitch_df$sig_cat == "Both"),
                sum(mitch_df$sig_cat == "One"),
                pw_cor))

# --- 8. Test save ---
ggsave(file.path(RPT_DIR, "test_panelD.pdf"), pD,
       width = 160, height = 150, units = "mm")
message("Panel D test saved")

# === 12. PANEL E — Reversal Classification (Diverging Lollipop) =============

message("Building Panel E: reversal classification (diverging lollipop)...")

# --- 1. Classify Aging DEPs by Training_Old behavior ---
rev_class <- dep_df %>%
  filter(sig_pi_Aging == 1) %>%
  dplyr::select(gene, logFC_A = logFC_Aging, logFC_TO = logFC_Training_Old,
         sig_TO = sig_pi_Training_Old) %>%
  filter(!is.na(logFC_A), !is.na(logFC_TO)) %>%
  mutate(
    opposite_dir = sign(logFC_A) != sign(logFC_TO),
    category = case_when(
      opposite_dir & sig_TO == 1 ~ "Fully reversed",
      opposite_dir & sig_TO == 0 ~ "Partially reversed",
      !opposite_dir & sig_TO == 1 ~ "Exacerbated",
      !opposite_dir & sig_TO == 0 ~ "Non-reversed"
    ),
    reversal_index = -logFC_A * logFC_TO / pmax(abs(logFC_A), 1e-6)
  ) %>%
  arrange(desc(reversal_index))

message(sprintf("Panel E: %d Aging DEPs classified", nrow(rev_class)))

# --- 2. If >40 proteins, keep top 35 by |reversal_index| for readability ---
n_total_rev <- nrow(rev_class)
if (n_total_rev > 40) {
  rev_plot_df <- rev_class %>% slice_head(n = 35)
  trunc_note <- sprintf("Top 35 of %d Aging DEPs by reversal index", n_total_rev)
  message(sprintf("  Truncating to top 35 of %d for plot", n_total_rev))
} else {
  rev_plot_df <- rev_class
  trunc_note <- NULL
}

# --- 3. Build diverging lollipop plot ---
rev_colors <- c("Fully reversed" = "#00897B", "Partially reversed" = "#80CBC4",
                "Non-reversed" = "grey60", "Exacerbated" = "#FF8F00")
cat_counts <- rev_class %>% count(category) %>% deframe()
cat_pcts   <- round(100 * cat_counts / sum(cat_counts), 1)

pE <- ggplot(rev_plot_df, aes(y = reorder(gene, reversal_index))) +
  geom_segment(aes(x = logFC_A, xend = logFC_TO, yend = reorder(gene, reversal_index),
                   color = category), linewidth = 0.4) +
  geom_point(aes(x = logFC_A), color = "#4CAF50", size = 1.0) +
  geom_point(aes(x = logFC_TO), color = "#5DA5DA", size = 1.0) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.3) +
  scale_color_manual(values = rev_colors, name = "Category") +
  labs(x = expression(log[2]*FC),
       y = NULL,
       title = "Reversal Classification",
       subtitle = paste(
         paste0(names(cat_counts), ": ", cat_counts, " (", cat_pcts, "%)"),
         collapse = "  |  ")) +
  THEME_PUB +
  theme(axis.text.y = element_text(size = 1.5),
        legend.position = "bottom",
        legend.key.size = unit(2, "mm"),
        legend.text = element_text(size = 5))

# Add truncation caption if needed
if (!is.null(trunc_note)) {
  pE <- pE + labs(caption = trunc_note) +
    theme(plot.caption = element_text(size = 5, color = "grey50", hjust = 0.5))
}

# --- 4. Export data and test save ---
write_csv(rev_class, file.path(DAT_DIR, "fig3_reversal_classification.csv"))
ggsave(file.path(RPT_DIR, "test_panelE.pdf"), pE,
       width = 120, height = 180, units = "mm")
message("Panel E test saved")

# === 13. PANEL F — Pathway Enrichment by Reversal Category ==================

message("Building Panel F: compareCluster enrichment by reversal category...")

# --- 1. Prepare gene lists by reversal category (categories with >= 5 genes) ---
gene_list_rev <- rev_class %>%
  group_by(category) %>%
  filter(n() >= 5) %>%
  ungroup() %>%
  split(.$category) %>%
  lapply(function(x) x$gene)

# Background universe: all measured proteins
all_genes <- unique(dep_df$gene)

message(sprintf("  Reversal categories with >= 5 genes: %s",
                paste(names(gene_list_rev), vapply(gene_list_rev, length, integer(1)),
                      sep = "=", collapse = ", ")))

if (length(gene_list_rev) < 1) {
  message("  WARNING: No reversal categories with >= 5 genes — skipping Panel F")
  pF <- ggplot() + annotate("text", x = 0.5, y = 0.5,
                             label = "No enrichment\n(empty gene sets)", size = 3) +
    theme_void()
} else {
  cc_res <- compareCluster(
    geneClusters = gene_list_rev,
    fun           = "enrichGO",
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    universe      = all_genes
  )

  cc_df <- as.data.frame(cc_res)
  message(sprintf("  compareCluster returned %d enriched terms", nrow(cc_df)))

  if (nrow(cc_df) == 0) {
    message("  WARNING: No significant GO:BP terms — Panel F will be a placeholder")
    pF <- ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "No significant GO:BP terms\n(p.adjust < 0.05)",
                               size = 3) +
      theme_void()
  } else {
    # --- 2. Build dot plot ---
    pF <- tryCatch({
      dotplot(cc_res, showCategory = 8, font.size = 6) +
        labs(title = "GO:BP Enrichment",
             subtitle = "By Reversal Category") +
        THEME_PUB +
        theme(axis.text.y = element_text(size = 5))
    }, error = function(e) {
      message("  dotplot() failed: ", e$message)
      message("  Building manual dot plot from compareCluster data frame...")

      # Select top terms per cluster
      plot_df <- cc_df %>%
        group_by(Cluster) %>%
        slice_min(order_by = p.adjust, n = 8, with_ties = FALSE) %>%
        ungroup() %>%
        mutate(
          GeneRatio_num = sapply(GeneRatio, function(x) {
            parts <- as.numeric(strsplit(x, "/")[[1]])
            parts[1] / parts[2]
          }),
          Description = str_wrap(Description, width = 45)
        )

      ggplot(plot_df, aes(x = Cluster, y = reorder(Description, GeneRatio_num))) +
        geom_point(aes(size = GeneRatio_num, color = p.adjust)) +
        scale_color_gradient(low = "#D6604D", high = "#4393C3",
                             name = "Adj. p-value") +
        scale_size_continuous(name = "Gene Ratio", range = c(1, 4)) +
        labs(title = "GO:BP Enrichment",
             subtitle = "By Reversal Category",
             x = NULL, y = NULL) +
        THEME_PUB +
        theme(axis.text.y = element_text(size = 5))
    })
  }

  # --- 3. Export data ---
  write_csv(cc_df, file.path(DAT_DIR, "fig3_reversal_enrichment.csv"))
}

# --- 4. Test save ---
ggsave(file.path(RPT_DIR, "test_panelF.pdf"), pF,
       width = 180, height = 160, units = "mm")
message("Panel F test saved")

# === 14. FINAL ASSEMBLY — 6-Panel Composite Figure =========================

message("Assembling Figure 3...")

# --- Strategy ---
# wrap_elements() for the volcano patchwork so it gets one tag.
# Render RRHO2 to a temp PNG and embed as a raster ggplot.
# Use pB_base (no marginals) for the composite.

# Panel A: wrap the 2-volcano patchwork as a single unit for tagging
pA_wrapped <- wrap_elements(full = pA)

# Panel C: render RRHO2 to temp PNG and read back as ggplot raster
tmp_rrho <- tempfile(fileext = ".png")
png(tmp_rrho, width = 1400, height = 1200, res = 300)
par(mar = c(2, 2, 2, 1))
RRHO2_heatmap(rrho_obj)
dev.off()
rrho_img <- png::readPNG(tmp_rrho)
pC_gg <- ggplot() +
  annotation_raster(rrho_img, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
  labs(title = "RRHO2 Reversal Map") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5))

# Compose: A=volcanos, C=RRHO2, B=scatter, D=mitch, E=lollipop, F=enrichment
fig3 <- (pA_wrapped | pC_gg) /
         (pB_base   | pD) /
         (pE        | pF) +
  plot_layout(
    widths  = c(0.6, 0.4),
    heights = c(0.35, 0.35, 0.30)
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 12)
    )
  )

# Save outputs
ggsave(file.path(RPT_DIR, "Figure_3.pdf"), fig3,
       width = 380, height = 500, units = "mm", limitsize = FALSE)
ggsave(file.path(RPT_DIR, "Figure_3.png"), fig3,
       width = 380, height = 500, units = "mm", dpi = 300, limitsize = FALSE)

message("Figure 3 saved to: ", RPT_DIR)
message("Figure 3 pipeline complete.")
