################################################################################
#   Figure 2 — Panel Export Pipeline
#   Exports individual panel PDFs + tool-ready CSVs for external assembly
#
#   Panels:
#     A — Volcano: Training Effect (Young)
#     B — Volcano: Training Effect (Old)
#     C — Concordance scatter (logFC x logFC)
#     D — RRHO2 threshold-free concordance map
#     E — fGSEA NES scatter (Hallmark + rrvgo-reduced GO:BP)
#     F — Interaction DEPs: Multi-Contrast Response & Pathway Enrichment
#         (merged dumbbell | sankey | bar panel)
#
#   Outputs per panel:
#     b_reports/panel_A_volcano.pdf         + c_data/panel_A/volcano_young.csv
#     b_reports/panel_B_volcano.pdf         + c_data/panel_B/volcano_old.csv
#     b_reports/panel_C_concordance.pdf     + c_data/panel_C/concordance.csv
#     b_reports/panel_D_rrho2.pdf           + c_data/panel_D/rrho2_*.csv
#     b_reports/panel_E_nes_bubble.pdf      + c_data/panel_E/nes_scatter.csv
#     b_reports/panel_F_interaction.pdf     + c_data/panel_F/*.csv
################################################################################

# === 1. PACKAGES ==============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(grid)
  library(RRHO2)
  library(msigdbr)
  library(fgsea)
  library(rrvgo)
  library(GOSemSim)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  # ggalluvial removed — manual Sankey construction below
})

# === 2. SETUP =================================================================

set.seed(42)
setwd("/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_YvO_2025")

RPT_DIR <- "04_Figures/F2/b_reports"
DAT_DIR <- "04_Figures/F2/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)
# Per-panel subdirectories for organized data output
for (pnl in c("panel_A", "panel_B", "panel_C", "panel_D", "panel_E", "panel_F", "shared"))
  dir.create(file.path(DAT_DIR, pnl), recursive = TRUE, showWarnings = FALSE)

# === 3. CANONICAL STYLE (from figure-style guide) ============================

CONTRAST_COLORS <- c(Aging = "#4CAF50", Training_Young = "#E05A4E",
                     Training_Old = "#5DA5DA", Interaction = "#9B7FBF")
AGE_COLORS <- c(Young = "#4393C3", Old = "#D6604D")
DIR_COLORS <- c(Up = "#D6604D", Down = "#4393C3")
KEY_TEXT  <- 2.2
KEY_TITLE <- 2.3

SIG_COLORS <- c(
  "Interaction"    = "#9B7FBF",
  "Sig Both"       = "#2E7D32",
  "Sig Young only" = "#E05A4E",
  "Sig Old only"   = "#5DA5DA",
  "NS"             = "grey55"
)

SIG_LABEL_FILL <- c(
  "Interaction"    = alpha("#9B7FBF", 0.75),
  "Sig Both"       = alpha("#2E7D32", 0.75),
  "Sig Young only" = alpha("#E05A4E", 0.75),
  "Sig Old only"   = alpha("#5DA5DA", 0.75),
  "NS"             = alpha("grey55",  0.75)
)
SIG_LABEL_TEXT <- c(
  "Interaction"    = "white",
  "Sig Both"       = "white",
  "Sig Young only" = "white",
  "Sig Old only"   = "white",
  "NS"             = "white"
)

THEME_PUB <- theme_bw(base_size = 8) +
  theme(
    plot.title       = element_text(face = "bold", size = 9),
    plot.subtitle    = element_text(size = 6.5, color = "grey30", face = "italic"),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 6.5),
    legend.key.size  = unit(3, "mm")
  )

# === 4. HELPERS ===============================================================

clean_pathway_name <- function(name) {
  name %>%
    str_remove("^HALLMARK_") %>%
    str_remove("^GOBP_") %>%
    str_remove("^GOCC_") %>%
    str_remove("^GOMF_") %>%
    str_remove("^REACTOME_") %>%
    str_remove("^KEGG_MEDICUS_") %>%
    str_remove("^KEGG_") %>%
    str_remove("^WP_") %>%
    str_replace_all("_", " ") %>%
    str_to_title() %>%
    str_replace("Mtorc1", "mTORC1") %>%
    str_replace("Myc ", "MYC ") %>%
    str_replace("E2f ", "E2F ") %>%
    str_replace("Dna ", "DNA ") %>%
    str_replace("Rna ", "RNA ") %>%
    str_replace("Tnfa ", "TNFa ") %>%
    str_replace("Uv ", "UV ") %>%
    str_replace("G2m ", "G2M ") %>%
    str_replace("Il6 ", "IL6 ") %>%
    str_replace("Il2 ", "IL2 ") %>%
    str_replace("Kras ", "KRAS ") %>%
    str_replace("P53 ", "p53 ") %>%
    str_replace("Tgf ", "TGF ") %>%
    str_replace("Nf Kb", "NF-kB") %>%
    str_replace("Atp ", "ATP ") %>%
    str_replace("Nadh ", "NADH ") %>%
    str_replace("Oxidative Phosphorylation", "OXPHOS") %>%
    str_replace("External Encapsulating Structure Or.*",
                "Extracellular Matrix Organization") %>%
    str_replace("Enzyme Linked Receptor Protein Signaling.*",
                "Receptor Protein Signaling") %>%
    str_trunc(45, ellipsis = "...")
}

sig_stars <- function(padj) {
  case_when(padj < 0.001 ~ "***",
            padj < 0.01  ~ "**",
            padj < 0.05  ~ "*",
            TRUE         ~ "")
}

classify_proteins <- function(pi_A, pi_B, pi_int, threshold = 0.05) {
  case_when(
    pi_int < threshold              ~ "Interaction",
    pi_A < threshold & pi_B < threshold ~ "Sig Both",
    pi_A < threshold                ~ "Sig Young only",
    pi_B < threshold                ~ "Sig Old only",
    TRUE                            ~ "NS"
  ) %>%
    factor(levels = c("Interaction", "Sig Both",
                       "Sig Young only", "Sig Old only", "NS"))
}

darken_color <- function(col, factor = 0.5) {
  rgb_vals <- col2rgb(col) / 255
  sapply(seq_along(col), function(i)
    rgb(rgb_vals[1, i] * factor, rgb_vals[2, i] * factor, rgb_vals[3, i] * factor))
}

# === 5. DATA LOADING ==========================================================

message("Loading data...")
dep_df <- read_csv("03_DEP/c_data/combined_results.csv", show_col_types = FALSE)
stopifnot(nrow(dep_df) > 2000)

# Imputation status (MAR/MNAR/Complete per protein)
imputation_df <- read_csv("02_Imputation/c_data/mar_mnar_classification.csv",
                          show_col_types = FALSE) %>%
  transmute(gene, imputed = classification != "Complete")
message(sprintf("  %d proteins with imputation status (%d imputed)",
                nrow(imputation_df), sum(imputation_df$imputed)))

fgsea_cache <- file.path(DAT_DIR, "shared", "fgsea_tstat_all_v2.csv")
if (!file.exists(fgsea_cache)) stop("fGSEA cache not found at ", fgsea_cache)
fgsea_all <- read_csv(fgsea_cache, show_col_types = FALSE)

message(sprintf("Loaded %d proteins, %d fGSEA results", nrow(dep_df), nrow(fgsea_all)))

# ==============================================================================
# PANELS A & B — Volcano Ring Composites (volcano + enrichment arcs)
# ==============================================================================

message("Panels A & B: volcano ring composites...")

# Source the ring utility (defines make_volcano_ring_pair)
source("04_Figures/F2/a_script/volcano_ring.R")

# Build the paired ring plots (saves PDFs + ring term CSVs)
make_volcano_ring_pair(
  de_df          = dep_df,
  go_df          = fgsea_all,
  contrast_young = "Training_Young",
  contrast_old   = "Training_Old",
  title_young    = "Training Effect (Young)",
  title_old      = "Training Effect (Old)",
  output_dir     = "04_Figures/F2",
  save_outputs   = TRUE
)

# Also export flat volcano CSVs for tool-ready data access
export_volcano_csv <- function(ctr, panel_dir, filename) {
  col_logFC <- paste0("logFC_", ctr)
  col_pval  <- paste0("P.Value_", ctr)
  col_pi    <- paste0("pi_score_", ctr)
  col_adjp  <- paste0("adj.P.Val_", ctr)

  dep_df %>%
    transmute(
      gene,
      log2_fold_change = round(.data[[col_logFC]], 4),
      neg_log10_pvalue = round(-log10(.data[[col_pval]]), 4),
      pi_score         = round(.data[[col_pi]], 6),
      adjusted_pvalue  = round(.data[[col_adjp]], 6),
      direction = case_when(
        .data[[col_pi]] < 0.05 & .data[[col_logFC]] > 0 ~ "Up",
        .data[[col_pi]] < 0.05 & .data[[col_logFC]] < 0 ~ "Down",
        TRUE ~ "NS"
      )
    ) %>%
    filter(!is.na(log2_fold_change), !is.na(neg_log10_pvalue)) %>%
    arrange(pi_score) %>%
    write_csv(file.path(DAT_DIR, panel_dir, filename))
}

export_volcano_csv("Training_Young", "panel_A", "volcano_young.csv")
export_volcano_csv("Training_Old",   "panel_B", "volcano_old.csv")
message("  Panels A & B saved")

# ==============================================================================
# PANEL C — Concordance Scatter
# ==============================================================================

message("Panel C: concordance scatter...")

scatter_df <- dep_df %>%
  transmute(
    gene,
    logFC_Training_Young = logFC_Training_Young,
    logFC_Training_Old   = logFC_Training_Old,
    pi_Young     = pi_score_Training_Young,
    pi_Old       = pi_score_Training_Old,
    pi_Interaction = pi_score_Interaction
  ) %>%
  filter(!is.na(logFC_Training_Young), !is.na(logFC_Training_Old)) %>%
  left_join(imputation_df, by = "gene") %>%
  mutate(
    imputed = replace_na(imputed, FALSE),
    significance = classify_proteins(pi_Young, pi_Old, pi_Interaction),
    quadrant = case_when(
      logFC_Training_Young > 0 & logFC_Training_Old > 0 ~ "Concordant Up",
      logFC_Training_Young < 0 & logFC_Training_Old < 0 ~ "Concordant Down",
      logFC_Training_Young > 0 & logFC_Training_Old < 0 ~ "Discordant (Young Up / Old Down)",
      TRUE ~ "Discordant (Young Down / Old Up)"
    ),
    border_col = ifelse(imputed, "black", "white"),
    bubble_alpha = case_when(
      significance == "NS"             ~ 0.30,
      significance == "Interaction"    ~ 0.50,
      significance == "Sig Both"       ~ 0.70,
      significance == "Sig Young only" ~ 0.80,
      significance == "Sig Old only"   ~ 0.80,
      TRUE ~ 0.30
    )
  )

# Correlation
cor_r   <- cor.test(scatter_df$logFC_Training_Young, scatter_df$logFC_Training_Old, method = "pearson")
cor_rho <- cor.test(scatter_df$logFC_Training_Young, scatter_df$logFC_Training_Old, method = "spearman")

# Sign concordance among proteins with |logFC| > 0.2 in at least one contrast
concordance_set <- scatter_df %>%
  filter(abs(logFC_Training_Young) > 0.2 | abs(logFC_Training_Old) > 0.2)
sign_concordance <- mean(sign(concordance_set$logFC_Training_Young) ==
                         sign(concordance_set$logFC_Training_Old)) * 100

xlim_range <- c(-2.5, 2)
ylim_range <- c(-1, 2)

# Quadrant counts
q_counts <- scatter_df %>%
  mutate(q = case_when(
    logFC_Training_Young > 0 & logFC_Training_Old > 0 ~ "Q1",
    logFC_Training_Young < 0 & logFC_Training_Old < 0 ~ "Q3",
    logFC_Training_Young > 0 & logFC_Training_Old < 0 ~ "Q4",
    TRUE ~ "Q2"
  )) %>% count(q) %>% deframe()

# Labels for significant proteins
label_df <- scatter_df %>%
  filter(significance != "NS") %>%
  group_by(significance) %>%
  arrange(desc(abs(logFC_Training_Young) + abs(logFC_Training_Old))) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    label_fill     = SIG_LABEL_FILL[as.character(significance)],
    label_text_col = SIG_LABEL_TEXT[as.character(significance)],
    nudge_y = case_when(
      significance == "Interaction"    ~ -0.03,
      significance == "Sig Young only" ~  0.03,
      significance == "Sig Both"       ~  0.04,
      significance == "Sig Old only"   ~ -0.04,
      TRUE ~ 0
    )
  )

# --- Sort: NS first (bottom layer), significant on top ---
plot_order <- scatter_df %>% arrange(desc(as.integer(significance)))

pC <- ggplot(plot_order, aes(x = logFC_Training_Young, y = logFC_Training_Old)) +
  # Quadrant shading
  annotate("rect", xmin = 0, xmax = Inf,  ymin = 0, ymax = Inf,
           fill = "#E88D6D", alpha = 0.18) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#E88D6D", alpha = 0.18) +
  annotate("rect", xmin = 0, xmax = Inf,  ymin = -Inf, ymax = 0,
           fill = "#7BAFD4", alpha = 0.18) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#7BAFD4", alpha = 0.18) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  # Bubble layer: fill = significance, border = imputation
  geom_point(aes(fill = significance),
             shape = 21,
             size = 1.8,
             color = plot_order$border_col,
             alpha = plot_order$bubble_alpha,
             stroke = 0.6) +
  scale_fill_manual(values = SIG_COLORS, name = "Significance",
                    guide = guide_legend(
                      order = 1,
                      override.aes = list(size = 3.5, alpha = 0.85,
                                          stroke = 0.6, color = "black"))) +
  # Gene labels (colored label boxes, matching Panel D style)
  geom_label_repel(data = label_df, aes(label = gene),
                   fill = label_df$label_fill, color = label_df$label_text_col,
                   nudge_y = label_df$nudge_y,
                   size = 2.2, fontface = "bold",
                   max.overlaps = 30,
                   segment.size = 0.2, segment.color = "grey50",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 0.5, force = 3, force_pull = 0.5,
                   label.padding = unit(1.5, "pt"),
                   label.r = unit(1, "pt"),
                   label.size = 0.15, seed = 42,
                   xlim = c(xlim_range[1] * 0.9, xlim_range[2] * 0.9),
                   ylim = c(ylim_range[1] * 0.9, ylim_range[2] * 0.9)) +
  # Quadrant labels (flush to panel corners)
  annotate("label", x = Inf, y = Inf,
           label = sprintf("Concordant Up\u2002n = %s", q_counts["Q1"]),
           hjust = 1, vjust = 1, size = 2.5, fontface = "bold",
           color = "#E88D6D", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = -Inf,
           label = sprintf("Concordant Down\u2002n = %s", q_counts["Q3"]),
           hjust = 0, vjust = 0, size = 2.5, fontface = "bold",
           color = "#E88D6D", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = Inf,
           label = sprintf("Discordant\u2002n = %s", q_counts["Q2"]),
           hjust = 0, vjust = 1, size = 2.5, fontface = "bold",
           color = "#7BAFD4", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = Inf, y = -Inf,
           label = sprintf("Discordant\u2002n = %s", q_counts["Q4"]),
           hjust = 1, vjust = 0, size = 2.5, fontface = "bold",
           color = "#7BAFD4", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  coord_cartesian(xlim = xlim_range, ylim = ylim_range, expand = FALSE) +
  labs(
    title = "Protein-Level Concordance of Training Response",
    subtitle = sprintf("logFC Training Young vs Old | %s proteins | r = %.2f, rho = %.2f, concordance = %.0f%% | border = imputation status",
                       format(nrow(scatter_df), big.mark = ","),
                       cor_r$estimate, cor_rho$estimate, sign_concordance),
    x = expression(log[2]*FC ~ "(Training Young)"),
    y = expression(log[2]*FC ~ "(Training Old)")
  ) +
  THEME_PUB +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 5.5),
        legend.title = element_text(size = 6, face = "bold"),
        legend.spacing.x = unit(4, "mm"))

# Imputation key strip (mirrors Panel D database key strip)
pC_imp_key_strip <- ggplot(tibble(x = c(1, 3), y = c(0, 0),
                                   label = c("Imputed", "Non-imputed")),
                            aes(x = x, y = y)) +
  annotate("text", x = 0.3, y = 0, label = "Border:",
           hjust = 0, size = 2.0, fontface = "bold", color = "grey30") +
  geom_point(shape = 21, size = 3.5, fill = "grey70",
             color = c("black", "white"), stroke = c(0.8, 1.2)) +
  geom_text(aes(label = label), hjust = -0.3, size = 1.8, color = "grey30") +
  scale_x_continuous(limits = c(-0.5, 5)) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

pC_combined <- pC / pC_imp_key_strip + plot_layout(heights = c(0.96, 0.04))

ggsave(file.path(RPT_DIR, "panel_C_concordance.pdf"), pC_combined,
       width = 200, height = 200, units = "mm", device = pdf)
ggsave(file.path(RPT_DIR, "panel_C_concordance.png"), pC_combined,
       width = 200, height = 200, units = "mm", dpi = 300)

# Clean CSV (now includes imputed column)
scatter_df %>%
  transmute(
    gene,
    logFC_Training_Young = round(logFC_Training_Young, 4),
    logFC_Training_Old   = round(logFC_Training_Old, 4),
    pi_score_Young       = round(pi_Young, 6),
    pi_score_Old         = round(pi_Old, 6),
    pi_score_Interaction = round(pi_Interaction, 6),
    significance         = as.character(significance),
    quadrant,
    imputed
  ) %>%
  arrange(significance, desc(abs(logFC_Training_Young) + abs(logFC_Training_Old))) %>%
  write_csv(file.path(DAT_DIR, "panel_C", "concordance.csv"))

message("  Panel C saved")

# ==============================================================================
# PANEL D — RRHO2 Concordance Map
# ==============================================================================

message("Panel D: RRHO2 concordance map...")

# Build ranked lists
rrho_list1 <- dep_df %>%
  dplyr::select(gene, t = t_Training_Young) %>%
  filter(!is.na(t)) %>% distinct(gene, .keep_all = TRUE) %>%
  as.data.frame()

rrho_list2 <- dep_df %>%
  dplyr::select(gene, t = t_Training_Old) %>%
  filter(!is.na(t)) %>% distinct(gene, .keep_all = TRUE) %>%
  as.data.frame()

shared_genes <- intersect(rrho_list1$gene, rrho_list2$gene)
rrho_list1 <- rrho_list1 %>% filter(gene %in% shared_genes)
rrho_list2 <- rrho_list2 %>% filter(gene %in% shared_genes)

rrho_obj <- RRHO2_initialize(
  list1 = rrho_list1, list2 = rrho_list2,
  labels = c("Training (Young)", "Training (Old)"),
  log10.ind = TRUE, method = "hyper", boundary = 0.04
)

hmat <- rrho_obj$hypermat
nr <- nrow(hmat); nc <- ncol(hmat)
mid_r <- floor(nr / 2); mid_c <- floor(nc / 2)

max_UU <- max(hmat[1:mid_r, 1:mid_c], na.rm = TRUE)
max_DD <- max(hmat[(mid_r+1):nr, (mid_c+1):nc], na.rm = TRUE)
max_UD <- max(hmat[1:mid_r, (mid_c+1):nc], na.rm = TRUE)
max_DU <- max(hmat[(mid_r+1):nr, 1:mid_c], na.rm = TRUE)

hmat_df <- expand.grid(row = 1:nr, col = 1:nc) %>%
  mutate(neg_log10_pvalue = as.vector(hmat))

pD <- ggplot(hmat_df, aes(x = row, y = col, fill = neg_log10_pvalue)) +
  geom_raster() +
  scale_fill_viridis_c(option = "viridis", name = expression(-log[10](P)),
                        guide = guide_colorbar(barwidth = unit(3, "cm"),
                                               barheight = unit(0.3, "cm"),
                                               title.position = "left",
                                               title.theme = element_text(size = 5.5, vjust = 0.8))) +
  geom_vline(xintercept = mid_r + 0.5, linetype = "dashed", color = "white", linewidth = 0.5) +
  geom_hline(yintercept = mid_c + 0.5, linetype = "dashed", color = "white", linewidth = 0.5) +
  annotate("text", x = mid_r * 0.5, y = mid_c * 0.5,
           label = sprintf("Concordant Up\nmax = %.1f", max_UU),
           color = "white", fontface = "bold", size = 2.5) +
  annotate("text", x = mid_r + (nr - mid_r) * 0.5, y = mid_c + (nc - mid_c) * 0.5,
           label = sprintf("Concordant Down\nmax = %.1f", max_DD),
           color = "white", fontface = "bold", size = 2.5) +
  annotate("text", x = mid_r * 0.5, y = mid_c + (nc - mid_c) * 0.5,
           label = sprintf("Discordant\nY Up / O Down\nmax = %.1f", max_UD),
           color = "white", fontface = "bold", size = 2.0) +
  annotate("text", x = mid_r + (nr - mid_r) * 0.5, y = mid_c * 0.5,
           label = sprintf("Discordant\nY Down / O Up\nmax = %.1f", max_DU),
           color = "white", fontface = "bold", size = 2.0) +
  # X-axis (Training Young) direction labels
  annotate("text", x = 1, y = -nc * 0.04,
           label = "<- Most upregulated", hjust = 0, size = 1.8, color = "grey30") +
  annotate("text", x = nr, y = -nc * 0.04,
           label = "Most downregulated ->", hjust = 1, size = 1.8, color = "grey30") +
  # Y-axis (Training Old) direction labels
  annotate("text", x = -nr * 0.04, y = 1, angle = 90,
           label = "<- Most upregulated", hjust = 0, size = 1.8, color = "grey30") +
  annotate("text", x = -nr * 0.04, y = nc, angle = 90,
           label = "Most downregulated ->", hjust = 1, size = 1.8, color = "grey30") +
  coord_cartesian(clip = "off") +
  labs(
    title = "Threshold-Free Concordance (RRHO2)",
    subtitle = sprintf("Hypergeometric overlap of ranked gene lists | %d shared genes", length(shared_genes)),
    x = "Training (Young) rank",
    y = "Training (Old) rank"
  ) +
  THEME_PUB +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom")

ggsave(file.path(RPT_DIR, "panel_D_rrho2.pdf"), pD,
       width = 180, height = 180, units = "mm", device = pdf)

# Clean CSV — matrix with metadata
rrho2_meta <- tibble(
  quadrant = c("Concordant_Up", "Concordant_Down",
               "Discordant_YoungUp_OldDown", "Discordant_YoungDown_OldUp"),
  max_neg_log10_pvalue = round(c(max_UU, max_DD, max_UD, max_DU), 2),
  matrix_rows = nr,
  matrix_cols = nc,
  n_shared_genes = length(shared_genes)
)
write_csv(rrho2_meta, file.path(DAT_DIR, "panel_D", "rrho2_summary.csv"))

# Also export the full matrix as a proper CSV with indices
rrho2_export <- as.data.frame(hmat)
colnames(rrho2_export) <- paste0("col_", 1:nc)
rrho2_export$row <- 1:nr
rrho2_export <- rrho2_export %>% dplyr::select(row, everything())
write_csv(rrho2_export, file.path(DAT_DIR, "panel_D", "rrho2_matrix.csv"))

message("  Panel D saved")

# ==============================================================================
# PANEL E — fGSEA NES Scatter (Hallmark + rrvgo-reduced GO:BP only)
# ==============================================================================

message("Panel E: NES scatter (Hallmark + GO:BP reduced)...")

# --- 1. Filter to Hallmark + GO:BP only ---
fgsea_hbp <- fgsea_all %>%
  filter(database %in% c("Hallmark", "GO:BP"),
         contrast %in% c("Training_Young", "Training_Old", "Interaction"))

# --- 2. Reduce GO:BP terms with rrvgo (threshold 0.85) ---
gobp_sig_names <- fgsea_hbp %>%
  filter(database == "GO:BP", padj < 0.05) %>%
  pull(pathway) %>% unique()

gobp_keep <- gobp_sig_names  # fallback: keep all

if (length(gobp_sig_names) > 5) {
  gobp_msigdb <- msigdbr(species = "Homo sapiens", collection = "C5",
                          subcollection = "GO:BP") %>%
    dplyr::select(gs_name, gs_exact_source) %>% distinct()
  name_to_goid <- setNames(gobp_msigdb$gs_exact_source, gobp_msigdb$gs_name)
  gobp_go_ids <- unique(na.omit(name_to_goid[gobp_sig_names]))

  if (length(gobp_go_ids) > 5) {
    tryCatch({
      hsGO <- GOSemSim::godata("org.Hs.eg.db", ont = "BP")
      simMatrix <- calculateSimMatrix(gobp_go_ids, orgdb = "org.Hs.eg.db",
                                       semdata = hsGO,
                                       ont = "BP", method = "Rel")
      # Build scores vector (mean -log10 padj across contrasts for each GO term)
      gobp_scores_df <- fgsea_hbp %>%
        filter(database == "GO:BP", padj < 0.05) %>%
        dplyr::select(pathway, padj) %>%
        left_join(gobp_msigdb %>% rename(pathway = gs_name), by = "pathway") %>%
        filter(gs_exact_source %in% gobp_go_ids) %>%
        group_by(gs_exact_source) %>%
        summarise(score = mean(-log10(padj)), .groups = "drop")
      scores_vec <- setNames(gobp_scores_df$score, gobp_scores_df$gs_exact_source)
      reducedTerms <- reduceSimMatrix(simMatrix, scores = scores_vec,
                                       threshold = 0.5,
                                       orgdb = "org.Hs.eg.db")
      goid_to_name <- setNames(gobp_msigdb$gs_name, gobp_msigdb$gs_exact_source)
      # Keep only parent (representative) terms, not all children
      parent_go_ids <- unique(reducedTerms$parent)
      gobp_keep <- unique(na.omit(goid_to_name[parent_go_ids]))
      message(sprintf("  rrvgo reduced GO:BP from %d to %d terms",
                      length(gobp_sig_names), length(gobp_keep)))
    }, error = function(e) {
      message("  rrvgo failed: ", e$message, " — keeping all GO:BP")
    })
  }
}

# Keep Hallmark pathways + reduced GO:BP pathways
fgsea_filtered <- fgsea_hbp %>%
  filter(database == "Hallmark" | pathway %in% gobp_keep)

# --- 3. Pivot wide and classify ---
fgsea_wide <- fgsea_filtered %>%
  dplyr::select(pathway, contrast, NES, padj, size, database) %>%
  pivot_wider(id_cols = c(pathway, database), names_from = contrast,
              values_from = c(NES, padj, size)) %>%
  filter(!is.na(NES_Training_Young), !is.na(NES_Training_Old)) %>%
  mutate(set_size = coalesce(size_Training_Young, size_Training_Old))

fgsea_sig <- fgsea_wide %>%
  filter(
    (!is.na(padj_Training_Young) & padj_Training_Young < 0.05) |
    (!is.na(padj_Training_Old)   & padj_Training_Old < 0.05) |
    (!is.na(padj_Interaction)    & padj_Interaction < 0.05)
  ) %>%
  mutate(
    sig_Y = !is.na(padj_Training_Young) & padj_Training_Young < 0.05,
    sig_O = !is.na(padj_Training_Old) & padj_Training_Old < 0.05,
    sig_I = !is.na(padj_Interaction) & padj_Interaction < 0.05,
    significance = case_when(
      sig_I         ~ "Interaction",
      sig_Y & sig_O ~ "Sig Both",
      sig_Y         ~ "Sig Young only",
      sig_O         ~ "Sig Old only",
      TRUE          ~ "NS"
    ) %>% factor(levels = names(SIG_COLORS)),
    pathway_label = clean_pathway_name(pathway)
  )

message(sprintf("  %d pathways after filtering (Hallmark: %d, GO:BP: %d)",
                nrow(fgsea_sig),
                sum(fgsea_sig$database == "Hallmark"),
                sum(fgsea_sig$database == "GO:BP")))

nes_cor <- cor.test(fgsea_sig$NES_Training_Young, fgsea_sig$NES_Training_Old)
nes_lim <- max(abs(c(fgsea_sig$NES_Training_Young, fgsea_sig$NES_Training_Old))) * 1.15

# Quadrant counts
nq1 <- sum(fgsea_sig$NES_Training_Young > 0 & fgsea_sig$NES_Training_Old > 0)
nq2 <- sum(fgsea_sig$NES_Training_Young < 0 & fgsea_sig$NES_Training_Old > 0)
nq3 <- sum(fgsea_sig$NES_Training_Young < 0 & fgsea_sig$NES_Training_Old < 0)
nq4 <- sum(fgsea_sig$NES_Training_Young > 0 & fgsea_sig$NES_Training_Old < 0)

# Label pathways with set_size >= 50
label_pw <- fgsea_sig %>%
  filter(set_size >= 50)

# Bubble border: Hallmark = solid black, GO:BP = transparent
fgsea_sig <- fgsea_sig %>%
  mutate(
    border_col = ifelse(database == "Hallmark", "black", "white"),
    bubble_alpha = case_when(
      significance == "NS"             ~ 0.45,
      significance == "Interaction"    ~ 0.60,
      significance == "Sig Both"       ~ 0.75,
      significance == "Sig Old only"   ~ 0.85,
      significance == "Sig Young only" ~ 0.85,
      TRUE ~ 0.60
    )
  )

label_pw <- label_pw %>%
  mutate(
    label_fill = SIG_LABEL_FILL[as.character(significance)],
    label_text_col = SIG_LABEL_TEXT[as.character(significance)]
  ) %>%
  # Condense verbose pathway labels for readability
  mutate(pathway_label = pathway_label %>%
    str_replace("Neg(ative)? Reg(ulation)? Of Programmed Cell Death", "Anti-Apoptosis") %>%
    str_replace("Pos(itive)? Reg(ulation)? Of Cellular Component.*", "Pos. Reg. Cell Comp. Org.") %>%
    str_replace("Positive Regulation Of ", "Pos. Reg. ") %>%
    str_replace("Negative Regulation Of ", "Neg. Reg. ") %>%
    str_replace("Regulation Of ", "Reg. ") %>%
    str_replace("Process Utilizing Autophagic Me.*", "Autophagy") %>%
    str_replace("Post Translational Protein Modi.*", "Protein PTMs") %>%
    str_replace("Proton Motive Force Driven.*", "PMF-Driven ATP Synthesis") %>%
    str_replace("ATP Synthesis Coupled Electron.*", "ETC / ATP Synthesis") %>%
    str_replace("Ribose Phosphate Biosynthetic.*", "Ribose-P Biosynthesis") %>%
    str_replace("Purine Containing Compound Bio.*", "Purine Biosynthesis") %>%
    str_replace("Modified Amino Acid Metabolic.*", "Modified AA Metabolism") %>%
    str_replace("Sulfur Compound Metabolic.*", "Sulfur Metabolism") %>%
    str_replace("Proteolysis Involved In Protei.*", "Protein Proteolysis") %>%
    str_replace("Muscle Structure Development", "Muscle Development") %>%
    str_replace("Mrna Metabolic.*", "mRNA Metabolism") %>%
    str_replace("Small Molecule Catabolic.*", "Small Molecule Catabolism") %>%
    str_replace("Lipid Catabolic Process", "Lipid Catabolism") %>%
    str_replace("Microtubule Based Process", "Microtubule Processes") %>%
    str_replace("Proton Transmembrane Transport", "H+ Transmembrane Transport") %>%
    str_replace(" Process$", "")
  ) %>%
  # Per-label nudge: bias same-colored labels vertically so they cluster,

  # while still using a single ggrepel layer for cross-group repulsion
  mutate(nudge_y = case_when(
    significance == "Interaction"    ~ -0.12,
    significance == "Sig Young only" ~  0.12,
    significance == "Sig Both"       ~  0.18,
    significance == "Sig Old only"   ~ -0.18,
    TRUE ~ 0
  )) %>%
  arrange(significance)

# Sort for plotting (less prominent drawn first)
plot_df <- fgsea_sig %>%
  mutate(draw_order = factor(significance,
    levels = c("NS", "Interaction", "Sig Both", "Sig Old only", "Sig Young only"))) %>%
  arrange(draw_order)

pE <- ggplot(plot_df, aes(x = NES_Training_Young, y = NES_Training_Old)) +
  # Quadrant shading
  annotate("rect", xmin = 0, xmax = Inf,  ymin = 0, ymax = Inf,
           fill = "#E88D6D", alpha = 0.18) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#E88D6D", alpha = 0.18) +
  annotate("rect", xmin = 0, xmax = Inf,  ymin = -Inf, ymax = 0,
           fill = "#7BAFD4", alpha = 0.18) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#7BAFD4", alpha = 0.18) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  # Bubbles: fill = significance, size = set_size, border = database
  geom_point(aes(fill = significance, size = set_size),
             shape = 21,
             color = plot_df$border_col,
             alpha = plot_df$bubble_alpha,
             stroke = 0.8) +
  scale_fill_manual(values = SIG_COLORS, name = "Significance",
                    guide = guide_legend(
                      order = 1,
                      override.aes = list(size = 3.5, alpha = 0.85,
                                          stroke = 0.8, color = "black"))) +
  scale_size_continuous(range = c(2, 8), name = "Set size",
                        breaks = c(20, 50, 100, 200),
                        guide = guide_legend(
                          order = 2,
                          override.aes = list(fill = "grey60",
                                              color = "black", alpha = 0.7))) +
  scale_x_continuous(expand = expansion(0, 0)) +
  scale_y_continuous(expand = expansion(0, 0)) +
  coord_cartesian(xlim = c(-3.5, 2.5), ylim = c(-2.5, 3)) +
  labs(
    title = "Pathway-Level Concordance (fGSEA)",
    subtitle = sprintf("Hallmark + GO:BP (rrvgo-reduced) | padj < 0.05 | %d pathways | r = %.2f, p %s",
                       nrow(fgsea_sig), nes_cor$estimate,
                       ifelse(nes_cor$p.value < 0.001, "< 0.001",
                              sprintf("= %.3f", nes_cor$p.value))),
    x = "NES (Training Young)",
    y = "NES (Training Old)"
  ) +
  THEME_PUB +
  theme(legend.position = "none")

# Pathway labels — single layer (all labels repel each other to prevent
# cross-color overlap) with per-significance nudge_y for same-color clustering
pE <- pE +
  geom_label_repel(data = label_pw, aes(label = pathway_label),
                   fill = label_pw$label_fill, color = label_pw$label_text_col,
                   nudge_y = label_pw$nudge_y,
                   size = 2.2, fontface = "bold",
                   max.overlaps = 30,
                   segment.size = 0.2, segment.color = "grey50",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 0.5, force = 3, force_pull = 0.5,
                   label.padding = unit(1.5, "pt"),
                   label.r = unit(1, "pt"),
                   label.size = 0.15, seed = 42,
                   xlim = c(-3.2, 2.2),
                   ylim = c(-2.2, 2.7))

# Quadrant count labels (on top of pathway labels)
pE <- pE +
  annotate("label", x = 2.5, y = 3.0,
           label = sprintf("Concordant Up\u2002n = %d", nq1),
           hjust = 1, vjust = 1, size = 2.5, fontface = "bold",
           color = "#E88D6D", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -3.5, y = -2.5,
           label = sprintf("Concordant Down\u2002n = %d", nq3),
           hjust = 0, vjust = 0, size = 2.5, fontface = "bold",
           color = "#E88D6D", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -3.5, y = 3.0,
           label = sprintf("Discordant\u2002n = %d", nq2),
           hjust = 0, vjust = 1, size = 2.5, fontface = "bold",
           color = "#7BAFD4", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = 2.5, y = -2.5,
           label = sprintf("Discordant\u2002n = %d", nq4),
           hjust = 1, vjust = 0, size = 2.5, fontface = "bold",
           color = "#7BAFD4", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt"))

# --- Hand-built legend: three columns side-by-side, items stacked vertically ---

sig_levels_e <- c("Interaction", "Sig Both", "Sig Young only", "Sig Old only")
ks_e <- 0.15

# Column 1: Significance
sig_key_df <- tibble(
  x = 0, y = rev(seq_along(sig_levels_e)) * ks_e,
  label = sig_levels_e,
  fill  = unname(SIG_COLORS[sig_levels_e])
)
# Column 2: Set size
size_breaks_e <- c(20, 50, 100)
size_range_e  <- c(2, 8)
size_key_df <- tibble(
  x = 3.5, y = rev(seq_along(size_breaks_e)) * ks_e,
  label = as.character(size_breaks_e),
  pt_size = scales::rescale(size_breaks_e, to = size_range_e, from = c(20, 200))
)
# Column 3: Database
db_key_df <- tibble(
  x = 6.0, y = c(2, 1) * ks_e,
  label = c("Hallmark", "GO:BP"),
  border = c("black", "white"),
  stroke = c(0.8, 1.2)
)

title_y_e <- (max(length(sig_levels_e), length(size_breaks_e), 2) + 1) * ks_e

pE_key <- ggplot() +
  # Significance column
  annotate("text", x = 0, y = title_y_e, label = "Significance",
           hjust = 0, size = 2.0, fontface = "bold", color = "grey30") +
  geom_point(data = sig_key_df, aes(x = x, y = y),
             shape = 21, size = 3.5, fill = sig_key_df$fill,
             color = "black", stroke = 0.8) +
  geom_text(data = sig_key_df, aes(x = x + 0.3, y = y, label = label),
            hjust = 0, size = 1.8, color = "grey30") +
  # Set size column
  annotate("text", x = 3.5, y = title_y_e, label = "Set size",
           hjust = 0, size = 2.0, fontface = "bold", color = "grey30") +
  geom_point(data = size_key_df, aes(x = x, y = y),
             shape = 21, size = size_key_df$pt_size, fill = "grey60",
             color = "black", alpha = 0.7) +
  geom_text(data = size_key_df, aes(x = x + 0.3, y = y, label = label),
            hjust = 0, size = 1.8, color = "grey30") +
  # Database column
  annotate("text", x = 6.0, y = title_y_e, label = "Database",
           hjust = 0, size = 2.0, fontface = "bold", color = "grey30") +
  geom_point(data = db_key_df, aes(x = x, y = y),
             shape = 21, size = 3.5, fill = "grey70",
             color = db_key_df$border, stroke = db_key_df$stroke) +
  geom_text(data = db_key_df, aes(x = x + 0.3, y = y, label = label),
            hjust = 0, size = 1.8, color = "grey30") +
  scale_x_continuous(limits = c(-0.3, 8.5)) +
  scale_y_continuous(limits = c(0, title_y_e + ks_e)) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

pE_combined <- pE / pE_key + plot_layout(heights = c(0.90, 0.10))

ggsave(file.path(RPT_DIR, "panel_E_nes_bubble.pdf"), pE_combined,
       width = 200, height = 200, units = "mm", device = pdf)
ggsave(file.path(RPT_DIR, "panel_E_nes_bubble.png"), pE_combined,
       width = 200, height = 200, units = "mm", dpi = 300)

# Clean CSV
fgsea_sig %>%
  transmute(
    pathway,
    pathway_label,
    database,
    NES_Training_Young = round(NES_Training_Young, 3),
    NES_Training_Old   = round(NES_Training_Old, 3),
    padj_Training_Young = signif(padj_Training_Young, 4),
    padj_Training_Old   = signif(padj_Training_Old, 4),
    padj_Interaction    = signif(padj_Interaction, 4),
    significance        = as.character(significance),
    set_size
  ) %>%
  arrange(significance, desc(abs(NES_Training_Young) + abs(NES_Training_Old))) %>%
  write_csv(file.path(DAT_DIR, "panel_E", "nes_scatter.csv"))

message("  Panel E saved")

# ==============================================================================
# SHARED: Interaction category definitions (used by Panel F)
# ==============================================================================

CATEGORY_COLORS <- c("Down Young / Up Old" = "#E05A4E",
                     "Attenuated"          = "#5DA5DA",
                     "Up Young / Down Old" = "#9B7FBF")
CATEGORY_ORDER  <- c("Down Young / Up Old", "Attenuated", "Up Young / Down Old")

# 4-pattern response classification (biological gradient for merged panel)
PATTERN_COLORS <- c(
  "Concordant"  = "#78909C",
  "Attenuated"  = "#5DA5DA",
  "Blunted"     = "#F4A460",
  "Reversed"    = "#E05A4E"
)
PATTERN_ORDER <- c("Concordant", "Attenuated", "Blunted", "Reversed")

# Classify all interaction DEPs into categories
int_class_all <- dep_df %>%
  filter(sig_pi_Interaction == TRUE) %>%
  transmute(
    gene,
    logFC_Y = logFC_Training_Young,
    logFC_O = logFC_Training_Old,
    logFC_Aging = logFC_Aging,
    divergence = abs(logFC_Training_Young - logFC_Training_Old),
    category = case_when(
      abs(logFC_O) < 0.15 ~ "Attenuated",
      logFC_Y > 0 & logFC_O < 0 ~ "Up Young / Down Old",
      logFC_Y < 0 & logFC_O > 0 ~ "Down Young / Up Old",
      sign(logFC_Y) == sign(logFC_O) & abs(logFC_O) < abs(logFC_Y) * 0.5 ~ "Attenuated",
      sign(logFC_Y) == sign(logFC_O) & abs(logFC_Y) < abs(logFC_O) * 0.5 ~ "Attenuated",
      logFC_Y < 0 & logFC_O <= 0 ~ "Down Young / Up Old",
      logFC_Y > 0 & logFC_O >= 0 ~ "Up Young / Down Old",
      TRUE ~ "Attenuated"
    )
  ) %>% filter(!is.na(logFC_Y), !is.na(logFC_O)) %>%
  mutate(
    response_pattern = case_when(
      logFC_O < -0.20                        ~ "Concordant",
      logFC_O >= -0.20 & logFC_O <= 0.15     ~ "Attenuated",
      logFC_O >  0.15  & logFC_O <= 0.30     ~ "Blunted",
      logFC_O >  0.30                         ~ "Reversed",
      TRUE                                    ~ "Attenuated"
    ) %>% factor(levels = PATTERN_ORDER)
  )

# Named lookups
gene_category_lookup <- setNames(int_class_all$category, int_class_all$gene)
gene_pattern_lookup  <- setNames(as.character(int_class_all$response_pattern),
                                  int_class_all$gene)

# ==============================================================================
# SHARED: ORA + Pathway Mapping (used by Panel F)
# ==============================================================================

message("Shared: ORA + pathway mapping for Panel F...")

int_class <- int_class_all
cat_counts <- int_class %>% count(category) %>% deframe()
int_class <- int_class %>%
  mutate(facet_label = sprintf("%s (n = %d)", category, cat_counts[category]))

all_int_genes   <- int_class$gene
all_genes       <- dep_df$gene

# --- ORA: Hallmark + GO:BP ---
hallmark_t2g <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  select(term = gs_name, gene = gene_symbol) %>% distinct()
gobp_t2g <- msigdbr(species = "Homo sapiens", collection = "C5",
                     subcollection = "GO:BP") %>%
  select(term = gs_name, gene = gene_symbol) %>% distinct()

run_ora <- function(t2g, db_name) {
  res <- tryCatch(
    enricher(gene = all_int_genes, universe = all_genes, TERM2GENE = t2g,
             pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
             minGSSize = 3, maxGSSize = 500),
    error = function(e) NULL)
  if (!is.null(res) && nrow(as.data.frame(res)) > 0)
    as.data.frame(res) %>% mutate(database = db_name) else tibble()
}

ora_combined <- bind_rows(run_ora(hallmark_t2g, "Hallmark"),
                          run_ora(gobp_t2g, "GO:BP")) %>%
  mutate(pathway_label = clean_pathway_name(ID))

# --- Select top pathways + greedy coverage ---
ora_sig <- ora_combined %>% filter(pvalue < 0.05) %>% arrange(pvalue)
ora_top <- bind_rows(
  ora_sig %>% filter(database == "Hallmark") %>% slice_head(n = 5),
  ora_sig %>% filter(database == "GO:BP")    %>% slice_head(n = 5)
) %>% group_by(pathway_label) %>% slice_min(pvalue, n = 1, with_ties = FALSE) %>%
  ungroup()
if (nrow(ora_top) < 4) ora_top <- ora_combined %>% arrange(pvalue) %>% slice_head(n = 8)

all_t2g    <- bind_rows(hallmark_t2g, gobp_t2g) %>% distinct()
mapped_now <- ora_top %>% select(geneID) %>% separate_rows(geneID, sep = "/") %>%
  filter(geneID %in% all_int_genes) %>% pull(geneID) %>% unique()
orphans    <- setdiff(all_int_genes, mapped_now)
remaining  <- ora_combined %>% filter(!ID %in% ora_top$ID) %>% arrange(pvalue)

# Greedy rescue: add pathways for unmapped genes (cap total at ~12)
max_total_pws <- 15
for (i in seq_len(15)) {
  if (length(orphans) == 0 || nrow(remaining) == 0) break
  if (nrow(ora_top) >= max_total_pws) break
  coverage <- remaining %>% rowwise() %>%
    mutate(hits = sum(str_split(geneID, "/")[[1]] %in% orphans)) %>%
    ungroup() %>% filter(hits >= 1) %>% arrange(desc(hits), pvalue)
  if (nrow(coverage) == 0) break
  best <- coverage %>% slice_head(n = 1)
  ora_top   <- bind_rows(ora_top, best)
  remaining <- remaining %>% filter(ID != best$ID)
  orphans   <- setdiff(orphans, str_split(best$geneID, "/")[[1]])
}
ora_top <- ora_top %>% group_by(pathway_label) %>%
  slice_min(pvalue, n = 1, with_ties = FALSE) %>% ungroup()

# --- Build ALL gene-pathway links (for ORA stats) ---
sankey_links_all <- ora_top %>%
  select(pathway = pathway_label, ID, geneID, pvalue, p.adjust, Count, database) %>%
  separate_rows(geneID, sep = "/") %>% rename(gene = geneID) %>%
  filter(gene %in% all_int_genes)

# Rescue orphans via membership in selected pathways
leftover <- setdiff(all_int_genes, unique(sankey_links_all$gene))
if (length(leftover) > 0) {
  rescue <- all_t2g %>% filter(gene %in% leftover, term %in% ora_top$ID) %>%
    left_join(ora_top %>% select(ID, pathway = pathway_label, pvalue, p.adjust,
                                  Count, database), by = c("term" = "ID")) %>%
    group_by(gene) %>% slice_min(pvalue, n = 1, with_ties = FALSE) %>% ungroup() %>%
    mutate(ID = term) %>% select(gene, pathway, ID, pvalue, p.adjust, Count, database)
  sankey_links_all <- bind_rows(sankey_links_all, rescue)
  leftover <- setdiff(leftover, rescue$gene)
}

# Final catch-all: orphans that matched no selected pathway -> "Other"
leftover <- setdiff(all_int_genes, unique(sankey_links_all$gene))
if (length(leftover) > 0) {
  other_links <- tibble(
    gene     = leftover,
    pathway  = "Other",
    ID       = "OTHER",
    pvalue   = 1,
    p.adjust = 1,
    Count    = length(leftover),
    database = "Other"
  )
  sankey_links_all <- bind_rows(sankey_links_all, other_links)

  # Add "Other" to ora_top so downstream code picks it up
  ora_top <- bind_rows(ora_top, tibble(
    ID = "OTHER", Description = "Other", GeneRatio = NA, BgRatio = NA,
    pvalue = 1, p.adjust = 1, qvalue = 1, geneID = paste(leftover, collapse = "/"),
    Count = length(leftover), database = "Other", pathway_label = "Other"
  ))
  message(sprintf("  %d orphan genes assigned to 'Other' pathway", length(leftover)))
}

sankey_links_all <- sankey_links_all %>%
  left_join(int_class %>% select(gene, category, logFC_Y, logFC_O, divergence),
            by = "gene")
mapped_class <- int_class %>% filter(gene %in% unique(sankey_links_all$gene))

# --- Force 1:1 mapping: each gene -> best pathway ---
links_1to1 <- sankey_links_all %>%
  group_by(gene) %>%
  slice_min(pvalue, n = 1, with_ties = FALSE) %>%
  ungroup()

# Redistribute from over-represented pathways
pw_counts <- links_1to1 %>% count(pathway, name = "n_1to1")
for (iter in seq_len(5)) {
  solo_pws <- pw_counts %>% filter(n_1to1 == 1) %>% pull(pathway)
  if (length(solo_pws) == 0) break
  solo_genes <- links_1to1 %>% filter(pathway %in% solo_pws) %>% pull(gene)
  multi_pws <- pw_counts %>% filter(n_1to1 >= 2) %>% pull(pathway)
  better <- sankey_links_all %>%
    filter(gene %in% solo_genes, pathway %in% multi_pws) %>%
    group_by(gene) %>% slice_min(pvalue, n = 1, with_ties = FALSE) %>% ungroup()
  if (nrow(better) == 0) break
  links_1to1 <- bind_rows(
    links_1to1 %>% filter(!gene %in% better$gene),
    better)
  pw_counts <- links_1to1 %>% count(pathway, name = "n_1to1")
}

# Drop pathways that lost all genes after 1:1 assignment
active_pws <- unique(links_1to1$pathway)
ora_top <- ora_top %>% filter(pathway_label %in% active_pws)
mapped_class <- int_class %>% filter(gene %in% unique(links_1to1$gene))
leftover <- setdiff(all_int_genes, unique(links_1to1$gene))

message(sprintf("  %d pathways, %d mapped (1:1), %d excluded",
                nrow(ora_top), nrow(links_1to1), length(leftover)))

# --- Pathway and gene ordering (used by Panel F) ---
pw_1to1_counts <- links_1to1 %>% count(pathway, name = "n_1to1")
pw_pvals <- ora_top %>% select(pathway_label, pvalue)
pw_order_df <- pw_1to1_counts %>%
  left_join(pw_pvals, by = c("pathway" = "pathway_label")) %>%
  arrange(desc(pvalue))  # least significant first -> level 1
pw_order <- pw_order_df$pathway
# Ensure "Other" is always last (remove from sorted position, re-append)
if ("Other" %in% unique(links_1to1$pathway)) {
  pw_order <- setdiff(pw_order, "Other")
  pw_order <- c(pw_order, "Other")
}

# Gene order: pattern-based (biological gradient), then pathway within group
gene_order <- links_1to1 %>%
  mutate(
    pattern     = gene_pattern_lookup[gene],
    pattern_idx = match(pattern, PATTERN_ORDER),
    pw_idx      = match(pathway, pw_order)
  ) %>%
  arrange(pattern_idx, pw_idx, logFC_O) %>%
  pull(gene)

message("  Shared ORA + pathway mapping complete")

# ==============================================================================
# PANEL F — Data Preparation (Interaction DEP multi-contrast data)
# ==============================================================================

message("Panel F data prep: interaction DEP multi-contrast data...")

# Select interaction DEPs and gather logFC across all 3 contrasts
int_deps <- dep_df %>%
  filter(sig_pi_Interaction == TRUE) %>%
  transmute(
    gene,
    logFC_Aging          = logFC_Aging,
    logFC_Training_Young = logFC_Training_Young,
    logFC_Training_Old   = logFC_Training_Old,
    pi_score_Interaction = pi_score_Interaction,
    # Compute interaction magnitude for ordering
    interaction_magnitude = abs(logFC_Training_Young - logFC_Training_Old)
  ) %>%
  filter(!is.na(logFC_Aging), !is.na(logFC_Training_Young), !is.na(logFC_Training_Old)) %>%
  arrange(desc(interaction_magnitude))

n_int <- nrow(int_deps)
message(sprintf("  %d interaction DEPs", n_int))

# Truncate to top 30 if too many for readability
max_show <- 30
if (n_int > max_show) {
  int_deps <- int_deps %>% slice_head(n = max_show)
  message(sprintf("  Showing top %d by interaction magnitude", max_show))
}

# Shared ordering: use pathway-grouped gene_order from Sankey;
# append any Panel F proteins not in Sankey at the end (by magnitude)
e_genes <- as.character(int_deps$gene)
shared_order <- c(gene_order, setdiff(e_genes, gene_order))
int_deps <- int_deps %>%
  mutate(gene = factor(gene, levels = rev(shared_order)))

# Assign each gene its interaction category color
int_deps <- int_deps %>%
  mutate(cat = gene_category_lookup[as.character(gene)])

# Pivot to long format for plotting
int_long <- int_deps %>%
  pivot_longer(
    cols = c(logFC_Aging, logFC_Training_Young, logFC_Training_Old),
    names_to = "contrast",
    values_to = "logFC",
    names_prefix = "logFC_"
  ) %>%
  mutate(
    contrast_label = case_when(
      contrast == "Aging"          ~ "Aging",
      contrast == "Training_Young" ~ "Training (Young)",
      contrast == "Training_Old"   ~ "Training (Old)"
    ),
    contrast_label = factor(contrast_label,
                            levels = c("Aging", "Training (Young)", "Training (Old)"))
  )

message("  Panel F data prep complete")

# ==============================================================================
# PANEL F — Shared Data Preparation (ORA, Sankey coordinates)
# ==============================================================================

message("Panel F data prep: Sankey + bar coordinates...")

# Save classification for reference
write_csv(int_class, file.path(DAT_DIR, "panel_F", "interaction_classification.csv"))

# ORA, pathway mapping, 1:1 assignment, gene_order, pw_order already computed
# in the shared section above.

# =========================================================================
# Sankey + Bar chart coordinate construction (used by merged Panel F)
# =========================================================================

# Pathway colors: muted pastels matching SRplot aesthetic
pw_palette <- c(
  "#F48FB1", "#FDAE91", "#E8E8A0", "#A8D8A8", "#8DD3C7",
  "#A2CEE5", "#B6C8E8", "#DEB4D4", "#C9A9A6", "#AED581",
  "#CE93D8", "#80DEEA", "#FFD54F", "#90A4AE", "#E57373",
  "#FFCC80", "#9FA8DA", "#B39DDB", "#80CBC4", "#81C784"
)
PW_COLORS <- setNames(pw_palette[seq_along(pw_order)], pw_order)
PW_COLORS["Other"] <- "#D0D0D0"  # neutral grey for catch-all

# Map each gene -> its interaction category color (from Panel F)
gene_cat <- setNames(mapped_class$category, mapped_class$gene)

# ---- Coordinate system: shared numeric y-space ----
n_genes <- nrow(links_1to1)
n_pw    <- length(pw_order)
Y_SPAN  <- n_genes

X_GENE <- 1.0;  X_PW <- 3.0;  BAR_W <- 0.06

gene_h   <- Y_SPAN / (n_genes * 1.15)
gene_gap <- (Y_SPAN - n_genes * gene_h) / max(n_genes - 1, 1)

pw_h   <- Y_SPAN / (n_pw * 1.4)
pw_gap <- (Y_SPAN - n_pw * pw_h) / max(n_pw - 1, 1)

# ---- Gene bar data frame (top-to-bottom in gene_order) ----
gene_bar_df <- tibble(
  gene = gene_order,
  idx  = seq_along(gene_order)
) %>% mutate(
  ymax = Y_SPAN - (idx - 1) * (gene_h + gene_gap),
  ymin = ymax - gene_h,
  y_ctr = (ymin + ymax) / 2,
  fill_col = PATTERN_COLORS[gene_pattern_lookup[gene]]
)

# ---- Pathway bar data frame (top-to-bottom in pw_order) ----
pw_bar_df <- tibble(
  pathway = pw_order,
  idx     = seq_along(pw_order)
) %>% mutate(
  ymax = Y_SPAN - (idx - 1) * (pw_h + pw_gap),
  ymin = ymax - pw_h,
  y_ctr = (ymin + ymax) / 2,
  fill_col   = PW_COLORS[pathway]
)

# ---- Stacking within pathway bars: gene slots ----
slot_df <- links_1to1 %>%
  mutate(pw_idx = match(pathway, pw_order)) %>%
  arrange(pw_idx, match(gene, gene_order)) %>%
  group_by(pathway) %>%
  mutate(
    k        = n(),
    slot_idx = row_number()
  ) %>%
  ungroup() %>%
  left_join(pw_bar_df %>% select(pathway, pw_ymin = ymin, pw_ymax = ymax),
            by = "pathway") %>%
  mutate(
    slot_h   = (pw_ymax - pw_ymin) / k,
    slot_ymax = pw_ymax - (slot_idx - 1) * slot_h,
    slot_ymin = slot_ymax - slot_h
  )

# ---- Sigmoid ribbon polygons ----
make_sigmoid_ribbon <- function(x0, x1, y0_top, y0_bot, y1_top, y1_bot,
                                n_pts = 50, ribbon_id) {
  t <- seq(0, 1, length.out = n_pts)
  blend <- (1 - cos(pi * t)) / 2
  tibble(
    x = c(x0 + (x1 - x0) * t, rev(x0 + (x1 - x0) * t)),
    y = c(y0_top + (y1_top - y0_top) * blend,
          rev(y0_bot + (y1_bot - y0_bot) * blend)),
    ribbon_id = ribbon_id
  )
}

ribbon_input <- slot_df %>%
  left_join(gene_bar_df %>% select(gene, gene_ymin = ymin, gene_ymax = ymax),
            by = "gene")

ribbon_polys <- pmap_dfr(ribbon_input, function(gene, pathway,
                                                  gene_ymin, gene_ymax,
                                                  slot_ymin, slot_ymax, ...) {
  make_sigmoid_ribbon(
    x0 = X_GENE + BAR_W / 2,  x1 = X_PW - BAR_W / 2,
    y0_top = gene_ymax, y0_bot = gene_ymin,
    y1_top = slot_ymax, y1_bot = slot_ymin,
    ribbon_id = paste(gene, pathway, sep = "->")
  ) %>% mutate(pathway = pathway, gene = gene)
})

# Ribbon fill: pathway color
ribbon_polys <- ribbon_polys %>%
  mutate(fill_col = PW_COLORS[pathway])

# Draw order: crossing ribbons behind shorter ones
ribbon_order <- ribbon_input %>%
  left_join(gene_bar_df %>% select(gene, gene_yctr = y_ctr), by = "gene") %>%
  mutate(
    slot_yctr  = (slot_ymin + slot_ymax) / 2,
    crossing_d = abs(gene_yctr - slot_yctr),
    ribbon_id  = paste(gene, pathway, sep = "->")
  ) %>%
  arrange(desc(crossing_d)) %>%
  pull(ribbon_id)

ribbon_polys <- ribbon_polys %>%
  mutate(ribbon_id = factor(ribbon_id, levels = ribbon_order))

# ---- Dot plot data (used by merged Panel F bar sub-panel + CSV export) ----
srplot_summary <- links_1to1 %>%
  group_by(pathway) %>%
  summarise(Count_1to1 = n(), .groups = "drop")

dot_df <- ora_top %>%
  left_join(srplot_summary, by = c("pathway_label" = "pathway")) %>%
  left_join(pw_bar_df %>% transmute(pathway_label = as.character(pathway),
                                     dot_y = y_ctr),
            by = "pathway_label") %>%
  mutate(gene_ratio  = Count_1to1 / length(all_int_genes),
         neg_log10_p = -log10(pvalue))

# --- Export CSVs (organized by panel) ---
dot_df %>%
  transmute(pathway = pathway_label, database, gene_count = Count_1to1,
            gene_ratio = round(gene_ratio, 4), pvalue = signif(pvalue, 4),
            p_adjust = signif(p.adjust, 4),
            genes = ora_top$geneID[match(pathway_label, ora_top$pathway_label)]) %>%
  arrange(pvalue) %>% write_csv(file.path(DAT_DIR, "panel_F", "sankey_dot.csv"))

links_1to1 %>%
  transmute(gene = as.character(gene), pathway = as.character(pathway), database,
            interaction_pattern = category, logFC_Training_Young = round(logFC_Y, 4),
            logFC_Training_Old = round(logFC_O, 4),
            pathway_pvalue = signif(pvalue, 4)) %>%
  arrange(interaction_pattern, pathway, gene) %>%
  write_csv(file.path(DAT_DIR, "panel_F", "sankey_links.csv"))

mapped_class %>%
  transmute(gene, interaction_pattern = category,
            logFC_Training_Young = round(logFC_Y, 4),
            logFC_Training_Old = round(logFC_O, 4),
            logFC_Aging = round(logFC_Aging, 4),
            divergence = round(divergence, 4)) %>%
  arrange(interaction_pattern, desc(divergence)) %>%
  write_csv(file.path(DAT_DIR, "panel_F", "interaction_patterns.csv"))

# --- SRplot-format export (bioinformatics.com.cn Sankey-Dot tool input) ---
# Use 1:1 assigned genes (matching the Sankey) — not raw ORA gene lists
srplot_genes <- links_1to1 %>%
  group_by(pathway) %>%
  summarise(geneID = paste(sort(as.character(gene)), collapse = "/"),
            Count  = n(), .groups = "drop")

srplot_df <- ora_top %>%
  filter(pathway_label %in% active_pws) %>%
  select(Description = pathway_label, pvalue) %>%
  left_join(srplot_genes, by = c("Description" = "pathway")) %>%
  mutate(GeneRatio = Count / length(all_int_genes)) %>%
  select(Description, GeneRatio, pvalue, geneID, Count) %>%
  filter(pvalue > 0) %>%   # SRplot requires no zero p-values
  arrange(pvalue)

write_csv(srplot_df, file.path(DAT_DIR, "panel_F", "srplot_input.csv"))
message(sprintf("  SRplot input: %d pathways exported", nrow(srplot_df)))

# Export interaction dot CSVs (Panel F data)
int_deps %>%
  transmute(
    gene         = as.character(gene),
    logFC_Aging          = round(logFC_Aging, 4),
    logFC_Training_Young = round(logFC_Training_Young, 4),
    logFC_Training_Old   = round(logFC_Training_Old, 4),
    pi_score_Interaction = round(pi_score_Interaction, 6),
    interaction_magnitude = round(interaction_magnitude, 4)
  ) %>%
  arrange(desc(interaction_magnitude)) %>%
  write_csv(file.path(DAT_DIR, "panel_F", "interaction_dot.csv"))

int_long %>%
  transmute(
    gene = as.character(gene),
    contrast = contrast_label,
    logFC = round(logFC, 4)
  ) %>%
  write_csv(file.path(DAT_DIR, "panel_F", "interaction_dot_long.csv"))

message("  Panel F data prep + CSV exports complete")

# ==============================================================================
# PANEL F — Interaction DEPs: Multi-Contrast Response & Pathway Enrichment
# ==============================================================================

message("Panel F: merged interaction panel...")

# --- Dumbbell data: all 36 genes in Sankey gene_order, y from gene_bar_df ---
dumbbell_df <- int_class %>%
  filter(gene %in% gene_order) %>%
  left_join(gene_bar_df %>% select(gene, y_ctr), by = "gene")

db_xrange <- range(c(dumbbell_df$logFC_Y, dumbbell_df$logFC_O,
                      dumbbell_df$logFC_Aging), na.rm = TRUE)
db_xpad <- diff(db_xrange) * 0.08

# --- Pattern group separators: y-positions between adjacent groups ---
dumbbell_df <- dumbbell_df %>%
  mutate(pattern = gene_pattern_lookup[gene])

pattern_genes <- split(gene_order, gene_pattern_lookup[gene_order])
sep_ys <- numeric(0)
for (i in seq_len(length(PATTERN_ORDER) - 1)) {
  grp_cur  <- pattern_genes[[PATTERN_ORDER[i]]]
  grp_next <- pattern_genes[[PATTERN_ORDER[i + 1]]]
  if (length(grp_cur) > 0 && length(grp_next) > 0) {
    y_last  <- gene_bar_df$ymin[gene_bar_df$gene == tail(grp_cur, 1)]
    y_first <- gene_bar_df$ymax[gene_bar_df$gene == head(grp_next, 1)]
    sep_ys  <- c(sep_ys, (y_last + y_first) / 2)
  }
}

# --- Pattern group midpoints for bracket labels ---
pattern_label_df <- tibble(pattern = PATTERN_ORDER) %>%
  rowwise() %>%
  mutate(
    genes_in = list(pattern_genes[[pattern]]),
    n_genes  = length(genes_in),
    mid_y    = if (n_genes > 0)
      mean(gene_bar_df$y_ctr[gene_bar_df$gene %in% genes_in]) else NA_real_
  ) %>%
  ungroup() %>%
  filter(n_genes > 0)

# --- Dumbbell sub-panel (left): logFC across contrasts per gene ---
pF_dumbbell <- ggplot(dumbbell_df) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  # Pattern group separators
  geom_hline(yintercept = sep_ys, linetype = "dotted", color = "grey50", linewidth = 0.4) +
  # Connecting line: Training Young to Training Old
  geom_segment(aes(x = logFC_Y, xend = logFC_O, y = y_ctr, yend = y_ctr),
               color = "grey65", linewidth = 0.4) +
  # Aging: open circle
  geom_point(aes(x = logFC_Aging, y = y_ctr), shape = 1, size = 1.8, stroke = 0.5,
             color = CONTRAST_COLORS["Aging"]) +
  # Training Young: filled
  geom_point(aes(x = logFC_Y, y = y_ctr), shape = 16, size = 1.8,
             color = CONTRAST_COLORS["Training_Young"]) +
  # Training Old: filled
  geom_point(aes(x = logFC_O, y = y_ctr), shape = 16, size = 1.8,
             color = CONTRAST_COLORS["Training_Old"]) +
  # Pattern group labels on the right edge
  annotate("text",
           x = db_xrange[2] + db_xpad * 0.9,
           y = pattern_label_df$mid_y,
           label = pattern_label_df$pattern,
           hjust = 0, size = 2.2, fontface = "bold.italic",
           color = PATTERN_COLORS[pattern_label_df$pattern]) +
  scale_y_continuous(limits = c(0, Y_SPAN), expand = expansion(mult = 0.02),
                     breaks = gene_bar_df$y_ctr, labels = gene_bar_df$gene) +
  scale_x_continuous(position = "top") +
  coord_cartesian(xlim = c(db_xrange[1] - db_xpad,
                            db_xrange[2] + db_xpad * 3.5),
                  clip = "off") +
  labs(x = expression(log[2]~FC), y = NULL) +
  THEME_PUB +
  theme(
    axis.text.y        = element_text(size = 6.5, face = "italic", hjust = 1),
    axis.ticks.y       = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    panel.border       = element_rect(color = "black", linewidth = 0.6),
    plot.margin        = margin(8, 0, 8, 2)
  )

# --- Sankey sub-panel (center): gene bars + ribbons + pathway bars + labels ---
# Filter out "Other" — orphan genes keep their bars but have no ribbon/pathway
ribbon_polys_ef <- ribbon_polys %>% filter(pathway != "Other")
pw_bar_ef       <- pw_bar_df %>% filter(pathway != "Other")

pF_sankey <- ggplot() +
  # Sigmoid ribbons (excluding "Other")
  geom_polygon(data = ribbon_polys_ef %>% arrange(ribbon_id),
               aes(x = x, y = y, group = ribbon_id, fill = fill_col),
               alpha = 0.32, color = NA) +
  scale_fill_identity() +
  # Gene bars (category-colored, ALL genes including orphans)
  geom_rect(data = gene_bar_df,
            aes(xmin = X_GENE - BAR_W / 2, xmax = X_GENE + BAR_W / 2,
                ymin = ymin, ymax = ymax),
            fill = gene_bar_df$fill_col, color = NA) +
  # Pathway bars (excluding "Other")
  geom_rect(data = pw_bar_ef,
            aes(xmin = X_PW - BAR_W / 2, xmax = X_PW + BAR_W / 2,
                ymin = ymin, ymax = ymax),
            fill = pw_bar_ef$fill_col, color = NA) +
  # Pathway labels (left of bars, excluding "Other")
  geom_text(data = pw_bar_ef,
            aes(x = X_PW - BAR_W / 2 - 0.05, y = y_ctr, label = pathway),
            hjust = 1, size = 2.8, fontface = "bold") +
  scale_y_continuous(limits = c(0, Y_SPAN), expand = expansion(mult = 0.02)) +
  coord_cartesian(xlim = c(X_GENE - 0.05, X_PW + 0.15), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(8, 0, 8, 0))

# --- Bar sub-panel (right): enrichment statistics (excluding "Other") ---
dot_df_ef <- dot_df %>% filter(pathway_label != "Other")

# Bar height matches pathway bar height for visual alignment
bar_half_h <- pw_h * 0.35

dot_df_ef <- dot_df_ef %>%
  mutate(
    bar_fill = PW_COLORS[pathway_label],
    ymin_bar = dot_y - bar_half_h,
    ymax_bar = dot_y + bar_half_h
  )

# Compute fraction of Y_SPAN occupied by pathway bars (for layout sizing)
dot_y_lo <- min(dot_df_ef$ymin_bar)
dot_y_hi <- max(dot_df_ef$ymax_bar)
dot_frac <- (dot_y_hi - dot_y_lo) / Y_SPAN

# Clip y-axis to just the pathway range (remove blank space below lowest pathway)
dot_y_pad <- pw_h * 0.8
dot_ylim  <- c(dot_y_lo - dot_y_pad, dot_y_hi + dot_y_pad)

pF_dot <- ggplot(dot_df_ef) +
  # Pathway-colored horizontal bars
  geom_rect(aes(xmin = 0, xmax = gene_ratio,
                ymin = ymin_bar, ymax = ymax_bar),
            fill = dot_df_ef$bar_fill, color = "grey30", linewidth = 0.2) +
  # Count label at end of bar
  geom_text(aes(x = gene_ratio + 0.005, y = dot_y,
                label = Count_1to1),
            hjust = 0, size = 2.5, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = 0)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  coord_cartesian(ylim = dot_ylim) +
  labs(x = "Gene Ratio", y = NULL) +
  THEME_PUB +
  theme(
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "black", linewidth = 0.6),
    legend.position    = "none",
    plot.margin        = margin(8, 5, 8, 0)
  )

# --- Vertical stacked keys: Contrast (left column) | Response (right column) ---
# Both keys are vertical lists, placed side by side, bottom-right under bar plot

ks <- 0.10  # vertical spacing between key items (tight)
# Both columns: items start at y = 3*ks down to 1*ks; titles at 4*ks
contrast_key_ef <- tibble(
  x     = 0,
  y     = c(3, 2, 1) * ks,
  label = c("Aging", "Training (Old)", "Training (Young)"),
  color = unname(CONTRAST_COLORS[c("Aging", "Training_Old", "Training_Young")]),
  shape = c(1, 16, 16)
)
pat_key_ef <- tibble(
  x     = 3.0,
  y     = c(3, 2, 1, 0) * ks,
  label = PATTERN_ORDER,
  color = unname(PATTERN_COLORS[PATTERN_ORDER])
)
title_y <- 4 * ks  # headers one step above top items

pF_key <- ggplot() +
  # Contrast column
  annotate("text", x = 0, y = title_y, label = "Contrast",
           hjust = 0, size = 2.8, fontface = "bold", color = "grey30") +
  geom_point(data = contrast_key_ef, aes(x = x, y = y),
             shape = contrast_key_ef$shape, size = 3.0,
             color = contrast_key_ef$color, stroke = 0.5) +
  geom_text(data = contrast_key_ef, aes(x = x + 0.3, y = y, label = label),
            hjust = 0, size = 2.5, fontface = "bold",
            color = contrast_key_ef$color) +
  # Response column
  annotate("text", x = 3.0, y = title_y, label = "Response",
           hjust = 0, size = 2.8, fontface = "bold", color = "grey30") +
  geom_point(data = pat_key_ef, aes(x = x, y = y),
             shape = 15, size = 3.0, color = pat_key_ef$color) +
  geom_text(data = pat_key_ef, aes(x = x + 0.3, y = y, label = label),
            hjust = 0, size = 2.5, fontface = "bold",
            color = pat_key_ef$color) +
  scale_x_continuous(limits = c(-0.3, 7)) +
  scale_y_continuous(limits = c(-0.3 * ks, title_y + 0.5 * ks)) +
  theme_void() +
  theme(plot.margin = margin(0, 5, 0, 0))

# --- Assemble with design: 4 panels on a 30-row x 100-col grid ---
# Key sits directly below bar panel (tight)
dot_grid_rows <- round(26 * dot_frac)
key_top <- dot_grid_rows + 1

ef_design <- c(
  patchwork::area(1,  1,  26, 28),              # dumbbell: full height, left 28%
  patchwork::area(1,  29, 26, 68),              # sankey:   full height, middle 40%
  patchwork::area(1,  69, dot_grid_rows, 100),  # bar:      top portion, right 32%
  patchwork::area(key_top, 69, 30, 100)         # key:      directly under bar
)

pF <- pF_dumbbell + pF_sankey + pF_dot + pF_key +
  plot_layout(design = ef_design) +
  plot_annotation(
    title    = "Interaction DEPs: Multi-Contrast Response & Pathway Enrichment",
    subtitle = sprintf("%d proteins with significant Age x Training interaction (pi < 0.05)",
                       length(all_int_genes)),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 7, color = "grey30", face = "italic")
    )
  )

ef_height <- max(200, nrow(links_1to1) * 7 + 60)

ggsave(file.path(RPT_DIR, "panel_F_interaction.pdf"), pF,
       width = 350, height = ef_height, units = "mm", device = pdf)
ggsave(file.path(RPT_DIR, "panel_F_interaction.png"), pF,
       width = 350, height = ef_height, units = "mm", dpi = 300)

message("  Panel F merged saved")

# ==============================================================================
# SUMMARY
# ==============================================================================
cat("\n", strrep("=", 61), "\nFigure 2 Panel Export Complete\n", strrep("=", 61), "\n")
cat("\nPDFs:\n")
for (f in c("panel_A_volcano", "panel_B_volcano", "panel_C_concordance",
            "panel_D_rrho2", "panel_E_nes_bubble", "panel_F_interaction"))
  cat(sprintf("  %s/%s.pdf\n", RPT_DIR, f))
cat("\nData (organized by panel):\n")
csv_map <- list(
  panel_A = "volcano_young.csv",
  panel_B = "volcano_old.csv",
  panel_C = "concordance.csv",
  panel_D = c("rrho2_summary.csv", "rrho2_matrix.csv"),
  panel_E = "nes_scatter.csv",
  panel_F = c("interaction_classification.csv", "interaction_dot.csv",
              "interaction_dot_long.csv", "sankey_dot.csv", "sankey_links.csv",
              "interaction_patterns.csv", "srplot_input.csv"),
  shared  = "fgsea_tstat_all_v2.csv"
)
for (pnl in names(csv_map))
  for (f in csv_map[[pnl]])
    cat(sprintf("  %s/%s/%s\n", DAT_DIR, pnl, f))
