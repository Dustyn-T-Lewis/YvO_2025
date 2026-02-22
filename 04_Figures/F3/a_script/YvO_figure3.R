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
  library(rrvgo)
  library(GOSemSim)
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

# --- Significance hierarchy colors (Panel B/D) ---
SIG_COLORS <- c(
  "Interaction"          = "#9B7FBF",
  "Sig Both"             = "#2E7D32",
  "Sig Aging only"       = "#4CAF50",
  "Sig Tr.Old only"      = "#5DA5DA",
  "NS"                   = "grey80"
)
SIG_SIZES  <- c(Interaction = 2.5, `Sig Both` = 2.0,
                `Sig Aging only` = 1.5, `Sig Tr.Old only` = 1.5, NS = 0.4)
SIG_ALPHAS <- c(Interaction = 0.90, `Sig Both` = 0.85,
                `Sig Aging only` = 0.70, `Sig Tr.Old only` = 0.70, NS = 0.20)

# --- Volcano contrast-specific coloring ---
VOLC_COLORS <- list(
  Aging        = c(Up = "#388E3C", Down = "#A5D6A7", NS = "grey80"),
  Training_Old = c(Up = "#2980B9", Down = "#85C1E9", NS = "grey80")
)

# --- Panel E reversal category colors ---
REVERSAL_CAT_COLORS <- c(
  "Fully Reversed"    = "#00897B",
  "Partially Reversed" = "#80CBC4",
  "Non-Reversed"       = "grey60",
  "Exacerbated"        = "#FF8F00"
)

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

classify_proteins <- function(pi_A, pi_B, pi_int,
                              label_A = "Young", label_B = "Old",
                              threshold = 0.05) {
  category <- case_when(
    pi_int < threshold                       ~ "Interaction",
    pi_A < threshold & pi_B < threshold      ~ "Sig Both",
    pi_A < threshold                         ~ paste0("Sig ", label_A, " only"),
    pi_B < threshold                         ~ paste0("Sig ", label_B, " only"),
    TRUE                                     ~ "NS"
  )
  factor(category, levels = c("Interaction", "Sig Both",
                               paste0("Sig ", label_A, " only"),
                               paste0("Sig ", label_B, " only"), "NS"))
}

quadrant_ora <- function(scatter_df, logFC_x_col, logFC_y_col,
                         pi_cols, threshold_strict = 0.05,
                         threshold_relax = 0.15, min_n = 20,
                         hallmark_t2g) {
  qdf <- scatter_df %>%
    mutate(
      quadrant = case_when(
        .data[[logFC_x_col]] > 0 & .data[[logFC_y_col]] > 0 ~ "Q1",
        .data[[logFC_x_col]] < 0 & .data[[logFC_y_col]] > 0 ~ "Q2",
        .data[[logFC_x_col]] < 0 & .data[[logFC_y_col]] < 0 ~ "Q3",
        TRUE ~ "Q4"
      ),
      relaxed_sig = Reduce(`|`, lapply(pi_cols, function(col) .data[[col]] < threshold_relax))
    )

  all_genes <- unique(scatter_df$gene)
  results <- list()

  for (q in c("Q1", "Q2", "Q3", "Q4")) {
    genes_q <- qdf %>% filter(quadrant == q, relaxed_sig) %>% pull(gene)

    if (length(genes_q) < min_n) {
      genes_q <- qdf %>%
        filter(quadrant == q,
               Reduce(`|`, lapply(pi_cols, function(col) .data[[col]] < 0.20))) %>%
        pull(gene)
    }

    if (length(genes_q) < min_n) {
      results[[q]] <- tibble(quadrant = q, Description = paste0("n = ", length(genes_q), "; underpowered"),
                             pvalue = NA_real_, p.adjust = NA_real_, n_genes = length(genes_q))
      next
    }

    ora <- tryCatch({
      enricher(gene = genes_q, TERM2GENE = hallmark_t2g,
               universe = all_genes, pvalueCutoff = 0.1, qvalueCutoff = 1)
    }, error = function(e) NULL)

    if (is.null(ora) || nrow(as.data.frame(ora)) == 0) {
      results[[q]] <- tibble(quadrant = q, Description = "No sig. Hallmark terms",
                             pvalue = NA_real_, p.adjust = NA_real_, n_genes = length(genes_q))
    } else {
      top3 <- as.data.frame(ora) %>%
        arrange(p.adjust) %>%
        slice_head(n = 3) %>%
        mutate(quadrant = q, Description = clean_pathway_name(Description),
               n_genes = length(genes_q))
      results[[q]] <- top3 %>% dplyr::select(quadrant, Description, pvalue, p.adjust, n_genes)
    }
  }

  bind_rows(results)
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

# --- Hallmark gene sets in TERM2GENE format for ORA ---
hallmark_t2g <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  as.data.frame()
message(sprintf("Loaded %d Hallmark gene-set mappings", nrow(hallmark_t2g)))

# === 8. PANEL A — Side-by-Side Volcanos (Aging | Training Old) ==============

message("Building Panel A: side-by-side volcanos (Aging | Training Old)...")

make_volcano <- function(ctr) {
  col_logFC  <- paste0("logFC_", ctr)
  col_pval   <- paste0("P.Value_", ctr)
  col_pi     <- paste0("pi_score_", ctr)

  volc_cols <- VOLC_COLORS[[ctr]]

  vdf <- dep_df %>%
    dplyr::select(gene,
           logFC    = all_of(col_logFC),
           P.Value  = all_of(col_pval),
           pi_score = all_of(col_pi)) %>%
    filter(!is.na(logFC), !is.na(P.Value)) %>%
    mutate(
      neg_log10p = -log10(P.Value),
      direction  = case_when(
        pi_score < 0.05 & logFC > 0 ~ "Up",
        pi_score < 0.05 & logFC < 0 ~ "Down",
        TRUE                         ~ "NS"
      )
    )

  n_up   <- sum(vdf$direction == "Up",   na.rm = TRUE)
  n_down <- sum(vdf$direction == "Down", na.rm = TRUE)

  dir_note_up <- ""
  dir_note_down <- ""
  if (n_up > 0 & n_down == 0) dir_note_up <- "\n(exclusively upregulated)"
  if (n_down > 0 & n_up == 0) dir_note_down <- "\n(exclusively downregulated)"

  top_genes <- vdf %>%
    filter(pi_score < 0.05) %>%
    arrange(pi_score) %>%
    slice_head(n = 6)

  pw_df <- fgsea_all %>%
    filter(contrast == ctr, database == "Hallmark", padj < 0.05)

  pw_up <- pw_df %>%
    filter(NES > 0) %>% arrange(desc(NES)) %>% slice_head(n = 4) %>%
    mutate(label = clean_pathway_name(pathway))

  pw_down <- pw_df %>%
    filter(NES < 0) %>% arrange(NES) %>% slice_head(n = 4) %>%
    mutate(label = clean_pathway_name(pathway))

  y_max_est <- max(vdf$neg_log10p, na.rm = TRUE)
  x_max_est <- max(abs(vdf$logFC), na.rm = TRUE)

  strip_title <- if (ctr == "Aging") "Aging" else paste0("Training (", str_extract(ctr, "Young|Old"), ")")

  p <- ggplot(vdf, aes(x = logFC, y = neg_log10p)) +
    geom_point(aes(color = direction), size = 0.5, alpha = 0.4) +
    scale_color_manual(values = volc_cols) +
    geom_text_repel(
      data = top_genes, aes(label = gene),
      size = KEY_TEXT, max.overlaps = 15,
      segment.size = 0.2, fontface = "italic",
      min.segment.length = 0
    ) +
    annotate("text",
             x = x_max_est * 0.95, y = y_max_est * 0.97,
             label = paste0("pi < 0.05: ", n_up, " up", dir_note_up),
             hjust = 1, vjust = 1, size = KEY_TITLE, color = volc_cols["Up"]) +
    annotate("text",
             x = -x_max_est * 0.95, y = y_max_est * 0.97,
             label = paste0("pi < 0.05: ", n_down, " down", dir_note_down),
             hjust = 0, vjust = 1, size = KEY_TITLE, color = volc_cols["Down"]) +
    labs(
      title    = strip_title,
      subtitle = "Colored points: pi < 0.05",
      x = expression(log[2]~fold~change),
      y = expression(-log[10]~italic(P))
    ) +
    THEME_PUB +
    theme(legend.position = "none")

  if (nrow(pw_up) > 0) {
    p <- p + annotate("label",
                       x = x_max_est * 0.95, y = y_max_est * 0.05,
                       label = paste(pw_up$label, collapse = "\n"),
                       hjust = 1, vjust = 0, size = 1.8,
                       fill = alpha("white", 0.85),
                       label.padding = unit(1.5, "pt"),
                       color = volc_cols["Up"])
  }
  if (nrow(pw_down) > 0) {
    p <- p + annotate("label",
                       x = -x_max_est * 0.95, y = y_max_est * 0.05,
                       label = paste(pw_down$label, collapse = "\n"),
                       hjust = 0, vjust = 0, size = 1.8,
                       fill = alpha("white", 0.85),
                       label.padding = unit(1.5, "pt"),
                       color = volc_cols["Down"])
  }

  return(p)
}

# --- Panel A assembly ---
volc_xlim <- max(abs(c(dep_df$logFC_Aging, dep_df$logFC_Training_Old)),
                 na.rm = TRUE) * 1.05
volc_ylim <- max(-log10(c(dep_df$P.Value_Aging, dep_df$P.Value_Training_Old)),
                 na.rm = TRUE)
volc_ylim <- min(volc_ylim, 15)

pA_left  <- make_volcano("Aging") +
  coord_cartesian(xlim = c(-volc_xlim, volc_xlim), ylim = c(0, volc_ylim))
pA_right <- make_volcano("Training_Old") +
  coord_cartesian(xlim = c(-volc_xlim, volc_xlim), ylim = c(0, volc_ylim)) +
  theme(axis.title.y = element_blank())

pA <- (pA_left | pA_right) + plot_layout(widths = c(1, 1))

ggsave(file.path(RPT_DIR, "test_panelA.pdf"), pA,
       width = 250, height = 120, units = "mm")
message("Panel A test saved")

# === 9. PANEL B — Reversal Scatter with 5-Category Classification ===========

message("Building Panel B: reversal scatter...")

# --- 1. Prepare data with 5-category hierarchy ---
scatter_df <- dep_df %>%
  transmute(gene,
            logFC_A  = logFC_Aging,
            logFC_TO = logFC_Training_Old,
            pi_A     = pi_score_Aging,
            pi_TO    = pi_score_Training_Old,
            pi_int   = pi_score_Interaction) %>%
  filter(!is.na(logFC_A), !is.na(logFC_TO)) %>%
  mutate(
    sig_cat = classify_proteins(pi_A, pi_TO, pi_int, "Aging", "Tr.Old"),
    reversed = (logFC_A > 0 & logFC_TO < 0) | (logFC_A < 0 & logFC_TO > 0),
    quadrant = case_when(
      logFC_A > 0 & logFC_TO > 0 ~ "Q1",
      logFC_A > 0 & logFC_TO < 0 ~ "Q2",
      logFC_A < 0 & logFC_TO > 0 ~ "Q4",
      TRUE ~ "Q3"
    )
  )

# --- 2. Compute correlation stats with CIs ---
cor_test_r   <- cor.test(scatter_df$logFC_A, scatter_df$logFC_TO, method = "pearson")
cor_test_rho <- cor.test(scatter_df$logFC_A, scatter_df$logFC_TO, method = "spearman")
cor_r    <- cor_test_r$estimate
cor_rho  <- cor_test_rho$estimate
cor_r_ci <- cor_test_r$conf.int

# Expected null correlation due to shared Old_Pre samples
r_null <- -0.50
message(sprintf("Panel B: observed r = %.3f [%.3f, %.3f], expected null r ~ %.2f (shared Old_Pre)",
                cor_r, cor_r_ci[1], cor_r_ci[2], r_null))

# Reversal ratio (among proteins with |logFC| > 0.2 in >= 1 contrast)
reversal_set <- scatter_df %>%
  filter(abs(logFC_A) > 0.2 | abs(logFC_TO) > 0.2)
reversal_ratio <- mean(reversal_set$reversed) * 100

# --- 3. Quadrant counts ---
q_counts <- scatter_df %>% count(quadrant) %>% deframe()

# --- 4. Per-quadrant Hallmark ORA ---
qora_results <- quadrant_ora(
  scatter_df, "logFC_A", "logFC_TO",
  pi_cols = c("pi_A", "pi_TO", "pi_int"),
  hallmark_t2g = hallmark_t2g
)

# --- 5. Axis range ---
axis_lim <- max(abs(c(scatter_df$logFC_A, scatter_df$logFC_TO)), na.rm = TRUE) * 1.15

# --- 6. Build scatter plot ---
scatter_ordered <- scatter_df %>%
  mutate(plot_order = as.integer(sig_cat)) %>%
  arrange(desc(plot_order))

pB_base <- ggplot(scatter_ordered, aes(x = logFC_A, y = logFC_TO)) +
  # Reversed quadrants (Q2, Q4) = teal
  annotate("rect", xmin = 0, xmax = axis_lim, ymin = -axis_lim, ymax = 0,
           fill = "#00897B", alpha = 0.06) +
  annotate("rect", xmin = -axis_lim, xmax = 0, ymin = 0, ymax = axis_lim,
           fill = "#00897B", alpha = 0.06) +
  # Exacerbated quadrants (Q1, Q3) = amber
  annotate("rect", xmin = 0, xmax = axis_lim, ymin = 0, ymax = axis_lim,
           fill = "#FF8F00", alpha = 0.06) +
  annotate("rect", xmin = -axis_lim, xmax = 0, ymin = -axis_lim, ymax = 0,
           fill = "#FF8F00", alpha = 0.06) +
  # Reference lines
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  # Points with 5-category encoding
  geom_point(aes(color = sig_cat, size = sig_cat, alpha = sig_cat)) +
  scale_color_manual(values = SIG_COLORS, name = "Significance") +
  scale_size_manual(values = SIG_SIZES, name = "Significance") +
  scale_alpha_manual(values = SIG_ALPHAS, name = "Significance") +
  # Stats annotation
  annotate("text", x = -axis_lim * 0.95, y = axis_lim * 0.95,
           label = sprintf("r = %.2f [%.2f, %.2f]\nrho = %.2f\nr(null) ~ %.2f\nReversal: %.0f%%",
                           cor_r, cor_r_ci[1], cor_r_ci[2], cor_rho, r_null, reversal_ratio),
           hjust = 0, vjust = 1, size = KEY_TITLE, fontface = "bold") +
  # Quadrant counts
  annotate("text", x = axis_lim, y = axis_lim,
           label = paste("Exac. n =", q_counts["Q1"]),
           hjust = 1.1, vjust = 1.5, size = 2.0, color = "#FF8F00") +
  annotate("text", x = axis_lim, y = -axis_lim,
           label = paste("Rev. n =", q_counts["Q2"]),
           hjust = 1.1, vjust = -0.5, size = 2.0, color = "#00897B") +
  annotate("text", x = -axis_lim, y = -axis_lim,
           label = paste("Exac. n =", q_counts["Q3"]),
           hjust = -0.1, vjust = -0.5, size = 2.0, color = "#FF8F00") +
  annotate("text", x = -axis_lim, y = axis_lim,
           label = paste("Rev. n =", q_counts["Q4"]),
           hjust = -0.1, vjust = 1.5, size = 2.0, color = "#00897B") +
  labs(
    title = "Protein Reversal Map",
    subtitle = sprintf("Dashed = perfect reversal | r(null) ~ %.2f (shared Old_Pre)", r_null),
    x = expression(log[2]*FC ~ "(Aging)"),
    y = expression(log[2]*FC ~ "(Training Old)")
  ) +
  coord_cartesian(xlim = c(-axis_lim, axis_lim),
                  ylim = c(-axis_lim, axis_lim)) +
  THEME_PUB +
  guides(color = guide_legend(override.aes = list(size = c(2.5, 2.0, 1.5, 1.5, 0.8),
                                                    alpha = c(0.9, 0.85, 0.7, 0.7, 0.3))),
         size = "none", alpha = "none") +
  theme(legend.position = "bottom",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 5))

# --- 7. Per-quadrant ORA text annotations ---
q_positions <- list(
  Q1 = c(x = axis_lim * 0.95, y = axis_lim * 0.60, hjust = 1, vjust = 1),
  Q2 = c(x = axis_lim * 0.95, y = -axis_lim * 0.60, hjust = 1, vjust = 0),
  Q3 = c(x = -axis_lim * 0.95, y = -axis_lim * 0.60, hjust = 0, vjust = 0),
  Q4 = c(x = -axis_lim * 0.95, y = axis_lim * 0.60, hjust = 0, vjust = 1)
)

for (q in names(q_positions)) {
  q_terms <- qora_results %>% filter(quadrant == q)
  if (nrow(q_terms) > 0) {
    label_text <- paste(q_terms$Description, collapse = "\n")
    pos <- q_positions[[q]]
    pB_base <- pB_base + annotate("label",
      x = as.numeric(pos["x"]), y = as.numeric(pos["y"]),
      label = label_text,
      hjust = as.numeric(pos["hjust"]), vjust = as.numeric(pos["vjust"]),
      size = 1.6, fill = alpha("white", 0.85),
      label.padding = unit(1.5, "pt"), color = "grey30")
  }
}

# --- 8. Protein labels per category ---
label_df <- scatter_df %>%
  filter(sig_cat != "NS") %>%
  group_by(sig_cat) %>%
  arrange(desc(abs(logFC_A) + abs(logFC_TO))) %>%
  slice_head(n = 6) %>%
  ungroup()

pB_base <- pB_base +
  geom_text_repel(
    data = label_df, aes(label = gene),
    size = KEY_TEXT, max.overlaps = 20,
    segment.size = 0.2, fontface = "italic",
    min.segment.length = 0
  )

# --- 9. Export and test save ---
write_csv(scatter_df, file.path(DAT_DIR, "fig3_reversal_scatter.csv"))
message(sprintf("Panel B: r = %.3f, rho = %.3f, reversal = %.1f%%",
                cor_r, cor_rho, reversal_ratio))

ggsave(file.path(RPT_DIR, "test_panelB.pdf"), pB_base,
       width = 170, height = 160, units = "mm")
message("Panel B test saved")

# === 10. PANEL C — RRHO2 Reversal Map (Annotated) ===========================

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

shared_genes <- intersect(rrho_list1$gene, rrho_list2$gene)
rrho_list1 <- rrho_list1 %>% dplyr::filter(gene %in% shared_genes)
rrho_list2 <- rrho_list2 %>% dplyr::filter(gene %in% shared_genes)

# --- 2. Run RRHO2 ---
rrho_obj <- RRHO2_initialize(
  list1 = rrho_list1, list2 = rrho_list2,
  labels = c("Aging", "Training (Old)"),
  log10.ind = TRUE, method = "hyper", boundary = 0.04
)

# --- 3. Extract max -log10(p) per quadrant ---
hmat <- rrho_obj$hypermat
nr <- nrow(hmat); nc <- ncol(hmat)
mid_r <- floor(nr / 2); mid_c <- floor(nc / 2)
max_UU <- max(hmat[1:mid_r, (mid_c+1):nc], na.rm = TRUE)
max_DD <- max(hmat[(mid_r+1):nr, 1:mid_c], na.rm = TRUE)
max_UD <- max(hmat[1:mid_r, 1:mid_c], na.rm = TRUE)
max_DU <- max(hmat[(mid_r+1):nr, (mid_c+1):nc], na.rm = TRUE)

# --- 4. Export ---
write.csv(hmat, file.path(DAT_DIR, "fig3_rrho2_matrix.csv"))

# --- 5. Render RRHO2 to PNG ---
tmp_rrho <- tempfile(fileext = ".png")
png(tmp_rrho, width = 1400, height = 1200, res = 300)
par(mar = c(2, 2, 2, 1))
RRHO2_heatmap(rrho_obj)
dev.off()
rrho_img <- png::readPNG(tmp_rrho)

# --- 6. Build annotated ggplot wrapper ---
pC_gg <- ggplot() +
  annotation_raster(rrho_img, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "white", linewidth = 0.5) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "white", linewidth = 0.5) +
  # Quadrant labels — F3 context: reversal vs exacerbation
  annotate("label", x = 0.75, y = 0.75, label = paste0("Reversal up-dn\nmax = ", round(max_UD, 1)),
           fill = alpha("white", 0.7), size = 2.0, fontface = "bold") +
  annotate("label", x = 0.25, y = 0.25, label = paste0("Reversal dn-up\nmax = ", round(max_DU, 1)),
           fill = alpha("white", 0.7), size = 2.0, fontface = "bold") +
  annotate("label", x = 0.25, y = 0.75, label = paste0("Exacerbation up-up\nmax = ", round(max_UU, 1)),
           fill = alpha("white", 0.7), size = 2.0, fontface = "bold") +
  annotate("label", x = 0.75, y = 0.25, label = paste0("Exacerbation dn-dn\nmax = ", round(max_DD, 1)),
           fill = alpha("white", 0.7), size = 2.0, fontface = "bold") +
  annotate("text", x = 0.95, y = 0.02, label = "Most upregulated ->",
           hjust = 1, size = 1.8, color = "grey30") +
  annotate("text", x = 0.05, y = 0.02, label = "<- Most downregulated",
           hjust = 0, size = 1.8, color = "grey30") +
  labs(title = "RRHO2 Reversal Map",
       subtitle = "-log10(p) hypergeometric overlap") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5),
        plot.subtitle = element_text(size = 6.5, color = "grey30", hjust = 0.5, face = "italic"))

# --- 7. Test save ---
pdf(file.path(RPT_DIR, "test_panelC.pdf"), width = 7, height = 5)
par(mar = c(2, 2, 2, 1))
RRHO2_heatmap(rrho_obj)
dev.off()

message(sprintf("RRHO2 quadrant max -log10(p): UU=%.1f DD=%.1f UD=%.1f DU=%.1f",
                max_UU, max_DD, max_UD, max_DU))
message("Panel C test saved")

# === 11. PANEL D — mitch 2D Pathway Enrichment (Reversal) ===================

message("Building Panel D: mitch 2D pathway enrichment...")

# --- 1. Build t-stat matrix ---
tstat_mat <- dep_df %>%
  dplyr::select(gene, Aging = t_Aging, Training_Old = t_Training_Old) %>%
  dplyr::filter(!is.na(Aging), !is.na(Training_Old)) %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# --- 2. Fetch gene sets ---
hallmark_gs <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol)
gobp_gs <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP") %>%
  dplyr::select(gs_name, gene_symbol)

gs_df <- bind_rows(hallmark_gs, gobp_gs)
genesets <- split(gs_df$gene_symbol, gs_df$gs_name)

# --- 3. Run mitch ---
mitch_res <- mitch_calc(x = tstat_mat, genesets = genesets,
                        priority = "effect", cores = 1, resrows = Inf)

# --- 4. Process results with 4-category classification ---
mitch_df <- mitch_res$enrichment_result %>%
  mutate(
    padj_A   = p.adjust(p.Aging, method = "BH"),
    padj_TO  = p.adjust(p.Training_Old, method = "BH"),
    sig_A    = padj_A < 0.05,
    sig_TO   = padj_TO < 0.05,
    # Reversal-like: joint test significant AND in reversal quadrant
    reversal_quadrant = (s.Aging > 0 & s.Training_Old < 0) |
                        (s.Aging < 0 & s.Training_Old > 0),
    sig_joint = p.adjustMANOVA < 0.05 & reversal_quadrant,
    sig_cat  = case_when(
      sig_joint             ~ "Interaction",
      sig_A & sig_TO        ~ "Sig Both",
      sig_A                 ~ "Sig Aging only",
      sig_TO                ~ "Sig Tr.Old only",
      TRUE                  ~ "NS"
    ),
    sig_cat = factor(sig_cat, levels = c("Interaction", "Sig Both",
                                          "Sig Aging only", "Sig Tr.Old only", "NS")),
    pathway_clean = clean_pathway_name(set)
  )

# --- 5. Identify pathways to label ---
label_keywords <- c("ribosom", "oxphos", "oxidative phosph", "mtorc1",
                     "extracellular matrix", "myogenes", "glycoly",
                     "proteasome", "translation", "unfolded protein",
                     "mitotic spindle", "muscle", "respiratory chain",
                     "fatty acid", "inflammatory")

label_df <- bind_rows(
  mitch_df %>% filter(sig_cat != "NS") %>%
    filter(str_detect(tolower(set), paste(label_keywords, collapse = "|"))),
  mitch_df %>% filter(sig_cat != "NS", s.Aging > 0, s.Training_Old < 0) %>%
    slice_max(abs(s.Aging) + abs(s.Training_Old), n = 3),
  mitch_df %>% filter(sig_cat != "NS", s.Aging < 0, s.Training_Old > 0) %>%
    slice_max(abs(s.Aging) + abs(s.Training_Old), n = 3),
  mitch_df %>% filter(sig_cat != "NS", s.Aging > 0, s.Training_Old > 0) %>%
    slice_max(abs(s.Aging) + abs(s.Training_Old), n = 3),
  mitch_df %>% filter(sig_cat != "NS", s.Aging < 0, s.Training_Old < 0) %>%
    slice_max(abs(s.Aging) + abs(s.Training_Old), n = 3)
) %>%
  distinct(set, .keep_all = TRUE) %>%
  slice_head(n = 15)

# --- 6. Correlations ---
pw_cor  <- cor(mitch_df$s.Aging, mitch_df$s.Training_Old)
pro_cor <- cor_r

# --- 7. Axis range ---
pw_lim <- max(abs(c(mitch_df$s.Aging, mitch_df$s.Training_Old)), na.rm = TRUE) * 1.1

# --- 8. Build scatter plot ---
pD <- ggplot(mitch_df %>% arrange(desc(as.integer(sig_cat))),
             aes(x = s.Aging, y = s.Training_Old)) +
  # Reversal quadrants (Q2, Q4) = teal
  annotate("rect", xmin = 0, xmax = pw_lim, ymin = -pw_lim, ymax = 0,
           fill = "#00897B", alpha = 0.04) +
  annotate("rect", xmin = -pw_lim, xmax = 0, ymin = 0, ymax = pw_lim,
           fill = "#00897B", alpha = 0.04) +
  # Exacerbation quadrants (Q1, Q3) = amber
  annotate("rect", xmin = 0, xmax = pw_lim, ymin = 0, ymax = pw_lim,
           fill = "#FF8F00", alpha = 0.04) +
  annotate("rect", xmin = -pw_lim, xmax = 0, ymin = -pw_lim, ymax = 0,
           fill = "#FF8F00", alpha = 0.04) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40", linewidth = 0.3) +
  geom_point(aes(color = sig_cat, size = setSize, alpha = sig_cat)) +
  scale_color_manual(values = SIG_COLORS, name = "Significance") +
  scale_alpha_manual(values = SIG_ALPHAS, guide = "none") +
  scale_size_continuous(range = c(0.3, 3.0), name = "Protein set size",
                        guide = guide_legend(override.aes = list(alpha = 0.8))) +
  geom_text_repel(data = label_df, aes(label = pathway_clean),
                  size = KEY_TEXT, max.overlaps = 25,
                  segment.size = 0.2, fontface = "italic",
                  min.segment.length = 0) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = sprintf("Pathway r = %.2f\nProtein r = %.2f", pw_cor, pro_cor),
           size = KEY_TITLE, fontface = "bold") +
  labs(x = "Enrichment Score (Aging)",
       y = "Enrichment Score Tr. (Old)",
       title = "Pathway-Level Reversal Map",
       subtitle = "Hallmark + GO:BP via mitch") +
  coord_cartesian(xlim = c(-pw_lim, pw_lim), ylim = c(-pw_lim, pw_lim)) +
  THEME_PUB +
  theme(legend.position = "bottom",
        legend.key.size = unit(2, "mm"),
        legend.text = element_text(size = 5))

# --- 9. Export ---
write_csv(mitch_df, file.path(DAT_DIR, "fig3_mitch_2d_results.csv"))
message(sprintf("mitch: %d pathways, Interaction=%d, Both=%d, Aging=%d, Tr.Old=%d, r = %.3f",
                nrow(mitch_df),
                sum(mitch_df$sig_cat == "Interaction"),
                sum(mitch_df$sig_cat == "Sig Both"),
                sum(mitch_df$sig_cat == "Sig Aging only"),
                sum(mitch_df$sig_cat == "Sig Tr.Old only"),
                pw_cor))

ggsave(file.path(RPT_DIR, "test_panelD.pdf"), pD,
       width = 170, height = 160, units = "mm")
message("Panel D test saved")

# === 12. PANEL E — Reversal Classification (Diverging Lollipop) =============

message("Building Panel E: reversal classification (diverging lollipop)...")

# --- 1. Classify Aging DEPs by Training_Old behavior ---
rev_class <- dep_df %>%
  filter(pi_score_Aging < 0.05) %>%
  dplyr::select(gene, logFC_A = logFC_Aging, logFC_TO = logFC_Training_Old,
         pi_TO = pi_score_Training_Old) %>%
  filter(!is.na(logFC_A), !is.na(logFC_TO)) %>%
  mutate(
    opposite_dir = sign(logFC_A) != sign(logFC_TO),
    category = case_when(
      opposite_dir & pi_TO < 0.05  ~ "Fully reversed",
      opposite_dir & pi_TO >= 0.05 ~ "Partially reversed",
      !opposite_dir & pi_TO < 0.05  ~ "Exacerbated",
      !opposite_dir & pi_TO >= 0.05 ~ "Non-reversed"
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
