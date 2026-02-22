# ============================================================================
# Figure 2: Concordance & Discordance of Young vs Old Training Responses
# ============================================================================
# Study: YvO — DIA-MS Proteomics of Skeletal Muscle
# Design: 2x2 mixed (Age x Time) with repeated measures
#
# Central question: To what extent do young and old adults remodel the same
# proteome in the same direction — and where do they diverge?
#
# Panels:
#   A — Side-by-side volcanos (Training Young | Training Old)
#   B — Concordance scatter (logFC x logFC)
#   C — RRHO2 threshold-free concordance map (Cahill et al. 2018)
#   D — mitch 2D pathway enrichment (Kaspi & Ziemann 2020)
#   E — Interaction DEP classification (diverging lollipop)
#   F — Pathway enrichment of concordant vs discordant sets (Yu et al. 2012)
#
# References:
#   Cahill et al. 2018, Scientific Reports 8:9588 (RRHO2)
#   Kaspi & Ziemann 2020, BMC Genomics 21:447 (mitch)
#   Yu et al. 2012, OMICS 16:284-287 (clusterProfiler)
#   Xiao et al. 2014 (pi-score methodology)
# ============================================================================

# ═══ 1. PACKAGES ═════════════════════════════════════════════════════════════

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

# ═══ 2. SEED ═════════════════════════════════════════════════════════════════

set.seed(42)

# ═══ 3. PATH RESOLUTION ═════════════════════════════════════════════════════

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

# ═══ 4. DIRECTORY CREATION ══════════════════════════════════════════════════

dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUP_DIR, recursive = TRUE, showWarnings = FALSE)

# ═══ 5. CANONICAL CONSTANTS ═════════════════════════════════════════════════

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
  "Sig Young only"       = "#E05A4E",
  "Sig Old only"         = "#5DA5DA",
  "NS"                   = "grey80"
)
SIG_SIZES  <- c(Interaction = 2.5, `Sig Both` = 2.0,
                `Sig Young only` = 1.5, `Sig Old only` = 1.5, NS = 0.4)
SIG_ALPHAS <- c(Interaction = 0.90, `Sig Both` = 0.85,
                `Sig Young only` = 0.70, `Sig Old only` = 0.70, NS = 0.20)

# --- Volcano contrast-specific coloring ---
VOLC_COLORS <- list(
  Training_Young = c(Up = "#C0392B", Down = "#E8A09A", NS = "grey80"),
  Training_Old   = c(Up = "#2980B9", Down = "#85C1E9", NS = "grey80")
)

# --- Panel E category colors ---
INTERACTION_CAT_COLORS <- c(
  "Attenuated"          = "#81C784",
  "Opposite Direction"  = "#C62828",
  "Old-Specific"        = "#5DA5DA",
  "Young-Specific"      = "#CE93D8"
)

# ═══ 6. HELPER FUNCTIONS ════════════════════════════════════════════════════

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
  # Define quadrants based on logFC signs
  qdf <- scatter_df %>%
    mutate(
      quadrant = case_when(
        .data[[logFC_x_col]] > 0 & .data[[logFC_y_col]] > 0 ~ "Q1",
        .data[[logFC_x_col]] < 0 & .data[[logFC_y_col]] > 0 ~ "Q2",
        .data[[logFC_x_col]] < 0 & .data[[logFC_y_col]] < 0 ~ "Q3",
        TRUE ~ "Q4"
      ),
      # Relaxed significance: any pi_score < threshold_relax
      relaxed_sig = if_any(all_of(pi_cols), ~ . < threshold_relax)
    )

  # For each quadrant, run Hallmark ORA on relaxed-sig proteins
  all_genes <- unique(scatter_df$gene)
  results <- list()

  for (q in c("Q1", "Q2", "Q3", "Q4")) {
    genes_q <- qdf %>% filter(quadrant == q, relaxed_sig) %>% pull(gene)

    if (length(genes_q) < min_n) {
      # Try even more relaxed threshold
      genes_q <- qdf %>%
        filter(quadrant == q,
               if_any(all_of(pi_cols), ~ . < 0.20)) %>%
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

# ═══ 7. DATA LOADING ════════════════════════════════════════════════════════

message("Loading data...")
dep_df <- read_csv(file.path(BASE_DIR, "03_DEP", "c_data", "combined_results.csv"),
                   show_col_types = FALSE)
fgsea_all <- read_csv(file.path(DAT_DIR, "fgsea_tstat_all_v2.csv"),
                      show_col_types = FALSE)

# Validation
stopifnot(nrow(dep_df) > 2000)
stopifnot(all(c("gene", "logFC_Training_Young", "t_Training_Young",
                "sig_pi_Interaction") %in% names(dep_df)))
message(sprintf("Loaded %d proteins, %d fGSEA results", nrow(dep_df), nrow(fgsea_all)))

# --- Hallmark gene sets in TERM2GENE format for ORA ---
hallmark_t2g <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  as.data.frame()
message(sprintf("Loaded %d Hallmark gene-set mappings", nrow(hallmark_t2g)))

# ═══ 8. PANEL A — Side-by-Side Volcanos with Inset Hallmark Text ═══════════

message("Building Panel A: side-by-side volcanos...")

make_volcano <- function(ctr) {
  # --- 1. Extract & rename columns for this contrast ---
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

  # --- 2. Summary counts ---
  n_up   <- sum(vdf$direction == "Up",   na.rm = TRUE)
  n_down <- sum(vdf$direction == "Down", na.rm = TRUE)

  # Direction note (e.g., "(exclusively upregulated)")
  dir_note_up <- ""
  dir_note_down <- ""
  if (n_up > 0 & n_down == 0) dir_note_up <- "\n(exclusively upregulated)"
  if (n_down > 0 & n_up == 0) dir_note_down <- "\n(exclusively downregulated)"

  # --- 3. Top 6 DEPs by |pi_score| among significant ---
  top_genes <- vdf %>%
    filter(pi_score < 0.05) %>%
    arrange(pi_score) %>%
    slice_head(n = 6)

  # --- 4. Pathway inset: Hallmark, padj < 0.05 ---
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

  # Strip title
  strip_title <- if (ctr == "Aging") "Aging" else paste0("Training (", str_extract(ctr, "Young|Old"), ")")

  # --- 5. Build ggplot ---
  p <- ggplot(vdf, aes(x = logFC, y = neg_log10p)) +
    # Points — contrast-specific coloring
    geom_point(aes(color = direction), size = 0.5, alpha = 0.4) +
    scale_color_manual(values = volc_cols) +
    # Gene labels
    geom_text_repel(
      data = top_genes, aes(label = gene),
      size = KEY_TEXT, max.overlaps = 15,
      segment.size = 0.2, fontface = "italic",
      min.segment.length = 0
    ) +
    # DEP count — upper-right (Up)
    annotate("text",
             x = x_max_est * 0.95, y = y_max_est * 0.97,
             label = paste0("pi < 0.05: ", n_up, " up", dir_note_up),
             hjust = 1, vjust = 1, size = KEY_TITLE, color = volc_cols["Up"]) +
    # DEP count — upper-left (Down)
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

  # --- 6. Pathway insets — bottom corners ---
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

# --- Panel A assembly: shared axis limits ---
volc_xlim <- max(abs(c(dep_df$logFC_Training_Young, dep_df$logFC_Training_Old)),
                 na.rm = TRUE) * 1.05
volc_ylim <- max(-log10(c(dep_df$P.Value_Training_Young, dep_df$P.Value_Training_Old)),
                 na.rm = TRUE)
volc_ylim <- min(volc_ylim, 15)   # cap at 15

pA_left  <- make_volcano("Training_Young") +
  coord_cartesian(xlim = c(-volc_xlim, volc_xlim), ylim = c(0, volc_ylim))
pA_right <- make_volcano("Training_Old") +
  coord_cartesian(xlim = c(-volc_xlim, volc_xlim), ylim = c(0, volc_ylim)) +
  theme(axis.title.y = element_blank())

pA <- (pA_left | pA_right) + plot_layout(widths = c(1, 1))

# --- Panel A test save ---
ggsave(file.path(RPT_DIR, "test_panelA.pdf"), pA,
       width = 250, height = 120, units = "mm")
message("Panel A test saved")

# ═══ 9. PANEL B — Concordance Scatter with 5-Category Classification ════════

message("Building Panel B: concordance scatter...")

# --- 1. Prepare data with 5-category hierarchy ---
scatter_df <- dep_df %>%
  transmute(gene,
            logFC_Y  = logFC_Training_Young,
            logFC_O  = logFC_Training_Old,
            pi_Y     = pi_score_Training_Young,
            pi_O     = pi_score_Training_Old,
            pi_int   = pi_score_Interaction) %>%
  filter(!is.na(logFC_Y), !is.na(logFC_O)) %>%
  mutate(
    sig_cat = classify_proteins(pi_Y, pi_O, pi_int, "Young", "Old"),
    quadrant = case_when(
      logFC_Y > 0 & logFC_O > 0 ~ "Q1",
      logFC_Y < 0 & logFC_O > 0 ~ "Q2",
      logFC_Y < 0 & logFC_O < 0 ~ "Q3",
      TRUE ~ "Q4"
    ),
    concordant = quadrant %in% c("Q1", "Q3")
  )

# --- 2. Compute correlation stats with CIs ---
cor_test_r   <- cor.test(scatter_df$logFC_Y, scatter_df$logFC_O, method = "pearson")
cor_test_rho <- cor.test(scatter_df$logFC_Y, scatter_df$logFC_O, method = "spearman")
cor_r    <- cor_test_r$estimate
cor_rho  <- cor_test_rho$estimate
cor_r_ci <- cor_test_r$conf.int

# Sign concordance ratio (among proteins with |logFC| > 0.2 in >= 1 contrast)
concordance_set <- scatter_df %>%
  filter(abs(logFC_Y) > 0.2 | abs(logFC_O) > 0.2)
sign_concordance <- mean(sign(concordance_set$logFC_Y) == sign(concordance_set$logFC_O)) * 100

# --- 3. Quadrant counts ---
q_counts <- scatter_df %>% count(quadrant) %>% deframe()

# --- 4. Per-quadrant Hallmark ORA ---
qora_results <- quadrant_ora(
  scatter_df, "logFC_Y", "logFC_O",
  pi_cols = c("pi_Y", "pi_O", "pi_int"),
  hallmark_t2g = hallmark_t2g
)

# --- 5. Axis range ---
axis_lim <- max(abs(c(scatter_df$logFC_Y, scatter_df$logFC_O)), na.rm = TRUE) * 1.15

# --- 6. Build scatter plot ---
# Order data so NS is plotted first, Interaction on top
scatter_ordered <- scatter_df %>%
  mutate(plot_order = as.integer(sig_cat)) %>%
  arrange(desc(plot_order))

pB_base <- ggplot(scatter_ordered, aes(x = logFC_Y, y = logFC_O)) +
  # Quadrant background shading
  annotate("rect", xmin = 0, xmax = axis_lim, ymin = 0, ymax = axis_lim,
           fill = "#00897B", alpha = 0.06) +
  annotate("rect", xmin = -axis_lim, xmax = 0, ymin = -axis_lim, ymax = 0,
           fill = "#00897B", alpha = 0.06) +
  annotate("rect", xmin = 0, xmax = axis_lim, ymin = -axis_lim, ymax = 0,
           fill = "#FF8F00", alpha = 0.06) +
  annotate("rect", xmin = -axis_lim, xmax = 0, ymin = 0, ymax = axis_lim,
           fill = "#FF8F00", alpha = 0.06) +
  # Reference lines
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  # Points with 5-category encoding
  geom_point(aes(color = sig_cat, size = sig_cat, alpha = sig_cat)) +
  scale_color_manual(values = SIG_COLORS, name = "Significance") +
  scale_size_manual(values = SIG_SIZES, name = "Significance") +
  scale_alpha_manual(values = SIG_ALPHAS, name = "Significance") +
  # Stats annotation
  annotate("text", x = -axis_lim * 0.95, y = axis_lim * 0.95,
           label = sprintf("r = %.2f [%.2f, %.2f]\nrho = %.2f\nConcordance: %.0f%%",
                           cor_r, cor_r_ci[1], cor_r_ci[2], cor_rho, sign_concordance),
           hjust = 0, vjust = 1, size = KEY_TITLE, fontface = "bold") +
  # Quadrant counts
  annotate("text", x = axis_lim, y = axis_lim,
           label = paste("n =", q_counts["Q1"]),
           hjust = 1.1, vjust = 1.5, size = 2.0, color = "grey40") +
  annotate("text", x = -axis_lim, y = axis_lim,
           label = paste("n =", q_counts["Q2"]),
           hjust = -0.1, vjust = 1.5, size = 2.0, color = "grey40") +
  annotate("text", x = -axis_lim, y = -axis_lim,
           label = paste("n =", q_counts["Q3"]),
           hjust = -0.1, vjust = -0.5, size = 2.0, color = "grey40") +
  annotate("text", x = axis_lim, y = -axis_lim,
           label = paste("n =", q_counts["Q4"]),
           hjust = 1.1, vjust = -0.5, size = 2.0, color = "grey40") +
  labs(
    title = "Protein Concordance",
    x = expression(log[2]*FC ~ "(Training Young)"),
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
# Position ORA labels in each quadrant corner
q_positions <- list(
  Q1 = c(x = axis_lim * 0.95, y = axis_lim * 0.60, hjust = 1, vjust = 1),
  Q2 = c(x = -axis_lim * 0.95, y = axis_lim * 0.60, hjust = 0, vjust = 1),
  Q3 = c(x = -axis_lim * 0.95, y = -axis_lim * 0.60, hjust = 0, vjust = 0),
  Q4 = c(x = axis_lim * 0.95, y = -axis_lim * 0.60, hjust = 1, vjust = 0)
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
  arrange(desc(abs(logFC_Y) + abs(logFC_O))) %>%
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
write_csv(scatter_df, file.path(DAT_DIR, "fig2_concordance_scatter.csv"))
message(sprintf("Panel B: r = %.3f, rho = %.3f, concordance = %.1f%%",
                cor_r, cor_rho, sign_concordance))

ggsave(file.path(RPT_DIR, "test_panelB.pdf"), pB_base,
       width = 170, height = 160, units = "mm")
message("Panel B test saved")

# ═══ 10. PANEL C — RRHO2 Threshold-Free Concordance Map (Annotated) ═════════

message("Building Panel C: RRHO2 concordance map...")

# --- 1. Prepare ranked lists ---
rrho_list1 <- dep_df %>%
  dplyr::select(gene, t = t_Training_Young) %>%
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
  labels = c("Training (Young)", "Training (Old)"),
  log10.ind = TRUE, method = "hyper", boundary = 0.04
)

# --- 3. Extract max -log10(p) per quadrant ---
hmat <- rrho_obj$hypermat
nr <- nrow(hmat); nc <- ncol(hmat)
mid_r <- floor(nr / 2); mid_c <- floor(nc / 2)
max_UU <- max(hmat[1:mid_r, (mid_c+1):nc], na.rm = TRUE)       # top-left = Concordant up-up
max_DD <- max(hmat[(mid_r+1):nr, 1:mid_c], na.rm = TRUE)       # bottom-right = Concordant dn-dn
max_UD <- max(hmat[1:mid_r, 1:mid_c], na.rm = TRUE)            # top-right = Discordant up-dn
max_DU <- max(hmat[(mid_r+1):nr, (mid_c+1):nc], na.rm = TRUE)  # bottom-left = Discordant dn-up

# --- 4. Export ---
write.csv(hmat, file.path(DAT_DIR, "fig2_rrho2_matrix.csv"))

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
  # White crosshairs at midpoint
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "white", linewidth = 0.5) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "white", linewidth = 0.5) +
  # Quadrant labels
  annotate("label", x = 0.25, y = 0.75, label = paste0("Concordant up-up\nmax = ", round(max_UU, 1)),
           fill = alpha("white", 0.7), size = 2.0, fontface = "bold") +
  annotate("label", x = 0.75, y = 0.25, label = paste0("Concordant dn-dn\nmax = ", round(max_DD, 1)),
           fill = alpha("white", 0.7), size = 2.0, fontface = "bold") +
  annotate("label", x = 0.75, y = 0.75, label = paste0("Discordant up-dn\nmax = ", round(max_UD, 1)),
           fill = alpha("white", 0.7), size = 2.0, fontface = "bold") +
  annotate("label", x = 0.25, y = 0.25, label = paste0("Discordant dn-up\nmax = ", round(max_DU, 1)),
           fill = alpha("white", 0.7), size = 2.0, fontface = "bold") +
  # Axis labels
  annotate("text", x = 0.95, y = 0.02, label = "Most upregulated ->",
           hjust = 1, size = 1.8, color = "grey30") +
  annotate("text", x = 0.05, y = 0.02, label = "<- Most downregulated",
           hjust = 0, size = 1.8, color = "grey30") +
  labs(title = "RRHO2 Concordance Map",
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

# ═══ 11. PANEL D — mitch 2D Pathway Enrichment ═══════════════════════════════

message("Building Panel D: mitch 2D pathway enrichment...")

# --- 1. Build t-stat matrix ---
tstat_mat <- dep_df %>%
  dplyr::select(gene, Training_Young = t_Training_Young, Training_Old = t_Training_Old) %>%
  dplyr::filter(!is.na(Training_Young), !is.na(Training_Old)) %>%
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

# --- 4. Process results with 4-category classification ---
mitch_df <- mitch_res$enrichment_result %>%
  mutate(
    padj_TY  = p.adjust(p.Training_Young, method = "BH"),
    padj_TO  = p.adjust(p.Training_Old, method = "BH"),
    sig_TY   = padj_TY < 0.05,
    sig_TO   = padj_TO < 0.05,
    # Interaction-like: joint test significant AND in discordant quadrant
    discordant_quadrant = (s.Training_Young > 0 & s.Training_Old < 0) |
                          (s.Training_Young < 0 & s.Training_Old > 0),
    sig_joint = p.adjustMANOVA < 0.05 & discordant_quadrant,
    sig_cat  = case_when(
      sig_joint             ~ "Interaction",
      sig_TY & sig_TO       ~ "Sig Both",
      sig_TY                ~ "Sig Young only",
      sig_TO                ~ "Sig Old only",
      TRUE                  ~ "NS"
    ),
    sig_cat = factor(sig_cat, levels = c("Interaction", "Sig Both",
                                          "Sig Young only", "Sig Old only", "NS")),
    pathway_clean = clean_pathway_name(set)
  )

# --- 5. Identify pathways to label (all quadrants) ---
label_keywords <- c("ribosom", "oxphos", "oxidative phosph", "mtorc1",
                     "extracellular matrix", "myogenes", "glycoly",
                     "proteasome", "translation", "unfolded protein",
                     "mitotic spindle", "muscle", "respiratory chain",
                     "fatty acid", "inflammatory")

# Label from significant pathways + keyword matches
label_df <- bind_rows(
  # Keyword matches among significant
  mitch_df %>% filter(sig_cat != "NS") %>%
    filter(str_detect(tolower(set), paste(label_keywords, collapse = "|"))),
  # Top 3 per quadrant by effect size
  mitch_df %>% filter(sig_cat != "NS", s.Training_Young > 0, s.Training_Old > 0) %>%
    slice_max(abs(s.Training_Young) + abs(s.Training_Old), n = 3),
  mitch_df %>% filter(sig_cat != "NS", s.Training_Young < 0, s.Training_Old < 0) %>%
    slice_max(abs(s.Training_Young) + abs(s.Training_Old), n = 3),
  mitch_df %>% filter(sig_cat != "NS", s.Training_Young > 0, s.Training_Old < 0) %>%
    slice_max(abs(s.Training_Young) + abs(s.Training_Old), n = 3),
  mitch_df %>% filter(sig_cat != "NS", s.Training_Young < 0, s.Training_Old > 0) %>%
    slice_max(abs(s.Training_Young) + abs(s.Training_Old), n = 3)
) %>%
  distinct(set, .keep_all = TRUE) %>%
  slice_head(n = 15)

# --- 6. Correlations ---
pw_cor  <- cor(mitch_df$s.Training_Young, mitch_df$s.Training_Old)
pro_cor <- cor_r  # from Panel B

# --- 7. Axis range for shading ---
pw_lim <- max(abs(c(mitch_df$s.Training_Young, mitch_df$s.Training_Old)), na.rm = TRUE) * 1.1

# --- 8. Build scatter plot ---
pD <- ggplot(mitch_df %>% arrange(desc(as.integer(sig_cat))),
             aes(x = s.Training_Young, y = s.Training_Old)) +
  # Quadrant background shading (concordant = teal, discordant = amber)
  annotate("rect", xmin = 0, xmax = pw_lim, ymin = 0, ymax = pw_lim,
           fill = "#00897B", alpha = 0.04) +
  annotate("rect", xmin = -pw_lim, xmax = 0, ymin = -pw_lim, ymax = 0,
           fill = "#00897B", alpha = 0.04) +
  annotate("rect", xmin = 0, xmax = pw_lim, ymin = -pw_lim, ymax = 0,
           fill = "#FF8F00", alpha = 0.04) +
  annotate("rect", xmin = -pw_lim, xmax = 0, ymin = 0, ymax = pw_lim,
           fill = "#FF8F00", alpha = 0.04) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40", linewidth = 0.3) +
  # Points with 4-category + size by gene-set size
  geom_point(aes(color = sig_cat, size = setSize, alpha = sig_cat)) +
  scale_color_manual(values = SIG_COLORS, name = "Significance") +
  scale_alpha_manual(values = SIG_ALPHAS, guide = "none") +
  scale_size_continuous(range = c(0.3, 3.0), name = "Protein set size",
                        guide = guide_legend(override.aes = list(alpha = 0.8))) +
  # Labels
  geom_text_repel(data = label_df, aes(label = pathway_clean),
                  size = KEY_TEXT, max.overlaps = 25,
                  segment.size = 0.2, fontface = "italic",
                  min.segment.length = 0) +
  # Correlation annotation
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = sprintf("Pathway r = %.2f\nProtein r = %.2f", pw_cor, pro_cor),
           size = KEY_TITLE, fontface = "bold") +
  labs(x = "Enrichment Score Tr. (Young)",
       y = "Enrichment Score Tr. (Old)",
       title = "Pathway-Level Enrichment Concordance",
       subtitle = "Hallmark + GO:BP via mitch") +
  coord_cartesian(xlim = c(-pw_lim, pw_lim), ylim = c(-pw_lim, pw_lim)) +
  THEME_PUB +
  theme(legend.position = "bottom",
        legend.key.size = unit(2, "mm"),
        legend.text = element_text(size = 5))

# --- 9. Export ---
write_csv(mitch_df, file.path(DAT_DIR, "fig2_mitch_2d_results.csv"))
message(sprintf("mitch: %d pathways, Interaction=%d, Both=%d, Young=%d, Old=%d, r = %.3f",
                nrow(mitch_df),
                sum(mitch_df$sig_cat == "Interaction"),
                sum(mitch_df$sig_cat == "Sig Both"),
                sum(mitch_df$sig_cat == "Sig Young only"),
                sum(mitch_df$sig_cat == "Sig Old only"),
                pw_cor))

ggsave(file.path(RPT_DIR, "test_panelD.pdf"), pD,
       width = 170, height = 160, units = "mm")
message("Panel D test saved")

# ═══ 12. PANEL E — Interaction DEP Classification (Stacked Bar + Schematics) ═

message("Building Panel E: interaction DEP classification...")

# --- 1. Classify interaction DEPs into 4 categories ---
int_df <- dep_df %>%
  filter(pi_score_Interaction < 0.05) %>%
  dplyr::select(gene,
         logFC_Y = logFC_Training_Young, logFC_O = logFC_Training_Old,
         pi_Y = pi_score_Training_Young, pi_O = pi_score_Training_Old) %>%
  filter(!is.na(logFC_Y), !is.na(logFC_O)) %>%
  mutate(
    same_dir = sign(logFC_Y) == sign(logFC_O),
    category = case_when(
      !same_dir                                        ~ "Opposite Direction",
      same_dir & abs(logFC_O) > 2 * abs(logFC_Y)      ~ "Old-Specific",
      same_dir & abs(logFC_Y) > 2 * abs(logFC_O)      ~ "Attenuated",
      same_dir & pi_Y < 0.05 & pi_O >= 0.05           ~ "Young-Specific",
      same_dir & pi_O < 0.05 & pi_Y >= 0.05           ~ "Old-Specific",
      TRUE                                              ~ "Attenuated"
    ),
    category = factor(category, levels = c("Attenuated", "Opposite Direction",
                                            "Old-Specific", "Young-Specific"))
  )

cat_counts <- int_df %>% count(category, .drop = FALSE) %>% deframe()
cat_pcts   <- round(100 * cat_counts / sum(cat_counts), 0)
n_total    <- nrow(int_df)

message(sprintf("Panel E: %d interaction DEPs -- %s",
                n_total,
                paste(names(cat_counts), cat_counts, sep = "=", collapse = ", ")))

# --- 2. ORA per category (Hallmark) ---
ora_labels <- list()
for (cat in levels(int_df$category)) {
  genes_cat <- int_df %>% filter(category == cat) %>% pull(gene)
  if (length(genes_cat) < 3) {
    ora_labels[[cat]] <- paste0("n = ", length(genes_cat))
    next
  }
  ora_res <- tryCatch({
    enricher(gene = genes_cat, TERM2GENE = hallmark_t2g,
             universe = unique(dep_df$gene), pvalueCutoff = 0.2, qvalueCutoff = 1)
  }, error = function(e) NULL)

  if (is.null(ora_res) || nrow(as.data.frame(ora_res)) == 0) {
    ora_labels[[cat]] <- paste0("n = ", length(genes_cat), "; no sig. terms")
  } else {
    top_terms <- as.data.frame(ora_res) %>%
      arrange(p.adjust) %>%
      slice_head(n = 5) %>%
      mutate(label = clean_pathway_name(Description))
    ora_labels[[cat]] <- paste(top_terms$label, collapse = "\n")
  }
}

# --- 3. Build stacked bar data ---
bar_df <- tibble(
  category = factor(names(cat_counts), levels = levels(int_df$category)),
  count    = as.integer(cat_counts),
  pct      = as.numeric(cat_pcts)
) %>%
  mutate(
    xmax = cumsum(pct),
    xmin = lag(xmax, default = 0),
    xmid = (xmin + xmax) / 2,
    label = paste0(pct, "%")
  )

# --- 4. Build stacked bar plot ---
pE_bar <- ggplot(bar_df) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = category),
            color = "white", linewidth = 0.3) +
  geom_text(aes(x = xmid, y = 0.5, label = label),
            size = 3.5, fontface = "bold", color = "white") +
  scale_fill_manual(values = INTERACTION_CAT_COLORS, guide = "none") +
  # Category names below bar
  geom_text(aes(x = xmid, y = -0.15, label = str_wrap(category, width = 12)),
            size = 2.0, lineheight = 0.85) +
  # ORA pathway labels above bar
  geom_text(aes(x = xmid, y = 1.15,
                label = sapply(category, function(c) ora_labels[[c]])),
            size = 1.5, lineheight = 0.85, vjust = 0, color = "grey20") +
  coord_cartesian(xlim = c(-2, 102), ylim = c(-0.4, 2.5), clip = "off") +
  labs(title = "Distribution of Interaction DEPs",
       subtitle = sprintf("n = %d interaction DEPs (pi < 0.05)", n_total)) +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 8, hjust = 0.5),
        plot.subtitle = element_text(size = 6, color = "grey30", hjust = 0.5, face = "italic"))

# --- 5. Build schematic line plots using real exemplar proteins ---
# Exemplar for Attenuated: largest |logFC_Y|/|logFC_O| ratio
exemplar_att <- int_df %>%
  filter(category == "Attenuated") %>%
  mutate(ratio = abs(logFC_Y) / pmax(abs(logFC_O), 0.01)) %>%
  slice_max(ratio, n = 1)

# Exemplar for Opposite Direction: largest |logFC_Y - logFC_O|
exemplar_opp <- int_df %>%
  filter(category == "Opposite Direction") %>%
  mutate(dev = abs(logFC_Y - logFC_O)) %>%
  slice_max(dev, n = 1)

make_schematic <- function(lfc_y, lfc_o, gene_name, title_text) {
  sdf <- tibble(
    time  = rep(c("Pre", "Post"), 2),
    group = rep(c("Young", "Old"), each = 2),
    value = c(0, lfc_y, 0, lfc_o)
  ) %>%
    mutate(time = factor(time, levels = c("Pre", "Post")))

  ggplot(sdf, aes(x = time, y = value, group = group, color = group, linetype = group)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    scale_color_manual(values = c(Young = "#E05A4E", Old = "#5DA5DA")) +
    scale_linetype_manual(values = c(Young = "solid", Old = "dashed")) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.2) +
    labs(title = title_text, subtitle = gene_name,
         x = NULL, y = "logFC") +
    THEME_PUB +
    theme(legend.position = "none",
          plot.title = element_text(size = 7, face = "bold"),
          plot.subtitle = element_text(size = 6, face = "italic"),
          axis.text = element_text(size = 5),
          axis.title.y = element_text(size = 5))
}

pE_sch1 <- make_schematic(
  exemplar_att$logFC_Y[1], exemplar_att$logFC_O[1],
  exemplar_att$gene[1], "Attenuated"
)
pE_sch2 <- if (nrow(exemplar_opp) > 0) {
  make_schematic(
    exemplar_opp$logFC_Y[1], exemplar_opp$logFC_O[1],
    exemplar_opp$gene[1], "Opposite Direction"
  )
} else {
  ggplot() + theme_void()  # placeholder if no opposite-direction DEPs
}

# --- 6. Assemble Panel E: bar (left 60%) | schematics (right 40%) ---
pE <- pE_bar | (pE_sch1 / pE_sch2) +
  plot_layout(widths = c(0.6, 0.4))

# --- 7. Export ---
write_csv(int_df, file.path(DAT_DIR, "fig2_interaction_classification.csv"))
ggsave(file.path(RPT_DIR, "test_panelE.pdf"), pE,
       width = 250, height = 100, units = "mm")
message("Panel E test saved")

# ═══ 13. PANEL F — GO:BP Enrichment: Concordant vs Discordant (rrvgo) ═══════

message("Building Panel F: GO:BP enrichment (rrvgo-reduced)...")

# --- 1. Prepare gene lists ---
concordant_genes <- scatter_df %>%
  filter(concordant, sig_cat %in% c("Sig Both", "Sig Young only", "Sig Old only")) %>%
  pull(gene)

discordant_genes <- int_df$gene  # interaction DEPs from Panel E

all_genes <- unique(dep_df$gene)
message(sprintf("  Concordant DEPs: %d | Discordant (interaction) DEPs: %d | Universe: %d",
                length(concordant_genes), length(discordant_genes), length(all_genes)))

# --- 2. Run enrichGO for each set ---
run_gobp_ora <- function(genes, label) {
  if (length(genes) < 5) return(tibble())
  res <- tryCatch({
    enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
             ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
             universe = all_genes)
  }, error = function(e) { message("  enrichGO error: ", e$message); NULL })

  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(tibble())

  # Apply rrvgo reduction
  res_df <- as.data.frame(res)
  sim_mat <- tryCatch({
    hsGO <- GOSemSim::godata(annoDb = "org.Hs.eg.db", ont = "BP")
    calculateSimMatrix(res_df$ID, orgdb = "org.Hs.eg.db", ont = "BP", semdata = hsGO,
                       method = "Rel")
  }, error = function(e) { message("  rrvgo sim matrix error: ", e$message); NULL })

  if (!is.null(sim_mat) && nrow(sim_mat) > 1) {
    scores <- setNames(-log10(res_df$p.adjust), res_df$ID)
    reduced <- reduceSimMatrix(sim_mat, scores = scores, threshold = 0.7,
                               orgdb = "org.Hs.eg.db")
    # Keep only parent terms
    parent_ids <- unique(reduced$parentTerm)
    res_df <- res_df %>% filter(ID %in% parent_ids)
  }

  res_df %>%
    arrange(p.adjust) %>%
    slice_head(n = 10) %>%
    mutate(Set = label, neg_log10_padj = -log10(p.adjust),
           Description = str_wrap(Description, width = 40))
}

conc_res <- run_gobp_ora(concordant_genes, "Concordant")
disc_res <- run_gobp_ora(discordant_genes, "Discordant")

enrich_df <- bind_rows(conc_res, disc_res)

message(sprintf("  After rrvgo: Concordant=%d, Discordant=%d terms",
                nrow(conc_res), nrow(disc_res)))

# --- 3. Build grouped bar chart ---
if (nrow(enrich_df) == 0) {
  pF <- ggplot() + annotate("text", x = 0.5, y = 0.5,
    label = "No significant GO:BP terms\nat FDR < 0.05", size = 3) + theme_void()
} else {
  # Add placeholder for empty categories
  if (nrow(conc_res) == 0) {
    enrich_df <- bind_rows(enrich_df,
      tibble(Description = paste0("No enriched terms\n(n = ", length(concordant_genes), ")"),
             neg_log10_padj = 0, Set = "Concordant"))
  }
  if (nrow(disc_res) == 0) {
    enrich_df <- bind_rows(enrich_df,
      tibble(Description = paste0("No enriched terms\n(n = ", length(discordant_genes), ")"),
             neg_log10_padj = 0, Set = "Discordant"))
  }

  pF <- ggplot(enrich_df, aes(x = neg_log10_padj,
                                y = reorder(Description, neg_log10_padj),
                                fill = Set)) +
    geom_col(position = position_dodge2(width = 0.8, preserve = "single"),
             width = 0.7) +
    scale_fill_manual(values = c(Concordant = "#00897B", Discordant = "#FF8F00"),
                      name = NULL) +
    labs(x = expression(-log[10]~(p[adj])),
         y = NULL,
         title = "GO:BP Enrichment",
         subtitle = sprintf("rrvgo-reduced (0.7) | Conc. n=%d, Disc. n=%d",
                            length(concordant_genes), length(discordant_genes))) +
    THEME_PUB +
    theme(axis.text.y = element_text(size = 5),
          legend.position = "bottom",
          legend.key.size = unit(3, "mm"),
          legend.text = element_text(size = 6))
}

# --- 4. Export ---
write_csv(enrich_df, file.path(DAT_DIR, "fig2_concordant_discordant_enrichment.csv"))
ggsave(file.path(RPT_DIR, "test_panelF.pdf"), pF,
       width = 180, height = 160, units = "mm")
message("Panel F test saved")

# ═══ 14. FINAL ASSEMBLY — 6-Panel Composite Figure ══════════════════════════

message("Assembling Figure 2...")

pA_wrapped <- wrap_elements(full = pA)

fig2 <- (pA_wrapped | pC_gg) /
         (pB_base   | pD) /
         (pE        | pF) +
  plot_layout(
    widths  = c(0.55, 0.45),
    heights = c(0.30, 0.38, 0.32)
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 12))
  )

ggsave(file.path(RPT_DIR, "Figure_2.pdf"), fig2,
       width = 380, height = 500, units = "mm", limitsize = FALSE)
ggsave(file.path(RPT_DIR, "Figure_2.png"), fig2,
       width = 380, height = 500, units = "mm", dpi = 300, limitsize = FALSE)

message("Figure 2 saved to: ", RPT_DIR)
