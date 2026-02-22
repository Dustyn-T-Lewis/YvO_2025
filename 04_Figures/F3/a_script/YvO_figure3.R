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

# --- Volcano unified direction coloring ---
VOLC_COLORS <- list(
  Aging        = c(Up = "#D6604D", Down = "#4393C3", NS = "grey80"),
  Training_Old = c(Up = "#D6604D", Down = "#4393C3", NS = "grey80")
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
      relaxed_sig = if_any(all_of(pi_cols), ~ . < threshold_relax)
    )

  all_genes <- unique(scatter_df$gene)
  results <- list()

  for (q in c("Q1", "Q2", "Q3", "Q4")) {
    genes_q <- qdf %>% filter(quadrant == q, relaxed_sig) %>% pull(gene)

    if (length(genes_q) < min_n) {
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

# --- NES bar inset builder ---
make_nes_inset <- function(pw_df, direction = "up") {
  if (is.null(pw_df) || nrow(pw_df) == 0) return(NULL)
  pw_df <- pw_df %>%
    mutate(stars = sig_stars(padj),
           label = clean_pathway_name(pathway))

  # Ensure minimum bar width for text containment
  min_width <- 0.8
  pw_df <- pw_df %>%
    mutate(bar_width = pmax(abs(NES), min_width))

  bar_fill <- if (direction == "up") "#D6604D" else "#4393C3"

  if (direction == "down") {
    pw_df <- pw_df %>% mutate(xmin = 0, xmax = bar_width)
  } else {
    pw_df <- pw_df %>% mutate(xmin = -bar_width, xmax = 0)
  }

  pw_df <- pw_df %>%
    mutate(y_pos = row_number(),
           bar_label = paste0(label, " ", stars))

  ggplot(pw_df) +
    geom_rect(aes(xmin = xmin, xmax = xmax,
                  ymin = y_pos - 0.4, ymax = y_pos + 0.4),
              fill = bar_fill, color = "black", linewidth = 0.15) +
    geom_text(aes(x = (xmin + xmax) / 2, y = y_pos, label = bar_label),
              color = "white", fontface = "bold", size = 1.8,
              hjust = 0.5) +
    scale_y_continuous(breaks = NULL) +
    theme_void() +
    theme(plot.background = element_blank(),
          plot.margin = margin(0, 0, 0, 0))
}

# --- ORA bar inset builder ---
make_ora_inset <- function(q_terms, direction, fill_color) {
  if (is.null(q_terms) || nrow(q_terms) == 0 || all(is.na(q_terms$p.adjust))) return(NULL)
  q_terms <- q_terms %>%
    filter(!is.na(p.adjust)) %>%
    mutate(neg_log10_padj = -log10(p.adjust),
           stars = sig_stars(p.adjust),
           label = paste0(Description, " ", stars))
  if (nrow(q_terms) == 0) return(NULL)
  min_width <- max(q_terms$neg_log10_padj) * 0.5
  q_terms <- q_terms %>% mutate(bar_width = pmax(neg_log10_padj, min_width))
  if (direction == "left") {
    q_terms <- q_terms %>% mutate(xmin = -bar_width, xmax = 0)
  } else {
    q_terms <- q_terms %>% mutate(xmin = 0, xmax = bar_width)
  }
  q_terms <- q_terms %>% mutate(y_pos = row_number())
  ggplot(q_terms) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = y_pos - 0.4, ymax = y_pos + 0.4),
              fill = fill_color, alpha = 0.7, color = "black", linewidth = 0.1) +
    geom_text(aes(x = (xmin + xmax) / 2, y = y_pos, label = label),
              color = "white", fontface = "bold", size = 1.5, hjust = 0.5) +
    theme_void() +
    theme(plot.background = element_blank(), plot.margin = margin(0, 0, 0, 0))
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
  strip_title <- if (ctr == "Aging") {
    "Age-Related Proteomic Changes"
  } else {
    "Training Response in Old Adults"
  }
  strip_subtitle <- "DEPs by pi-score < 0.05 | Hallmark pathway enrichment inset"

  # --- Pi-score boundary curve ---
  pi_threshold <- -log10(0.05)
  curve_logFC <- seq(0.05, x_max_est, length.out = 200)
  pi_curve <- data.frame(
    logFC = c(-rev(curve_logFC), curve_logFC),
    neg_log10p = c(rev(pi_threshold / curve_logFC), pi_threshold / curve_logFC)
  )
  pi_curve <- pi_curve %>% filter(neg_log10p <= y_max_est * 1.1)

  # --- 5. Build ggplot ---
  p <- ggplot(vdf, aes(x = logFC, y = neg_log10p)) +
    # Points — unified direction coloring
    geom_point(aes(color = direction), size = 1.2, alpha = 0.65) +
    scale_color_manual(values = volc_cols) +
    # Pi-score boundary curve
    geom_line(data = pi_curve, aes(x = logFC, y = neg_log10p),
              color = "grey40", linetype = "dashed", linewidth = 0.3,
              inherit.aes = FALSE) +
    # Gene labels
    geom_text_repel(
      data = top_genes, aes(label = gene),
      size = KEY_TEXT, max.overlaps = 15,
      segment.size = 0.2, fontface = "italic",
      min.segment.length = 0
    ) +
    # DEP count — upper-right (Up) — boxed
    annotate("label",
             x = x_max_est * 0.98, y = y_max_est * 0.99,
             label = paste0(n_up, " Up", dir_note_up),
             hjust = 1, vjust = 1, size = KEY_TITLE,
             color = "#D6604D",
             fill = alpha("#D6604D", 0.12),
             label.padding = unit(2.5, "pt"),
             fontface = "bold") +
    # DEP count — upper-left (Down) — boxed
    annotate("label",
             x = -x_max_est * 0.98, y = y_max_est * 0.99,
             label = paste0(n_down, " Down", dir_note_down),
             hjust = 0, vjust = 1, size = KEY_TITLE,
             color = "#4393C3",
             fill = alpha("#4393C3", 0.12),
             label.padding = unit(2.5, "pt"),
             fontface = "bold") +
    labs(
      title    = strip_title,
      subtitle = strip_subtitle,
      x = expression(log[2]~fold~change),
      y = expression(-log[10]~italic(P))
    ) +
    THEME_PUB +
    theme(legend.position = "none")

  # --- 6. Pathway NES bar insets — bottom corners via inset_element() ---
  inset_up   <- make_nes_inset(pw_up, direction = "up")
  inset_down <- make_nes_inset(pw_down, direction = "down")

  if (!is.null(inset_up)) {
    p <- p + inset_element(inset_up,
                            left = 0.58, right = 1.0,
                            bottom = 0.0, top = 0.25)
  }
  if (!is.null(inset_down)) {
    p <- p + inset_element(inset_down,
                            left = 0.0, right = 0.42,
                            bottom = 0.0, top = 0.25)
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
axis_lim <- quantile(abs(c(scatter_df$logFC_A, scatter_df$logFC_TO)), 0.98, na.rm = TRUE) * 1.1

# --- 6. Build scatter plot ---
scatter_ordered <- scatter_df %>%
  mutate(plot_order = as.integer(sig_cat)) %>%
  arrange(desc(plot_order))

pB_base <- ggplot(scatter_ordered, aes(x = logFC_A, y = logFC_TO)) +
  # Reversed quadrants (Q2, Q4) = teal
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "#00897B", alpha = 0.18) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#00897B", alpha = 0.18) +
  # Exacerbated quadrants (Q1, Q3) = amber
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "#FF8F00", alpha = 0.18) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#FF8F00", alpha = 0.18) +
  # Reference lines
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  # NS layer (small, faded)
  geom_point(data = . %>% filter(sig_cat == "NS"),
             aes(color = sig_cat), size = 0.4, alpha = 0.15) +
  # Significant layer (prominent)
  geom_point(data = . %>% filter(sig_cat != "NS"),
             aes(color = sig_cat), size = 1.5, alpha = 0.85) +
  scale_color_manual(values = SIG_COLORS, name = "Significance") +
  # Boxed quadrant labels
  annotate("label", x = axis_lim * 0.95, y = axis_lim * 0.95,
           label = paste0("Exacerbated Up-Up\nn = ", q_counts["Q1"]),
           hjust = 1, vjust = 1, size = 2.0, fontface = "bold",
           color = "white", fill = alpha("#FF8F00", 0.7),
           label.padding = unit(2, "pt"), label.size = 0) +
  annotate("label", x = axis_lim * 0.95, y = -axis_lim * 0.95,
           label = paste0("Reversed Up-Down\nn = ", q_counts["Q2"]),
           hjust = 1, vjust = 0, size = 2.0, fontface = "bold",
           color = "white", fill = alpha("#00897B", 0.7),
           label.padding = unit(2, "pt"), label.size = 0) +
  annotate("label", x = -axis_lim * 0.95, y = -axis_lim * 0.95,
           label = paste0("Exacerbated Dn-Dn\nn = ", q_counts["Q3"]),
           hjust = 0, vjust = 0, size = 2.0, fontface = "bold",
           color = "white", fill = alpha("#FF8F00", 0.7),
           label.padding = unit(2, "pt"), label.size = 0) +
  annotate("label", x = -axis_lim * 0.95, y = axis_lim * 0.95,
           label = paste0("Reversed Dn-Up\nn = ", q_counts["Q4"]),
           hjust = 0, vjust = 1, size = 2.0, fontface = "bold",
           color = "white", fill = alpha("#00897B", 0.7),
           label.padding = unit(2, "pt"), label.size = 0) +
  labs(
    title = "Protein-Level Reversal of Age-Related Changes",
    subtitle = sprintf("r = %.2f [%.2f, %.2f] | rho = %.2f | r(null) ~ %.2f | Reversal: %.0f%%",
                       cor_r, cor_r_ci[1], cor_r_ci[2], cor_rho, r_null, reversal_ratio),
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

# --- 7. Build ORA bar insets for each quadrant ---
q_insets <- list(
  Q1 = make_ora_inset(qora_results %>% filter(quadrant == "Q1"), "left", "#FF8F00"),
  Q2 = make_ora_inset(qora_results %>% filter(quadrant == "Q2"), "left", "#00897B"),
  Q3 = make_ora_inset(qora_results %>% filter(quadrant == "Q3"), "right", "#FF8F00"),
  Q4 = make_ora_inset(qora_results %>% filter(quadrant == "Q4"), "right", "#00897B")
)

# Attach insets to quadrant corners
if (!is.null(q_insets$Q1))
  pB_base <- pB_base + inset_element(q_insets$Q1, left = 0.55, right = 1.0, bottom = 0.75, top = 1.0)
if (!is.null(q_insets$Q2))
  pB_base <- pB_base + inset_element(q_insets$Q2, left = 0.55, right = 1.0, bottom = 0.0, top = 0.25)
if (!is.null(q_insets$Q3))
  pB_base <- pB_base + inset_element(q_insets$Q3, left = 0.0, right = 0.45, bottom = 0.0, top = 0.25)
if (!is.null(q_insets$Q4))
  pB_base <- pB_base + inset_element(q_insets$Q4, left = 0.0, right = 0.45, bottom = 0.75, top = 1.0)

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

# --- 5. Render as ggplot with viridis colorscale ---
hmat_df <- expand.grid(row = 1:nr, col = 1:nc) %>%
  mutate(value = as.vector(hmat))

# Quadrant label placement follows F2 convention:
# In ggplot with row on y-axis: row 1 = bottom, col 1 = left
# hmat indices map to ggplot coordinates directly.
#
# Existing quadrant extraction:
#   max_UU = hmat[1:mid_r, (mid_c+1):nc]  -> low rows, high cols -> bottom-right in ggplot
#   max_DD = hmat[(mid_r+1):nr, 1:mid_c]  -> high rows, low cols -> top-left in ggplot
#   max_UD = hmat[1:mid_r, 1:mid_c]       -> low rows, low cols  -> bottom-left in ggplot
#   max_DU = hmat[(mid_r+1):nr, (mid_c+1):nc] -> high rows, high cols -> top-right in ggplot
#
# Biological meaning (list1=Aging t-stats, list2=Training_Old t-stats):
#   Bottom-left (max_UD area): Aging Down x Tr.Old Down = Exacerbated Down-Down
#   Top-right (max_DU area):   Aging Up x Tr.Old Up     = Exacerbated Up-Up
#   Top-left (max_DD area):    Aging Up x Tr.Old Down   = Reversed
#   Bottom-right (max_UU area): Aging Down x Tr.Old Up  = Reversed

pC_gg <- ggplot(hmat_df, aes(x = col, y = row, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c(option = "viridis", name = expression(-log[10](P)),
                        guide = guide_colorbar(barwidth = unit(3, "cm"),
                                               barheight = unit(0.3, "cm"),
                                               title.position = "left",
                                               title.theme = element_text(size = 5.5, vjust = 0.8))) +
  # White crosshair lines at midpoint
  geom_hline(yintercept = mid_r + 0.5, linetype = "dashed", color = "white", linewidth = 0.5) +
  geom_vline(xintercept = mid_c + 0.5, linetype = "dashed", color = "white", linewidth = 0.5) +
  # Quadrant annotations — bold white, reversal terminology
  # Top-right (high rows, high cols) = max_DU area = Exacerbated Up-Up
  annotate("text", x = mid_c + (nc - mid_c) * 0.5, y = mid_r + (nr - mid_r) * 0.5,
           label = sprintf("Exacerbated\nAging Up / Tr. Up\nmax = %.1f", max_DU),
           color = "white", fontface = "bold", size = 2.0) +
  # Bottom-left (low rows, low cols) = max_UD area = Exacerbated Down-Down
  annotate("text", x = mid_c * 0.5, y = mid_r * 0.5,
           label = sprintf("Exacerbated\nAging Dn / Tr. Dn\nmax = %.1f", max_UD),
           color = "white", fontface = "bold", size = 2.0) +
  # Top-left (high rows, low cols) = max_DD area = Reversed (Aging Up, Tr. Down)
  annotate("text", x = mid_c * 0.5, y = mid_r + (nr - mid_r) * 0.5,
           label = sprintf("Reversed\nAging Up / Tr. Down\nmax = %.1f", max_DD),
           color = "white", fontface = "bold", size = 2.0) +
  # Bottom-right (low rows, high cols) = max_UU area = Reversed (Aging Dn, Tr. Up)
  annotate("text", x = mid_c + (nc - mid_c) * 0.5, y = mid_r * 0.5,
           label = sprintf("Reversed\nAging Dn / Tr. Up\nmax = %.1f", max_UU),
           color = "white", fontface = "bold", size = 2.0) +
  # Axis labels
  labs(title = "Threshold-Free Reversal of Age-Related Changes",
       subtitle = "RRHO2 hypergeometric overlap, -log10(P)",
       x = "Aging rank",
       y = "Training (Old) rank") +
  # Direction annotations on axes
  annotate("text", x = 1, y = -nr * 0.04,
           label = "<- Most downregulated",
           hjust = 0, size = 1.8, color = "grey30") +
  annotate("text", x = nc, y = -nr * 0.04,
           label = "Most upregulated ->",
           hjust = 1, size = 1.8, color = "grey30") +
  annotate("text", x = -nc * 0.04, y = 1, angle = 90,
           label = "<- Most downregulated",
           hjust = 0, size = 1.8, color = "grey30") +
  annotate("text", x = -nc * 0.04, y = nr, angle = 90,
           label = "Most upregulated ->",
           hjust = 1, size = 1.8, color = "grey30") +
  coord_cartesian(clip = "off") +
  THEME_PUB +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom")

# --- 6. Test save ---
ggsave(file.path(RPT_DIR, "test_panelC.pdf"), pC_gg,
       width = 170, height = 150, units = "mm")

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

mitch_df <- mitch_df %>%
  mutate(database = ifelse(str_starts(set, "HALLMARK_"), "Hallmark", "GO:BP"))

# --- 5. Per-quadrant counts split by database ---
q_counts_d <- mitch_df %>%
  filter(sig_cat != "NS") %>%
  mutate(quadrant = case_when(
    s.Aging > 0 & s.Training_Old > 0 ~ "Q1",
    s.Aging > 0 & s.Training_Old < 0 ~ "Q2",
    s.Aging < 0 & s.Training_Old < 0 ~ "Q3",
    TRUE ~ "Q4"
  )) %>%
  count(quadrant, database) %>%
  pivot_wider(names_from = database, values_from = n, values_fill = 0)

make_q_label <- function(q_name, label_name, counts_df) {
  row <- counts_df %>% filter(quadrant == q_name)
  h_count <- if (nrow(row) > 0 && "Hallmark" %in% names(row)) row$Hallmark else 0
  bp_count <- if (nrow(row) > 0 && "GO:BP" %in% names(row)) row$`GO:BP` else 0
  paste0(label_name, "\nHallmark: ", h_count, " / GO:BP: ", bp_count)
}

# --- 6. Tissue-relevance filter ---
exclude_keywords <- c("sperm", "zona pellucida", "spermat", "oocyte",
                       "embryonic", "placenta", "retina", "olfactory",
                       "photoreceptor", "taste", "pollen", "fertiliz")

filter_relevant <- function(df) {
  df %>% filter(!str_detect(tolower(set), paste(exclude_keywords, collapse = "|")))
}

# --- 7. Identify pathways to label (set-size-prioritized) ---
label_df <- mitch_df %>%
  filter(sig_cat != "NS") %>%
  filter_relevant() %>%
  group_by(sig_cat) %>%
  arrange(desc(setSize)) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  distinct(set, .keep_all = TRUE)

# --- 8. Correlations ---
pw_cor  <- cor(mitch_df$s.Aging, mitch_df$s.Training_Old)
pro_cor <- cor_r

# --- 9. Axis range ---
pw_lim <- max(abs(c(mitch_df$s.Aging, mitch_df$s.Training_Old)), na.rm = TRUE) * 1.1

# --- 10. Build scatter plot ---
pD <- ggplot(mitch_df %>% arrange(desc(as.integer(sig_cat))),
             aes(x = s.Aging, y = s.Training_Old)) +
  # Reversed quadrants (Q2: Aging>0 & Tr.Old<0, Q4: Aging<0 & Tr.Old>0) = teal
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "#00897B", alpha = 0.18) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#00897B", alpha = 0.18) +
  # Exacerbated quadrants (Q1: both>0, Q3: both<0) = amber
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "#FF8F00", alpha = 0.18) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#FF8F00", alpha = 0.18) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40", linewidth = 0.3) +
  geom_point(aes(color = sig_cat, size = setSize, alpha = sig_cat)) +
  scale_color_manual(values = SIG_COLORS, name = "Significance") +
  scale_alpha_manual(values = SIG_ALPHAS, guide = "none") +
  scale_size_continuous(range = c(0.8, 5.0),
                        breaks = c(50, 150, 300, 500),
                        name = "Set size",
                        guide = guide_legend(override.aes = list(alpha = 0.8))) +
  geom_label_repel(data = label_df, aes(label = pathway_clean, fill = sig_cat),
                   color = "white", fontface = "bold", size = KEY_TEXT,
                   max.overlaps = 25, segment.size = 0.2,
                   label.padding = unit(1.5, "pt"), label.size = 0,
                   min.segment.length = 0, show.legend = FALSE) +
  scale_fill_manual(values = SIG_COLORS, guide = "none") +
  # Quadrant labels with counts
  annotate("label", x = pw_lim * 0.95, y = pw_lim * 0.95,
           label = make_q_label("Q1", "Exacerbated Up", q_counts_d),
           hjust = 1, vjust = 1, size = 1.8, fontface = "bold",
           color = "white", fill = alpha("#FF8F00", 0.7),
           label.padding = unit(2, "pt"), label.size = 0) +
  annotate("label", x = pw_lim * 0.95, y = -pw_lim * 0.95,
           label = make_q_label("Q2", "Reversed", q_counts_d),
           hjust = 1, vjust = 0, size = 1.8, fontface = "bold",
           color = "white", fill = alpha("#00897B", 0.7),
           label.padding = unit(2, "pt"), label.size = 0) +
  annotate("label", x = -pw_lim * 0.95, y = -pw_lim * 0.95,
           label = make_q_label("Q3", "Exacerbated Down", q_counts_d),
           hjust = 0, vjust = 0, size = 1.8, fontface = "bold",
           color = "white", fill = alpha("#FF8F00", 0.7),
           label.padding = unit(2, "pt"), label.size = 0) +
  annotate("label", x = -pw_lim * 0.95, y = pw_lim * 0.95,
           label = make_q_label("Q4", "Reversed", q_counts_d),
           hjust = 0, vjust = 1, size = 1.8, fontface = "bold",
           color = "white", fill = alpha("#00897B", 0.7),
           label.padding = unit(2, "pt"), label.size = 0) +
  labs(
    title = "Pathway-Level Reversal of Age-Related Changes",
    subtitle = sprintf("Hallmark + GO:BP via mitch | Pathway r = %.2f, Protein r = %.2f", pw_cor, pro_cor),
    x = "Enrichment Score (Aging)",
    y = "Enrichment Score (Training Old)"
  ) +
  coord_cartesian(xlim = c(-pw_lim, pw_lim), ylim = c(-pw_lim, pw_lim)) +
  THEME_PUB +
  theme(legend.position = "bottom",
        legend.key.size = unit(2, "mm"),
        legend.text = element_text(size = 5))

# --- 11. Export ---
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

# === 12. PANEL E — Reversal Classification (Dumbbell Plot) ================

message("Building Panel E: reversal classification (dumbbell)...")

# --- 1. Classify Aging DEPs by Training_Old behavior ---
rev_class <- dep_df %>%
  filter(pi_score_Aging < 0.05) %>%
  dplyr::select(gene, logFC_A = logFC_Aging, logFC_TO = logFC_Training_Old,
         pi_TO = pi_score_Training_Old) %>%
  filter(!is.na(logFC_A), !is.na(logFC_TO)) %>%
  mutate(
    opposite_dir = sign(logFC_A) != sign(logFC_TO),
    category = case_when(
      opposite_dir & pi_TO < 0.05          ~ "Fully Reversed",
      opposite_dir & abs(logFC_TO) > 0.1   ~ "Partially Reversed",
      !opposite_dir & pi_TO < 0.05         ~ "Exacerbated",
      TRUE                                  ~ "Non-Reversed"
    ),
    category = factor(category, levels = c("Fully Reversed", "Partially Reversed",
                                            "Non-Reversed", "Exacerbated"))
  )

n_total_original <- nrow(rev_class)  # Store original count
cat_counts <- rev_class %>% count(category) %>% deframe()
cat_pcts   <- round(100 * cat_counts / sum(cat_counts), 1)

message(sprintf("Panel E: %d Aging DEPs classified", n_total_original))

# --- 2. Dynamic protein count handling ---
# If too many proteins for readable dumbbell, truncate per category
max_per_cat <- 15
if (n_total_original > 60) {
  rev_class <- rev_class %>%
    group_by(category) %>%
    arrange(category, desc(abs(logFC_TO))) %>%
    slice_head(n = max_per_cat) %>%
    ungroup()
  message(sprintf("  Truncated to top %d per category (%d shown of %d total)",
                  max_per_cat, nrow(rev_class), n_total_original))
}

# --- 3. Facet labels with counts and percentages ---
rev_class <- rev_class %>%
  mutate(
    facet_label = sprintf("%s (n = %d, %.0f%%)",
                          as.character(category),
                          cat_counts[as.character(category)],
                          cat_pcts[as.character(category)])
  )

facet_levels <- levels(rev_class$category) %>%
  sapply(function(cat) {
    sprintf("%s (n = %d, %.0f%%)", cat, cat_counts[cat], cat_pcts[cat])
  })
rev_class$facet_label <- factor(rev_class$facet_label, levels = facet_levels)

# Sort by reversal magnitude within each category
rev_class <- rev_class %>%
  arrange(category, abs(logFC_TO)) %>%
  mutate(gene = factor(gene, levels = gene))

# --- 4. Build dumbbell plot ---
pE <- ggplot(rev_class) +
  # Background shading per facet
  geom_rect(aes(fill = category),
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
            alpha = 0.12, show.legend = FALSE) +
  scale_fill_manual(values = REVERSAL_CAT_COLORS) +
  # Connecting segments
  geom_segment(aes(x = logFC_A, xend = logFC_TO,
                   y = gene, yend = gene),
               color = "grey40", linewidth = 0.4) +
  # Aging dots (circle) — green
  geom_point(aes(x = logFC_A, y = gene),
             color = "#4CAF50", size = 2.0, shape = 16) +
  # Training_Old dots (triangle) — blue
  geom_point(aes(x = logFC_TO, y = gene),
             color = "#5DA5DA", size = 2.0, shape = 17) +
  # Reference line
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  # Facet by category
  facet_wrap(~ facet_label, scales = "free_y", ncol = 1) +
  labs(
    title = "Protein-Level Reversal by Training in Old Adults",
    subtitle = sprintf("%d Aging DEPs classified by training reversal pattern, pi < 0.05", n_total_original),
    x = expression(log[2]~fold~change),
    y = NULL
  ) +
  THEME_PUB +
  theme(axis.text.y = element_text(size = 5),
        strip.text = element_text(size = 6, face = "bold"))

# --- 5. Export ---
write_csv(rev_class, file.path(DAT_DIR, "fig3_reversal_classification.csv"))
ggsave(file.path(RPT_DIR, "test_panelE.pdf"), pE,
       width = 180, height = 180, units = "mm")
message(sprintf("Panel E: %d proteins shown (%s)",
                nrow(rev_class),
                paste(names(cat_counts), cat_counts, sep = "=", collapse = ", ")))
message("Panel E test saved")

# === 13. PANEL F — Pathway Reversal Waterfall ================================

message("Building Panel F: pathway reversal waterfall...")

# --- 1. Get fgsea results for both contrasts ---
hallmark_fgsea <- fgsea_all %>%
  filter(database == "Hallmark") %>%
  dplyr::select(pathway, contrast, NES, padj)

# --- 2. Get GO:BP fgsea results and reduce with rrvgo ---
gobp_fgsea <- fgsea_all %>%
  filter(database == "GOBP") %>%
  dplyr::select(pathway, contrast, NES, padj)

# Reduce GO:BP terms with rrvgo (threshold = 0.85 for detailed panels)
gobp_sig_ids <- gobp_fgsea %>%
  filter(padj < 0.05) %>%
  pull(pathway) %>% unique() %>%
  str_remove("^GOBP_") %>%
  str_replace_all("_", " ")

if (length(gobp_sig_ids) > 5) {
  gobp_go_ids <- tryCatch({
    go_terms <- AnnotationDbi::select(org.Hs.eg.db,
                                       keys = gobp_sig_ids,
                                       keytype = "GOALL",
                                       columns = "GOALL")
    unique(go_terms$GOALL)
  }, error = function(e) {
    message("rrvgo GO ID mapping failed, keeping all GO:BP terms")
    NULL
  })

  if (!is.null(gobp_go_ids) && length(gobp_go_ids) > 5) {
    tryCatch({
      simMatrix <- calculateSimMatrix(gobp_go_ids, orgdb = "org.Hs.eg.db",
                                       ont = "BP", method = "Rel")
      reducedTerms <- reduceSimMatrix(simMatrix, threshold = 0.85)
      gobp_fgsea <- gobp_fgsea %>%
        mutate(go_clean = str_remove(pathway, "^GOBP_") %>% str_replace_all("_", " ")) %>%
        filter(go_clean %in% reducedTerms$term |
               pathway %in% paste0("GOBP_", str_replace_all(reducedTerms$term, " ", "_")))
    }, error = function(e) {
      message("rrvgo reduction failed: ", e$message, " — keeping all GO:BP terms")
    })
  }
}

# Fallback: if too many GO:BP terms remain, take top by significance
if (nrow(gobp_fgsea) > 50) {
  top_gobp <- gobp_fgsea %>%
    filter(padj < 0.05) %>%
    group_by(pathway) %>%
    summarise(min_padj = min(padj), .groups = "drop") %>%
    arrange(min_padj) %>%
    slice_head(n = 20) %>%
    pull(pathway)
  gobp_fgsea <- gobp_fgsea %>% filter(pathway %in% top_gobp)
}

# --- 3. Combine and pivot wide ---
combined_fgsea <- bind_rows(
  hallmark_fgsea %>% mutate(db = "Hallmark"),
  gobp_fgsea %>% mutate(db = "GO:BP")
)

fw <- combined_fgsea %>%
  pivot_wider(id_cols = c(pathway, db), names_from = contrast,
              values_from = c(NES, padj)) %>%
  filter(!is.na(NES_Aging) | !is.na(NES_Training_Old)) %>%
  mutate(
    NES_Aging        = replace_na(NES_Aging, 0),
    NES_Training_Old = replace_na(NES_Training_Old, 0),
    delta_NES = NES_Aging - NES_Training_Old,
    # Classify
    sig_A = !is.na(padj_Aging) & padj_Aging < 0.05,
    sig_O = !is.na(padj_Training_Old) & padj_Training_Old < 0.05,
    sig_I = !is.na(padj_Interaction) & padj_Interaction < 0.25,
    opposite_sign = sign(NES_Aging) != sign(NES_Training_Old),
    pw_cat = case_when(
      sig_I                              ~ "Reversed",
      sig_A & sig_O & opposite_sign      ~ "Reversed",
      sig_A & sig_O & !opposite_sign     ~ "Exacerbated",
      sig_A | sig_O                      ~ "Attenuated",
      TRUE                               ~ "NS"
    ),
    pw_cat = factor(pw_cat, levels = c("Reversed", "Exacerbated", "Attenuated", "NS")),
    pathway_clean = clean_pathway_name(pathway),
    int_stars = ifelse(!is.na(padj_Interaction) & padj_Interaction < 0.05,
                       sig_stars(padj_Interaction), "")
  ) %>%
  filter(pw_cat != "NS") %>%
  arrange(desc(abs(delta_NES))) %>%
  slice_head(n = 25) %>%
  mutate(pathway_clean = factor(pathway_clean, levels = rev(pathway_clean)))

# --- 4. Category colors ---
PF_COLORS <- c(Reversed = "#00897B", Exacerbated = "#FF8F00", Attenuated = "#9E9E9E")

# --- 5. Build waterfall ---
pF <- ggplot(fw, aes(x = delta_NES, y = pathway_clean)) +
  # Delta NES bars
  geom_col(aes(fill = pw_cat), width = 0.7, color = "black", linewidth = 0.15) +
  scale_fill_manual(values = PF_COLORS, name = "Category") +
  # Individual NES dots: Aging (green circle), Training_Old (blue triangle)
  geom_point(aes(x = NES_Aging), color = "#4CAF50",
             size = 1.5, shape = 16, alpha = 0.8) +
  geom_point(aes(x = NES_Training_Old), color = "#5DA5DA",
             size = 1.5, shape = 17, alpha = 0.8) +
  # Significance stars at bar tips
  geom_text(aes(x = delta_NES + sign(delta_NES) * 0.05, label = int_stars),
            hjust = ifelse(fw$delta_NES > 0, 0, 1),
            size = 2.2, fontface = "bold", color = "grey20") +
  # Reference line
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  labs(
    title = "Pathway-Level Reversal Magnitude",
    subtitle = expression(Delta*"NES (Aging - Training Old), Hallmark + GO:BP"),
    x = expression(Delta~NES),
    y = NULL
  ) +
  THEME_PUB +
  theme(axis.text.y = element_text(size = 5),
        legend.position = "bottom",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6))

# --- 6. Export ---
write_csv(fw, file.path(DAT_DIR, "fig3_reversal_waterfall.csv"))
ggsave(file.path(RPT_DIR, "test_panelF.pdf"), pF,
       width = 180, height = 160, units = "mm")
message(sprintf("Panel F: %d pathways — Reversed=%d, Exacerbated=%d, Attenuated=%d",
                nrow(fw),
                sum(fw$pw_cat == "Reversed"),
                sum(fw$pw_cat == "Exacerbated"),
                sum(fw$pw_cat == "Attenuated")))
message("Panel F test saved")

# === 14. FINAL ASSEMBLY — 6-Panel Composite Figure =========================

message("Assembling Figure 3...")

pA_wrapped <- wrap_elements(full = pA)

fig3 <- (pA_wrapped | pC_gg) /
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

# --- Save with fallback ---
fig3_pdf <- file.path(RPT_DIR, "Figure_3.pdf")
fig3_png <- file.path(RPT_DIR, "Figure_3.png")

tryCatch({
  ggsave(fig3_pdf, fig3,
         width = 380, height = 500, units = "mm", limitsize = FALSE)
  message("Figure 3 PDF saved via ggsave")
}, error = function(e) {
  message("ggsave PDF failed: ", e$message, " — using pdf() fallback")
  pdf(fig3_pdf, width = 380/25.4, height = 500/25.4)
  print(fig3)
  dev.off()
  message("Figure 3 PDF saved via pdf() device")
})

tryCatch({
  ggsave(fig3_png, fig3,
         width = 380, height = 500, units = "mm", dpi = 300, limitsize = FALSE)
  message("Figure 3 PNG saved via ggsave")
}, error = function(e) {
  message("ggsave PNG failed: ", e$message, " — using png() fallback")
  png(fig3_png, width = 380, height = 500, units = "mm", res = 300)
  print(fig3)
  dev.off()
  message("Figure 3 PNG saved via png() device")
})

message("Figure 3 saved to: ", RPT_DIR)
message("Figure 3 pipeline complete.")
