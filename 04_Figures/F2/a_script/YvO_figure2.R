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

# --- Volcano unified direction coloring ---
VOLC_COLORS <- list(
  Training_Young = c(Up = "#D6604D", Down = "#4393C3", NS = "grey80"),
  Training_Old   = c(Up = "#D6604D", Down = "#4393C3", NS = "grey80")
)

# --- Panel E category colors ---
INTERACTION_CAT_COLORS <- c(
  "Attenuated"          = "#81C784",
  "Up Young / Down Old" = "#C62828",
  "Down Young / Up Old" = "#1565C0"
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
    # --- Canonical replacements (from figure-style guide) ---
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
    # --- Additional biological abbreviations ---
    str_replace("Ifn", "IFN") %>%
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

# ═══ 8. PANEL A — Side-by-Side Volcanos with NES Bar Insets ════════════════

message("Building Panel A: side-by-side volcanos...")

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
  strip_title <- if (ctr == "Training_Young") {
    "Training Response in Young Adults"
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

# --- Panel A assembly: shared x-axis, per-contrast y-axis ---
volc_xlim <- max(abs(c(dep_df$logFC_Training_Young, dep_df$logFC_Training_Old)),
                 na.rm = TRUE) * 1.05

volc_ylim_young <- max(-log10(dep_df$P.Value_Training_Young), na.rm = TRUE)
volc_ylim_young <- min(volc_ylim_young, 15)

volc_ylim_old <- max(-log10(dep_df$P.Value_Training_Old), na.rm = TRUE)
volc_ylim_old <- min(volc_ylim_old * 1.1, 15)   # tighter for Training_Old

pA_left  <- make_volcano("Training_Young") +
  coord_cartesian(xlim = c(-volc_xlim, volc_xlim), ylim = c(0, volc_ylim_young))
pA_right <- make_volcano("Training_Old") +
  coord_cartesian(xlim = c(-volc_xlim, volc_xlim), ylim = c(0, volc_ylim_old)) +
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
axis_lim <- quantile(abs(c(scatter_df$logFC_Y, scatter_df$logFC_O)), 0.98, na.rm = TRUE) * 1.1

# --- 6. Build scatter plot ---
# Order data so NS is plotted first, Interaction on top
scatter_ordered <- scatter_df %>%
  mutate(plot_order = as.integer(sig_cat)) %>%
  arrange(desc(plot_order))

pB_base <- ggplot(scatter_ordered, aes(x = logFC_Y, y = logFC_O)) +
  # Quadrant background shading (canonical: salmon=concordant, blue=discordant)
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "#E88D6D", alpha = 0.18) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#E88D6D", alpha = 0.18) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "#7BAFD4", alpha = 0.18) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#7BAFD4", alpha = 0.18) +
  # Reference lines
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  # Points — NS layer (small, faded)
  geom_point(data = . %>% filter(sig_cat == "NS"),
             aes(color = sig_cat), size = 0.4, alpha = 0.15) +
  # Points — significant layer (uniform, prominent)
  geom_point(data = . %>% filter(sig_cat != "NS"),
             aes(color = sig_cat), size = 1.5, alpha = 0.85) +
  scale_color_manual(values = SIG_COLORS, name = "Significance") +
  # Quadrant labels with counts (boxed, white bold)
  annotate("label", x = axis_lim * 0.95, y = axis_lim * 0.95,
           label = paste0("Concordant Up\nn = ", q_counts["Q1"]),
           hjust = 1, vjust = 1, size = 2.0, fontface = "bold",
           color = "white", fill = alpha("#E88D6D", 0.7),
           label.padding = unit(2, "pt"), label.size = 0) +
  annotate("label", x = -axis_lim * 0.95, y = -axis_lim * 0.95,
           label = paste0("Concordant Down\nn = ", q_counts["Q3"]),
           hjust = 0, vjust = 0, size = 2.0, fontface = "bold",
           color = "white", fill = alpha("#E88D6D", 0.7),
           label.padding = unit(2, "pt"), label.size = 0) +
  annotate("label", x = -axis_lim * 0.95, y = axis_lim * 0.95,
           label = paste0("Discordant\nn = ", q_counts["Q2"]),
           hjust = 0, vjust = 1, size = 2.0, fontface = "bold",
           color = "white", fill = alpha("#7BAFD4", 0.7),
           label.padding = unit(2, "pt"), label.size = 0) +
  annotate("label", x = axis_lim * 0.95, y = -axis_lim * 0.95,
           label = paste0("Discordant\nn = ", q_counts["Q4"]),
           hjust = 1, vjust = 0, size = 2.0, fontface = "bold",
           color = "white", fill = alpha("#7BAFD4", 0.7),
           label.padding = unit(2, "pt"), label.size = 0) +
  labs(
    title = "Protein-Level Concordance of Training Response",
    subtitle = sprintf("r = %.2f [%.2f, %.2f] | rho = %.2f | Sign concordance: %.0f%%",
                       cor_r, cor_r_ci[1], cor_r_ci[2], cor_rho, sign_concordance),
    x = expression(log[2]*FC ~ "(Training Young)"),
    y = expression(log[2]*FC ~ "(Training Old)")
  ) +
  coord_cartesian(xlim = c(-axis_lim, axis_lim),
                  ylim = c(-axis_lim, axis_lim)) +
  THEME_PUB +
  guides(color = guide_legend(override.aes = list(size = c(2.5, 2.0, 1.5, 1.5, 0.8),
                                                    alpha = c(0.9, 0.85, 0.7, 0.7, 0.3)))) +
  theme(legend.position = "bottom",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 5))

# --- 7. Protein labels per category (must add before insets) ---
label_df <- scatter_df %>%
  filter(sig_cat != "NS") %>%
  group_by(sig_cat) %>%
  arrange(desc(abs(logFC_Y) + abs(logFC_O))) %>%
  slice_head(n = 6) %>%
  ungroup()

pB_base <- pB_base +
  geom_label_repel(
    data = label_df, aes(label = gene, fill = sig_cat),
    color = "white", fontface = "bold", size = KEY_TEXT,
    max.overlaps = 20, segment.size = 0.2,
    label.padding = unit(1.5, "pt"), label.size = 0,
    min.segment.length = 0, show.legend = FALSE
  ) +
  scale_fill_manual(values = SIG_COLORS, guide = "none")

# --- 8. Per-quadrant ORA bar insets ---
make_ora_inset <- function(q_terms, direction, fill_color) {
  if (is.null(q_terms) || nrow(q_terms) == 0 ||
      all(is.na(q_terms$p.adjust))) return(NULL)
  q_terms <- q_terms %>%
    filter(!is.na(p.adjust)) %>%
    mutate(neg_log10_padj = -log10(p.adjust),
           stars = sig_stars(p.adjust),
           label = paste0(Description, " ", stars))
  if (nrow(q_terms) == 0) return(NULL)

  min_width <- max(q_terms$neg_log10_padj) * 0.5
  q_terms <- q_terms %>%
    mutate(bar_width = pmax(neg_log10_padj, min_width))

  if (direction == "left") {
    q_terms <- q_terms %>% mutate(xmin = -bar_width, xmax = 0)
  } else {
    q_terms <- q_terms %>% mutate(xmin = 0, xmax = bar_width)
  }

  q_terms <- q_terms %>% mutate(y_pos = row_number())

  ggplot(q_terms) +
    geom_rect(aes(xmin = xmin, xmax = xmax,
                  ymin = y_pos - 0.4, ymax = y_pos + 0.4),
              fill = fill_color, alpha = 0.7, color = "black", linewidth = 0.1) +
    geom_text(aes(x = (xmin + xmax) / 2, y = y_pos, label = label),
              color = "white", fontface = "bold", size = 1.5, hjust = 0.5) +
    theme_void() +
    theme(plot.background = element_blank(),
          plot.margin = margin(0, 0, 0, 0))
}

# Build insets for each quadrant
q_insets <- list(
  Q1 = make_ora_inset(qora_results %>% filter(quadrant == "Q1"), "left", "#E88D6D"),
  Q2 = make_ora_inset(qora_results %>% filter(quadrant == "Q2"), "right", "#7BAFD4"),
  Q3 = make_ora_inset(qora_results %>% filter(quadrant == "Q3"), "right", "#E88D6D"),
  Q4 = make_ora_inset(qora_results %>% filter(quadrant == "Q4"), "left", "#7BAFD4")
)

# Attach insets to quadrant corners (converts pB_base to patchwork)
if (!is.null(q_insets$Q1))
  pB_base <- pB_base + inset_element(q_insets$Q1,
    left = 0.55, right = 1.0, bottom = 0.75, top = 1.0)
if (!is.null(q_insets$Q2))
  pB_base <- pB_base + inset_element(q_insets$Q2,
    left = 0.0, right = 0.45, bottom = 0.75, top = 1.0)
if (!is.null(q_insets$Q3))
  pB_base <- pB_base + inset_element(q_insets$Q3,
    left = 0.0, right = 0.45, bottom = 0.0, top = 0.25)
if (!is.null(q_insets$Q4))
  pB_base <- pB_base + inset_element(q_insets$Q4,
    left = 0.55, right = 1.0, bottom = 0.0, top = 0.25)

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
max_UU <- max(hmat[1:mid_r, 1:mid_c], na.rm = TRUE)             # bottom-left in image() = Concordant up-up
max_DD <- max(hmat[(mid_r+1):nr, (mid_c+1):nc], na.rm = TRUE)   # top-right in image() = Concordant dn-dn
max_UD <- max(hmat[1:mid_r, (mid_c+1):nc], na.rm = TRUE)        # top-left in image() = Discordant Y↑O↓
max_DU <- max(hmat[(mid_r+1):nr, 1:mid_c], na.rm = TRUE)        # bottom-right in image() = Discordant Y↓O↑

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
  # Quadrant labels — corrected positions matching image() display
  # Concordant up-up → bottom-left (small rows, small cols)
  annotate("label", x = 0.25, y = 0.25,
           label = paste0("Concordant\nBoth Up\nmax = ", round(max_UU, 1)),
           fill = alpha("white", 0.7), size = 2.0, fontface = "bold") +
  # Concordant dn-dn → top-right (large rows, large cols)
  annotate("label", x = 0.75, y = 0.75,
           label = paste0("Concordant\nBoth Down\nmax = ", round(max_DD, 1)),
           fill = alpha("white", 0.7), size = 2.0, fontface = "bold") +
  # Discordant Y↑O↓ → top-left (small rows, large cols)
  annotate("label", x = 0.25, y = 0.75,
           label = paste0("Discordant\nY Up / O Dn\nmax = ", round(max_UD, 1)),
           fill = alpha("white", 0.7), size = 2.0, fontface = "bold") +
  # Discordant Y↓O↑ → bottom-right (large rows, small cols)
  annotate("label", x = 0.75, y = 0.25,
           label = paste0("Discordant\nY Dn / O Up\nmax = ", round(max_DU, 1)),
           fill = alpha("white", 0.7), size = 2.0, fontface = "bold") +
  # Axis labels with direction arrows
  annotate("text", x = 0.5, y = -0.02,
           label = "Training (Young) rank", size = 2.2, color = "grey20") +
  annotate("text", x = 0.95, y = -0.06, label = "Most upregulated ->",
           hjust = 1, size = 1.8, color = "grey30") +
  annotate("text", x = 0.05, y = -0.06, label = "<- Most downregulated",
           hjust = 0, size = 1.8, color = "grey30") +
  annotate("text", x = -0.02, y = 0.5, angle = 90,
           label = "Training (Old) rank", size = 2.2, color = "grey20") +
  annotate("text", x = -0.06, y = 0.95, angle = 90,
           label = "Most upregulated ->",
           hjust = 1, size = 1.8, color = "grey30") +
  annotate("text", x = -0.06, y = 0.05, angle = 90,
           label = "<- Most downregulated",
           hjust = 0, size = 1.8, color = "grey30") +
  labs(title = "RRHO2 Concordance Map",
       subtitle = "-log10(p) hypergeometric overlap") +
  coord_cartesian(xlim = c(-0.1, 1.05), ylim = c(-0.1, 1.05), clip = "off") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5),
        plot.subtitle = element_text(size = 6.5, color = "grey30", hjust = 0.5, face = "italic"))

# --- 7. Test save ---
pdf(file.path(RPT_DIR, "test_panelC.pdf"),
    width = 170 / 25.4, height = 130 / 25.4)
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
# Exclude GO terms clearly irrelevant to skeletal muscle proteomics
exclude_keywords <- c("sperm", "zona pellucida", "spermat", "oocyte",
                       "embryonic", "placenta", "retina", "olfactory",
                       "photoreceptor", "taste", "pollen", "fertiliz")

label_keywords <- c("ribosom", "oxphos", "oxidative phosph", "mtorc1",
                     "extracellular matrix", "myogenes", "glycoly",
                     "proteasome", "translation", "unfolded protein",
                     "mitotic spindle", "muscle", "respiratory chain",
                     "fatty acid", "inflammatory")

# Helper: filter out tissue-irrelevant GO terms
filter_relevant <- function(df) {
  df %>% filter(!str_detect(tolower(set), paste(exclude_keywords, collapse = "|")))
}

# Label from significant pathways + keyword matches
label_df <- bind_rows(
  # Keyword matches among significant
  mitch_df %>% filter(sig_cat != "NS") %>%
    filter(str_detect(tolower(set), paste(label_keywords, collapse = "|"))),
  # Top 3 per quadrant by effect size (excluding irrelevant terms)
  mitch_df %>% filter(sig_cat != "NS", s.Training_Young > 0, s.Training_Old > 0) %>%
    filter_relevant() %>%
    slice_max(abs(s.Training_Young) + abs(s.Training_Old), n = 3),
  mitch_df %>% filter(sig_cat != "NS", s.Training_Young < 0, s.Training_Old < 0) %>%
    filter_relevant() %>%
    slice_max(abs(s.Training_Young) + abs(s.Training_Old), n = 3),
  mitch_df %>% filter(sig_cat != "NS", s.Training_Young > 0, s.Training_Old < 0) %>%
    filter_relevant() %>%
    slice_max(abs(s.Training_Young) + abs(s.Training_Old), n = 3),
  mitch_df %>% filter(sig_cat != "NS", s.Training_Young < 0, s.Training_Old > 0) %>%
    filter_relevant() %>%
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
  # Quadrant background shading (canonical: salmon=concordant, blue=discordant)
  annotate("rect", xmin = 0, xmax = pw_lim, ymin = 0, ymax = pw_lim,
           fill = "#E88D6D", alpha = 0.18) +
  annotate("rect", xmin = -pw_lim, xmax = 0, ymin = -pw_lim, ymax = 0,
           fill = "#E88D6D", alpha = 0.18) +
  annotate("rect", xmin = 0, xmax = pw_lim, ymin = -pw_lim, ymax = 0,
           fill = "#7BAFD4", alpha = 0.18) +
  annotate("rect", xmin = -pw_lim, xmax = 0, ymin = 0, ymax = pw_lim,
           fill = "#7BAFD4", alpha = 0.18) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40", linewidth = 0.3) +
  # Points with 4-category + size by gene-set size
  geom_point(aes(color = sig_cat, size = setSize, alpha = sig_cat)) +
  scale_color_manual(values = SIG_COLORS, name = "Significance") +
  scale_alpha_manual(values = SIG_ALPHAS, guide = "none") +
  scale_size_continuous(range = c(0.8, 5.0),
                        breaks = seq(50, 500, by = 50),
                        name = "Protein set size",
                        guide = guide_legend(override.aes = list(alpha = 0.8))) +
  # Labels — colored box labels
  geom_label_repel(data = label_df, aes(label = pathway_clean, fill = sig_cat),
                   color = "white", fontface = "bold", size = KEY_TEXT,
                   max.overlaps = 25, segment.size = 0.2,
                   label.padding = unit(1.5, "pt"), label.size = 0,
                   min.segment.length = 0, show.legend = FALSE) +
  scale_fill_manual(values = SIG_COLORS, guide = "none") +
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

# ═══ 12. PANEL E — Interaction DEP Classification (Vertical Bar + Trajectories)

message("Building Panel E: interaction DEP classification...")

# --- 1. Classify interaction DEPs into 3 categories ---
int_df <- dep_df %>%
  filter(pi_score_Interaction < 0.05) %>%
  dplyr::select(gene,
         logFC_Y = logFC_Training_Young, logFC_O = logFC_Training_Old,
         pi_Y = pi_score_Training_Young, pi_O = pi_score_Training_Old) %>%
  filter(!is.na(logFC_Y), !is.na(logFC_O)) %>%
  mutate(
    same_dir = sign(logFC_Y) == sign(logFC_O),
    category = case_when(
      !same_dir & logFC_Y > 0 ~ "Up Young / Down Old",
      !same_dir & logFC_Y < 0 ~ "Down Young / Up Old",
      TRUE                     ~ "Attenuated"
    ),
    category = factor(category,
      levels = c("Attenuated", "Up Young / Down Old", "Down Young / Up Old"))
  )

cat_counts <- int_df %>% count(category, .drop = FALSE) %>% deframe()
cat_pcts   <- round(100 * cat_counts / sum(cat_counts), 0)
n_total    <- nrow(int_df)

message(sprintf("Panel E: %d interaction DEPs -- %s",
                n_total,
                paste(names(cat_counts), cat_counts, sep = "=", collapse = ", ")))

# --- 2. Build stacked bar data ---
bar_df <- tibble(
  category = factor(names(cat_counts), levels = levels(int_df$category)),
  count    = as.integer(cat_counts),
  pct      = as.numeric(cat_pcts)
)

# --- 3. Build vertical stacked bar plot ---
pE_bar <- ggplot(bar_df, aes(x = 1, y = count, fill = category)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.3) +
  geom_text(aes(label = paste0(count, "\n(", pct, "%)")),
            position = position_stack(vjust = 0.5),
            size = 2.5, fontface = "bold", color = "white") +
  scale_fill_manual(values = INTERACTION_CAT_COLORS) +
  coord_cartesian(xlim = c(0.3, 1.7)) +
  labs(title = sprintf("Interaction DEPs\n(n = %d)", n_total), y = "Count") +
  THEME_PUB +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 5),
        legend.key.size = unit(3, "mm"))

# --- 4. Pathway trajectory divergence plots ---
make_trajectory_plot <- function(pathway_name, int_df, hallmark_t2g, dep_df) {
  gs_key <- paste0("HALLMARK_", toupper(str_replace_all(pathway_name, " ", "_")))
  pathway_genes <- hallmark_t2g %>%
    filter(gs_name == gs_key) %>%
    pull(gene_symbol) %>% unique()

  traj_df <- dep_df %>%
    filter(gene %in% pathway_genes) %>%
    dplyr::select(gene, logFC_Y = logFC_Training_Young, logFC_O = logFC_Training_Old) %>%
    filter(!is.na(logFC_Y), !is.na(logFC_O)) %>%
    mutate(is_int_dep = gene %in% int_df$gene)

  if (nrow(traj_df) == 0) return(ggplot() + theme_void())

  # Build long format: pre=0, post=logFC for each group
  traj_long <- traj_df %>%
    pivot_longer(cols = c(logFC_Y, logFC_O),
                 names_to = "group_raw", values_to = "post_val") %>%
    mutate(
      group = ifelse(group_raw == "logFC_Y", "Young", "Old"),
      pre_val = 0
    ) %>%
    pivot_longer(cols = c(pre_val, post_val),
                 names_to = "time", values_to = "value") %>%
    mutate(time = factor(ifelse(time == "pre_val", "Pre", "Post"),
                         levels = c("Pre", "Post")))

  ggplot(traj_long, aes(x = time, y = value,
                         group = interaction(gene, group),
                         color = group)) +
    geom_line(aes(alpha = is_int_dep, linewidth = is_int_dep)) +
    geom_point(size = 0.8) +
    scale_color_manual(values = c(Young = "#E05A4E", Old = "#5DA5DA")) +
    scale_alpha_manual(values = c(`TRUE` = 0.8, `FALSE` = 0.25), guide = "none") +
    scale_linewidth_manual(values = c(`TRUE` = 0.8, `FALSE` = 0.3), guide = "none") +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.2) +
    labs(title = pathway_name,
         subtitle = sprintf("%d proteins (%d interaction DEPs)",
                            length(unique(traj_df$gene)),
                            sum(traj_df$is_int_dep)),
         x = NULL, y = expression(log[2]~FC)) +
    THEME_PUB +
    theme(legend.position = "none",
          plot.title = element_text(size = 7, face = "bold"),
          plot.subtitle = element_text(size = 5.5, face = "italic"),
          axis.text = element_text(size = 5),
          axis.title.y = element_text(size = 5))
}

pE_myogen <- make_trajectory_plot("Myogenesis", int_df, hallmark_t2g, dep_df)
pE_emt    <- make_trajectory_plot("Epithelial Mesenchymal Transition", int_df, hallmark_t2g, dep_df) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 6),
        legend.key.size = unit(3, "mm")) +
  labs(color = "Age Group")

# --- 5. Assemble Panel E: bar (left 35%) | trajectories (right 65%) ---
pE <- (pE_bar | (pE_myogen / pE_emt)) +
  plot_layout(widths = c(0.35, 0.65))

# --- 6. Export ---
write_csv(int_df, file.path(DAT_DIR, "fig2_interaction_classification.csv"))
ggsave(file.path(RPT_DIR, "test_panelE.pdf"), pE,
       width = 250, height = 120, units = "mm")
message("Panel E test saved")

# ═══ 13. PANEL F — Dumbbell Dot Plot: Concordant/Attenuated/Discordant ═══════

message("Building Panel F: pathway response dumbbell plot...")

# --- 1. Extract Hallmark fgsea results per contrast ---
hallmark_fgsea <- fgsea_all %>%
  filter(database == "Hallmark") %>%
  dplyr::select(pathway, contrast, NES, padj)

# --- 2. Pivot to wide format (one row per pathway) ---
hw <- hallmark_fgsea %>%
  pivot_wider(id_cols = pathway, names_from = contrast,
              values_from = c(NES, padj))

# --- 3. Classify pathways ---
hw <- hw %>%
  mutate(
    sig_Y = !is.na(padj_Training_Young) & padj_Training_Young < 0.05,
    sig_O = !is.na(padj_Training_Old) & padj_Training_Old < 0.05,
    sig_I = !is.na(padj_Interaction) & padj_Interaction < 0.25,
    same_sign = sign(NES_Training_Young) == sign(NES_Training_Old),
    pw_cat = case_when(
      sig_I                           ~ "Discordant",
      sig_Y & sig_O & same_sign       ~ "Concordant",
      sig_Y & sig_O & !same_sign      ~ "Discordant",
      sig_Y | sig_O                   ~ "Attenuated",
      TRUE                            ~ "NS"
    ),
    pw_cat = factor(pw_cat, levels = c("Concordant", "Attenuated", "Discordant", "NS")),
    pathway_clean = clean_pathway_name(pathway),
    max_nes = pmax(abs(NES_Training_Young), abs(NES_Training_Old), na.rm = TRUE)
  ) %>%
  filter(pw_cat != "NS")

# Replace NA NES values with 0 for plotting
hw <- hw %>%
  mutate(
    NES_Training_Young = replace_na(NES_Training_Young, 0),
    NES_Training_Old   = replace_na(NES_Training_Old, 0)
  )

message(sprintf("  Panel F: Concordant=%d, Attenuated=%d, Discordant=%d pathways",
                sum(hw$pw_cat == "Concordant"),
                sum(hw$pw_cat == "Attenuated"),
                sum(hw$pw_cat == "Discordant")))

# --- 4. Build dumbbell dot plot ---
PF_COLORS <- c(Concordant = "#00897B", Attenuated = "#FFA726", Discordant = "#C62828")

pF <- ggplot(hw, aes(y = reorder(pathway_clean, max_nes))) +
  # Connecting segments
  geom_segment(aes(x = NES_Training_Young, xend = NES_Training_Old,
                   yend = reorder(pathway_clean, max_nes),
                   color = pw_cat), linewidth = 0.4) +
  # Circle = Young
  geom_point(aes(x = NES_Training_Young, color = pw_cat, size = max_nes),
             alpha = 0.8, shape = 16) +
  # Triangle = Old
  geom_point(aes(x = NES_Training_Old, color = pw_cat, size = max_nes),
             alpha = 0.8, shape = 17) +
  scale_color_manual(values = PF_COLORS, name = "Category") +
  scale_size_continuous(range = c(1.5, 4), guide = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  labs(x = "NES", y = NULL,
       title = "Pathway Response by Age Group",
       subtitle = "Circle = Young, Triangle = Old") +
  THEME_PUB +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 5),
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6))

# --- 5. Export ---
write_csv(hw, file.path(DAT_DIR, "fig2_concordant_discordant_enrichment.csv"))
ggsave(file.path(RPT_DIR, "test_panelF.pdf"), pF,
       width = 180, height = 160, units = "mm")
message("Panel F test saved")

# ═══ 14. FINAL ASSEMBLY — 6-Panel Composite Figure ══════════════════════════

message("Assembling Figure 2...")

pA_wrapped <- wrap_elements(full = pA)
pE_wrapped <- wrap_elements(full = pE)

fig2 <- (pA_wrapped | pC_gg) /
         (pB_base   | pD) /
         (pE_wrapped | pF) +
  plot_layout(
    widths  = c(0.55, 0.45),
    heights = c(0.28, 0.37, 0.35)
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 12))
  )

# --- Save with fallback ---
fig2_pdf <- file.path(RPT_DIR, "Figure_2.pdf")
fig2_png <- file.path(RPT_DIR, "Figure_2.png")

tryCatch({
  ggsave(fig2_pdf, fig2,
         width = 380, height = 520, units = "mm", limitsize = FALSE)
  message("Figure 2 PDF saved via ggsave")
}, error = function(e) {
  message("ggsave PDF failed: ", e$message, " — using pdf() fallback")
  pdf(fig2_pdf, width = 380/25.4, height = 520/25.4)
  print(fig2)
  dev.off()
  message("Figure 2 PDF saved via pdf() device")
})

tryCatch({
  ggsave(fig2_png, fig2,
         width = 380, height = 520, units = "mm", dpi = 300, limitsize = FALSE)
  message("Figure 2 PNG saved via ggsave")
}, error = function(e) {
  message("ggsave PNG failed: ", e$message, " — using png() fallback")
  png(fig2_png, width = 380, height = 520, units = "mm", res = 300)
  print(fig2)
  dev.off()
  message("Figure 2 PNG saved via png() device")
})

message("Figure 2 saved to: ", RPT_DIR)
