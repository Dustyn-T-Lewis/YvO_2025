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
#   D — fGSEA NES scatter (pathway-level concordance)
#
# References:
#   Cahill et al. 2018, Scientific Reports 8:9588 (RRHO2)
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
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(fgsea)
  library(ggExtra)
  library(png)
})

source("04_Figures/F2/a_script/volcano_ring.R")

# ═══ 2. SEED ═════════════════════════════════════════════════════════════════

set.seed(42)

# ═══ 3. PATH RESOLUTION ═════════════════════════════════════════════════════

setwd("/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_YvO_2025")
RPT_DIR <- "04_Figures/F2/b_reports"
DAT_DIR <- "04_Figures/F2/c_data"
SUP_DIR <- file.path(RPT_DIR, "supplementary")

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
  "NS"                   = "grey55"
)
SIG_SIZES  <- c(Interaction = 2.5, `Sig Both` = 2.0,
                `Sig Young only` = 1.5, `Sig Old only` = 1.5, NS = 0.4)
SIG_ALPHAS <- c(Interaction = 0.90, `Sig Both` = 0.85,
                `Sig Young only` = 0.70, `Sig Old only` = 0.70, NS = 0.30)

# --- Volcano unified direction coloring ---
VOLC_COLORS <- list(
  Training_Young = c(Up = "#D6604D", Down = "#4393C3", NS = "grey80"),
  Training_Old   = c(Up = "#D6604D", Down = "#4393C3", NS = "grey80")
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
dep_df <- read_csv("03_DEP/c_data/combined_results.csv",
                   show_col_types = FALSE)
# Validation
stopifnot(nrow(dep_df) > 2000)
stopifnot(all(c("gene", "logFC_Training_Young", "t_Training_Young",
                "sig_pi_Interaction") %in% names(dep_df)))

# --- fGSEA computation (or cache read) ---
fgsea_cache <- file.path(DAT_DIR, "fgsea_tstat_all_v2.csv")

if (file.exists(fgsea_cache)) {
  message("Reading cached fGSEA results: ", fgsea_cache)
  fgsea_all <- read_csv(fgsea_cache, show_col_types = FALSE)
} else {
  message("Computing fGSEA across all contrasts and gene-set databases...")

  # 1. Gene-set collections from MSigDB
  gs_collections <- list(
    Hallmark = msigdbr(species = "Homo sapiens", collection = "H"),
    GOBP     = msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"),
    GOCC     = msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:CC"),
    GOMF     = msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:MF"),
    Reactome = msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
  )

  # Database label mapping (GOBP -> GO:BP, etc.)
  db_labels <- c(Hallmark = "Hallmark", GOBP = "GO:BP", GOCC = "GO:CC",
                 GOMF = "GO:MF", Reactome = "Reactome")

  # Convert to named-list format for fgsea
  gs_lists <- lapply(gs_collections, function(db) {
    split(db$gene_symbol, db$gs_name)
  })

  # 2. Contrasts and their t-statistic columns
  contrasts <- c("Training_Young", "Training_Old", "Aging", "Interaction")

  # 3. Run fgsea for each contrast x database combination
  fgsea_results <- list()

  for (ctr in contrasts) {
    t_col <- paste0("t_", ctr)
    if (!t_col %in% names(dep_df)) {
      message("  Skipping contrast ", ctr, " — column ", t_col, " not found")
      next
    }

    # Build named ranked vector (gene -> t-stat), drop NAs and duplicates
    rank_df <- dep_df %>%
      dplyr::select(gene, t_val = all_of(t_col)) %>%
      dplyr::filter(!is.na(t_val)) %>%
      dplyr::distinct(gene, .keep_all = TRUE)
    ranks <- setNames(rank_df$t_val, rank_df$gene)

    for (db_name in names(gs_lists)) {
      message(sprintf("  fgsea: %s x %s (%d gene sets)",
                      ctr, db_name, length(gs_lists[[db_name]])))
      res <- fgsea::fgseaMultilevel(
        pathways  = gs_lists[[db_name]],
        stats     = ranks,
        minSize   = 15,
        maxSize   = 200,
        nPermSimple = 10000,
        eps       = 0
      )
      if (nrow(res) > 0) {
        res$contrast <- ctr
        res$database <- db_labels[[db_name]]
        # Collapse leadingEdge list column to semicolon-delimited string
        res$leadingEdge <- sapply(res$leadingEdge, paste, collapse = ";")
        fgsea_results[[paste0(ctr, "_", db_name)]] <- res
      }
    }
  }

  fgsea_all <- bind_rows(fgsea_results)
  write_csv(fgsea_all, fgsea_cache)
  message(sprintf("fGSEA complete: %d results written to %s", nrow(fgsea_all), fgsea_cache))
}

message(sprintf("Loaded %d proteins, %d fGSEA results", nrow(dep_df), nrow(fgsea_all)))

# --- Hallmark gene sets in TERM2GENE format for ORA ---
hallmark_t2g <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  as.data.frame()
message(sprintf("Loaded %d Hallmark gene-set mappings", nrow(hallmark_t2g)))

# ═══ 8. PANEL A — Side-by-Side Volcanos with NES Bar Insets ════════════════

message("Building Panel A: side-by-side volcanos...")

# --- NES bar inset builder ---
make_nes_inset <- function(pw_df, direction = "up", max_bars = 5) {
  if (is.null(pw_df) || nrow(pw_df) == 0) return(NULL)
  pw_df <- pw_df %>%
    slice_head(n = max_bars) %>%
    mutate(stars = sig_stars(padj),
           label = clean_pathway_name(pathway))

  # Width proportional to |NES|, with minimum for visibility
  max_nes <- max(abs(pw_df$NES), na.rm = TRUE)
  pw_df <- pw_df %>%
    mutate(bar_width = pmax(abs(NES), 0.3))

  x_extent <- max_nes + 0.3

  if (direction == "down") {
    # Down bars left-anchored
    pw_df <- pw_df %>% mutate(xmin = 0, xmax = bar_width)
  } else {
    # Up bars right-anchored
    pw_df <- pw_df %>% mutate(xmin = x_extent - bar_width, xmax = x_extent)
  }

  pw_df <- pw_df %>% mutate(y_pos = row_number())

  p <- ggplot(pw_df) +
    geom_rect(aes(xmin = xmin, xmax = xmax,
                  ymin = y_pos - 0.4, ymax = y_pos + 0.4,
                  fill = NES),
              color = "black", linewidth = 0.15) +
    scale_fill_gradient2(low = "#4393C3", mid = "grey95", high = "#D6604D",
                         midpoint = 0, guide = "none") +
    scale_y_continuous(limits = c(0.5, max_bars + 0.5), breaks = NULL)

  if (direction == "down") {
    # Text outside bars (right side), stars at bar tip
    p <- p +
      geom_text(aes(x = xmax + 0.05, y = y_pos, label = label),
                size = 1.5, color = "black", fontface = "plain", hjust = 0) +
      geom_text(aes(x = xmax + 0.02, y = y_pos, label = stars),
                size = 1.8, color = "grey20", fontface = "bold", hjust = 0)
  } else {
    # Text outside bars (left side), stars at bar tip
    p <- p +
      geom_text(aes(x = xmin - 0.05, y = y_pos, label = label),
                size = 1.5, color = "black", fontface = "plain", hjust = 1) +
      geom_text(aes(x = xmin - 0.02, y = y_pos, label = stars),
                size = 1.8, color = "grey20", fontface = "bold", hjust = 1)
  }

  p +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.background = element_blank(),
          plot.margin = margin(1, 1, 1, 1, "mm"))
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

  # FDR count
  adj_col <- paste0("adj.P.Val_", ctr)
  n_fdr <- if (adj_col %in% names(dep_df)) sum(dep_df[[adj_col]] < 0.05, na.rm = TRUE) else NA

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
    filter(NES > 0) %>% arrange(desc(NES)) %>% slice_head(n = 5) %>%
    mutate(label = clean_pathway_name(pathway))

  pw_down <- pw_df %>%
    filter(NES < 0) %>% arrange(NES) %>% slice_head(n = 5) %>%
    mutate(label = clean_pathway_name(pathway))

  y_max_est <- max(vdf$neg_log10p, na.rm = TRUE)
  x_max_est <- max(abs(vdf$logFC), na.rm = TRUE)

  # Strip title
  strip_title <- if (ctr == "Training_Young") {
    "A  Training Response in Young Adults"
  } else {
    "Training Response in Old Adults"
  }
  strip_subtitle <- sprintf("Pi < 0.05 | FDR < 0.05: %s", ifelse(is.na(n_fdr), "N/A", as.character(n_fdr)))

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

# ═══ 8. PANEL A — Volcano-in-Ring Composite ═══════════════════════════════════

message("Building Panel A: volcano-in-ring composites...")

pA <- make_volcano_ring_pair(dep_df, fgsea_all)

# Standalone Panel A export handled by make_volcano_ring_pair()

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

# Interaction DEP breakdown for annotation
int_proteins <- scatter_df %>% filter(sig_cat == "Interaction")
n_int_total <- nrow(int_proteins)
n_attenuated <- sum(sign(int_proteins$logFC_Y) == sign(int_proteins$logFC_O))
n_up_y_down_o <- sum(int_proteins$logFC_Y > 0 & int_proteins$logFC_O < 0)
n_down_y_up_o <- sum(int_proteins$logFC_Y < 0 & int_proteins$logFC_O > 0)
int_annotation <- sprintf("Interaction DEPs: %d\nAttenuated: %d\nUp Y / Down O: %d\nDown Y / Up O: %d",
                           n_int_total, n_attenuated, n_up_y_down_o, n_down_y_up_o)

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
             aes(color = sig_cat), size = 0.4, alpha = 0.4) +
  # Points — significant layer (uniform, prominent)
  geom_point(data = . %>% filter(sig_cat != "NS"),
             aes(color = sig_cat), size = 1.5, alpha = 0.85) +
  scale_color_manual(values = SIG_COLORS, name = "Significance") +
  # Quadrant labels with counts (boxed, white background)
  annotate("label", x = axis_lim * 0.98, y = axis_lim * 0.98,
           label = paste0("Concordant Up\nn = ", q_counts["Q1"]),
           hjust = 1, vjust = 1, size = 2.0, fontface = "bold",
           color = "#E88D6D", fill = alpha("white", 0.85),
           label.padding = unit(2, "pt")) +
  annotate("label", x = -axis_lim * 0.98, y = -axis_lim * 0.98,
           label = paste0("Concordant Down\nn = ", q_counts["Q3"]),
           hjust = 0, vjust = 0, size = 2.0, fontface = "bold",
           color = "#E88D6D", fill = alpha("white", 0.85),
           label.padding = unit(2, "pt")) +
  annotate("label", x = -axis_lim * 0.98, y = axis_lim * 0.98,
           label = paste0("Discordant\nn = ", q_counts["Q2"]),
           hjust = 0, vjust = 1, size = 2.0, fontface = "bold",
           color = "#7BAFD4", fill = alpha("white", 0.85),
           label.padding = unit(2, "pt")) +
  annotate("label", x = axis_lim * 0.98, y = -axis_lim * 0.98,
           label = paste0("Discordant\nn = ", q_counts["Q4"]),
           hjust = 1, vjust = 0, size = 2.0, fontface = "bold",
           color = "#7BAFD4", fill = alpha("white", 0.85),
           label.padding = unit(2, "pt")) +
  labs(
    title = "B  Protein-Level Concordance of Training Response",
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

# --- 7. Interaction annotation + protein labels ---
pB_base <- pB_base +
  annotate("label", x = -axis_lim * 0.95, y = -axis_lim * 0.45,
           label = int_annotation, hjust = 0, vjust = 1, size = 1.8,
           fontface = "bold", color = "#9B7FBF",
           fill = alpha("white", 0.90), label.padding = unit(2, "pt"))

label_df <- scatter_df %>%
  filter(sig_cat != "NS") %>%
  group_by(sig_cat) %>%
  arrange(desc(abs(logFC_Y) + abs(logFC_O))) %>%
  slice_head(n = 5) %>%
  ungroup()

pB_base <- pB_base +
  geom_text_repel(
    data = label_df, aes(label = gene, color = sig_cat),
    size = 2.0, max.overlaps = 20, segment.size = 0.2,
    min.segment.length = 0, show.legend = FALSE
  )

# --- 8. Consolidated ORA inset (bottom-right) ---
make_ora_inset_consolidated <- function(ora_df) {
  if (is.null(ora_df) || nrow(ora_df) == 0 ||
      all(is.na(ora_df$p.adjust))) return(NULL)
  ora_sig <- ora_df %>%
    filter(!is.na(p.adjust), p.adjust < 0.1) %>%
    mutate(neg_log10_padj = -log10(p.adjust),
           stars = sig_stars(p.adjust),
           label = paste0(Description, " ", stars))
  if (nrow(ora_sig) == 0) return(NULL)
  # Top 3 per direction (based on quadrant)
  ora_top <- ora_sig %>%
    arrange(p.adjust) %>%
    slice_head(n = 6) %>%
    mutate(y_pos = row_number(),
           bar_width = pmax(neg_log10_padj, max(neg_log10_padj) * 0.3),
           xmin = 0, xmax = bar_width)
  ggplot(ora_top) +
    geom_rect(aes(xmin = xmin, xmax = xmax,
                  ymin = y_pos - 0.4, ymax = y_pos + 0.4),
              fill = "#E88D6D", alpha = 0.7, color = "black", linewidth = 0.1) +
    geom_text(aes(x = (xmin + xmax) / 2, y = y_pos, label = label),
              color = "white", fontface = "bold", size = 0.8, hjust = 0.5) +
    theme_void() +
    theme(plot.background = element_blank(),
          plot.margin = margin(0, 0, 0, 0))
}

ora_inset <- make_ora_inset_consolidated(qora_results)
if (!is.null(ora_inset))
  pB_base <- pB_base + inset_element(ora_inset,
    left = 0.65, right = 0.98, bottom = 0.02, top = 0.25)

# --- 9. Marginal density distributions ---
density_top <- ggplot(scatter_ordered, aes(x = logFC_Y, fill = sig_cat)) +
  geom_density(alpha = 0.4, linewidth = 0.3) +
  scale_fill_manual(values = SIG_COLORS, guide = "none") +
  coord_cartesian(xlim = c(-axis_lim, axis_lim)) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

density_right <- ggplot(scatter_ordered, aes(x = logFC_O, fill = sig_cat)) +
  geom_density(alpha = 0.4, linewidth = 0.3) +
  scale_fill_manual(values = SIG_COLORS, guide = "none") +
  coord_flip(xlim = c(-axis_lim, axis_lim)) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

empty_corner <- plot_spacer()

pB_compound <- (density_top + empty_corner + pB_base + density_right) +
  plot_layout(widths = c(0.85, 0.15), heights = c(0.15, 0.85))
pB_base <- wrap_elements(full = pB_compound)

# --- 10. Export ---
write_csv(scatter_df, file.path(DAT_DIR, "fig2_concordance_scatter.csv"))
message(sprintf("Panel B: r = %.3f, rho = %.3f, concordance = %.1f%%",
                cor_r, cor_rho, sign_concordance))

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

# --- 5. Render as ggplot with viridis colorscale ---
hmat_df <- expand.grid(row = 1:nr, col = 1:nc) %>%
  mutate(value = as.vector(hmat))

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
  # Quadrant annotations — bold white, with max -log10(p)
  annotate("text", x = mid_c * 0.5, y = mid_r * 0.5,
           label = sprintf("Concordant Up\nmax = %.1f", max_UU),
           color = "white", fontface = "bold", size = 2.5) +
  annotate("text", x = mid_c + (nc - mid_c) * 0.5, y = mid_r + (nr - mid_r) * 0.5,
           label = sprintf("Concordant Down\nmax = %.1f", max_DD),
           color = "white", fontface = "bold", size = 2.5) +
  annotate("text", x = mid_c * 0.5, y = mid_r + (nr - mid_r) * 0.5,
           label = sprintf("Discordant\nY Up / O Down\nmax = %.1f", max_UD),
           color = "white", fontface = "bold", size = 2.0) +
  annotate("text", x = mid_c + (nc - mid_c) * 0.5, y = mid_r * 0.5,
           label = sprintf("Discordant\nY Down / O Up\nmax = %.1f", max_DU),
           color = "white", fontface = "bold", size = 2.0) +
  # Axis labels
  labs(title = "C  Threshold-Free Concordance of Training Response",
       subtitle = "RRHO2 hypergeometric overlap, -log10(P)",
       x = "Training (Young) rank",
       y = "Training (Old) rank") +
  # Direction annotations on axes
  annotate("text", x = 1, y = -nr * 0.04,
           label = "<- Most upregulated",
           hjust = 0, size = 1.8, color = "grey30") +
  annotate("text", x = nc, y = -nr * 0.04,
           label = "Most downregulated ->",
           hjust = 1, size = 1.8, color = "grey30") +
  annotate("text", x = -nc * 0.04, y = 1, angle = 90,
           label = "<- Most upregulated",
           hjust = 0, size = 1.8, color = "grey30") +
  annotate("text", x = -nc * 0.04, y = nr, angle = 90,
           label = "Most downregulated ->",
           hjust = 1, size = 1.8, color = "grey30") +
  coord_cartesian(clip = "off") +
  THEME_PUB +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom")

message(sprintf("RRHO2 quadrant max -log10(p): UU=%.1f DD=%.1f UD=%.1f DU=%.1f",
                max_UU, max_DD, max_UD, max_DU))

# ═══ NEW PANEL D — fGSEA NES Scatter ═════════════════════════════════════════

message("Building Panel D: fGSEA NES scatter...")

# --- 1. Read fGSEA results ---
fgsea_wide <- fgsea_all %>%
  filter(contrast %in% c("Training_Young", "Training_Old", "Interaction")) %>%
  dplyr::select(pathway, contrast, NES, padj, size) %>%
  pivot_wider(id_cols = c(pathway), names_from = contrast,
              values_from = c(NES, padj, size)) %>%
  filter(!is.na(NES_Training_Young), !is.na(NES_Training_Old))

# Use size from Training_Young (or whichever is available)
fgsea_wide <- fgsea_wide %>%
  mutate(set_size = coalesce(size_Training_Young, size_Training_Old))

# --- 2. Filter: significant in at least one contrast OR significant interaction ---
fgsea_sig <- fgsea_wide %>%
  filter(
    (!is.na(padj_Training_Young) & padj_Training_Young < 0.05) |
    (!is.na(padj_Training_Old) & padj_Training_Old < 0.05) |
    (!is.na(padj_Interaction) & padj_Interaction < 0.05)
  )

# --- 3. Classify by significance hierarchy (same as protein scatter) ---
fgsea_sig <- fgsea_sig %>%
  mutate(
    sig_Y = !is.na(padj_Training_Young) & padj_Training_Young < 0.05,
    sig_O = !is.na(padj_Training_Old) & padj_Training_Old < 0.05,
    sig_I = !is.na(padj_Interaction) & padj_Interaction < 0.05,
    sig_cat = case_when(
      sig_I                  ~ "Interaction",
      sig_Y & sig_O          ~ "Sig Both",
      sig_Y                  ~ "Sig Young only",
      sig_O                  ~ "Sig Old only",
      TRUE                   ~ "NS"
    ),
    sig_cat = factor(sig_cat, levels = names(SIG_COLORS)),
    pathway_clean = clean_pathway_name(pathway)
  )

# --- 4. Correlation ---
nes_cor <- cor.test(fgsea_sig$NES_Training_Young, fgsea_sig$NES_Training_Old)
nes_r <- nes_cor$estimate
nes_pval <- nes_cor$p.value

# --- 5. Axis limits ---
nes_lim <- max(abs(c(fgsea_sig$NES_Training_Young, fgsea_sig$NES_Training_Old)),
               na.rm = TRUE) * 1.15

# --- 6. Quadrant counts ---
nq1 <- sum(fgsea_sig$NES_Training_Young > 0 & fgsea_sig$NES_Training_Old > 0)
nq2 <- sum(fgsea_sig$NES_Training_Young < 0 & fgsea_sig$NES_Training_Old > 0)
nq3 <- sum(fgsea_sig$NES_Training_Young < 0 & fgsea_sig$NES_Training_Old < 0)
nq4 <- sum(fgsea_sig$NES_Training_Young > 0 & fgsea_sig$NES_Training_Old < 0)

# --- 7. Top pathway labels (top 4 per category, excluding NS) ---
label_pw <- fgsea_sig %>%
  filter(sig_cat != "NS") %>%
  group_by(sig_cat) %>%
  arrange(desc(abs(NES_Training_Young) + abs(NES_Training_Old))) %>%
  slice_head(n = 4) %>%
  ungroup()

# --- 8. Build scatter ---
pD <- ggplot(fgsea_sig %>% arrange(desc(as.integer(sig_cat))),
             aes(x = NES_Training_Young, y = NES_Training_Old)) +
  # Quadrant shading (same as Panel B)
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
  # Points: circles, colored by significance category, sized by gene set size
  geom_point(aes(color = sig_cat, size = set_size), shape = 16, alpha = 0.75) +
  scale_color_manual(values = SIG_COLORS, name = "Significance") +
  scale_size_continuous(range = c(1.0, 5.0), name = "Set size",
                        breaks = c(20, 50, 100, 200),
                        guide = guide_legend(override.aes = list(alpha = 0.8))) +
  # Pathway labels
  geom_text_repel(data = label_pw, aes(label = pathway_clean, color = sig_cat),
                  size = 2.0, max.overlaps = 20, segment.size = 0.2,
                  min.segment.length = 0, show.legend = FALSE) +
  # Quadrant labels
  annotate("label", x = nes_lim * 0.95, y = nes_lim * 0.95,
           label = paste0("Concordant Up\nn = ", nq1),
           hjust = 1, vjust = 1, size = 2.0, fontface = "bold",
           color = "#E88D6D", fill = alpha("white", 0.85),
           label.padding = unit(2, "pt")) +
  annotate("label", x = -nes_lim * 0.95, y = -nes_lim * 0.95,
           label = paste0("Concordant Down\nn = ", nq3),
           hjust = 0, vjust = 0, size = 2.0, fontface = "bold",
           color = "#E88D6D", fill = alpha("white", 0.85),
           label.padding = unit(2, "pt")) +
  annotate("label", x = -nes_lim * 0.95, y = nes_lim * 0.95,
           label = paste0("Discordant\nn = ", nq2),
           hjust = 0, vjust = 1, size = 2.0, fontface = "bold",
           color = "#7BAFD4", fill = alpha("white", 0.85),
           label.padding = unit(2, "pt")) +
  annotate("label", x = nes_lim * 0.95, y = -nes_lim * 0.95,
           label = paste0("Discordant\nn = ", nq4),
           hjust = 1, vjust = 0, size = 2.0, fontface = "bold",
           color = "#7BAFD4", fill = alpha("white", 0.85),
           label.padding = unit(2, "pt")) +
  # Correlation annotation
  annotate("text", x = -nes_lim * 0.95, y = -nes_lim * 0.75,
           label = sprintf("r = %.2f, p %s", nes_r,
                           ifelse(nes_pval < 0.001, "< 0.001",
                                  sprintf("= %.3f", nes_pval))),
           hjust = 0, size = 2.2, fontface = "italic", color = "grey30") +
  labs(
    title = "D  Pathway-Level Concordance (fGSEA)",
    subtitle = sprintf("NES: Training Young vs Old | padj < 0.05 in >= 1 contrast | %d pathways", nrow(fgsea_sig)),
    x = "NES (Training Young)",
    y = "NES (Training Old)"
  ) +
  coord_fixed(ratio = 1, xlim = c(-nes_lim, nes_lim),
              ylim = c(-nes_lim, nes_lim)) +
  THEME_PUB +
  theme(legend.position = "bottom",
        legend.key.size = unit(2, "mm"),
        legend.text = element_text(size = 5))

write_csv(fgsea_sig, file.path(DAT_DIR, "fig2_fgsea_nes_scatter.csv"))
message(sprintf("Panel D: %d pathways, r = %.3f", nrow(fgsea_sig), nes_r))

# ═══ FINAL ASSEMBLY — 4-Panel Composite Figure ═══════════════════════════════

message("Assembling Figure 2...")

pA_wrapped <- wrap_elements(full = pA)

fig2 <- (pA_wrapped | pC_gg) /
         (pB_base   | pD) +
  plot_layout(widths = c(0.55, 0.45), heights = c(0.40, 0.60))

fig2_pdf <- file.path(RPT_DIR, "Figure_2.pdf")
fig2_png <- file.path(RPT_DIR, "Figure_2.png")

ggsave(fig2_pdf, fig2,
       width = 380, height = 400, units = "mm", limitsize = FALSE)
ggsave(fig2_png, fig2,
       width = 380, height = 400, units = "mm", dpi = 300, limitsize = FALSE)

message("Figure 2 saved to: ", RPT_DIR)
