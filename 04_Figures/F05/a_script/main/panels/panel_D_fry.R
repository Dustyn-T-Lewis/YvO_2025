# F05 Panel D — fry Rotation Test: Aging Reversal Barcode
# Tests whether aging-significant protein sets collectively reverse with
# training in Old subjects, using limma's fry gene-set rotation test.
#
# fry is an exact rotation-based gene-set test (analogue of GSEA but uses
# rotation instead of permutation, accounting for inter-gene correlation).
#
# Circularity caveat: Aging (Old_Pre - Young_Pre) and Training_Old (Old_Post -
# Old_Pre) share Old_Pre with opposite signs -> structural negative correlation.
# Pi < 0.05 threshold reduces noise-driven bias vs nominal P < 0.05.
#
# References:
#   Wu & Smyth 2010, Bioinformatics — ROAST/fry
#   Melov et al. 2007, PLoS ONE — original reversal test
# ---------------------------------------------------------------------------
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

library(tidyverse)
library(limma)
library(fgsea)
library(patchwork)

set.seed(42)

RPT_PNG <- "04_Figures/F05/b_reports/main/png/panels"
RPT_PDF <- "04_Figures/F05/b_reports/main/pdf/panels"
RPT_SUP_PNG <- "04_Figures/F05/b_reports/supp/png/panels"
RPT_SUP_PDF <- "04_Figures/F05/b_reports/supp/pdf/panels"
DAT <- "04_Figures/F05/c_data"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUP_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUP_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_D_fry"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()
PE_W <- 270

# -- Step 1: Load data --------------------------------------------------------

dal      <- readRDS("02_Imputation/c_data/01_DAList_imputed.rds")
dep_df   <- read_csv("03_DEP/c_data/03_combined_results.csv",
                      show_col_types = FALSE)
imp_csv  <- read_csv("02_Imputation/c_data/01_imputed.csv",
                      show_col_types = FALSE)

meta <- dal$metadata
sample_cols <- meta$Col_ID

# -- Step 2: Build imputed matrix ----------------------------------------------

mat_imp <- imp_csv %>%
  select(uniprot_id, all_of(sample_cols)) %>%
  column_to_rownames("uniprot_id") %>%
  as.matrix()

# -- Step 3: Design matrix + duplicateCorrelation -----------------------------

meta$Group_Time <- factor(meta$Group_Time,
  levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"))
design <- model.matrix(~ 0 + Group_Time, data = meta)
colnames(design) <- gsub("^Group_Time", "", colnames(design))

block_id <- sub("_(Pre|Post)$", "", meta$Col_ID)

corfit_imp <- duplicateCorrelation(mat_imp, design, block = block_id)
cor_imp <- corfit_imp$consensus.correlation
message(sprintf("Within-subject cor: %.4f", cor_imp))

cm <- makeContrasts(
  Training_Old = Old_Post - Old_Pre,
  levels = design
)

# -- Step 4: Circularity diagnostic --------------------------------------------

circ_r <- cor(dep_df$t_Aging, dep_df$t_Training_Old, use = "complete.obs")
message(sprintf("Circularity: r(t_Aging, t_TO) = %.3f", circ_r))

# -- Step 5: Define aging gene sets --------------------------------------------

imp_ids <- rownames(mat_imp)

sig_pi <- dep_df %>%
  filter(pi_score_Aging < 0.05, uniprot_id %in% imp_ids)

sets_pi <- list(
  up       = match(sig_pi$uniprot_id[sig_pi$logFC_Aging > 0], imp_ids),
  down     = match(sig_pi$uniprot_id[sig_pi$logFC_Aging < 0], imp_ids),
  up_ids   = sig_pi$uniprot_id[sig_pi$logFC_Aging > 0],
  down_ids = sig_pi$uniprot_id[sig_pi$logFC_Aging < 0]
)

message(sprintf("Gene sets (Pi < 0.05): up = %d, down = %d",
                length(sets_pi$up), length(sets_pi$down)))

# -- Step 6: Run fry -----------------------------------------------------------

run_fry_set <- function(idx, set_name) {
  if (length(idx) < 3) return(tibble(set = set_name, n = length(idx),
                                      direction = NA_character_,
                                      PValue = NA_real_, PValue.Mixed = NA_real_))
  res <- fry(mat_imp, index = idx, design = design,
             contrast = cm[, "Training_Old"], block = block_id,
             correlation = cor_imp)
  tibble(set = set_name, n = length(idx), direction = res$Direction[1],
         PValue = res$PValue[1], PValue.Mixed = res$PValue.Mixed[1])
}

fry_up <- run_fry_set(sets_pi$up, "aging_up") %>%
  mutate(expected = "Down", consistent = direction == expected)
fry_dn <- run_fry_set(sets_pi$down, "aging_down") %>%
  mutate(expected = "Up", consistent = direction == expected)
fry_all <- bind_rows(fry_up, fry_dn) %>%
  mutate(cor_within = cor_imp, circularity_r = circ_r)

write_csv(fry_all, file.path(DAT, "panel_D_fry", "fry_results_all.csv"))

# -- Step 7: Driving proteins + leading-edge ORA -------------------------------

driving_up <- dep_df %>%
  filter(uniprot_id %in% sets_pi$up_ids, uniprot_id %in% imp_ids,
         t_Training_Old < 0) %>%
  transmute(gene, uniprot_id, set = "aging_up",
            t_aging = t_Aging, t_training_old = t_Training_Old,
            logFC_Aging, logFC_Training_Old, pi_score_Aging)

driving_dn <- dep_df %>%
  filter(uniprot_id %in% sets_pi$down_ids, uniprot_id %in% imp_ids,
         t_Training_Old > 0) %>%
  transmute(gene, uniprot_id, set = "aging_down",
            t_aging = t_Aging, t_training_old = t_Training_Old,
            logFC_Aging, logFC_Training_Old, pi_score_Aging)

driving_df <- bind_rows(driving_up, driving_dn)
write_csv(driving_df, file.path(DAT, "panel_D_fry", "driving_proteins.csv"))

# Leading-edge ORA on driving proteins
pw_collection <- build_pathway_collection(min_size = 10, max_size = 500,
                                           include_goslim = TRUE,
                                           exclude_variants = TRUE)
all_genes <- dep_df$gene[dep_df$uniprot_id %in% imp_ids]

ora_leading_up <- if (nrow(driving_up) >= 5) {
  tryCatch(
    run_ora_deduplicated(genes = driving_up$gene, universe = all_genes,
                          pathways = pw_collection, jaccard_cutoff = 0.5,
                          min_size = 10, max_size = 500, padj_cutoff = 0.1) %>%
      mutate(pathway_label = clean_pathway_name(pathway)) %>%
      slice_head(n = 5),
    error = function(e) tibble())
} else tibble()

ora_leading_dn <- if (nrow(driving_dn) >= 5) {
  tryCatch(
    run_ora_deduplicated(genes = driving_dn$gene, universe = all_genes,
                          pathways = pw_collection, jaccard_cutoff = 0.5,
                          min_size = 10, max_size = 500, padj_cutoff = 0.1) %>%
      mutate(pathway_label = clean_pathway_name(pathway)) %>%
      slice_head(n = 5),
    error = function(e) tibble())
} else tibble()

# -- Step 8: Barcode data -----------------------------------------------------

t_rank <- dep_df %>%
  filter(uniprot_id %in% imp_ids, !is.na(t_Training_Old)) %>%
  arrange(desc(t_Training_Old)) %>%
  mutate(rank = row_number(),
         in_up   = uniprot_id %in% sets_pi$up_ids,
         in_down = uniprot_id %in% sets_pi$down_ids)

running_es <- function(t_vals, in_set) {
  n <- length(t_vals); n_h <- sum(in_set)
  if (n_h == 0) return(rep(0, n))
  hit_w <- ifelse(in_set, abs(t_vals), 0)
  miss_w <- 1 / (n - n_h)
  cumsum(ifelse(in_set, hit_w / sum(hit_w), -miss_w))
}

t_rank$es_up   <- running_es(t_rank$t_Training_Old, t_rank$in_up)
t_rank$es_down <- running_es(t_rank$t_Training_Old, t_rank$in_down)

n_all <- nrow(t_rank)
txt_s <- scale_text(BASE_STAT, PE_W)

# -- Step 9: Barcode visualization — both Up and Down -------------------------

make_barcode <- function(t_df, in_col, es_col, fry_row, title, color) {
  marks <- t_df %>% filter(.data[[in_col]])
  is_sig <- !is.na(fry_row$PValue) && fry_row$PValue < 0.05
  line_color <- if (is_sig) color else scales::alpha(color, 0.4)

  p_label <- sprintf("fry %s, %s (n = %d)%s",
                      fry_row$direction, fmt_p(fry_row$PValue),
                      fry_row$n,
                      if (fry_row$consistent) "" else " \u2717")
  p_color <- if (fry_row$consistent) "grey20" else "#DC2626"

  p_es <- ggplot(t_df, aes(x = rank, y = .data[[es_col]])) +
    geom_area(fill = scales::alpha(line_color, 0.15), color = NA) +
    geom_line(color = line_color, linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60",
               linewidth = 0.3) +
    labs(title = title, subtitle = p_label, y = "ES") +
    scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    FIG_THEME +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title    = element_text(size = 9, face = "bold", color = "grey15",
                                       margin = margin(b = 0.5, unit = "mm")),
          plot.subtitle = element_text(size = 8, face = "bold", color = p_color,
                                       margin = margin(b = 1, unit = "mm")),
          plot.margin   = margin(3, 4, 0, 4, "mm"))

  p_bc <- ggplot(marks, aes(x = rank, xend = rank, y = 0, yend = 1)) +
    geom_segment(color = line_color, linewidth = 0.3, alpha = 0.7) +
    scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    FIG_THEME +
    theme(axis.text = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey97"),
          plot.margin = margin(0, 4, 0, 4, "mm"))

  list(es = p_es, bc = p_bc)
}

up_title <- sprintf("Aging-Up DEPs (\u03a0 < 0.05, n = %d) \u2192 Tr.(O) ranked t",
                     length(sets_pi$up))
dn_title <- sprintf("Aging-Down DEPs (\u03a0 < 0.05, n = %d) \u2192 Tr.(O) ranked t%s",
                     length(sets_pi$down),
                     if (fry_dn$PValue > 0.05) "  (n.s.)" else "")

# F05 reversal: both curves green (Aging contrast color)
FIG_COLOR <- unname(CONTRAST_COLORS["Aging"])  # #4CAF50

p1 <- make_barcode(t_rank, "in_up", "es_up", fry_up, up_title, FIG_COLOR)
p2 <- make_barcode(t_rank, "in_down", "es_down", fry_dn, dn_title, FIG_COLOR)

T_STAT_COLOR <- unname(CONTRAST_COLORS["Training_Old"])  # blue — matches F04 & F03
p_t <- ggplot(t_rank, aes(x = rank, y = t_Training_Old)) +
  geom_area(fill = scales::alpha(T_STAT_COLOR, 0.20),
            color = T_STAT_COLOR, linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  annotate("text", x = n_all * 0.98, y = Inf,
           label = sprintf("Circularity: r(t_Aging, t_TO) = %.3f", circ_r),
           hjust = 1, vjust = 1.5, size = txt_s * 0.75, color = "grey50") +
  labs(x = sprintf("Protein rank by t(Training Old)  [n = %d]", n_all),
       y = "t-stat") +
  scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
  FIG_THEME +
  theme(plot.margin = margin(2, 4, 4, 4, "mm"))

# -- Flanking ORA bar builder (compact, for right column) ---------------------

# F05-specific explicit label abbreviations
F05_LABEL_MAP <- c(
  "Amino Acid Catabolic Process"   = "AA Catabolic Process",
  "Fatty Acid Catabolic Process"   = "FA Catabolic Process",
  "Amino Acid Metabolic Process"   = "AA Metabolic Process",
  "Establishment Or Maintenance Of Cell Polarity" = "Cell Polarity",
  "Generation Of Precursor Metabolites And Energy" = "Precursor Metab. & Energy",
  "Protein Localization To Plasma Membrane" = "Plasma Membrane Protein Loc.",
  "Organic Acid Catabolic Process" = "Organic Acid Catabolism",
  "Membraneless Organelle Assembly" = "Membraneless Org. Assembly"
)

shorten_ora_label <- function(x, max_chars = 28, explicit_map = NULL) {
  # Strip "Reference " prefix FIRST so explicit-map keys don't need both variants
  x <- gsub("Reference ", "", x)
  if (!is.null(explicit_map)) {
    idx <- x %in% names(explicit_map)
    x[idx] <- explicit_map[x[idx]]
  }
  # Generic abbreviations
  x <- gsub("Regulation Of ", "Reg. ", x)
  x <- gsub("Negative Regulation Of ", "Neg. Reg. ", x)
  x <- gsub("Positive Regulation Of ", "Pos. Reg. ", x)
  x <- gsub("Epithelial Mesenchymal Transition", "EMT", x)
  x <- gsub(" Involved In ", " in ", x)
  x <- gsub("Regulated Microtubule Minus End Directed Transport", "MT Transport", x)
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 1), "\u2026"), x)
}

make_flanking_ora <- function(ora_df, set_label, bar_color,
                               label_map = NULL) {
  if (is.null(ora_df) || nrow(ora_df) == 0) {
    return(ggplot() + theme_void() +
             annotate("text", x = 0.5, y = 0.5, label = "No sig. pathways",
                      size = 3, color = "grey60"))
  }
  bars <- ora_df %>%
    slice_head(n = 5) %>%
    mutate(neg_log_padj = -log10(pmax(padj, 1e-20)),
           significant  = padj < 0.05,
           bar_fill     = ifelse(significant,
                                 scales::alpha(bar_color, 0.85),
                                 scales::alpha(bar_color, 0.30)),
           short_label  = shorten_ora_label(pathway_label,
                                             explicit_map = label_map),
           # Vertical asterisk stack — newline-separation forces vertical rendering
           star_raw = sig_stars(padj),
           star = gsub("\\*", "*\n", star_raw) %>% sub("\n$", "", .),
           y = rev(row_number()),
           bar_h = 0.7)

  x_max <- max(bars$neg_log_padj, na.rm = TRUE)
  x_display_max <- x_max * 1.15

  # Adaptive label placement (mirrors Panel A): inside-centered when bar is
  # wide enough to contain text, outside-right when bar is too short. Prevents
  # leading-character clipping on short bars where centered text would extend
  # left of x=0.
  bars <- bars %>%
    mutate(label_inside = neg_log_padj >= x_max * 0.55,
           label_x      = ifelse(label_inside,
                                 neg_log_padj * 0.5,
                                 neg_log_padj + x_max * 0.03),
           label_hjust  = ifelse(label_inside, 0.5, 0),
           label_color  = ifelse(label_inside,
                                 ifelse(significant, "white", "grey15"),
                                 "grey20"),
           text_size    = case_when(
             nchar(short_label) > 25 ~ 2.6,
             nchar(short_label) > 20 ~ 3.0,
             TRUE                    ~ 3.4))

  ggplot(bars, aes(y = y)) +
    geom_rect(aes(xmin = 0, xmax = neg_log_padj,
                  ymin = y - bar_h / 2, ymax = y + bar_h / 2),
              fill = bars$bar_fill, color = NA) +
    geom_text(aes(x = label_x, y = y, label = short_label),
              hjust = bars$label_hjust, size = bars$text_size, fontface = "bold",
              color = bars$label_color, lineheight = 0.85) +
    geom_text(aes(x = neg_log_padj + x_max * 0.03, label = star),
              hjust = 0, vjust = 0.5, size = 3.6, fontface = "bold",
              color = "black", lineheight = 0.5) +
    labs(title = set_label, x = expression(-log[10](p[adj])), y = NULL) +
    scale_x_continuous(limits = c(0, x_display_max),
                       breaks = scales::pretty_breaks(n = 3),
                       expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(limits = c(0.3, nrow(bars) + 0.7), expand = c(0, 0)) +
    theme_minimal(base_size = 9) +
    theme(panel.grid    = element_blank(),
          axis.text.y   = element_blank(),
          axis.ticks.y  = element_blank(),
          axis.title.y  = element_blank(),
          axis.text.x   = element_text(size = 8.5),
          axis.title.x  = element_text(size = 10),
          axis.line.x   = element_line(color = "grey40", linewidth = 0.3),
          plot.title    = element_text(face = "bold", size = 10, hjust = 0.5),
          plot.margin   = margin(3, 2, 0, 0, "mm"))
}

p_ora_flank_up <- make_flanking_ora(ora_leading_up, "Reversed (Up\u2192Down)",
                                     unname(DIR_COLORS["Up"]),
                                     label_map = F05_LABEL_MAP)
p_ora_flank_dn <- make_flanking_ora(ora_leading_dn, "Reversed (Down\u2192Up)",
                                     unname(DIR_COLORS["Down"]),
                                     label_map = F05_LABEL_MAP)

# -- Composite: fry (left) + flanking ORA bars (right) ------------------------
# Design: 5 rows (es1, bc1, es2, bc2, tstat) x 2 cols (fry, ora)
# ORA bars span the es+bc rows for each direction
fry_design <- c(
  area(1, 1, 1, 1),  # p1$es  (row 1, col 1)
  area(2, 1, 2, 1),  # p1$bc  (row 2, col 1)
  area(3, 1, 3, 1),  # p2$es  (row 3, col 1)
  area(4, 1, 4, 1),  # p2$bc  (row 4, col 1)
  area(5, 1, 5, 1),  # p_t    (row 5, col 1)
  area(1, 2, 2, 2),  # ora_up (rows 1-2, col 2)
  area(3, 2, 4, 2)   # ora_dn (rows 3-4, col 2)
)

pD_fry <- p1$es + p1$bc + p2$es + p2$bc + p_t +
  p_ora_flank_up + p_ora_flank_dn +
  plot_layout(design = fry_design,
              heights = c(2.5, 0.4, 2.5, 0.4, 1.2),
              widths  = c(3, 1.8)) +
  plot_annotation(
    title = "fry Gene-Set Rotation Test: Aging Reversal",
    subtitle = sprintf("Rotation-based set test (exact GSEA analogue) | Circularity r = %.3f | dupCor = %.3f | n = %d proteins",
                        circ_r, cor_imp, n_all),
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0,
                                            margin = margin(l = 12, unit = "mm")),
                  plot.subtitle = element_text(size = 10, color = "grey30", hjust = 0,
                                               margin = margin(l = 12, unit = "mm")),
                  plot.title.position = "panel")
  )

ggsave(file.path(RPT_PNG, "MAIN_panel_D_fry.png"), pD_fry,
       width = PE_W + 80, height = 175, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "MAIN_panel_D_fry.pdf"), pD_fry,
       width = PE_W + 80, height = 175, units = "mm", device = pdf_device)

# -- Step 10: Leading-edge ORA bar chart (SUPP) --------------------------------

# ORA composite-style bar builder (pathway names inside bars, sig stars at end)
make_ora_bars <- function(ora_df, set_label, bar_color, show_xaxis = FALSE) {
  if (nrow(ora_df) == 0) return(NULL)

  bars <- ora_df %>%
    slice_head(n = 3) %>%
    mutate(neg_log_padj = -log10(pmax(padj, 1e-20)),
           star = sig_stars(padj),
           y = rev(row_number()),
           bar_h = 0.7,
           text_size = case_when(
             nchar(pathway_label) > 35 ~ 2.8,
             nchar(pathway_label) > 25 ~ 3.4,
             TRUE ~ 4.0
           ))

  x_max <- max(bars$neg_log_padj, na.rm = TRUE)
  x_display_max <- x_max * 1.15

  ggplot(bars, aes(y = y)) +
    geom_rect(aes(xmin = 0, xmax = neg_log_padj,
                  ymin = y - bar_h / 2, ymax = y + bar_h / 2),
              fill = bar_color, color = NA) +
    geom_text(aes(x = neg_log_padj / 2, y = y, label = pathway_label),
              hjust = 0.5, size = bars$text_size, fontface = "bold",
              color = "white", lineheight = 0.85) +
    geom_text(aes(x = neg_log_padj + x_max * 0.03, label = star),
              hjust = 0, vjust = 0.5, size = 3.5, fontface = "bold",
              color = "black") +
    annotate("segment", x = 0, xend = x_display_max, y = -Inf, yend = -Inf,
             color = "grey40", linewidth = 0.3) +
    labs(title = set_label,
         x = if (show_xaxis) expression(-log[10](p[adj])) else NULL,
         y = NULL) +
    scale_x_continuous(limits = c(0, x_display_max),
                       breaks = scales::pretty_breaks(n = 3),
                       expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(limits = c(0.5, nrow(bars) + 0.5), expand = c(0, 0)) +
    theme_minimal(base_size = 9) +
    theme(panel.grid    = element_blank(),
          axis.text.y   = element_blank(),
          axis.ticks.y  = element_blank(),
          axis.title.y  = element_blank(),
          axis.text.x   = if (show_xaxis) element_text(size = 8.5, face = "bold")
                          else element_blank(),
          axis.title.x  = if (show_xaxis) element_text(size = 9.5, face = "bold")
                          else element_blank(),
          axis.line.x   = element_blank(),
          axis.ticks.x  = if (show_xaxis) element_line(color = "grey40", linewidth = 0.3)
                          else element_blank(),
          plot.title    = element_text(face = "bold", size = 10, hjust = 0.5),
          plot.margin   = margin(2, 6, 2, 2, "mm"))
}

p_ora_up <- make_ora_bars(ora_leading_up, "Aging-Up (Reversed)", unname(DIR_COLORS["Up"]),
                           show_xaxis = FALSE)
p_ora_dn <- make_ora_bars(ora_leading_dn, "Aging-Down (Reversed)", unname(DIR_COLORS["Down"]),
                           show_xaxis = TRUE)

if (!is.null(p_ora_up) || !is.null(p_ora_dn)) {
  ora_panels <- Filter(Negate(is.null), list(p_ora_up, p_ora_dn))
  p_ora <- wrap_plots(ora_panels, ncol = 1) +
    plot_annotation(
      title = "Leading-Edge ORA: fry Driving Proteins (Reversal)",
      subtitle = "Hypergeometric ORA on reversal-driving proteins | top 3 per set",
      theme = theme(plot.title    = element_text(size = 12, face = "bold", hjust = 0.5),
                    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey30")))

  ggsave(file.path(RPT_SUP_PNG, "SUPP_panel_D_fry_ora.png"), p_ora,
         width = 160, height = 100, units = "mm", dpi = 300)
  ggsave(file.path(RPT_SUP_PDF, "SUPP_panel_D_fry_ora.pdf"), p_ora,
         width = 160, height = 100, units = "mm", device = pdf_device)
  message("F05 Panel D ORA (supp) done")
}

message("F05 Panel D (fry) done")
