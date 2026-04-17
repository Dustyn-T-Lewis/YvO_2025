# Pilot: Swapped fry Barcode — Training_Old DEPs in Aging-ranked t-statistics
# Complement to panel_F_fry.R which tests Aging DEPs in Training_Old ranking.
# Here we test: do Training_Old DEPs collectively shift in the Aging ranking?
#
# Circularity caveat: same as panel_F — Aging and Training_Old share Old_Pre
# with opposite signs -> structural negative correlation.
# ---------------------------------------------------------------------------
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(tidyverse)
library(limma)
library(patchwork)

set.seed(42)

RPT <- "04_Figures/F05/b_reports/supp/panels"
DAT <- "04_Figures/F05/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()
PE_W <- 220

# -- Load data ---------------------------------------------------------------

dal      <- readRDS("02_Imputation/c_data/01_DAList_imputed.rds")
dep_df   <- read_csv("03_DEP/c_data/03_combined_results.csv",
                      show_col_types = FALSE)
imp_csv  <- read_csv("02_Imputation/c_data/01_imputed.csv",
                      show_col_types = FALSE)

meta <- dal$metadata
sample_cols <- meta$Col_ID

# -- Build imputed matrix ----------------------------------------------------

mat_imp <- imp_csv %>%
  select(uniprot_id, all_of(sample_cols)) %>%
  column_to_rownames("uniprot_id") %>%
  as.matrix()

n_imp <- nrow(mat_imp)
imp_ids <- rownames(mat_imp)
message(sprintf("Imputed matrix: %d proteins x %d samples", n_imp, ncol(mat_imp)))

# -- Design matrix + duplicateCorrelation ------------------------------------

meta$Group_Time <- factor(meta$Group_Time,
  levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"))
design <- model.matrix(~ 0 + Group_Time, data = meta)
colnames(design) <- gsub("^Group_Time", "", colnames(design))

block_id <- sub("_(Pre|Post)$", "", meta$Col_ID)

corfit_imp <- duplicateCorrelation(mat_imp, design, block = block_id)
cor_imp <- corfit_imp$consensus.correlation
message(sprintf("Within-subject cor: %.4f", cor_imp))

# SWAPPED: contrast is Aging (not Training_Old)
cm <- makeContrasts(
  Aging = Old_Pre - Young_Pre,
  levels = design
)

# -- Define Training_Old gene sets (SWAPPED) ---------------------------------

define_sets_TO <- function(dep, ids, use_pi) {
  sig <- if (use_pi) dep %>% filter(pi_score_Training_Old < 0.05)
         else dep %>% filter(P.Value_Training_Old < 0.05)
  sig <- sig %>% filter(uniprot_id %in% ids)
  list(
    up   = match(sig$uniprot_id[sig$logFC_Training_Old > 0], ids),
    down = match(sig$uniprot_id[sig$logFC_Training_Old < 0], ids),
    up_ids   = sig$uniprot_id[sig$logFC_Training_Old > 0],
    down_ids = sig$uniprot_id[sig$logFC_Training_Old < 0]
  )
}

sets_pi <- define_sets_TO(dep_df, imp_ids, TRUE)
sets_p  <- define_sets_TO(dep_df, imp_ids, FALSE)

message(sprintf("Training_Old gene sets: Pi up=%d dn=%d | P up=%d dn=%d",
                length(sets_pi$up), length(sets_pi$down),
                length(sets_p$up),  length(sets_p$down)))

# -- Circularity diagnostic --------------------------------------------------

circ_r <- cor(dep_df$t_Aging, dep_df$t_Training_Old, use = "complete.obs")
message(sprintf("Circularity: r(t_Aging, t_TO) = %.3f", circ_r))

# -- Run fry with Aging contrast ---------------------------------------------

run_fry <- function(mat, sets, design, cm, block, cor_val, config) {
  map_dfr(c("up", "down"), function(dir) {
    idx <- sets[[dir]]
    if (length(idx) < 3) return(tibble(config = config, set = paste0("to_", dir),
                                        n = length(idx), direction = NA_character_,
                                        PValue = NA_real_, PValue.Mixed = NA_real_))
    res <- fry(mat, index = idx, design = design,
               contrast = cm[, "Aging"], block = block, correlation = cor_val)
    tibble(config = config, set = paste0("to_", dir),
           n = length(idx), direction = res$Direction[1],
           PValue = res$PValue[1], PValue.Mixed = res$PValue.Mixed[1])
  })
}

fry_all <- bind_rows(
  run_fry(mat_imp, sets_pi, design, cm, block_id, cor_imp, "Imp_Pi"),
  run_fry(mat_imp, sets_p,  design, cm, block_id, cor_imp, "Imp_P")
) %>%
  mutate(
    expected = ifelse(set == "to_up", "Down", "Up"),
    consistent = direction == expected,
    cor_within = cor_imp,
    circularity_r = circ_r
  )

message("\n=== fry results (swapped: TO DEPs in Aging ranking) ===")
print(fry_all %>% select(config, set, n, direction, PValue, expected, consistent))

# -- Barcode data: rank by t_Aging (SWAPPED) ---------------------------------
# Use P-value sets for barcode (Pi sets too sparse: 8/10)

t_rank <- dep_df %>%
  filter(uniprot_id %in% imp_ids, !is.na(t_Aging)) %>%
  arrange(desc(t_Aging)) %>%
  mutate(rank = row_number(),
         in_up   = uniprot_id %in% sets_p$up_ids,
         in_down = uniprot_id %in% sets_p$down_ids)

running_es <- function(t_vals, in_set) {
  n <- length(t_vals); n_h <- sum(in_set)
  if (n_h == 0) return(rep(0, n))
  hit_w <- ifelse(in_set, abs(t_vals), 0)
  miss_w <- 1 / (n - n_h)
  cumsum(ifelse(in_set, hit_w / sum(hit_w), -miss_w))
}

t_rank$es_up   <- running_es(t_rank$t_Aging, t_rank$in_up)
t_rank$es_down <- running_es(t_rank$t_Aging, t_rank$in_down)

# -- Barcode visualization ---------------------------------------------------

fry_up <- fry_all %>% filter(config == "Imp_P", set == "to_up")
fry_dn <- fry_all %>% filter(config == "Imp_P", set == "to_down")

txt_s <- scale_text(BASE_STAT, PE_W)
n_all <- nrow(t_rank)

make_barcode <- function(t_df, in_col, es_col, fry_row, title, color) {
  marks <- t_df %>% filter(.data[[in_col]])

  p_es <- ggplot(t_df, aes(x = rank, y = .data[[es_col]])) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = color, alpha = 0.04) +
    geom_area(fill = scales::alpha(color, 0.2), color = NA) +
    geom_line(color = color, linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60",
               linewidth = 0.3) +
    annotate("text", x = n_all * 0.98, y = Inf,
             label = sprintf("fry %s, %s (n = %d)",
                              fry_row$direction, fmt_p(fry_row$PValue),
                              fry_row$n),
             hjust = 1, vjust = 1.3, size = txt_s * 1.15, fontface = "bold",
             color = ifelse(fry_row$consistent, "grey20", "#DC2626")) +
    labs(y = "ES", title = title) +
    scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    FIG_THEME +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(4, 4, 0, 4, "mm"),
          plot.title = element_text(size = 10, face = "bold"))

  p_bc <- ggplot(marks, aes(x = rank, xend = rank, y = 0, yend = 1)) +
    geom_segment(color = color, linewidth = 0.3, alpha = 0.7) +
    scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    FIG_THEME +
    theme(axis.text = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey97"),
          plot.margin = margin(0, 4, 0, 4, "mm"))

  list(es = p_es, bc = p_bc)
}

# Training_Old sets colored blue; grey if n.s.
up_color <- ifelse(fry_up$PValue > 0.05, "grey65", CONTRAST_COLORS["Training_Old"])
up_title <- ifelse(fry_up$PValue > 0.05,
                   "Tr.(O)-Up DEPs (P < 0.05) \u2192 Aging ranked t-statistics (n.s.)",
                   "Tr.(O)-Up DEPs (P < 0.05) \u2192 Aging ranked t-statistics")
p1 <- make_barcode(t_rank, "in_up", "es_up", fry_up, up_title, up_color)

dn_color <- ifelse(fry_dn$PValue > 0.05, "grey65", CONTRAST_COLORS["Training_Old"])
dn_title <- ifelse(fry_dn$PValue > 0.05,
                   "Tr.(O)-Down DEPs \u2192 Aging ranked t-statistics (n.s.)",
                   "Tr.(O)-Down DEPs \u2192 Aging ranked t-statistics")
p2 <- make_barcode(t_rank, "in_down", "es_down", fry_dn, dn_title, dn_color)

# t-stat area: green for Aging ranking
p_t <- ggplot(t_rank, aes(x = rank, y = t_Aging)) +
  geom_area(fill = scales::alpha(CONTRAST_COLORS["Aging"], 0.25),
            color = CONTRAST_COLORS["Aging"], linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  annotate("text", x = n_all * 0.98, y = Inf,
           label = sprintf("Circularity: r = %.3f", circ_r),
           hjust = 1, vjust = 1.5, size = txt_s * 0.75, color = "grey50") +
  labs(x = sprintf("Protein rank by t(Aging)  [%d imputed proteins]", n_all),
       y = "t-stat") +
  scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
  FIG_THEME +
  theme(plot.margin = margin(2, 4, 4, 4, "mm"))

pE_fry <- p1$es / p1$bc / p_t +
  plot_layout(heights = c(3, 0.5, 1.5)) +
  plot_annotation(
    title = "PILOT: fry Rotation Test (Swapped) \u2014 Tr.(O) DEP sets in Aging ranked t-statistics",
    subtitle = sprintf("Tr.(O) Up set (P < 0.05, n = %d) \u2192 Aging | fry p = %s | Circularity: r = %.3f | dupCor = %.3f",
                        length(sets_p$up), fmt_p(fry_up$PValue), circ_r, cor_imp),
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 9))
  )

ggsave(file.path(RPT, "SUPP_pilot_fry_swapped.png"), pE_fry,
       width = PE_W, height = 130, units = "mm", dpi = 300)

message("\nPilot saved: ", file.path(RPT, "SUPP_pilot_fry_swapped.png"))
message("Pilot fry (swapped) done")
