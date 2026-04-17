# F04 Supplementary: Directional Asymmetry Analysis
# Q1: Is concordant-down overlap stronger than concordant-up?
# Q2: Are specific biological processes asymmetrically conserved?
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

library(tidyverse)
library(patchwork)

RPT <- "04_Figures/F04/b_reports/supp"
DAT <- "04_Figures/F04/c_data/supp"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

rr <- dep_df %>%
  transmute(gene,
            t_y = t_Training_Young, t_o = t_Training_Old,
            lfc_y = logFC_Training_Young, lfc_o = logFC_Training_Old,
            pi_y = pi_score_Training_Young, pi_o = pi_score_Training_Old,
            sig_y = sig_pi_Training_Young, sig_o = sig_pi_Training_Old) %>%
  filter(!is.na(t_y), !is.na(t_o))

n_total <- nrow(rr)

# ── 1. Protein-level directional overlap ─────────────────────────────────────

# Define directional gene sets by sign of t-statistic
y_up   <- rr$gene[rr$t_y > 0]
y_down <- rr$gene[rr$t_y < 0]
o_up   <- rr$gene[rr$t_o > 0]
o_down <- rr$gene[rr$t_o < 0]

# Concordant overlap
conc_up   <- intersect(y_up, o_up)      # both increase with training
conc_down <- intersect(y_down, o_down)   # both decrease with training

# Expected by chance
exp_conc_up   <- length(y_up) * length(o_up) / n_total
exp_conc_down <- length(y_down) * length(o_down) / n_total

# Hypergeometric enrichment (one-sided, greater than expected)
p_conc_up <- phyper(length(conc_up) - 1, length(o_up),
                    n_total - length(o_up), length(y_up), lower.tail = FALSE)
p_conc_down <- phyper(length(conc_down) - 1, length(o_down),
                      n_total - length(o_down), length(y_down), lower.tail = FALSE)

# Enrichment ratios
or_up   <- length(conc_up) / exp_conc_up
or_down <- length(conc_down) / exp_conc_down

# Fisher's exact test: is concordance asymmetric between up and down?
# 2x2 table: concordant vs discordant × up vs down direction
fisher_tbl <- matrix(c(
  length(conc_up),                                     # conc up
  length(y_up) - length(conc_up),                      # disc: Y up, O down
  length(conc_down),                                   # conc down
  length(y_down) - length(conc_down)                   # disc: Y down, O up
), nrow = 2, byrow = TRUE,
dimnames = list(c("Up_direction", "Down_direction"),
                c("Concordant", "Discordant")))

fisher_res <- fisher.test(fisher_tbl)

cat("=== Protein-level directional asymmetry (F04) ===\n")
cat(sprintf("Concordant Up:   %d / %d expected (OR = %.3f, p = %.2e)\n",
            length(conc_up), round(exp_conc_up), or_up, p_conc_up))
cat(sprintf("Concordant Down: %d / %d expected (OR = %.3f, p = %.2e)\n",
            length(conc_down), round(exp_conc_down), or_down, p_conc_down))
cat(sprintf("Fisher asymmetry test: OR = %.3f, p = %.4f\n",
            fisher_res$estimate, fisher_res$p.value))
cat(sprintf("  Interpretation: %s\n",
            ifelse(fisher_res$p.value < 0.05,
                   ifelse(fisher_res$estimate > 1,
                          "Up-direction MORE concordant than down",
                          "Down-direction MORE concordant than up"),
                   "No significant directional asymmetry")))

# ── 2. Magnitude-weighted asymmetry ──────────────────────────────────────────

# For significant proteins only (Pi < 0.05 in either contrast)
sig_either <- rr %>% filter(abs(sig_y) > 0 | abs(sig_o) > 0)
sig_conc_up   <- sig_either %>% filter(t_y > 0, t_o > 0)
sig_conc_down <- sig_either %>% filter(t_y < 0, t_o < 0)
sig_disc      <- sig_either %>% filter(sign(t_y) != sign(t_o))

cat(sprintf("\nSignificant proteins (Pi < 0.05 in either):\n"))
cat(sprintf("  Concordant Up:   %d (mean |t_y| = %.2f, mean |t_o| = %.2f)\n",
            nrow(sig_conc_up), mean(abs(sig_conc_up$t_y)), mean(abs(sig_conc_up$t_o))))
cat(sprintf("  Concordant Down: %d (mean |t_y| = %.2f, mean |t_o| = %.2f)\n",
            nrow(sig_conc_down), mean(abs(sig_conc_down$t_y)), mean(abs(sig_conc_down$t_o))))
cat(sprintf("  Discordant:      %d\n", nrow(sig_disc)))

# Wilcoxon: is |t_old| larger for concordant-up vs concordant-down proteins?
if (nrow(sig_conc_up) > 5 && nrow(sig_conc_down) > 5) {
  w_test <- wilcox.test(abs(sig_conc_up$t_o), abs(sig_conc_down$t_o))
  cat(sprintf("  Wilcoxon |t_Old|: up vs down p = %.4f (median %.2f vs %.2f)\n",
              w_test$p.value, median(abs(sig_conc_up$t_o)), median(abs(sig_conc_down$t_o))))
}

# ── 3. Pathway-level directional asymmetry ───────────────────────────────────

# Load fGSEA results
fgsea_all <- read_csv("04_Figures/F04/c_data/shared/fgsea_tstat_all_v2.csv",
                      show_col_types = FALSE)

fgsea_wide <- fgsea_all %>%
  filter(contrast %in% c("Training_Young", "Training_Old")) %>%
  select(pathway, contrast, NES, padj) %>%
  pivot_wider(names_from = contrast, values_from = c(NES, padj),
              names_glue = "{.value}_{contrast}") %>%
  filter(!is.na(NES_Training_Young), !is.na(NES_Training_Old))

# Classify pathways by direction
pw_class <- fgsea_wide %>%
  mutate(
    sig_y = padj_Training_Young < 0.05,
    sig_o = padj_Training_Old < 0.05,
    conc_up   = NES_Training_Young > 0 & NES_Training_Old > 0,
    conc_down = NES_Training_Young < 0 & NES_Training_Old < 0,
    concordant = conc_up | conc_down,
    sig_both  = sig_y & sig_o
  )

cat("\n=== Pathway-level directional asymmetry ===\n")
cat(sprintf("Total pathways: %d\n", nrow(pw_class)))
cat(sprintf("Concordant Up:   %d (%d sig both)\n",
            sum(pw_class$conc_up), sum(pw_class$conc_up & pw_class$sig_both)))
cat(sprintf("Concordant Down: %d (%d sig both)\n",
            sum(pw_class$conc_down), sum(pw_class$conc_down & pw_class$sig_both)))

# Compare NES magnitude between concordant-up and concordant-down pathways
cu <- pw_class %>% filter(conc_up)
cd <- pw_class %>% filter(conc_down)

cat(sprintf("Mean |NES_Young| — Up: %.2f, Down: %.2f\n",
            mean(abs(cu$NES_Training_Young)), mean(abs(cd$NES_Training_Young))))
cat(sprintf("Mean |NES_Old|   — Up: %.2f, Down: %.2f\n",
            mean(abs(cu$NES_Training_Old)), mean(abs(cd$NES_Training_Old))))

# Blunting ratio: |NES_Old / NES_Young| (closer to 1 = more preserved)
cu_ratio <- abs(cu$NES_Training_Old / cu$NES_Training_Young)
cd_ratio <- abs(cd$NES_Training_Old / cd$NES_Training_Young)

w_ratio <- wilcox.test(cu_ratio, cd_ratio)
cat(sprintf("\nBlunting ratio |NES_Old/NES_Young| — Up: %.2f, Down: %.2f (Wilcoxon p = %.4f)\n",
            median(cu_ratio, na.rm = TRUE), median(cd_ratio, na.rm = TRUE), w_ratio$p.value))
cat(sprintf("  Interpretation: %s direction is more preserved across age\n",
            ifelse(median(cu_ratio, na.rm = TRUE) > median(cd_ratio, na.rm = TRUE),
                   "Up", "Down")))

# ── 4. Visualization: Asymmetry summary panel ────────────────────────────────

# Panel A: Protein-level enrichment comparison
protein_df <- tibble(
  direction = c("Concordant Up\n(both increase)", "Concordant Down\n(both decrease)"),
  observed  = c(length(conc_up), length(conc_down)),
  expected  = c(exp_conc_up, exp_conc_down),
  OR        = c(or_up, or_down),
  p_val     = c(p_conc_up, p_conc_down),
  n_hotspot = c(213, 344)  # from RRHO2
)

pA <- ggplot(protein_df, aes(x = direction)) +
  geom_col(aes(y = observed, fill = direction), width = 0.6, alpha = 0.85) +
  geom_point(aes(y = expected), shape = 4, size = 4, stroke = 1.2) +
  geom_text(aes(y = observed + 15,
                label = sprintf("OR = %.2f\np = %.1e", OR, p_val)),
            size = 2.8, lineheight = 0.85) +
  scale_fill_manual(values = c("Concordant Up\n(both increase)" = unname(DIR_COLORS["Up"]),
                                "Concordant Down\n(both decrease)" = unname(DIR_COLORS["Down"])),
                    guide = "none") +
  labs(title = "Protein-Level Overlap",
       subtitle = sprintf("Fisher asymmetry: OR = %.2f, p = %.3f",
                          fisher_res$estimate, fisher_res$p.value),
       y = "Proteins with concordant sign", x = NULL) +
  FIG_THEME

# Panel B: Pathway-level blunting ratio
ratio_df <- bind_rows(
  tibble(direction = "Up", ratio = cu_ratio),
  tibble(direction = "Down", ratio = cd_ratio)
) %>% filter(is.finite(ratio), ratio < 5)  # cap outliers

pB <- ggplot(ratio_df, aes(x = direction, y = ratio, fill = direction)) +
  geom_boxplot(width = 0.5, alpha = 0.7, outlier.size = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("Up" = unname(DIR_COLORS["Up"]), "Down" = unname(DIR_COLORS["Down"])),
                    guide = "none") +
  annotate("text", x = 1.5, y = max(ratio_df$ratio) * 0.95,
           label = sprintf("Wilcoxon p = %.3f", w_ratio$p.value),
           size = 3) +
  labs(title = "Pathway Preservation Ratio",
       subtitle = "|NES_Old / NES_Young| (1 = fully preserved)",
       y = "Blunting ratio", x = NULL) +
  FIG_THEME

# Panel C: Scatter of NES Young vs Old colored by direction
pC <- ggplot(pw_class, aes(x = NES_Training_Young, y = NES_Training_Old)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
           fill = DIR_COLORS["up"], alpha = 0.08) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = DIR_COLORS["down"], alpha = 0.08) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = case_when(
    conc_up & sig_both ~ "Conc. Up (sig both)",
    conc_down & sig_both ~ "Conc. Down (sig both)",
    concordant ~ "Concordant (NS)",
    TRUE ~ "Discordant"
  ), size = -log10(pmin(padj_Training_Young, padj_Training_Old))),
  alpha = 0.7) +
  scale_color_manual(values = c(
    "Conc. Up (sig both)" = unname(DIR_COLORS["Up"]),
    "Conc. Down (sig both)" = unname(DIR_COLORS["Down"]),
    "Concordant (NS)" = "grey60",
    "Discordant" = "grey80"
  ), name = NULL) +
  scale_size_continuous(range = c(1, 4), guide = "none") +
  labs(title = "Pathway NES: Young vs Old Training",
       subtitle = "Identity line = fully preserved; deviation = blunting",
       x = "NES (Training Young)", y = "NES (Training Old)") +
  coord_fixed() +
  FIG_THEME +
  theme(legend.position = "bottom", legend.text = element_text(size = 7))

p_composite <- pA + pB + pC + plot_layout(widths = c(1, 1, 1.5)) +
  plot_annotation(
    title = "Directional Asymmetry of Training Concordance",
    subtitle = sprintf("2,138 proteins | limma + dupCor | %d pathways (GO Slim + Hallmark)",
                        nrow(pw_class)),
    theme = theme(plot.title = element_text(size = 12, face = "bold"),
                  plot.subtitle = element_text(size = 9, color = "grey30"))
  )

ggsave(file.path(RPT, "supp_directional_asymmetry.pdf"), p_composite,
       width = 320, height = 130, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "supp_directional_asymmetry.png"), p_composite,
       width = 320, height = 130, units = "mm", dpi = 300)

# ── 5. Export audit data ─────────────────────────────────────────────────────

write_csv(protein_df, file.path(DAT, "asymmetry_protein_level.csv"))
write_csv(
  pw_class %>% select(pathway, NES_Training_Young, NES_Training_Old,
                      padj_Training_Young, padj_Training_Old,
                      conc_up, conc_down, concordant, sig_both),
  file.path(DAT, "asymmetry_pathway_level.csv")
)

message("F04 directional asymmetry analysis done")
