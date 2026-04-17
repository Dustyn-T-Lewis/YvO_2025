# F05 Supplementary: Directional Asymmetry of Reversal
# Q1: Is reversal of aging-up proteins (by training-down) stronger than
#     reversal of aging-down proteins (by training-up)?
# Q2: Which biological direction is more reversible?
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

library(tidyverse)
library(patchwork)

RPT <- "04_Figures/F05/b_reports/supp"
DAT <- "04_Figures/F05/c_data/supp"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

rr <- dep_df %>%
  transmute(gene,
            t_aging = t_Aging, t_old = t_Training_Old,
            lfc_aging = logFC_Aging, lfc_old = logFC_Training_Old,
            sig_aging = sig_pi_Aging, sig_old = sig_pi_Training_Old) %>%
  filter(!is.na(t_aging), !is.na(t_old))

n_total <- nrow(rr)

# ── 1. Protein-level reversal asymmetry ──────────────────────────────────────

# Reversed = opposite signs (Aging up → Training down, or Aging down → Training up)
aging_up   <- rr$gene[rr$t_aging > 0]
aging_down <- rr$gene[rr$t_aging < 0]
tr_up      <- rr$gene[rr$t_old > 0]
tr_down    <- rr$gene[rr$t_old < 0]

# Reversal directions
rev_up_down <- intersect(aging_up, tr_down)    # aging increases, training decreases
rev_down_up <- intersect(aging_down, tr_up)    # aging decreases, training restores up

# Expected by chance
exp_up_down <- length(aging_up) * length(tr_down) / n_total
exp_down_up <- length(aging_down) * length(tr_up) / n_total

# Note: Because Aging and Training_Old share Old_Pre (opposite signs),
# there's a structural negative correlation (r ~ -0.33).
# This means reversal overlap is EXPECTED to be enriched even under the null.
# The asymmetry test (comparing up→down vs down→up) is still valid because
# the structural correlation affects both directions equally.

# Hypergeometric (excess over chance)
p_up_down <- phyper(length(rev_up_down) - 1, length(tr_down),
                    n_total - length(tr_down), length(aging_up), lower.tail = FALSE)
p_down_up <- phyper(length(rev_down_up) - 1, length(tr_up),
                    n_total - length(tr_up), length(aging_down), lower.tail = FALSE)

or_up_down <- length(rev_up_down) / exp_up_down
or_down_up <- length(rev_down_up) / exp_down_up

# Fisher's exact: asymmetry between reversal directions
fisher_tbl <- matrix(c(
  length(rev_up_down),                              # reversed: age up → tr down
  length(aging_up) - length(rev_up_down),           # not reversed: age up, tr up
  length(rev_down_up),                              # reversed: age down → tr up
  length(aging_down) - length(rev_down_up)          # not reversed: age down, tr down
), nrow = 2, byrow = TRUE,
dimnames = list(c("Aging_Up", "Aging_Down"), c("Reversed", "Not_reversed")))

fisher_res <- fisher.test(fisher_tbl)

cat("=== Protein-level reversal asymmetry (F05) ===\n")
cat(sprintf("Aging Up → Training Down:  %d / %.0f expected (OR = %.3f, p = %.2e)\n",
            length(rev_up_down), exp_up_down, or_up_down, p_up_down))
cat(sprintf("Aging Down → Training Up:  %d / %.0f expected (OR = %.3f, p = %.2e)\n",
            length(rev_down_up), exp_down_up, or_down_up, p_down_up))
cat(sprintf("Fisher asymmetry: OR = %.3f, p = %.4f\n",
            fisher_res$estimate, fisher_res$p.value))
cat("Note: structural correlation from shared Old_Pre inflates both; asymmetry unaffected\n")

# ── 2. Magnitude-weighted: which reversal direction has stronger effects? ────

sig_aging <- rr %>% filter(abs(sig_aging) > 0)
sig_rev_up   <- sig_aging %>% filter(t_aging > 0, t_old < 0)  # age up, tr down
sig_rev_down <- sig_aging %>% filter(t_aging < 0, t_old > 0)  # age down, tr up
sig_exac     <- sig_aging %>% filter(sign(t_aging) == sign(t_old))

cat(sprintf("\nAging-significant proteins (Pi < 0.05):\n"))
cat(sprintf("  Reversed (Age Up → Tr Down):   %d\n", nrow(sig_rev_up)))
cat(sprintf("  Reversed (Age Down → Tr Up):   %d\n", nrow(sig_rev_down)))
cat(sprintf("  Exacerbated (same direction):  %d\n", nrow(sig_exac)))
cat(sprintf("  Reversal %%: %.1f%%\n",
            100 * (nrow(sig_rev_up) + nrow(sig_rev_down)) /
            (nrow(sig_rev_up) + nrow(sig_rev_down) + nrow(sig_exac))))

# Compare the magnitude of reversal (|lfc_old| for reversed proteins)
if (nrow(sig_rev_up) > 5 && nrow(sig_rev_down) > 5) {
  # Reversal magnitude = |logFC_Training_Old / logFC_Aging| (capped at 2)
  sig_rev_up$rev_mag   <- pmin(abs(sig_rev_up$lfc_old / sig_rev_up$lfc_aging), 2)
  sig_rev_down$rev_mag <- pmin(abs(sig_rev_down$lfc_old / sig_rev_down$lfc_aging), 2)

  w_mag <- wilcox.test(sig_rev_up$rev_mag, sig_rev_down$rev_mag)
  cat(sprintf("\nReversal magnitude |lfc_old/lfc_aging|:\n"))
  cat(sprintf("  Age Up → Tr Down:  median = %.3f\n", median(sig_rev_up$rev_mag)))
  cat(sprintf("  Age Down → Tr Up:  median = %.3f\n", median(sig_rev_down$rev_mag)))
  cat(sprintf("  Wilcoxon p = %.4f\n", w_mag$p.value))
}

# ── 3. Pathway-level reversal asymmetry ──────────────────────────────────────

fgsea_all <- read_csv("04_Figures/F05/c_data/shared/fgsea_tstat_all_v2.csv",
                      show_col_types = FALSE)

fgsea_wide <- fgsea_all %>%
  filter(contrast %in% c("Aging", "Training_Old")) %>%
  select(pathway, contrast, NES, padj) %>%
  pivot_wider(names_from = contrast, values_from = c(NES, padj),
              names_glue = "{.value}_{contrast}") %>%
  filter(!is.na(NES_Aging), !is.na(NES_Training_Old))

pw_class <- fgsea_wide %>%
  mutate(
    sig_a = padj_Aging < 0.05,
    sig_o = padj_Training_Old < 0.05,
    # Reversal = opposite NES signs
    rev_up_down = NES_Aging > 0 & NES_Training_Old < 0,  # aging up, training reverses down
    rev_down_up = NES_Aging < 0 & NES_Training_Old > 0,  # aging down, training restores up
    reversed    = rev_up_down | rev_down_up,
    sig_both    = sig_a & sig_o,
    reversal_ratio = ifelse(abs(NES_Aging) > 0.1,
                            -NES_Training_Old / NES_Aging, NA_real_)
  )

cat("\n=== Pathway-level reversal asymmetry ===\n")
cat(sprintf("Reversed (Age Up → Tr Down): %d (%d sig both)\n",
            sum(pw_class$rev_up_down), sum(pw_class$rev_up_down & pw_class$sig_both)))
cat(sprintf("Reversed (Age Down → Tr Up): %d (%d sig both)\n",
            sum(pw_class$rev_down_up), sum(pw_class$rev_down_up & pw_class$sig_both)))

rud <- pw_class %>% filter(rev_up_down)
rdu <- pw_class %>% filter(rev_down_up)

if (nrow(rud) > 3 && nrow(rdu) > 3) {
  rud_ratio <- rud$reversal_ratio[is.finite(rud$reversal_ratio)]
  rdu_ratio <- rdu$reversal_ratio[is.finite(rdu$reversal_ratio)]
  w_pw <- wilcox.test(rud_ratio, rdu_ratio)
  cat(sprintf("Reversal ratio (-NES_Old/NES_Aging):\n"))
  cat(sprintf("  Age Up → Tr Down:  median = %.3f\n", median(rud_ratio)))
  cat(sprintf("  Age Down → Tr Up:  median = %.3f\n", median(rdu_ratio)))
  cat(sprintf("  Wilcoxon p = %.4f\n", w_pw$p.value))
}

# ── 4. Visualization ─────────────────────────────────────────────────────────

# Panel A: Protein-level reversal enrichment
prot_df <- tibble(
  direction = c("Aging Up\nTraining Down", "Aging Down\nTraining Up"),
  observed  = c(length(rev_up_down), length(rev_down_up)),
  expected  = c(exp_up_down, exp_down_up),
  OR        = c(or_up_down, or_down_up),
  p_val     = c(p_up_down, p_down_up),
  n_hotspot = c(615, 500)  # from RRHO2
)

pA <- ggplot(prot_df, aes(x = direction)) +
  geom_col(aes(y = observed, fill = direction), width = 0.6, alpha = 0.85) +
  geom_point(aes(y = expected), shape = 4, size = 4, stroke = 1.2) +
  geom_text(aes(y = observed + 15,
                label = sprintf("OR = %.2f\np = %.1e", OR, p_val)),
            size = 2.8, lineheight = 0.85) +
  scale_fill_manual(values = c("Aging Up\nTraining Down" = "#2563EB",
                                "Aging Down\nTraining Up" = "#2563EB"),
                    guide = "none") +
  labs(title = "Protein-Level Reversal Overlap",
       subtitle = sprintf("Fisher asymmetry: OR = %.2f, p = %.3f",
                          fisher_res$estimate, fisher_res$p.value),
       y = "Proteins with reversed sign", x = NULL) +
  FIG_THEME

# Panel B: Reversal magnitude for aging-sig proteins
if (exists("sig_rev_up") && nrow(sig_rev_up) > 5) {
  mag_df <- bind_rows(
    sig_rev_up %>% transmute(direction = "Age Up\nTr Down", rev_mag),
    sig_rev_down %>% transmute(direction = "Age Down\nTr Up", rev_mag)
  )

  pB <- ggplot(mag_df, aes(x = direction, y = rev_mag, fill = direction)) +
    geom_boxplot(width = 0.5, alpha = 0.7, outlier.size = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
    scale_fill_manual(values = c("Age Up\nTr Down" = "#DC2626",
                                  "Age Down\nTr Up" = "#2563EB"),
                      guide = "none") +
    annotate("text", x = 1.5, y = 1.9,
             label = sprintf("Wilcoxon p = %.3f", w_mag$p.value), size = 3) +
    annotate("text", x = 0.5, y = 1, label = "full reversal",
             hjust = 0, size = 2.5, color = "grey50", fontface = "italic") +
    labs(title = "Reversal Magnitude",
         subtitle = "Aging-sig. proteins (Pi < 0.05) | capped at 2",
         y = "|logFC_Old / logFC_Aging|", x = NULL) +
    FIG_THEME
} else {
  pB <- ggplot() + theme_void() + labs(title = "Insufficient data")
}

# Panel C: NES scatter colored by reversal direction
pC <- ggplot(pw_class, aes(x = NES_Aging, y = NES_Training_Old)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "#DCEEFF", alpha = 0.3) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#DCEEFF", alpha = 0.3) +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = case_when(
    rev_up_down & sig_both ~ "Rev: Age Up -> Tr Down (sig)",
    rev_down_up & sig_both ~ "Rev: Age Down -> Tr Up (sig)",
    reversed ~ "Reversed (NS)",
    TRUE ~ "Not reversed"
  ), size = -log10(pmin(padj_Aging, padj_Training_Old))),
  alpha = 0.7) +
  scale_color_manual(values = c(
    "Rev: Age Up -> Tr Down (sig)" = "#DC2626",
    "Rev: Age Down -> Tr Up (sig)" = "#2563EB",
    "Reversed (NS)" = "grey60",
    "Not reversed" = "grey80"
  ), name = NULL) +
  scale_size_continuous(range = c(1, 4), guide = "none") +
  labs(title = "Pathway NES: Aging vs Training Old",
       subtitle = "Anti-diagonal = full reversal | blue shading = reversed quadrants",
       x = "NES (Aging)", y = "NES (Training Old)") +
  coord_fixed() +
  FIG_THEME +
  theme(legend.position = "bottom", legend.text = element_text(size = 7))

p_composite <- pA + pB + pC + plot_layout(widths = c(1, 1, 1.5)) +
  plot_annotation(
    title = "Directional Asymmetry of Aging Reversal",
    subtitle = sprintf("2,138 proteins | %d pathways | structural r(t_Aging, t_Old) ~ -0.33",
                        nrow(pw_class)),
    theme = theme(plot.title = element_text(size = 12, face = "bold"),
                  plot.subtitle = element_text(size = 9, color = "grey30"))
  )

ggsave(file.path(RPT, "supp_directional_asymmetry.pdf"), p_composite,
       width = 320, height = 130, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "supp_directional_asymmetry.png"), p_composite,
       width = 320, height = 130, units = "mm", dpi = 300)

# Export
write_csv(prot_df, file.path(DAT, "asymmetry_protein_level.csv"))
write_csv(
  pw_class %>% select(pathway, NES_Aging, NES_Training_Old,
                      padj_Aging, padj_Training_Old,
                      rev_up_down, rev_down_up, reversed, sig_both, reversal_ratio),
  file.path(DAT, "asymmetry_pathway_level.csv")
)

message("F05 directional asymmetry analysis done")
