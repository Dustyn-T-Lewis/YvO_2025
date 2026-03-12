################################################################################
#   Figure 1 — Panel D: DEPs per Contrast (Pseudo-log Stacked Bar)
#
#   Requires from setup: dep_df, CONTRASTS, CONTRAST_COLORS, THEME_PUB,
#                         RPT_DIR, DAT_DIR
#   Outputs: pD (ggplot object)
#   Side-effects: creates sig_sets, dir_map (used by Panel E)
################################################################################
#
# STAT AUDIT (2026-02-27)
# ---------------------------------------------------------------------------
# 1. Test appropriateness:
#    - This is a descriptive panel showing proportions of DEPs at three
#      significance thresholds across four contrasts. No formal test is
#      applied — this is appropriate for a summary display.            PASS
#
# 2. Assumption checking:
#    - N/A — descriptive.                                              PASS
#
# 3. Multiple comparison correction:
#    - The underlying q-values and pi-scores already incorporate FDR
#      correction from the limma pipeline.                             PASS
#
# 4. Effect sizes:
#    - Proportions (% of proteome) reported.                           PASS
#
# 5. Sample size adequacy:
#    - N/A — reporting counts from the full proteome.                  PASS
#
# 6. Confidence intervals:
#    - No CI on proportion of proteome altered per contrast.           ISSUE
#      FIX: Add binomial (Clopper-Pearson) 95% CI on the pi < 0.05
#      fraction for each contrast.
#
# 7. Reproducibility:
#    - Deterministic.                                                  PASS
# ---------------------------------------------------------------------------

if (!exists("meta")) source("04_Figures/F1/a_script/YvO_F1_setup.R")

SET_LABELS <- c(Aging = "Aging", Training_Young = "Training (Young)",
                Training_Old = "Training (Old)", Interaction = "Interaction")
all_genes <- unique(dep_df$gene[!is.na(dep_df$gene)])

# Pi-score significant sets (shared with Panel E)
sig_sets <- list()
dir_map  <- list()

for (ctr in CONTRASTS) {
  pi_vals  <- dep_df[[paste0("pi_score_", ctr)]]
  lfc_vals <- dep_df[[paste0("logFC_", ctr)]]
  is_sig   <- !is.na(pi_vals) & pi_vals < 0.05
  sig_sets[[ctr]] <- dep_df$gene[is_sig]
  dir_map[[ctr]]  <- setNames(ifelse(lfc_vals[is_sig] > 0, "Up", "Down"),
                               dep_df$gene[is_sig])
}

# --- AUDIT FIX: Clopper-Pearson 95% CI on pi < 0.05 fraction per contrast ---
n_total <- length(all_genes)
pi_ci <- data.frame(
  contrast = CONTRASTS,
  n_sig    = sapply(sig_sets, length),
  n_total  = n_total,
  pct      = 100 * sapply(sig_sets, length) / n_total,
  ci_lo    = sapply(sig_sets, function(s)
    100 * binom.test(length(s), n_total)$conf.int[1]),
  ci_hi    = sapply(sig_sets, function(s)
    100 * binom.test(length(s), n_total)$conf.int[2])
)
cat("DEP fractions with 95% Clopper-Pearson CI:\n"); print(pi_ci)

# Compute FDR and Pi totals for subtitle
pi_total  <- sum(sapply(sig_sets, length))
fdr_total <- sum(sapply(CONTRASTS, function(ctr) {
  fdr_col <- paste0("adj.P.Val_", ctr)
  if (fdr_col %in% names(dep_df)) sum(dep_df[[fdr_col]] < 0.05, na.rm = TRUE) else 0
}))

frac_list <- lapply(CONTRASTS, function(ctr) {
  fdr_col <- paste0("adj.P.Val_", ctr)
  p_col   <- paste0("P.Value_", ctr)
  tibble(
    contrast  = SET_LABELS[ctr],
    threshold = c("p < 0.05", "q < 0.05", "\u03A0 < 0.05"),
    n = c(sum(!is.na(dep_df[[p_col]])   & dep_df[[p_col]]   < 0.05),
          sum(!is.na(dep_df[[fdr_col]]) & dep_df[[fdr_col]] < 0.05),
          length(sig_sets[[ctr]]))
  )
})
frac_df <- bind_rows(frac_list) |>
  mutate(
    contrast  = factor(contrast,
                       levels = rev(c("Aging", "Training (Young)",
                                      "Training (Old)", "Interaction"))),
    threshold = factor(threshold, levels = c("p < 0.05", "q < 0.05", "\u03A0 < 0.05")),
    pct       = 100 * n / length(all_genes),
    fill_key  = paste(contrast, threshold, sep = "___")
  ) |>
  filter(n > 1)

SET_DISPLAY_COLORS <- c("Aging"            = unname(CONTRAST_COLORS["Aging"]),
                        "Training (Young)" = unname(CONTRAST_COLORS["Training_Young"]),
                        "Training (Old)"   = unname(CONTRAST_COLORS["Training_Old"]),
                        "Interaction"      = unname(CONTRAST_COLORS["Interaction"]))

FRAC_FILL <- c()
for (cname in names(SET_DISPLAY_COLORS)) {
  col <- unname(SET_DISPLAY_COLORS[cname])
  FRAC_FILL[paste(cname, "p < 0.05",      sep = "___")] <- adjustcolor(col, alpha.f = 0.25)
  FRAC_FILL[paste(cname, "q < 0.05",      sep = "___")] <- adjustcolor(col, alpha.f = 0.55)
  FRAC_FILL[paste(cname, "\u03A0 < 0.05", sep = "___")] <- col
}

THRESH_LABEL <- c("p < 0.05" = "p <= 0.05", "q < 0.05" = "FDR <= 0.05",
                  "\u03A0 < 0.05" = "Pi <= 0.05")

label_df <- frac_df |>
  group_by(contrast) |> arrange(contrast, threshold) |>
  mutate(label     = THRESH_LABEL[as.character(threshold)],
         next_pct  = lead(pct, default = 0),
         seg_width = pct - next_pct,
         label_y   = (next_pct + pct) / 2,
         text_col  = if_else(threshold == "p < 0.05", "grey20", "white")) |>
  filter(seg_width > 1.5) |>   # skip labels for segments too narrow to read
  ungroup()

pD <- ggplot(frac_df, aes(x = contrast, y = pct, fill = fill_key)) +
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Aging"], alpha = 0.20,
           color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Training_Young"], alpha = 0.20,
           color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Training_Old"], alpha = 0.20,
           color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Interaction"], alpha = 0.20,
           color = "grey70", linewidth = 0.2) +
  geom_col(position = "identity", width = 0.75, color = "black", linewidth = 0.3) +
  geom_text(data = label_df,
            aes(x = contrast, y = label_y, label = label, color = I(text_col)),
            inherit.aes = FALSE, hjust = 0.5, size = 2.2, fontface = "bold",
            parse = TRUE) +
  scale_fill_manual(values = FRAC_FILL) +
  # pseudo-log transform: behaves like log at large values but linear near zero,
  # avoiding the arbitrary power-law exponent (sigma=1 sets the linear region)
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = exp(1)),
                     expand = expansion(mult = c(0, 0.05)),
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  coord_flip() +
  labs(title = "D  DEPs per Contrast",
       subtitle = sprintf("Fraction of %s filtered proteins",
                          format(length(all_genes), big.mark = ",")),
       x = NULL, y = "% of proteome (pseudo-log scale)") +
  THEME_PUB + theme(legend.position = "none",
                    axis.text.y = element_text(face = "bold", size = 7))

# --- AUDIT: Export binomial CI table ---
write.csv(pi_ci,
          file.path(DAT_DIR, "audit_panel_D_dep_fraction_ci.csv"), row.names = FALSE)

cat("Panel D done\n")

ggsave(file.path(RPT_DIR, "panel_D_dep_counts.pdf"), pD,
       width = 170, height = 70, units = "mm", device = pdf)
ggsave(file.path(RPT_DIR, "panel_D_dep_counts.png"), pD,
       width = 170, height = 70, units = "mm", dpi = 300)
