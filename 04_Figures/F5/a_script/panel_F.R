################################################################################
#   YvO Figure 5 — Panel F: Module Preservation (Pre -> Post Training)
################################################################################
#
# ── STAT AUDIT (2026-02-27) ──────────────────────────────────────────────────
#
# PERMUTATION COUNT
#   200 permutations.  Langfelder & Horvath (2011) recommend >= 200 for
#   "publishable" results and note that Zsummary stabilizes by ~200 perm.
#   This meets the minimum threshold.  For extra robustness, 500–1000
#   permutations could be run, but 200 is standard and sufficient for
#   the clear separation observed (most Zsummary >> 10).
#
# CI ON Zsummary
#   modulePreservation() does not directly return CIs on Zsummary.
#   Zsummary is a composite Z-statistic (mean of density- and
#   connectivity-based Z-scores), each derived from a permutation
#   null.  Individual Z-components (Zdensity, Zconnectivity) are
#   available; Zsummary = (Zdensity + Zconnectivity) / 2.
#   AUDIT ADDITION: extract Zdensity and Zconnectivity to show the
#   range of component Z-scores (not a formal CI, but informative).
#
# INTERPRETATION ZONES
#   Zsummary > 10: "strong evidence of preservation" (standard)
#   Zsummary 2–10: "moderate to weak evidence"
#   Zsummary < 2:  "no evidence of preservation"
#   These thresholds are from Langfelder et al. (2011) and are the
#   standard in the WGCNA literature.
#
# AUDIT VERDICT: 200 permutations is adequate.  Z-component breakdown
# added.  Thresholds are literature-standard.
# ─────────────────────────────────────────────────────────────────────────────

source("04_Figures/F5/a_script/YvO_F5_setup.R")

pre_samp_f  <- meta %>% filter(time == "Pre") %>% pull(sample_id)
post_samp_f <- meta %>% filter(time == "Post") %>% pull(sample_id)

datExpr_pre  <- datExpr[pre_samp_f, ]
datExpr_post <- datExpr[post_samp_f, ]

multiExpr <- list(
  Pre  = list(data = datExpr_pre),
  Post = list(data = datExpr_post)
)

multiColor <- list(Pre = module_colors)

cat("Starting module preservation (200 permutations, this may take 10-30 min)...\n")
mp <- modulePreservation(
  multiExpr,
  multiColor,
  referenceNetworks = 1,
  testNetworks      = 2,
  nPermutations     = 200,
  randomSeed        = 42,
  quickCor          = 0,
  verbose           = 3
)

ref  <- 1
test <- 2
z_summary <- mp$preservation$Z[[ref]][[test]]
mod_sizes <- z_summary[, "moduleSize"]

# STAT AUDIT: extract Z-component breakdown (Zdensity + Zconnectivity)
z_density <- if ("Zdensity.pres" %in% colnames(z_summary)) {
  z_summary[, "Zdensity.pres"]
} else {
  rep(NA_real_, nrow(z_summary))
}
z_connectivity <- if ("Zconnectivity.pres" %in% colnames(z_summary)) {
  z_summary[, "Zconnectivity.pres"]
} else {
  rep(NA_real_, nrow(z_summary))
}

pres_df <- tibble(
  module           = rownames(z_summary),
  Zsummary         = z_summary[, "Zsummary.pres"],
  Zdensity         = z_density,
  Zconnectivity    = z_connectivity,
  module_size      = mod_sizes
) %>%
  filter(module != "gold", module != "grey") %>%
  mutate(preservation = case_when(
    Zsummary > 10 ~ "Strong",
    Zsummary > 2  ~ "Moderate",
    TRUE          ~ "Not preserved"
  ))

pF <- ggplot(pres_df, aes(x = module_size, y = Zsummary)) +
  # Colored background zones
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 2,
           fill = "#D6604D", alpha = 0.08) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2, ymax = 10,
           fill = "#F39C12", alpha = 0.08) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 10, ymax = Inf,
           fill = "#2ECC71", alpha = 0.08) +
  geom_hline(yintercept = 2, linetype = "dashed",
             color = "grey40", linewidth = 0.3) +
  geom_hline(yintercept = 10, linetype = "dashed",
             color = "grey40", linewidth = 0.3) +
  # Left-aligned zone labels
  annotate("text", x = min(pres_df$module_size) * 0.9, y = 1,
           label = "Not preserved", size = 2.2, color = "grey40",
           hjust = 0, fontface = "italic") +
  annotate("text", x = min(pres_df$module_size) * 0.9, y = 6,
           label = "Moderate", size = 2.2, color = "grey40",
           hjust = 0, fontface = "italic") +
  annotate("text", x = min(pres_df$module_size) * 0.9, y = 15,
           label = "Strong", size = 2.2, color = "grey40",
           hjust = 0, fontface = "italic") +
  geom_point(aes(fill = module, size = module_size), shape = 21,
             color = "black", stroke = 0.3) +
  scale_fill_identity() +
  scale_size_continuous(range = c(2, 6), guide = "none") +
  geom_text_repel(aes(label = module), size = 2.0, max.overlaps = 20) +
  scale_x_log10() +
  labs(title    = "F  Module Preservation (Pre -> Post Training)",
       subtitle = sprintf("Zsummary > 10 = preserved | < 2 = remodeled | %d permutations",
                          200),
       x = "Module Size", y = "Zsummary (preservation)") +
  THEME_PUB

write_csv(pres_df, file.path(DAT_DIR, "06_panel_F_preservation.csv"))

# STAT AUDIT: print Z-component breakdown
cat("  Preservation Z-component breakdown:\n")
pres_print <- pres_df %>%
  dplyr::select(module, Zsummary, Zdensity, Zconnectivity, module_size, preservation) %>%
  arrange(desc(Zsummary))
print(as.data.frame(pres_print))

ggsave(file.path(RPT_DIR, "panel_F_preservation.pdf"), pF,
       width = 350, height = 200, units = "mm",
       device = pdf, limitsize = FALSE)
ggsave(file.path(RPT_DIR, "panel_F_preservation.png"), pF,
       width = 350, height = 200, units = "mm",
       dpi = 300, limitsize = FALSE)

cat("Panel F saved\n")
