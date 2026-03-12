################################################################################
#   Figure 5 — Panel F: Module Preservation (Pre -> Post Training)
#   Self-contained panel script for 04_Figures_v2
#   Generates: panel_F_preservation.pdf/png, c_data/06_panel_F_preservation.csv
################################################################################
#
# ── STAT AUDIT (2026-02-27) ──────────────────────────────────────────────────
# 200 permutations (standard minimum for publishable Zsummary).
# Zsummary > 10: strong preservation; 2-10: moderate; < 2: none.
# Z-component breakdown (Zdensity + Zconnectivity) exported.
# ─────────────────────────────────────────────────────────────────────────────

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures_v2/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(WGCNA)
})

RPT <- "04_Figures_v2/F5/b_reports"
DAT <- "04_Figures_v2/F5/c_data"

pdf_device <- get_pdf_device()

# --- Load data ---
meta           <- read_csv(file.path(DAT, "meta.csv"), show_col_types = FALSE)
datExpr        <- readRDS(file.path(DAT, "datExpr.rds"))
module_colors  <- readRDS(file.path(DAT, "module_colors.rds"))

# --- Biological labels helper ---
mod_bio_labels_df  <- read_csv(file.path(DAT, "mod_bio_labels.csv"), show_col_types = FALSE)
mod_bio_labels_vec <- setNames(mod_bio_labels_df$bio_label, mod_bio_labels_df$module_color)
mod_display_label  <- function(color) {
  lbl <- mod_bio_labels_vec[color]
  lbl[is.na(lbl)] <- stringr::str_to_title(color[is.na(lbl)])
  paste0(lbl, " (", stringr::str_to_title(color), ")")
}

message("Panel F: module preservation...")

PF_W <- 180  # panel width mm
PF_H <- 120  # panel height mm

# --- Text sizes ---
txt_bar   <- scale_text(BASE_COUNT, PF_W)
txt_annot <- scale_text(BASE_GENE, PF_W)

# --- Split data by timepoint ---
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

# Extract Z-component breakdown
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

# --- Horizontal bar chart (compact) ---
pres_df <- pres_df %>%
  mutate(bio_label = mod_display_label(module)) %>%
  arrange(Zsummary) %>%
  mutate(bio_label = factor(bio_label, levels = bio_label))

pF <- ggplot(pres_df, aes(x = Zsummary, y = bio_label)) +
  geom_col(fill = pres_df$module, color = "black", linewidth = 0.3, width = 0.7) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = 2,  linetype = "dotted", color = "grey40", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.0f", Zsummary), x = Zsummary - 0.5),
            hjust = 1, size = txt_bar, fontface = "bold", color = "white") +
  annotate("text", x = 10.5, y = 0.5, label = "Strong (>10)",
           size = txt_annot, color = "grey40", hjust = 0, fontface = "italic") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title    = "Module Preservation (Pre \u2192 Post Training)",
       subtitle = sprintf("All %d modules Zsummary > 10 | %d permutations",
                          nrow(pres_df), 200),
       x = "Zsummary", y = NULL) +
  FIG_THEME +
  theme(panel.grid.major.y = element_blank(),
        legend.position    = "none")

write_csv(pres_df, file.path(DAT, "06_panel_F_preservation.csv"))

# Print Z-component breakdown
cat("  Preservation Z-component breakdown:\n")
pres_print <- pres_df %>%
  dplyr::select(module, Zsummary, Zdensity, Zconnectivity, module_size, preservation) %>%
  arrange(desc(Zsummary))
print(as.data.frame(pres_print))

ggsave(file.path(RPT, "panel_F_preservation.pdf"), pF,
       width = PF_W, height = PF_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT, "panel_F_preservation.png"), pF,
       width = PF_W, height = PF_H, units = "mm",
       dpi = 300, limitsize = FALSE)

message("  Panel F saved")
