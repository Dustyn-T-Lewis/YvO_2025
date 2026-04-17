# F04 SUPP — Convergent Pathways: RRHO intersection fry ORA (Concordance)
# Identifies high-confidence concordance pathways supported by both threshold-free
# (RRHO) and design-aware (fry) methods.
#
# Ordering: reads F04/c_data/panel_E/rrho2_ora_concordant.csv and
# F04/c_data/panel_F_fry/driving_ora.csv — both are absorbed into
# F04_supplementary.xlsx and removed by cleanup_after_workbook() in the main
# stitcher. Run this script BEFORE F04/a_script/90_stitch_figure.R to produce
# the SUPP_panel_F_convergent.png output; running after produces an empty-CSV
# soft-skip (exit 0, message only).
# ---------------------------------------------------------------------------
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

library(tidyverse)

RPT <- "04_Figures/F04/b_reports/supp/panels"
DAT <- "04_Figures/F04/c_data"
dir.create(RPT, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(DAT, "panel_F_fry"), showWarnings = FALSE, recursive = TRUE)

pdf_device <- get_pdf_device()
PF_W <- 200

# -- Load RRHO ORA and fry ORA ------------------------------------------------

rrho_ora_path <- file.path(DAT, "panel_E", "rrho2_ora_concordant.csv")
fry_ora_path  <- file.path(DAT, "panel_F_fry", "driving_ora.csv")

if (!file.exists(rrho_ora_path)) {
  message("panel_F_convergent: RRHO ORA not found -- run panel_E.R first")
  write_csv(tibble(pathway = character(), rrho_padj = numeric(), fry_padj = numeric(),
                   rrho_overlap = integer(), fry_overlap = integer(),
                   quadrant = character(), database = character()),
            file.path(DAT, "panel_F_fry", "convergent_pathways.csv"))
  quit(save = "no", status = 0)
}

if (!file.exists(fry_ora_path)) {
  message("panel_F_convergent: fry ORA not found -- run panel_F_ORA.R first")
  write_csv(tibble(note = "fry ORA not available"),
            file.path(DAT, "panel_F_fry", "convergent_pathways.csv"))
  quit(save = "no", status = 0)
}

rrho_ora <- read_csv(rrho_ora_path, show_col_types = FALSE)
fry_ora  <- read_csv(fry_ora_path,  show_col_types = FALSE)

# Guard: skip if either is empty or has no real pathways
if (nrow(rrho_ora) == 0 || nrow(fry_ora) == 0 || "note" %in% names(fry_ora)) {
  message("panel_F_convergent: one or both ORA results empty -- skipping")
  write_csv(tibble(pathway = character(), rrho_padj = numeric(), fry_padj = numeric(),
                   rrho_overlap = integer(), fry_overlap = integer(),
                   quadrant = character(), database = character()),
            file.path(DAT, "panel_F_fry", "convergent_pathways.csv"))
  quit(save = "no", status = 0)
}

# -- Find shared pathways -----------------------------------------------------

normalize_pw <- function(x) tolower(str_trim(x))

rrho_clean <- rrho_ora %>%
  mutate(pw_key = normalize_pw(pathway)) %>%
  select(pathway, pw_key, rrho_padj = padj, rrho_overlap = overlap, quadrant) %>%
  distinct(pw_key, .keep_all = TRUE)

fry_clean <- fry_ora %>%
  mutate(pw_key = normalize_pw(pathway)) %>%
  select(pathway_fry = pathway, pw_key, fry_padj = padj,
         fry_overlap = overlap, database) %>%
  distinct(pw_key, .keep_all = TRUE)

convergent <- inner_join(rrho_clean, fry_clean, by = "pw_key") %>%
  mutate(pathway_label = clean_pathway_name(pathway)) %>%
  arrange(pmin(rrho_padj, fry_padj)) %>%
  # Deduplicate by pathway_label (same name from different DBs → keep best padj)
  group_by(pathway_label) %>%
  slice_min(pmin(rrho_padj, fry_padj), n = 1, with_ties = FALSE) %>%
  ungroup()

if (nrow(convergent) == 0) {
  message("panel_F_convergent: no shared pathways between RRHO and fry ORA")
  write_csv(tibble(pathway = character(), rrho_padj = numeric(), fry_padj = numeric(),
                   rrho_overlap = integer(), fry_overlap = integer(),
                   quadrant = character(), database = character()),
            file.path(DAT, "panel_F_fry", "convergent_pathways.csv"))
  quit(save = "no", status = 0)
}

write_csv(convergent %>% select(pathway = pathway_label, rrho_padj, fry_padj,
                                  rrho_overlap, fry_overlap, quadrant, database),
          file.path(DAT, "panel_F_fry", "convergent_pathways.csv"))
message(sprintf("Convergent pathways: %d", nrow(convergent)))

# -- Paired comparison chart ---------------------------------------------------

conv_long <- convergent %>%
  select(pathway_label, rrho_padj, fry_padj) %>%
  pivot_longer(cols = c(rrho_padj, fry_padj),
               names_to = "method", values_to = "padj") %>%
  mutate(
    method = recode(method, rrho_padj = "RRHO", fry_padj = "fry"),
    neg_log10_padj = -log10(padj),
    pathway_label = fct_reorder(pathway_label, neg_log10_padj, .fun = max)
  )

txt_s <- scale_text(BASE_STAT, PF_W)

conv_labels <- conv_long %>%
  group_by(pathway_label) %>%
  summarise(max_x = max(neg_log10_padj), .groups = "drop")

pConv <- ggplot(conv_long, aes(x = neg_log10_padj, y = pathway_label,
                                 fill = method)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(data = conv_labels,
            aes(x = 0.05, y = pathway_label, label = pathway_label),
            hjust = 0, size = txt_s * 0.85, inherit.aes = FALSE,
            color = "white", fontface = "bold") +
  scale_fill_manual(values = c(RRHO = "#2563EB", fry = "#D97706"), name = "Method") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = expression("Convergent Pathways: RRHO" ~ intersect ~ "fry (Concordance)"),
    subtitle = sprintf("%d pathways significant in both methods", nrow(convergent)),
    x = expression(-log[10](p[adj])),
    y = NULL
  ) +
  FIG_THEME +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.background = element_blank(),
    plot.background  = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom"
  )

conv_h <- max(100, nrow(convergent) * 10 + 60)
ggsave(file.path(RPT, "SUPP_panel_F_convergent.png"), pConv,
       width = PF_W, height = conv_h, units = "mm", dpi = 300)

message("F04 SUPP convergent pathways panel done")
