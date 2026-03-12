# Figure 6 — Panel D: kME x GS Scatter (Hub Proteins)
# Module membership vs |gene significance| for best predictive module
#
# STAT AUDIT (2026-02-27):
# 1. kME: Pearson r (expression ~ eigengene, all samples). CI via cor.test().
# 2. GS: |Pearson r| (baseline expression ~ delta VL). CI via cor.test().
# 3. Hub thresholds: kME>0.7, |GS|>0.3 (Langfelder & Horvath 2008).
#    Sensitivity at 0.6/0.2 and 0.8/0.4 annotated.
# 4. kME-GS association: Spearman rho with p-value.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
})

RPT <- "04_Figures/F6/b_reports"
DAT <- "04_Figures/F6/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

pred_cor  <- read_csv(file.path(DAT, "pred_cor.csv"), show_col_types = FALSE)
module_df <- read_csv(file.path(DAT, "module_df.csv"), show_col_types = FALSE)
kME_all   <- readRDS(file.path(DAT, "kME_all.rds"))
gs_vl     <- readRDS(file.path(DAT, "gs_vl.rds"))
shared    <- readRDS(file.path(DAT, "shared_objects.rds"))
common_subj <- shared$common_subj

message("Panel D: kME x GS scatter (hub proteins)...")

PD_W <- 250  # panel width mm
PD_H <- 250  # panel height mm

txt_gene  <- scale_text(BASE_GENE, PD_W)
txt_annot <- scale_text(BASE_STAT, PD_W)
best_mod_name  <- pred_cor %>% arrange(desc(max_r)) %>%
  slice(1) %>% pull(module)
best_mod_color <- gsub("^ME", "", best_mod_name)
kme_col <- paste0("kME", best_mod_color)

mod_proteins <- module_df$uniprot_id[module_df$module_color == best_mod_color]

mod_proteins <- intersect(mod_proteins, rownames(kME_all))
mod_proteins <- intersect(mod_proteins, rownames(gs_vl))

kme_gs_df <- tibble(
  uniprot_id = mod_proteins,
  kME = kME_all[mod_proteins, kme_col],
  GS  = abs(gs_vl[mod_proteins, 1])
) %>%
  left_join(module_df %>% dplyr::select(uniprot_id, gene),
            by = "uniprot_id") %>%
  filter(!is.na(kME), !is.na(GS))

kme_thresh <- 0.7
gs_thresh  <- 0.3
n_candidates <- sum(kme_gs_df$kME > kme_thresh & kme_gs_df$GS > gs_thresh)

kme_gs_df <- kme_gs_df %>%
  mutate(is_hub = kME > kme_thresh & GS > gs_thresh,
         label  = ifelse(is_hub | rank(-kME * GS) <= 15, gene, NA_character_))

sens_df <- tibble(
  kme_cut = c(0.6, 0.7, 0.8),
  gs_cut  = c(0.2, 0.3, 0.4),
  n_hubs  = c(
    sum(kme_gs_df$kME > 0.6 & kme_gs_df$GS > 0.2),
    sum(kme_gs_df$kME > 0.7 & kme_gs_df$GS > 0.3),
    sum(kme_gs_df$kME > 0.8 & kme_gs_df$GS > 0.4)
  )
)
kme_gs_cor <- cor.test(kme_gs_df$kME, kme_gs_df$GS, method = "spearman",
                        exact = FALSE)
message(sprintf("  rho(kME,|GS|)=%.3f, p=%.3g; hubs: %s",
                kme_gs_cor$estimate, kme_gs_cor$p.value,
                paste(sens_df$n_hubs, collapse="/")))

pD <- ggplot(kme_gs_df, aes(x = kME, y = GS)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = best_mod_color, alpha = 0.05) +
  geom_vline(xintercept = kme_thresh, linetype = "dashed",
             color = "grey40", linewidth = 0.3) +
  geom_hline(yintercept = gs_thresh, linetype = "dashed",
             color = "grey40", linewidth = 0.3) +
  geom_point(aes(color = is_hub, alpha = is_hub), size = 1.5) +
  scale_color_manual(values = c("FALSE" = scales::alpha(best_mod_color, 0.4),
                                 "TRUE" = best_mod_color),
                     guide = "none") +
  scale_alpha_manual(values = c("FALSE" = 0.4, "TRUE" = 0.8), guide = "none") +
  geom_text_repel(aes(label = label), size = txt_gene, max.overlaps = 20,
                  segment.size = 0.2, na.rm = TRUE) +
  annotate("text",
           x = max(kme_gs_df$kME, na.rm = TRUE),
           y = max(kme_gs_df$GS, na.rm = TRUE) * 0.95,
           label = sprintf("%d candidates (kME>%.1f & |GS|>%.1f)\nSensitivity: %d (0.6/0.2), %d (0.8/0.4)\nrho(kME, |GS|) = %.2f, p = %.2g",
                           n_candidates, kme_thresh, gs_thresh,
                           sens_df$n_hubs[1], sens_df$n_hubs[3],
                           kme_gs_cor$estimate, kme_gs_cor$p.value),
           hjust = 1, size = txt_annot, fontface = "bold", color = "grey25") +
  labs(title = sprintf("Hub Proteins (%s Module)",
                       str_to_title(best_mod_color)),
       subtitle = "Module Membership vs |Gene Significance| for delta VL",
       x = sprintf("Module Membership (kME, %s)", best_mod_color),
       y = "|Gene Significance| (delta VL)") +
  FIG_THEME +
  theme(legend.position = "none")

ggsave(file.path(RPT, "panel_D_hub_scatter.pdf"), pD,
       width = PD_W, height = PD_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT, "panel_D_hub_scatter.png"), pD,
       width = PD_W, height = PD_H, units = "mm", dpi = 300)

write_csv(kme_gs_df, file.path(DAT, "04_panel_D_kme_gs.csv"))
write_csv(sens_df,   file.path(DAT, "04_panel_D_hub_sensitivity.csv"))

message("  Panel D saved")
