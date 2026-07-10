#!/usr/bin/env Rscript
# Sourced by 02_supp_panels.R — expects style.R already loaded.
# F05 Supplementary Panel A: ORA Dedup Sensitivity
# Diagnostic for main Panel A — shows enrichment count stability across Jaccard cutoffs.
# Runs hypergeometric ORA on the two Reversed quadrants at cutoffs 0.3, 0.5, 0.7, 1.0.
# Grouped bar: x = quadrant, fill = cutoff, y = # enriched pathways (FDR < 0.05).

source("04_Figures/shared/pathway_utils.R")

pacman::p_load(tidyverse, fgsea)

BASE    <- "04_Figures/F05"
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
DAT     <- file.path(BASE, "c_data", "panel_supp")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT,     recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

dep <- read_csv("03_DEP/c_data/03_combined_results.csv",
                show_col_types = FALSE)

scatter_df <- dep |>
  transmute(gene, logFC_Aging, logFC_TO = logFC_Training_Old) |>
  filter(!is.na(logFC_Aging), !is.na(logFC_TO)) |>
  mutate(quadrant = case_when(
    logFC_Aging > 0 & logFC_TO < 0 ~ "Reversed Up",
    logFC_Aging < 0 & logFC_TO > 0 ~ "Reversed Down",
    logFC_Aging > 0 & logFC_TO > 0 ~ "Exacerbated Up",
    TRUE                            ~ "Exacerbated Down"
  ))

universe <- scatter_df$gene
pw_collection <- build_pathway_collection(min_size = 15, max_size = 500,
                                           include_goslim = FALSE,
                                           exclude_variants = TRUE)

# Sweep Jaccard cutoffs for Reversed quadrants
cutoffs  <- c(0.3, 0.5, 0.7, 1.0)
quad_set <- c("Reversed Up", "Reversed Down")

results <- list()
for (quad in quad_set) {
  genes_q <- scatter_df$gene[scatter_df$quadrant == quad]
  if (length(genes_q) < 5) next
  for (jc in cutoffs) {
    ora_res <- tryCatch(
      run_ora_deduplicated(
        genes         = genes_q,
        universe      = universe,
        pathways      = pw_collection,
        jaccard_cutoff = jc,
        min_size      = 15,
        max_size      = 500,
        padj_cutoff   = 1
      ),
      error = function(e) tibble()
    )
    n_sig <- if (nrow(ora_res) > 0) sum(ora_res$padj < 0.05) else 0L
    n_tot <- nrow(ora_res)
    results[[length(results) + 1]] <- tibble(
      quadrant       = quad,
      jaccard_cutoff = jc,
      n_enriched     = n_sig,
      n_total_tested = n_tot
    )
  }
}
sens_df <- bind_rows(results) |>
  mutate(cutoff_label = factor(sprintf("J = %.1f", jaccard_cutoff)))

write_csv(sens_df, file.path(DAT, "SUPP_ora_dedup_sensitivity.csv"))

pS_ora_dedup <- ggplot(sens_df, aes(x = quadrant, y = n_enriched, fill = cutoff_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "grey30", linewidth = 0.3) +
  scale_fill_brewer(palette = "Blues", name = "Jaccard cutoff") +
  labs(title    = "ORA Dedup Sensitivity (Reversal Quadrants)",
       subtitle = sprintf("Hypergeometric ORA | %d universe | Reversed quadrants only",
                          length(universe)),
       x = NULL, y = "Enriched pathways (FDR < 0.05)") +
  FIG_THEME +
  theme(legend.position = "right",
        axis.text.x = element_text(size = FIG_AXIS_TEXT, face = "bold"))

PW <- 89; PH <- 70
ggsave(file.path(RPT_PNG, "SUPP_ora_dedup.png"), pS_ora_dedup,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_ora_dedup.pdf"), pS_ora_dedup,
       width = PW, height = PH, units = "mm", device = pdf_device)

message("F05 SUPP Panel A (ORA dedup sensitivity) saved")

# Expose for composite
pS_ora_title    <- "ORA Dedup Sensitivity (Reversal Quadrants)"
pS_ora_subtitle <- sprintf("Hypergeometric ORA | %d universe | Reversed quadrants only",
                           length(universe))
pS_ora_dedup    <- strip_for_composite(pS_ora_dedup)
