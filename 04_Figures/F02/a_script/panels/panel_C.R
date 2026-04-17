# F02 — Panel C: DEPs per Contrast (Pseudo-log Stacked Bar)
# Adapted from 04_Figures/F03/a_script/panel_A.R
# Outputs: pC (ggplot object), panel_C_dep_counts_MAIN.pdf/.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F02/a_script/style.R")

library(ggplot2)
library(dplyr)
library(readr)
library(tibble)


DEP_FILE <- "03_DEP/c_data/03_combined_results.csv"
RPT      <- "04_Figures/F02/b_reports/panels"
DAT      <- "04_Figures/F02/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

CONTRASTS <- c("Aging", "Training_Young", "Training_Old", "Interaction")
dep_df    <- read_csv(DEP_FILE, show_col_types = FALSE)
pdf_device <- get_pdf_device()
PA_W <- 165
PA_H <- 80

SET_LABELS <- c(Aging = "Aging", Training_Young = "Tr.(Y)",
                Training_Old = "Tr.(O)", Interaction = "Tr.(O)\u2013Tr.(Y)")
all_genes <- unique(dep_df$gene[!is.na(dep_df$gene)])

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

n_total <- length(all_genes)
pi_ci <- data.frame(
  contrast = CONTRASTS,
  n_sig    = sapply(sig_sets, length),
  n_total  = n_total,
  pct      = 100 * sapply(sig_sets, length) / n_total,
  ci_lo    = sapply(sig_sets, function(s) 100 * binom.test(length(s), n_total)$conf.int[1]),
  ci_hi    = sapply(sig_sets, function(s) 100 * binom.test(length(s), n_total)$conf.int[2])
)

pi_total  <- sum(sapply(sig_sets, length))
fdr_total <- sum(sapply(CONTRASTS, function(ctr) {
  fdr_col <- paste0("adj.P.Val_", ctr)
  if (fdr_col %in% names(dep_df)) sum(dep_df[[fdr_col]] < 0.05, na.rm = TRUE) else 0
}))
p_total <- sum(sapply(CONTRASTS, function(ctr) {
  p_col <- paste0("P.Value_", ctr)
  if (p_col %in% names(dep_df)) sum(dep_df[[p_col]] < 0.05, na.rm = TRUE) else 0
}))

frac_list <- lapply(CONTRASTS, function(ctr) {
  fdr_col <- paste0("adj.P.Val_", ctr)
  p_col   <- paste0("P.Value_", ctr)
  tibble(
    contrast  = SET_LABELS[ctr],
    threshold = c("p < 0.05", "q < 0.05", "Pi < 0.05"),
    n = c(sum(!is.na(dep_df[[p_col]])   & dep_df[[p_col]]   < 0.05),
          sum(!is.na(dep_df[[fdr_col]]) & dep_df[[fdr_col]] < 0.05),
          length(sig_sets[[ctr]]))
  )
})
frac_df <- bind_rows(frac_list) |>
  mutate(
    contrast  = factor(contrast,
                       levels = rev(c("Aging", "Tr.(Y)",
                                      "Tr.(O)", "Tr.(O)\u2013Tr.(Y)"))),
    threshold = factor(threshold, levels = c("p < 0.05", "q < 0.05", "Pi < 0.05")),
    pct       = 100 * n / length(all_genes),
    fill_key  = paste(contrast, threshold, sep = "___")
  ) |>
  filter(n > 1)

SET_DISPLAY_COLORS <- c("Aging"  = unname(CONTRAST_COLORS["Aging"]),
                        "Tr.(Y)" = unname(CONTRAST_COLORS["Training_Young"]),
                        "Tr.(O)" = unname(CONTRAST_COLORS["Training_Old"]),
                        "Tr.(O)\u2013Tr.(Y)" = unname(CONTRAST_COLORS["Interaction"]))

FRAC_FILL <- c()
for (cname in names(SET_DISPLAY_COLORS)) {
  col <- unname(SET_DISPLAY_COLORS[cname])
  FRAC_FILL[paste(cname, "p < 0.05",      sep = "___")] <- adjustcolor(col, alpha.f = 0.15)
  FRAC_FILL[paste(cname, "q < 0.05",      sep = "___")] <- adjustcolor(col, alpha.f = 0.40)
  FRAC_FILL[paste(cname, "Pi < 0.05", sep = "___")] <- col
}

THRESH_LABEL <- c("p < 0.05" = "p", "q < 0.05" = "FDR",
                  "Pi < 0.05" = "\u03A0")

label_df <- frac_df |>
  group_by(contrast) |> arrange(contrast, threshold) |>
  mutate(label     = THRESH_LABEL[as.character(threshold)],
         next_pct  = lead(pct, default = 0),
         seg_width = pct - next_pct,
         label_y   = (next_pct + pct) / 2,
         text_col  = if_else(threshold == "p < 0.05", "grey20", "white")) |>
  filter(seg_width > 0.3) |>   # skip labels for segments too narrow to read
  ungroup()

# Pi count labels at bar endpoints
pi_label_df <- frac_df |>
  filter(threshold == "Pi < 0.05") |>
  mutate(label = as.character(n))

pC <- ggplot(frac_df, aes(x = contrast, y = pct, fill = fill_key)) +
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
            inherit.aes = FALSE, hjust = 0.5,
            size = 2.2, fontface = "bold") +
  scale_fill_manual(values = FRAC_FILL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)),
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  coord_flip() +
  labs(title = "DEPs per Contrast",
       subtitle = sprintf("%s proteins | \u03A0 %d | FDR %d | p %d",
                          format(length(all_genes), big.mark = ","),
                          pi_total, fdr_total, p_total),
       x = NULL, y = "% of proteome",
       tag = "C") +
  FIG_THEME + theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE,
                                                face = "bold.italic", color = "grey40"),
                    legend.position = "none",
                    axis.text.y = element_text(face = "bold",
                                               size = FIG_AXIS_TEXT - 0.5))

write.csv(pi_ci, file.path(DAT, "audit_panel_C_dep_fraction_ci.csv"), row.names = FALSE)

ggsave(file.path(RPT, "MAIN_panel_C_dep_counts.png"), pC,
       width = PA_W, height = PA_H, units = "mm", dpi = 300)
