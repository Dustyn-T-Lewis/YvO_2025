# Sourced by 02_supp_panels.R — expects style.R + pathway_utils.R already loaded.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

BASE <- "04_Figures_v2/F06"

RPT_SUPP_PNG <- file.path(BASE, "b_reports", "supp", "png", "modules")
RPT_SUPP_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "modules")
DAT <- file.path(BASE, "c_data")
dir.create(RPT_SUPP_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PDF, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

message("Panel B: WGCNA per-module triptych...")

stopifnot(
  "WGCNA module assignments missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "wgcna/wgcna_module_assignments.csv")),
  "group_z.rds missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "group_z.rds")),
  "MEs.rds missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "MEs.rds")),
  "meta.csv missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "meta.csv")),
  "mod_bio_labels.csv missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "mod_bio_labels.csv")),
  "LMM contrast check missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "wgcna/wgcna_lmm_contrast_check.csv")),
  "LMM stratified check missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "wgcna/wgcna_lmm_stratified_check.csv"))
)

mod_assign_raw <- read_csv(file.path(DAT, "wgcna/wgcna_module_assignments.csv"))
mod_size_order <- mod_assign_raw |>
  filter(module_color != "grey") |>
  count(module_color, sort = TRUE) |>
  pull(module_color)
KEY_MODULES <- mod_size_order

message("  Generating triptych data for: ", paste(KEY_MODULES, collapse = ", "))

mod_assign <- mod_assign_raw
ann <- read_csv(file.path(DAT, "imp_annotations.csv"))
group_z <- readRDS(file.path(DAT, "group_z.rds"))
MEs <- readRDS(file.path(DAT, "MEs.rds"))
meta <- read_csv(file.path(DAT, "meta.csv"))

gene_map <- setNames(ann$uniprot_id, ann$gene)
z_long <- as.data.frame(group_z) |>
  tibble::rownames_to_column("gene") |>
  pivot_longer(-gene, names_to = "group", values_to = "z") |>
  mutate(uniprot_id = gene_map[gene]) |>
  inner_join(mod_assign |> select(uniprot_id, module_color), by = "uniprot_id") |>
  filter(module_color %in% KEY_MODULES) |>
  select(uniprot_id, gene, group, z, module = module_color)
write_csv(z_long, file.path(DAT, "03_panel_B_heatmap_zscores.csv"))
message(
  "    z-scores: ", n_distinct(z_long$gene), " genes across ",
  n_distinct(z_long$module), " modules"
)

me_long <- MEs |>
  tibble::rownames_to_column("sample_id") |>
  pivot_longer(starts_with("ME"), names_to = "me_col", values_to = "eigengene") |>
  mutate(module = gsub("^ME", "", me_col)) |>
  filter(module %in% KEY_MODULES) |>
  inner_join(meta |> select(sample_id, subject, age, time, group), by = "sample_id") |>
  mutate(group_label = paste0(substr(age, 1, 1), "-", time)) |>
  select(sample_id, group, subject, age, time, eigengene, module, group_label)
write_csv(me_long, file.path(DAT, "03_panel_B_eigengene_data.csv"))
message(
  "    eigengene: ", n_distinct(me_long$sample_id), " samples x ",
  n_distinct(me_long$module), " modules"
)

pw_collection <- build_pathway_collection(min_size = 15, max_size = 500)
universe <- unique(mod_assign$gene)

enrich_list <- lapply(KEY_MODULES, function(mod) {
  mod_genes <- mod_assign |>
    filter(module_color == mod) |>
    pull(gene) |>
    unique()
  if (length(mod_genes) < 5) {
    return(NULL)
  }
  res <- tryCatch(
    run_ora_deduplicated(mod_genes, universe, pw_collection,
      min_size = 15, max_size = 500, padj_cutoff = 1.0
    ),
    error = function(e) {
      message("    ORA failed for ", mod, ": ", e$message)
      NULL
    }
  )
  if (is.null(res) || nrow(res) == 0) {
    return(NULL)
  }
  res$module <- mod
  res$Description <- clean_pathway_name(res$pathway)
  res
})
enrich_all <- bind_rows(enrich_list)
write_csv(enrich_all, file.path(DAT, "03_panel_B_triptych_enrichment.csv"))
message(
  "    enrichment: ", nrow(enrich_all), " pathway hits across ",
  n_distinct(enrich_all$module), " modules"
)

z_scores <- z_long
me_data <- me_long
enrich <- enrich_all
mod_bio <- read_csv(file.path(DAT, "mod_bio_labels.csv"))
lmm_audit <- read_csv(file.path(DAT, "wgcna/wgcna_lmm_contrast_check.csv"))

if (!"display_label" %in% colnames(mod_bio)) {
  mod_bio <- mod_bio |>
    mutate(display_label = paste0(bio_label, " (", str_to_title(module_color), ")"))
}
mod_labels <- setNames(mod_bio$display_label, mod_bio$module_color)

lmm_interp <- lmm_audit |>
  filter(module %in% paste0("ME", KEY_MODULES)) |>
  group_by(module) |>
  arrange(p_bh) |>
  slice(1) |>
  ungroup() |>
  mutate(
    mod_color = gsub("^ME", "", module),
    interp = sprintf("%s, r = %+.2f", gsub("_", " ", contrast), r_equiv)
  )
interp_map <- setNames(lmm_interp$interp, lmm_interp$mod_color)

PB_W <- 280
PB_H <- 300

txt_heat <- scale_text(BASE_GENE, PB_W) * 0.7
txt_axis <- scale_text(BASE_STAT, PB_W) * 1.0
txt_title <- scale_text(BASE_GENE, PB_W) * 1.3
txt_bar <- scale_text(BASE_GENE, PB_W) * 0.95
txt_sig <- scale_text(BASE_GENE, PB_W) * 0.85

group_order <- c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post")
group_labels <- c(
  Young_Pre = "Y-Pre", Young_Post = "Y-Post",
  Old_Pre = "O-Pre", Old_Post = "O-Post"
)
lmm_stats <- lmm_audit |>
  filter(contrast %in% c("Training_Young", "Training_Old", "Aging")) |>
  mutate(module = gsub("^ME", "", module)) |>
  dplyr::select(module, contrast, p_bh)

triptych_row <- function(mod, show_xlab = FALSE) {
  label <- mod_labels[mod]
  n_mod <- z_scores |>
    filter(module == mod) |>
    distinct(gene) |>
    nrow()
  title_txt <- paste0(label, " (n=", n_mod, ")")

  z_mod <- z_scores |>
    filter(module == mod) |>
    mutate(group = factor(group, levels = group_order))

  gene_order <- z_mod |>
    filter(group == "Young_Pre") |>
    arrange(z) |>
    pull(gene)
  z_mod$gene <- factor(z_mod$gene, levels = gene_order)

  p_heat <- ggplot(z_mod, aes(x = group, y = gene, fill = z)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "#4393C3", mid = "white", high = "#D6604D",
      midpoint = 0, limits = c(-2, 2), oob = scales::squish,
      guide = "none"
    ) +
    scale_x_discrete(labels = group_labels, position = "bottom") +
    labs(title = title_txt, y = NULL, x = NULL) +
    FIG_THEME +
    theme(
      plot.title = element_text(size = txt_title, face = "bold"),
      axis.text.x = if (show_xlab) {
        element_text(size = txt_axis * 0.85, angle = 45, hjust = 1)
      } else {
        element_blank()
      },
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(2, 1, 2, 2)
    )

  me_mod <- me_data |>
    filter(module == mod) |>
    mutate(
      time = factor(time, levels = c("Pre", "Post")),
      age = factor(age, levels = c("Young", "Old"))
    )

  p_ty <- lmm_stats |>
    filter(module == mod, contrast == "Training_Young") |>
    pull(p_bh)
  p_to <- lmm_stats |>
    filter(module == mod, contrast == "Training_Old") |>
    pull(p_bh)
  p_ag <- lmm_stats |>
    filter(module == mod, contrast == "Aging") |>
    pull(p_bh)

  fmt_sig <- function(p) {
    if (length(p) == 0 || is.na(p)) {
      return("ns")
    }
    if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
  }

  me_means <- me_mod |>
    group_by(age, time) |>
    summarise(mean_me = mean(eigengene, na.rm = TRUE), .groups = "drop")
  n_young_pairs <- n_distinct(me_mod$subject[me_mod$age == "Young"])
  n_old_pairs <- n_distinct(me_mod$subject[me_mod$age == "Old"])
  n_all <- nrow(me_mod)

  p_eigen <- ggplot(me_mod, aes(x = time, y = eigengene)) +
    geom_line(aes(group = subject, color = age), alpha = 0.2, linewidth = 0.3) +
    stat_summary(aes(group = age, color = age), fun = mean, geom = "line", linewidth = 1.2) +
    stat_summary(aes(group = age, color = age), fun = mean, geom = "point", size = 2.5) +
    scale_color_manual(values = AGE_COLORS, guide = "none") +
    annotate("text",
      x = 1.5, y = max(me_mod$eigengene) * 0.95,
      label = paste0(fmt_sig(p_ty), " (n=", n_young_pairs, ")"),
      size = txt_sig, fontface = "bold", color = AGE_COLORS["Young"]
    ) +
    annotate("text",
      x = 1.5, y = min(me_mod$eigengene) * 0.95,
      label = paste0(fmt_sig(p_to), " (n=", n_old_pairs, ")"),
      size = txt_sig, fontface = "bold", color = AGE_COLORS["Old"]
    ) +
    annotate("text",
      x = 0.65, y = mean(c(max(me_mod$eigengene), min(me_mod$eigengene))),
      label = paste0(fmt_sig(p_ag), " (", n_all, ")"),
      size = txt_sig, fontface = "bold", color = "grey40", angle = 90
    ) +
    labs(y = "Eigengene", x = NULL) +
    FIG_THEME +
    theme(
      axis.text.x = if (show_xlab) element_text(size = txt_axis * 0.85) else element_blank(),
      axis.text.y = element_text(size = txt_axis * 0.75),
      axis.title.y = element_text(size = txt_axis * 0.8),
      panel.border = element_blank(),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
      plot.margin = margin_auto(2, 1)
    )

  bar_data <- enrich |>
    filter(module == mod, padj < 0.05) |>
    arrange(padj) |>
    head(5) |>
    mutate(
      neg_log10_p = -log10(padj),
      clean_name  = stringr::str_trunc(clean_pathway_name(Description), 40, ellipsis = "\u2026"),
      db_fill     = DB_COLORS[database]
    ) |>
    mutate(clean_name = make.unique(clean_name, sep = " ")) |>
    mutate(clean_name = factor(clean_name, levels = rev(clean_name)))

  if (nrow(bar_data) == 0) {
    p_bars <- ggplot() +
      annotate("text",
        x = 0.5, y = 0.5, label = "No sig.\nenrichment",
        size = txt_bar, color = "grey50"
      ) +
      theme_void() +
      theme(plot.margin = margin(2, 2, 2, 1))
  } else {
    p_bars <- ggplot(bar_data, aes(x = neg_log10_p, y = clean_name)) +
      geom_col(aes(fill = db_fill), color = "black", linewidth = 0.3, width = 0.7) +
      geom_text(aes(label = clean_name, x = 0.3),
        hjust = 0, size = txt_bar,
        fontface = "bold",
        color = ifelse(bar_data$db_fill %in% c("#AA336A", "#1565C0", "#00796B"),
          "white", "grey20"
        )
      ) +
      scale_fill_identity() +
      scale_x_continuous(
        expand = expansion(mult = c(0, 0)),
        breaks = scales::breaks_pretty(n = 3),
        name = if (show_xlab) expression(-log[10](p[adj])) else NULL
      ) +
      scale_y_discrete(labels = NULL) +
      labs(y = NULL) +
      FIG_THEME +
      theme(
        axis.text.x  = if (show_xlab) element_text(size = txt_axis * 0.75) else element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x  = element_line(color = "black", linewidth = 0.3),
        panel.grid   = element_blank(),
        plot.margin  = margin(2, 6, 2, 1)
      )
  }

  p_heat + p_eigen + p_bars + plot_layout(widths = c(3, 2, 3))
}

z_legend <- ggplot(
  data.frame(z = seq(-2, 2, length.out = 100)),
  aes(x = z, y = 1, fill = z)
) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#4393C3", mid = "white", high = "#D6604D",
    midpoint = 0, limits = c(-2, 2),
    name = "Z-score", guide = guide_colorbar(
      barwidth = unit(40, "mm"), barheight = unit(3, "mm")
    )
  ) +
  theme_void() +
  theme_sub_legend(position = "bottom", text = element_text(size = txt_axis * 0.7))

supp_letters <- setNames(LETTERS[seq_along(KEY_MODULES)], KEY_MODULES)
pathway_slug <- setNames(
  gsub("[/ ]+", "_", tolower(mod_bio$bio_label)),
  mod_bio$module_color
)

SINGLE_W <- 280
SINGLE_H <- 110 # single row + legend

for (i in seq_along(KEY_MODULES)) {
  mod <- KEY_MODULES[i]
  letter <- supp_letters[mod]
  row <- triptych_row(mod, show_xlab = TRUE)

  interp_text <- if (mod %in% names(interp_map)) interp_map[mod] else ""

  single <- row / wrap_elements(z_legend) +
    plot_layout(heights = c(1, 0.1)) +
    plot_annotation(
      title = mod_labels[mod],
      subtitle = interp_text,
      caption = "* FDR < .05  ** FDR < .01  *** FDR < .001 (LMM, BH-corrected)",
      theme = theme_sub_plot(
        title = element_text(face = "bold", size = txt_axis * 1.5),
        subtitle = element_text(size = txt_axis, color = "grey30", face = "italic"),
        caption = element_text(size = txt_axis * 0.7, color = "grey40", hjust = 0)
      )
    )

  slug <- pathway_slug[mod]
  if (is.na(slug)) slug <- mod
  fname <- sprintf("SUPP_triptych_%s_%s", mod, slug)
  ggsave(file.path(RPT_SUPP_PNG, paste0(fname, ".png")), single,
    width = SINGLE_W, height = SINGLE_H, units = "mm", dpi = 300
  )
  ggsave(file.path(RPT_SUPP_PDF, paste0(fname, ".pdf")), single,
    width = SINGLE_W, height = SINGLE_H, units = "mm", device = pdf_device
  )
  message(sprintf("  Saved %s (%s: %s)", fname, mod, interp_text))
}

message("  Panel B/C/D (individual triptychs) saved")
