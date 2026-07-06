# Sourced by 02_supp_panels.R — expects style.R already loaded.

source("04_Figures_v2/shared/figure_supplement_helpers.R") # read_sheet_df

library(tidyverse)

BASE <- "04_Figures_v2/F06"

RPT_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
DAT_OUT <- file.path(BASE, "c_data", "supp")
DAT <- file.path(BASE, "c_data")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_OUT, recursive = TRUE, showWarnings = FALSE)

hpa <- read.delim("00_input/HPA_skeletal_muscle_annotations.tsv",
  stringsAsFactors = FALSE
)
F06_SUPP <- file.path(DAT, "F06_supplementary.xlsx")
stopifnot(
  "F06 stitcher must run first: missing F06_supplementary.xlsx" =
    file.exists(F06_SUPP)
)
mods <- read_sheet_df(F06_SUPP, "WGCNA_module_assignments")

bio_labels <- tryCatch(
  read_sheet_df(F06_SUPP, "WGCNA_mod_bio_labels"),
  error = function(e) NULL
)

hpa_sub <- hpa |>
  filter(Subcellular.main.location != "", Gene %in% mods$gene) |>
  select(gene = Gene, location = Subcellular.main.location)

hpa_long <- hpa_sub |>
  separate_rows(location, sep = ",\\s*") |>
  filter(location != "")

COMPARTMENTS <- c(
  "Mitochondria",
  "Cytosol",
  "Nucleoplasm",
  "Plasma membrane",
  "Endoplasmic reticulum",
  "Golgi apparatus",
  "Vesicles",
  "Cytoskeleton" # aggregates: Microtubules, Actin filaments, Intermediate filaments
)

hpa_long <- hpa_long |>
  mutate(compartment = case_when(
    location %in% c(
      "Microtubules", "Actin filaments",
      "Intermediate filaments", "Cytoskeleton"
    ) ~ "Cytoskeleton",
    location %in% COMPARTMENTS ~ location,
    TRUE ~ NA_character_
  )) |>
  filter(!is.na(compartment))

compartment_sets <- hpa_long |>
  distinct(gene, compartment) |>
  group_by(compartment) |>
  summarise(genes = list(gene), n = n(), .groups = "drop")

message(sprintf("Compartment gene sets built (%d compartments):", nrow(compartment_sets)))
for (i in seq_len(nrow(compartment_sets))) {
  message(sprintf("  %s: %d genes", compartment_sets$compartment[i], compartment_sets$n[i]))
}

universe <- unique(hpa_long$gene)
n_universe <- length(universe)
message(sprintf("Universe: %d genes (WGCNA genes with HPA location)", n_universe))

module_colors <- setdiff(unique(mods$module_color), "grey")

results <- lapply(module_colors, function(mc) {
  mod_genes <- mods$gene[mods$module_color == mc]
  mod_genes_in_universe <- intersect(mod_genes, universe)
  n_mod <- length(mod_genes_in_universe)

  comp_rows <- lapply(seq_len(nrow(compartment_sets)), function(i) {
    comp <- compartment_sets$compartment[i]
    comp_genes <- compartment_sets$genes[[i]]
    n_comp <- length(comp_genes)

    in_both <- length(intersect(mod_genes_in_universe, comp_genes))
    in_mod <- n_mod - in_both
    in_comp <- n_comp - in_both
    in_neither <- n_universe - in_both - in_mod - in_comp

    ct <- matrix(c(in_both, in_mod, in_comp, in_neither),
      nrow = 2,
      dimnames = list(
        c("in_compartment", "not_in_compartment"),
        c("in_module", "not_in_module")
      )
    )

    ft <- fisher.test(ct, alternative = "greater")

    tibble(
      module = mc,
      compartment = comp,
      odds_ratio = ft$estimate,
      p_raw = ft$p.value,
      n_overlap = in_both,
      n_module = n_mod,
      n_compartment = n_comp,
      n_universe = n_universe
    )
  })
  bind_rows(comp_rows)
})

enrich <- bind_rows(results) |>
  mutate(p_bh = p.adjust(p_raw, method = "BH")) |>
  arrange(p_bh)

write_csv(enrich, file.path(DAT_OUT, "a03_compartment_enrichment.csv"))
message(sprintf(
  "\nSignificant enrichments (FDR < 0.05): %d / %d",
  sum(enrich$p_bh < 0.05), nrow(enrich)
))

validate <- function(module_label, expected_compartment) {
  row <- enrich |> filter(module == module_label, compartment == expected_compartment)
  if (nrow(row) == 0) {
    message(sprintf("  CHECK MISSING: %s -> %s", module_label, expected_compartment))
  } else {
    status <- if (row$p_bh < 0.05) "PASS" else "WEAK"
    message(sprintf(
      "  %s: %s -> %s (OR=%.1f, p_bh=%.4f, n=%d)",
      status, module_label, expected_compartment,
      row$odds_ratio, row$p_bh, row$n_overlap
    ))
  }
}

message("\nBiological validation:")
validate("green", "Mitochondria")
validate("yellow", "Cytoskeleton")
validate("brown", "Cytosol")
validate("black", "Nucleoplasm")
validate("turquoise", "Mitochondria")

if (!is.null(bio_labels)) {
  if (!"display_label" %in% colnames(bio_labels)) {
    bio_labels <- bio_labels |>
      mutate(display_label = ifelse(is.na(bio_label),
        stringr::str_to_title(module_color),
        paste0(bio_label, " (", stringr::str_to_title(module_color), ")")
      ))
  }
  label_map <- setNames(bio_labels$display_label, bio_labels$module_color)
  enrich <- enrich |>
    mutate(module_label = ifelse(module %in% names(label_map),
      label_map[module], module
    ))
} else {
  enrich <- enrich |> mutate(module_label = module)
}

mod_order <- enrich |>
  group_by(module_label) |>
  summarise(total_sig = -sum(log10(pmax(p_bh, 1e-10))), .groups = "drop") |>
  arrange(desc(total_sig)) |>
  pull(module_label)

enrich <- enrich |>
  mutate(
    module_label = factor(module_label, levels = rev(mod_order)),
    compartment = factor(compartment, levels = COMPARTMENTS),
    neg_log10_p = -log10(pmax(p_bh, 1e-10)),
    sig_label = case_when(
      p_bh < 0.001 ~ "***",
      p_bh < 0.01 ~ "**",
      p_bh < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

p <- ggplot(enrich, aes(compartment, module_label, fill = neg_log10_p)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sig_label), size = 3.5, vjust = 0.75) +
  scale_fill_gradient(
    low = "grey95", high = "#1B5E20",
    limits = c(0, max(enrich$neg_log10_p, na.rm = TRUE)),
    name = expression(-log[10](p[BH]))
  ) +
  labs(
    x = NULL, y = NULL,
    title = "Module Subcellular Compartment Enrichment",
    subtitle = sprintf(
      "One-sided Fisher exact test, BH-corrected (%d tests)  |  * p<0.05  ** p<0.01  *** p<0.001",
      nrow(enrich)
    )
  ) +
  FIG_THEME +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, size = 8),
    axis.text.y = element_text(size = 7.5),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10)
  )

pdf_device <- get_pdf_device()

ggsave(file.path(RPT_PNG, "SUPP_compartment_enrichment.png"), p,
  width = 180, height = 130, units = "mm", dpi = 300
)
ggsave(file.path(RPT_PDF, "SUPP_compartment_enrichment.pdf"), p,
  width = 180, height = 130, units = "mm", device = pdf_device
)

message("  Compartment enrichment complete.")
