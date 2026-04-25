# WGCNA F06 Supplementary — Subcellular Compartment Enrichment
#
# Tests whether WGCNA modules are enriched for proteins from specific
# subcellular compartments using Human Protein Atlas (HPA) annotations.
#
# Method follows Johnson et al. 2020 (PMID 32284590, Nat Med): one-sided
# Fisher exact test per module x compartment, BH correction across all tests.
# That paper uses cell-type markers from brain scRNA-seq; here we use
# subcellular localization as a proxy since muscle single-cell proteomics
# references are not yet in our dataset.
#
# Inputs:  00_input/HPA_skeletal_muscle_annotations.tsv
#          04_Figures/F06/c_data/wgcna/wgcna_module_assignments.csv
# Outputs: c_data/supp/a03_compartment_enrichment.csv
#          b_reports/supp/01_QC/SUPP_compartment_enrichment.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/figure_supplement_helpers.R")  # read_sheet_df

library(tidyverse)

RPT_PNG <- "04_Figures/F06/b_reports/supp/png/panels"
RPT_PDF <- "04_Figures/F06/b_reports/supp/pdf/panels"
DAT_OUT <- "04_Figures/F06/c_data/supp"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_OUT, recursive = TRUE, showWarnings = FALSE)

# --- Load data
hpa <- read.delim("00_input/HPA_skeletal_muscle_annotations.tsv",
                   stringsAsFactors = FALSE)
F06_SUPP <- "04_Figures/F06/c_data/F06_supplementary.xlsx"
stopifnot("F06 stitcher must run first: missing F06_supplementary.xlsx" =
  file.exists(F06_SUPP))
mods <- read_sheet_df(F06_SUPP, "WGCNA_module_assignments")

# Load biological labels for display (from F06 Excel sheet)
bio_labels <- tryCatch(
  read_sheet_df(F06_SUPP, "WGCNA_mod_bio_labels"),
  error = function(e) NULL
)

# --- Parse HPA subcellular locations
# Each protein may have multiple comma-separated locations
hpa_sub <- hpa %>%
  filter(Subcellular.main.location != "", Gene %in% mods$gene) %>%
  select(gene = Gene, location = Subcellular.main.location)

# Expand to one row per gene x location
hpa_long <- hpa_sub %>%
  separate_rows(location, sep = ",\\s*") %>%
  filter(location != "")

# Define major compartments of biological interest for skeletal muscle
COMPARTMENTS <- c(
  "Mitochondria",
  "Cytosol",
  "Nucleoplasm",
  "Plasma membrane",
  "Endoplasmic reticulum",
  "Golgi apparatus",
  "Vesicles",
  "Cytoskeleton"  # aggregate: Microtubules, Actin filaments, Intermediate filaments
)

# Map cytoskeletal subtypes to "Cytoskeleton"
hpa_long <- hpa_long %>%
  mutate(compartment = case_when(
    location %in% c("Microtubules", "Actin filaments",
                     "Intermediate filaments", "Cytoskeleton") ~ "Cytoskeleton",
    location %in% COMPARTMENTS ~ location,
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(compartment))

# Build compartment gene sets
compartment_sets <- hpa_long %>%
  distinct(gene, compartment) %>%
  group_by(compartment) %>%
  summarise(genes = list(gene), n = n(), .groups = "drop")

cat(sprintf("Compartment gene sets built (%d compartments):\n", nrow(compartment_sets)))
for (i in seq_len(nrow(compartment_sets))) {
  cat(sprintf("  %s: %d genes\n", compartment_sets$compartment[i], compartment_sets$n[i]))
}

# --- Universe: all WGCNA genes with HPA location
universe <- unique(hpa_long$gene)
n_universe <- length(universe)
cat(sprintf("Universe: %d genes (WGCNA genes with HPA location)\n", n_universe))

# --- Fisher exact test: module x compartment
module_colors <- setdiff(unique(mods$module_color), "grey")

results <- list()
for (mc in module_colors) {
  mod_genes <- mods$gene[mods$module_color == mc]
  mod_genes_in_universe <- intersect(mod_genes, universe)
  n_mod <- length(mod_genes_in_universe)

  for (i in seq_len(nrow(compartment_sets))) {
    comp <- compartment_sets$compartment[i]
    comp_genes <- compartment_sets$genes[[i]]
    n_comp <- length(comp_genes)

    # 2x2 contingency: module membership x compartment membership
    in_both  <- length(intersect(mod_genes_in_universe, comp_genes))
    in_mod   <- n_mod - in_both
    in_comp  <- n_comp - in_both
    in_neither <- n_universe - in_both - in_mod - in_comp

    ct <- matrix(c(in_both, in_mod, in_comp, in_neither), nrow = 2,
                 dimnames = list(c("in_compartment", "not_in_compartment"),
                                 c("in_module", "not_in_module")))

    # One-sided Fisher exact (enrichment only, following Johnson 2020)
    ft <- fisher.test(ct, alternative = "greater")

    results[[length(results) + 1]] <- tibble(
      module       = mc,
      compartment  = comp,
      odds_ratio   = ft$estimate,
      p_raw        = ft$p.value,
      n_overlap    = in_both,
      n_module     = n_mod,
      n_compartment = n_comp,
      n_universe   = n_universe
    )
  }
}

enrich <- bind_rows(results) %>%
  mutate(p_bh = p.adjust(p_raw, method = "BH")) %>%
  arrange(p_bh)

write_csv(enrich, file.path(DAT_OUT, "a03_compartment_enrichment.csv"))
cat(sprintf("\nSignificant enrichments (FDR < 0.05): %d / %d\n",
            sum(enrich$p_bh < 0.05), nrow(enrich)))

# --- Validation checks
validate <- function(module_label, expected_compartment) {
  row <- enrich %>% filter(module == module_label, compartment == expected_compartment)
  if (nrow(row) == 0) {
    cat(sprintf("  CHECK MISSING: %s -> %s\n", module_label, expected_compartment))
  } else {
    status <- if (row$p_bh < 0.05) "PASS" else "WEAK"
    cat(sprintf("  %s: %s -> %s (OR=%.1f, p_bh=%.4f, n=%d)\n",
                status, module_label, expected_compartment,
                row$odds_ratio, row$p_bh, row$n_overlap))
  }
}

cat("\nBiological validation:\n")
validate("green",     "Mitochondria")
validate("yellow",    "Cytoskeleton")
validate("brown",     "Cytosol")
validate("black",     "Nucleoplasm")
validate("turquoise", "Mitochondria")

# --- Heatmap visualization
# Add biological labels if available
if (!is.null(bio_labels)) {
  if (!"display_label" %in% colnames(bio_labels)) {
    bio_labels <- bio_labels %>%
      mutate(display_label = ifelse(is.na(bio_label),
                                    stringr::str_to_title(module_color),
                                    paste0(bio_label, " (", stringr::str_to_title(module_color), ")")))
  }
  label_map <- setNames(bio_labels$display_label, bio_labels$module_color)
  enrich <- enrich %>%
    mutate(module_label = ifelse(module %in% names(label_map),
                                 label_map[module], module))
} else {
  enrich <- enrich %>% mutate(module_label = module)
}

# Order modules by total enrichment signal
mod_order <- enrich %>%
  group_by(module_label) %>%
  summarise(total_sig = -sum(log10(pmax(p_bh, 1e-10))), .groups = "drop") %>%
  arrange(desc(total_sig)) %>%
  pull(module_label)

enrich <- enrich %>%
  mutate(
    module_label = factor(module_label, levels = rev(mod_order)),
    compartment  = factor(compartment, levels = COMPARTMENTS),
    neg_log10_p  = -log10(pmax(p_bh, 1e-10)),
    sig_label    = case_when(
      p_bh < 0.001 ~ "***",
      p_bh < 0.01  ~ "**",
      p_bh < 0.05  ~ "*",
      TRUE          ~ ""
    )
  )

p <- ggplot(enrich, aes(compartment, module_label, fill = neg_log10_p)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sig_label), size = 3.5, vjust = 0.75) +
  scale_fill_gradient(low = "grey95", high = "#1B5E20",
                      limits = c(0, max(enrich$neg_log10_p, na.rm = TRUE)),
                      name = expression(-log[10](p[BH]))) +
  labs(x = NULL, y = NULL,
       title = "Module Subcellular Compartment Enrichment",
       subtitle = sprintf("One-sided Fisher exact test, BH-corrected (%d tests)  |  * p<0.05  ** p<0.01  *** p<0.001",
                          nrow(enrich))) +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8),
        axis.text.y = element_text(size = 7.5),
        plot.title    = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10))

pdf_device <- get_pdf_device()

ggsave(file.path(RPT_PNG, "SUPP_compartment_enrichment.png"), p,
       width = 180, height = 130, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_compartment_enrichment.pdf"), p,
       width = 180, height = 130, units = "mm", device = pdf_device)

cat("  Compartment enrichment complete.\n")
