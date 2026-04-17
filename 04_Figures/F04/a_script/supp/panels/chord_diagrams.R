# F04 SUPP: Quadrant-Based Chord Diagrams (Training Concordance)
# Pathway-first approach: top enriched pathways per quadrant -> leading-edge DEPs
# Protein arc colored by WGCNA module, chords colored by logFC, pathway arcs by NES
# Outputs: b_reports/supp/panels/SUPP_panel_G_chord{,_<quadrant>}.png (2x2: ConcUp, DiscYupOdn, ConcDown, DiscYdnOup)

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")
source("04_Figures/shared/volcano_ring.R")  # for clean_ring_label()
source("04_Figures/shared/figure_supplement_helpers.R")  # read_sheet_df

library(tidyverse)
library(circlize)
library(msigdbr)
library(fgsea)
library(patchwork)
library(grid)
library(png)

PW <- 190; PH <- 190
RPT <- "04_Figures/F04/b_reports/supp/panels"
DAT <- "04_Figures/F04/c_data"
dir.create(file.path(DAT, "panel_G"), recursive = TRUE, showWarnings = FALSE)
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# --- Load data ---
F06_SUPP <- "04_Figures/F06/c_data/F06_supplementary.xlsx"
stopifnot("F06 must run first: missing F06_supplementary.xlsx" = file.exists(F06_SUPP))
dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
wgcna_mod <- read_sheet_df(F06_SUPP, "WGCNA_module_assignments") %>%
  select(gene, module_color)

# fGSEA cache (shared, single source)
fgsea_cache <- "04_Figures/shared/fgsea_tstat_all_v2.csv"
stopifnot("fGSEA cache missing: frozen cache — see 04_Figures/shared/README.md" =
  file.exists(fgsea_cache))
fgsea_all <- read_csv(fgsea_cache, show_col_types = FALSE)

# --- Classify proteins into quadrants ---
quad_df <- dep_df %>%
  filter(!is.na(logFC_Training_Young), !is.na(logFC_Training_Old)) %>%
  mutate(
    sig_TY  = !is.na(pi_score_Training_Young) & pi_score_Training_Young < 0.05,
    sig_TO  = !is.na(pi_score_Training_Old)   & pi_score_Training_Old   < 0.05,
    sig_any = sig_TY | sig_TO,
    quadrant = case_when(
      logFC_Training_Young > 0 & logFC_Training_Old > 0 ~ "Concordant Up",
      logFC_Training_Young < 0 & logFC_Training_Old < 0 ~ "Concordant Down",
      logFC_Training_Young > 0 & logFC_Training_Old < 0 ~ "Discordant Y+/O-",
      logFC_Training_Young < 0 & logFC_Training_Old > 0 ~ "Discordant Y-/O+",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(quadrant)) %>%
  left_join(wgcna_mod, by = "gene")

# Fill missing module with grey
quad_df$module_color[is.na(quad_df$module_color)] <- "grey"

QUADRANTS <- c("Concordant Up", "Concordant Down",
               "Discordant Y+/O-", "Discordant Y-/O+")

# Quadrant-specific colors for the title background
QUAD_TITLE_COLORS <- c(
  "Concordant Up"    = "#E57373",
  "Concordant Down"  = "#64B5F6",
  "Discordant Y+/O-" = "#FFB74D",
  "Discordant Y-/O+" = "#81C784"
)

# --- Build pathway collection for ORA ---
pw_collection <- build_pathway_collection(min_size = 15, max_size = 500,
                                           include_goslim = FALSE,
                                           exclude_variants = TRUE)

# --- Pathway-first chord data for each quadrant ---
MAX_PATHWAYS <- 6
MAX_PROTEINS <- 35

build_chord_data <- function(quad_name) {
  # Get DEPs in this quadrant
  qd <- quad_df %>% filter(quadrant == quad_name, sig_any)

  if (nrow(qd) < 5) {
    message(sprintf("  %s: only %d DEPs, skipping", quad_name, nrow(qd)))
    return(NULL)
  }

  # Run ORA on this quadrant's genes
  ora_res <- tryCatch(
    run_ora_deduplicated(genes = qd$gene, universe = quad_df$gene,
                          pathways = pw_collection, jaccard_cutoff = 0.5,
                          min_size = 15, max_size = 500, padj_cutoff = 0.05),
    error = function(e) { message("  ORA error: ", e$message); tibble() }
  )

  if (nrow(ora_res) == 0) {
    message(sprintf("  %s: no significant pathways", quad_name))
    return(NULL)
  }

  # Take top pathways
  top_pw <- ora_res %>%
    arrange(padj) %>%
    slice_head(n = MAX_PATHWAYS)

  # Extract overlap genes (leading-edge equivalent for ORA)
  # overlapGenes is the intersection of pathway genes and input genes
  chord_links <- top_pw %>%
    rowwise() %>%
    mutate(genes = list(intersect(
      unlist(strsplit(as.character(overlapGenes), ";")),
      qd$gene
    ))) %>%
    ungroup() %>%
    select(pathway, padj, genes) %>%
    unnest(genes) %>%
    rename(gene = genes)

  # Limit proteins shown
  gene_counts <- chord_links %>% count(gene, sort = TRUE)
  if (nrow(gene_counts) > MAX_PROTEINS) {
    keep_genes <- gene_counts$gene[1:MAX_PROTEINS]
    chord_links <- chord_links %>% filter(gene %in% keep_genes)
  }

  # Clean pathway names
  chord_links <- chord_links %>%
    mutate(pathway_label = clean_ring_label(pathway))

  # Add protein metadata
  chord_links <- chord_links %>%
    left_join(qd %>% select(gene, logFC_Training_Young, logFC_Training_Old,
                             module_color),
              by = "gene")

  # Add pathway metadata (NES from fGSEA if available)
  nes_lookup <- fgsea_all %>%
    filter(contrast == "Training_Young") %>%
    select(pathway, NES) %>%
    distinct()
  chord_links <- chord_links %>%
    left_join(nes_lookup, by = "pathway")

  chord_links
}

# --- Render a single chord diagram to a temp PNG ---
render_chord_png <- function(chord_data, quad_name, filepath) {
  if (is.null(chord_data) || nrow(chord_data) == 0) return(NULL)

  # Unique proteins and pathways
  proteins <- unique(chord_data$gene)
  pathways <- unique(chord_data$pathway_label)

  # Sector sizes: proteins get equal width, pathways proportional to gene count
  pw_sizes <- chord_data %>% count(pathway_label) %>% tibble::deframe()
  prot_sizes <- setNames(rep(1, length(proteins)), proteins)

  all_sectors <- c(prot_sizes, pw_sizes)
  sector_order <- c(proteins, rev(pathways))

  # Module colors for protein sectors
  mod_colors <- chord_data %>%
    select(gene, module_color) %>%
    distinct() %>%
    tibble::deframe()

  # NES colors for pathway sectors
  nes_df <- chord_data %>% select(pathway_label, NES) %>% distinct()
  nes_color_map <- setNames(
    sapply(nes_df$NES, function(v) {
      if (is.na(v)) return("grey70")
      v <- pmax(-3, pmin(3, v))
      if (v >= 0) {
        scales::seq_gradient_pal("white", "#D6604D")(v / 3)
      } else {
        scales::seq_gradient_pal("#4393C3", "white")((v + 3) / 3)
      }
    }),
    nes_df$pathway_label
  )

  # Sector colors
  sector_colors <- c(
    setNames(mod_colors[proteins], proteins),
    nes_color_map[pathways]
  )
  # Fix NAs
  sector_colors[is.na(sector_colors)] <- "grey70"

  # logFC color for chords
  lfc_vals <- chord_data %>%
    select(gene, logFC_Training_Young) %>%
    distinct() %>%
    tibble::deframe()

  # Render
  png(filepath, width = PW, height = PH, units = "mm", res = 300)

  circos.clear()
  circos.par(
    gap.after = c(rep(1, length(proteins) - 1), 8,
                  rep(2, length(pathways) - 1), 8),
    start.degree = 90,
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.01, 0.01)
  )

  circos.initialize(
    factors = factor(sector_order, levels = sector_order),
    xlim = cbind(rep(0, length(sector_order)),
                 all_sectors[sector_order])
  )

  # Track 1: sector labels
  circos.track(
    ylim = c(0, 1), track.height = 0.08,
    bg.border = NA,
    panel.fun = function(x, y) {
      sector_name <- CELL_META$sector.index
      # Only label pathways (proteins labeled in track 2)
      if (sector_name %in% pathways) {
        circos.text(CELL_META$xcenter, 0.5,
                    sector_name, facing = "clockwise",
                    niceFacing = TRUE, adj = c(0, 0.5),
                    cex = 0.45, font = 2)
      }
    }
  )

  # Track 2: colored blocks + protein labels
  circos.track(
    ylim = c(0, 1), track.height = 0.06,
    bg.border = "grey50", bg.lwd = 0.3,
    panel.fun = function(x, y) {
      sector_name <- CELL_META$sector.index
      col <- sector_colors[sector_name]
      circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1,
                  col = col, border = NA)
      if (sector_name %in% proteins) {
        circos.text(CELL_META$xcenter, 0.5,
                    sector_name, facing = "clockwise",
                    niceFacing = TRUE, adj = c(0, 0.5),
                    cex = 0.35)
      }
    }
  )

  # Draw chords
  for (i in seq_len(nrow(chord_data))) {
    gene_i <- chord_data$gene[i]
    pw_i   <- chord_data$pathway_label[i]
    lfc_i  <- chord_data$logFC_Training_Young[i]

    # Chord color from logFC
    lfc_clamped <- pmax(-1, pmin(1, lfc_i))
    if (lfc_clamped >= 0) {
      chord_col <- scales::seq_gradient_pal("#FFFFFF", "#D6604D")(lfc_clamped)
    } else {
      chord_col <- scales::seq_gradient_pal("#4393C3", "#FFFFFF")((lfc_clamped + 1))
    }
    chord_col <- scales::alpha(chord_col, 0.4)

    # Find position of this gene in the pathway sector
    pw_genes_in_sector <- chord_data %>%
      filter(pathway_label == pw_i) %>%
      pull(gene)
    gene_idx <- which(pw_genes_in_sector == gene_i)

    tryCatch(
      circos.link(gene_i, c(0, 1),
                  pw_i, c(gene_idx - 1, gene_idx),
                  col = chord_col, border = NA),
      error = function(e) NULL
    )
  }

  circos.clear()
  dev.off()

  filepath
}

# --- Generate all 4 chord diagrams ---
chord_files <- list()
for (quad in QUADRANTS) {
  message(sprintf("Building chord: %s", quad))
  cd <- build_chord_data(quad)

  if (!is.null(cd)) {
    write_csv(cd, file.path(DAT, "panel_G",
                             sprintf("chord_%s.csv",
                                     tolower(gsub("[/ +]", "_", quad)))))
  }

  fname <- file.path(RPT, sprintf("SUPP_panel_G_chord_%s.png",
                                    tolower(gsub("[/ +]", "_", quad))))
  result <- render_chord_png(cd, quad, fname)
  chord_files[[quad]] <- result
}

# --- Composite 2x2 ---
valid_files <- chord_files[!sapply(chord_files, is.null)]

if (length(valid_files) >= 2) {
  # Read PNGs and compose with patchwork
  plots <- lapply(names(valid_files), function(qn) {
    img <- readPNG(valid_files[[qn]])
    g <- rasterGrob(img, interpolate = TRUE)
    wrap_elements(g) + ggtitle(qn) +
      theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5))
  })

  # Pad to 4 if needed
  while (length(plots) < 4) {
    plots <- c(plots, list(plot_spacer()))
  }

  composite <- (plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]])

  COMP_W <- PW * 2 + 10
  COMP_H <- PH * 2 + 20

  ggsave(file.path(RPT, "SUPP_panel_G_chord.png"), composite,
         width = COMP_W, height = COMP_H, units = "mm", dpi = 300)
}

message("F04 SUPP chord diagrams done")
