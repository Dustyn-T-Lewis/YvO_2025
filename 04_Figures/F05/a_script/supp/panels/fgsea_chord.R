# F05 Panel F (SUPP): Merged fGSEA Leading-Edge Chord Diagrams
# One chord per contrast: top 3 UP + top 3 DOWN pathways merged, sorted by NES.
# Pathway arcs colored by NES gradient (blue → white → red), width ∝ gene count.
# Union of leading-edge genes, all labeled, capped ~100.
# ---------------------------------------------------------------------------
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

library(tidyverse)
library(circlize)

set.seed(42)

RPT_PNG <- "04_Figures/F05/b_reports/supp/png/panels"
RPT_PDF <- "04_Figures/F05/b_reports/supp/pdf/panels"
DAT <- "04_Figures/F05/c_data"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_F"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# -- Constants ----------------------------------------------------------------

CONTRASTS <- c("Aging", "Training_Old")
TOP_N_PW <- 3   # per direction (3 UP + 3 DOWN = 6 total)
MAX_GENES <- 100
PANEL_W   <- 210  # mm — full mode (chord + legend)
PANEL_H   <- 180
PANEL_SQ  <- 180  # mm — chordonly mode (square, for composites)
DPI       <- 300

CONTRAST_COL <- c(
  Aging          = "#D6604D",
  Training_Young = "#4393C3",
  Training_Old   = "#F57C00",
  Interaction    = "#7B1FA2"
)
MULTI_COL <- "#388E3C"

RIBBON_ALPHA <- 0.25

PANEL_TITLES <- c(
  Aging          = "Aging",
  Training_Young = "Training (Young)",
  Training_Old   = "Training (Old)",
  Interaction    = "Interaction"
)

SHORT_NAMES <- c(
  Aging          = "Aging",
  Training_Young = "TY",
  Training_Old   = "TO",
  Interaction    = "Inter"
)

# -- Pathway name cleaner (self-contained) -----------------------------------

clean_pw_name <- function(x, max_chars = 30) {
  x <- gsub("^(HALLMARK_|REACTOME_|KEGG_MEDICUS_|KEGG_|GOBP_|GOSLIM_)", "", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  abbrev <- c(
    "Actomyosin Structure Organization" = "Actomyosin Struct. Org.",
    "Degradation of the Extracellular Matrix" = "ECM Degradation",
    "Extracellular Matrix Organization" = "ECM Organization",
    "Assembly of Collagen Fibrils and Other Multimeric Structures" = "Collagen Fibril Assembly",
    "Collagen Biosynthesis and Modifying Enzymes" = "Collagen Biosynthesis",
    "Integrin Cell Surface Interactions" = "Integrin Cell Surface",
    "Non Integrin Membrane Ecm Interactions" = "Non-Integrin ECM",
    "Reference Translation Initiation" = "Translation Initiation",
    "Cellular Response to Acid Chemical" = "Acid Chemical Response",
    "Striated Muscle Contraction" = "Striated Muscle Contr.",
    "Muscle Cell Development" = "Muscle Cell Dev.",
    "Muscle Cell Differentiation" = "Muscle Cell Diff.",
    "Cytoplasmic Translation" = "Cytoplasmic Transl.",
    "Ribonucleoprotein Complex Biogenesis" = "RNP Complex Biogenesis",
    "Proton Motive Force Driven Atp Synthesis" = "PMF-Driven ATP Synth.",
    "Respiratory Electron Transport" = "Respiratory ETC",
    "Oxidative Phosphorylation" = "OXPHOS",
    "Fatty Acid Catabolic Process" = "FA Catabolism",
    "Organic Acid Catabolic Process" = "Organic Acid Catabolism",
    "Membraneless Organelle Assembly" = "Membraneless Org. Assembly"
  )
  x <- ifelse(x %in% names(abbrev), abbrev[x], x)
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3), "..."), x)
}

# -- Color helpers -----------------------------------------------------------

nes_color <- function(nes) {
  # Continuous NES gradient: blue (-3) → white (0) → red (+3)
  nes_clamped <- pmax(-3, pmin(3, nes))
  colorRamp2(c(-3, 0, 3), c("#4393C3", "white", "#D6604D"))(nes_clamped)
}

logfc_color <- function(lfc) {
  lfc_clamped <- pmin(pmax(lfc, -2), 2)
  colorRamp2(c(-2, 0, 2), c("#4393C3", "white", "#D6604D"))(lfc_clamped)
}

contrast_sq_color <- function(n_ctr, memberships) {
  if (is.na(n_ctr) || n_ctr == 0) return("grey80")
  if (n_ctr >= 3) return(MULTI_COL)
  ctrs <- strsplit(as.character(memberships), ",")[[1]]
  if (length(ctrs) == 0) return("grey80")
  CONTRAST_COL[ctrs[1]]
}

padj_fill <- function(p) {
  nlog <- -log10(max(p, 1e-20))
  frac <- pmin(1, pmax(0, (nlog - 1) / 14))
  v <- 1 - (0.25 + 0.55 * frac)
  rgb(v, v, v)
}

# -- Step 1: Load data -------------------------------------------------------

fgsea_cache <- read_csv("04_Figures/shared/fgsea_f1_panel_F.csv",
                         show_col_types = FALSE)
dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv",
                    show_col_types = FALSE)

message(sprintf("fGSEA cache: %d rows, contrasts: %s",
                nrow(fgsea_cache),
                paste(unique(fgsea_cache$contrast), collapse = ", ")))

# -- Step 2: DEP membership (for contrast squares on proteins) ---------------

pi_cols    <- paste0("pi_score_", CONTRASTS)
logfc_cols <- paste0("logFC_", CONTRASTS)

dep_membership <- dep_df %>%
  select(gene, all_of(pi_cols)) %>%
  mutate(across(all_of(pi_cols), ~ !is.na(.) & . < 0.05))
names(dep_membership)[-1] <- CONTRASTS

primary_df <- dep_df %>%
  select(gene, all_of(pi_cols)) %>%
  pivot_longer(-gene, names_to = "contrast_col", values_to = "pi_val") %>%
  mutate(contrast = gsub("pi_score_", "", contrast_col)) %>%
  filter(!is.na(pi_val)) %>%
  group_by(gene) %>%
  slice_min(pi_val, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(gene, primary_contrast = contrast)

membership_str <- dep_membership %>%
  rowwise() %>%
  mutate(
    n_contrasts = sum(c_across(all_of(CONTRASTS))),
    memberships = paste(CONTRASTS[c_across(all_of(CONTRASTS))], collapse = ",")
  ) %>%
  ungroup() %>%
  select(gene, n_contrasts, memberships)

gene_master <- dep_df %>%
  select(gene, all_of(logfc_cols)) %>%
  left_join(primary_df, by = "gene") %>%
  left_join(membership_str, by = "gene") %>%
  mutate(n_contrasts = replace_na(n_contrasts, 0),
         memberships = replace_na(memberships, ""))

# -- Step 3: Build pathway collection for Jaccard dedup ----------------------

pw_collection <- build_pathway_collection(min_size = 15, max_size = 500,
                                           include_goslim = FALSE,
                                           exclude_variants = TRUE)

# -- Step 4: Process fGSEA: merge UP + DOWN per contrast ---------------------

process_contrast_merged <- function(fgsea_sub, contrast_name, logfc_col) {
  sig <- fgsea_sub %>% filter(padj < 0.05)
  if (nrow(sig) == 0) {
    message(sprintf("  %s: 0 sig pathways — skipping", contrast_name))
    return(NULL)
  }

  sig_up   <- sig %>% filter(NES > 0)
  sig_down <- sig %>% filter(NES < 0)

  message(sprintf("  %s: %d UP, %d DOWN sig pathways before dedup",
                  contrast_name, nrow(sig_up), nrow(sig_down)))

  # Dedup each direction separately
  dedup_direction <- function(sig_dir) {
    if (nrow(sig_dir) < 2) return(sig_dir)
    gobp_rows  <- sig_dir %>% filter(database == "GO:BP")
    other_rows <- sig_dir %>% filter(database != "GO:BP")

    if (nrow(gobp_rows) >= 2) {
      gobp_reduced <- tryCatch({
        reduce_gobp_rrvgo(gobp_rows, threshold = 0.7)
      }, error = function(e) {
        message("  rrvgo failed: ", e$message)
        gobp_rows
      })
    } else {
      gobp_reduced <- gobp_rows
    }

    if (nrow(other_rows) >= 2) {
      other_dedup <- deduplicate_enrichment_flat(other_rows, pw_collection,
                                                  jaccard_cutoff = 0.5)
    } else {
      other_dedup <- other_rows
    }

    bind_rows(gobp_reduced, other_dedup)
  }

  up_dedup   <- dedup_direction(sig_up) %>% arrange(desc(NES))
  down_dedup <- dedup_direction(sig_down) %>% arrange(NES)

  # Top N per direction
  top_up   <- up_dedup %>% slice_head(n = TOP_N_PW)
  top_down <- down_dedup %>% slice_head(n = TOP_N_PW)

  # Merge: positive NES first (descending), then negative NES (descending = least negative first)
  top_merged <- bind_rows(top_up, top_down) %>%
    arrange(desc(NES)) %>%
    mutate(pathway_label = clean_pw_name(pathway))

  if (nrow(top_merged) == 0) {
    message(sprintf("  %s: 0 pathways after dedup — skipping", contrast_name))
    return(NULL)
  }

  message(sprintf("  %s: %d merged pathways: %s",
                  contrast_name, nrow(top_merged),
                  paste(sprintf("%s (%.2f)", top_merged$pathway_label, top_merged$NES),
                        collapse = " | ")))

  # Extract leading-edge proteins from all pathways
  le_links <- list()
  for (i in seq_len(nrow(top_merged))) {
    le_str <- top_merged$leadingEdge[i]
    if (is.na(le_str) || le_str == "") next
    genes <- strsplit(as.character(le_str), ";\\s*")[[1]]
    genes <- genes[genes != "" & !is.na(genes)]
    if (length(genes) > 0) {
      le_links[[i]] <- tibble(
        gene = genes,
        pathway_label = top_merged$pathway_label[i],
        pathway = top_merged$pathway[i]
      )
    }
  }
  if (length(le_links) == 0) return(NULL)
  link_df <- bind_rows(le_links)

  # Merge with gene master for logFC and membership
  chord_df <- link_df %>%
    left_join(gene_master %>% select(gene, all_of(logfc_col),
                                      primary_contrast, n_contrasts, memberships),
              by = "gene") %>%
    rename(logFC = !!logfc_col) %>%
    filter(!is.na(logFC))

  if (nrow(chord_df) == 0) return(NULL)

  # Cap genes at MAX_GENES: prioritize genes in multiple pathways, then |logFC|
  unique_genes <- chord_df %>%
    group_by(gene) %>%
    summarise(n_pw = n_distinct(pathway_label),
              abs_lfc = max(abs(logFC)), .groups = "drop") %>%
    arrange(desc(n_pw), desc(abs_lfc))

  if (nrow(unique_genes) > MAX_GENES) {
    keep_genes <- unique_genes$gene[seq_len(MAX_GENES)]
    chord_df <- chord_df %>% filter(gene %in% keep_genes)
    message(sprintf("  %s: capped at %d genes (from %d)",
                    contrast_name, MAX_GENES, nrow(unique_genes)))
  }

  pw_meta <- top_merged %>%
    select(pathway_label, padj, NES, size) %>%
    rename(n_genes = size)

  list(chord = chord_df, pw = pw_meta)
}

# Process each contrast
chord_data_list <- list()
for (ctr in CONTRASTS) {
  fgsea_sub <- fgsea_cache %>% filter(contrast == ctr)
  lfc_col <- paste0("logFC_", ctr)
  chord_data_list[[ctr]] <- process_contrast_merged(fgsea_sub, ctr, lfc_col)
}

n_panels <- sum(!vapply(chord_data_list, is.null, logical(1)))
message(sprintf("\n%d merged panels to render (of %d possible)",
                n_panels, length(chord_data_list)))

# -- Save chord data CSVs ---------------------------------------------------

for (ctr in names(chord_data_list)) {
  cd <- chord_data_list[[ctr]]
  if (!is.null(cd)) {
    write_csv(cd$chord %>% select(-any_of("pathway")),
              file.path(DAT, "panel_F", paste0("fgsea_chord_", ctr, ".csv")))
    write_csv(cd$pw,
              file.path(DAT, "panel_F", paste0("fgsea_chord_", ctr, "_pathways.csv")))
  }
}
message("Chord data saved to c_data/panel_F/")

# -- Step 5: Render merged chord diagrams ------------------------------------

render_merged_chord <- function(chord_data, pw_data, contrast_name,
                                 filepath_base, mode = "full") {
  display_name <- PANEL_TITLES[contrast_name]
  if (is.null(chord_data) || nrow(chord_data) == 0) {
    message(sprintf("  %s: empty — skipping", contrast_name))
    return(invisible(NULL))
  }

  # --- Pathways sorted by NES descending (most positive → most negative) ---
  pathways <- pw_data %>% arrange(desc(NES)) %>%
    pull(pathway_label) %>% unique()
  n_pw <- length(pathways)
  if (n_pw == 0) return(invisible(NULL))

  pw_padj <- setNames(pw_data$padj[match(pathways, pw_data$pathway_label)], pathways)
  pw_nes  <- setNames(pw_data$NES[match(pathways, pw_data$pathway_label)], pathways)

  # NES-based fill for each pathway arc
  pw_nes_col <- setNames(vapply(pw_nes, function(n) nes_color(n), character(1)),
                          pathways)

  # --- Links: protein <-> pathway membership ---
  link_df <- chord_data %>%
    filter(pathway_label %in% pathways) %>%
    distinct(gene, pathway_label)

  # --- Order proteins by primary pathway (highest NES), then |logFC| desc ---
  gene_primary_pw <- link_df %>%
    left_join(pw_data %>% select(pathway_label, NES), by = "pathway_label") %>%
    group_by(gene) %>%
    slice_max(NES, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(gene, primary_pw = pathway_label)

  lfc_lookup <- chord_data %>% distinct(gene, logFC) %>% tibble::deframe()

  proteins <- gene_primary_pw %>%
    mutate(pw_rank = match(primary_pw, pathways),
           abs_lfc = abs(lfc_lookup[gene])) %>%
    arrange(pw_rank, desc(abs_lfc)) %>%
    pull(gene) %>% unique()

  n_prot <- length(proteins)
  if (n_prot == 0) return(invisible(NULL))

  # --- Contrast membership lookup ---
  membership_lookup <- chord_data %>%
    distinct(gene, n_contrasts, memberships)
  nctr_lookup <- setNames(membership_lookup$n_contrasts, membership_lookup$gene)
  memb_lookup <- setNames(membership_lookup$memberships, membership_lookup$gene)

  # --- Sector sizes: proteins = 1, pathways = leading-edge count (min 8) ---
  pw_le_count <- link_df %>% count(pathway_label) %>% tibble::deframe()
  pw_sizes <- pmax(pw_le_count[pathways], 8)
  pw_sizes[is.na(pw_sizes)] <- 8

  all_sectors <- c(proteins, pathways)
  sector_sizes <- setNames(c(rep(1, n_prot), pw_sizes), all_sectors)

  # --- Gaps: per-pathway spacing, larger transition gaps ---
  prot_gap <- max(0.2, 12 / n_prot)
  gap_vec <- c(rep(prot_gap, n_prot - 1), 15,
               rep(3, max(n_pw - 1, 0)), 15)

  # --- Dynamic font size based on protein count ---
  gene_cex <- if (n_prot > 80) 0.35
              else if (n_prot > 60) 0.40
              else if (n_prot > 40) 0.45
              else 0.50
  pw_cex <- if (n_pw > 8) 0.45
            else if (n_pw > 5) 0.50
            else 0.55

  # --- Render ---
  for (ext in c("png")) {
    outfile <- paste0(filepath_base, ".", ext)
    pw <- if (mode == "chordonly") PANEL_SQ else PANEL_W
    ph <- if (mode == "chordonly") PANEL_SQ else PANEL_H
    png(outfile, width = pw, height = ph, units = "mm", res = DPI,
        bg = "white")

    if (mode == "chordonly") {
      par(mar = c(0, 0, 4, 0))
    } else {
      # Layout: 72% chord, 28% legend
      layout(matrix(c(1, 2), ncol = 2), widths = c(0.72, 0.28))
      par(mar = c(0, 0, 3, 0))
    }

    circos.clear()
    circos.par(
      gap.after = gap_vec,
      start.degree = 90, clock.wise = TRUE,
      cell.padding = c(0, 0, 0, 0),
      track.margin = c(0.005, 0.005)
    )
    circos.initialize(
      factors = factor(all_sectors, levels = all_sectors),
      xlim = cbind(rep(0, length(all_sectors)), sector_sizes)
    )

    # TRACK A (outermost, thin): contrast identity | padj grey
    circos.track(ylim = c(0, 1), track.height = 0.020, bg.border = NA,
      panel.fun = function(x, y) {
        s <- CELL_META$sector.index
        if (s %in% proteins) {
          col <- contrast_sq_color(nctr_lookup[s], memb_lookup[s])
          circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1,
                      col = col, border = NA)
        } else if (s %in% pathways) {
          col <- padj_fill(pw_padj[s])
          circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1,
                      col = col, border = NA)
        }
      })

    # TRACK 1 (labels): gene symbol | pathway name on NES-colored bg
    circos.track(ylim = c(0, 1), track.height = 0.14, bg.border = NA,
      panel.fun = function(x, y) {
        s <- CELL_META$sector.index
        if (s %in% proteins) {
          circos.text(CELL_META$xcenter, 0.5, s,
                      cex = gene_cex, facing = "clockwise", niceFacing = TRUE,
                      font = 1)
        } else if (s %in% pathways) {
          bg_col <- pw_nes_col[s]
          circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1,
                      col = bg_col, border = NA)
          # White text on saturated colors, black on near-white
          abs_nes <- abs(pw_nes[s])
          txt_col <- if (abs_nes > 0.8) "white" else "grey20"
          circos.text(CELL_META$xcenter, 0.5, s,
                      cex = pw_cex, font = 2, col = txt_col,
                      facing = "bending.inside", niceFacing = TRUE)
        }
      })

    # TRACK 2 (inner): logFC bar (proteins) | NES color (pathways)
    circos.track(ylim = c(0, 1), track.height = 0.065, bg.border = NA,
      panel.fun = function(x, y) {
        s <- CELL_META$sector.index
        if (s %in% proteins) {
          col <- logfc_color(lfc_lookup[s])
          circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1,
                      col = col, border = NA)
        } else if (s %in% pathways) {
          circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1,
                      col = pw_nes_col[s], border = NA)
        }
      })

    # RIBBONS: protein Track 2 → pathway Track 2
    for (i in seq_len(nrow(link_df))) {
      gene_i <- link_df$gene[i]
      pw_i   <- link_df$pathway_label[i]
      if (!gene_i %in% proteins) next

      pw_genes <- link_df %>% filter(pathway_label == pw_i) %>% pull(gene)
      pw_genes <- pw_genes[pw_genes %in% proteins]
      idx <- which(pw_genes == gene_i)
      if (length(idx) == 0) next

      pw_xlim <- get.cell.meta.data("cell.xlim", sector.index = pw_i,
                                     track.index = 3)
      slot_w <- (pw_xlim[2] - pw_xlim[1]) / length(pw_genes)

      # Ribbon colored by pathway NES
      ribbon_col <- adjustcolor(pw_nes_col[pw_i], alpha.f = RIBBON_ALPHA)

      tryCatch(
        circos.link(
          gene_i, c(0, 1),
          pw_i, c(pw_xlim[1] + (idx - 1) * slot_w,
                  pw_xlim[1] + idx * slot_w),
          col = ribbon_col, border = NA,
          rou1 = get.cell.meta.data("cell.bottom.radius",
                                     sector.index = gene_i, track.index = 3),
          rou2 = get.cell.meta.data("cell.bottom.radius",
                                     sector.index = pw_i, track.index = 3)
        ),
        error = function(e) NULL
      )
    }

    if (mode == "chordonly") {
      title(main = display_name, cex.main = 1.6, font.main = 2, line = 1.8)
      mtext(sprintf("%d proteins \u00b7 %d pathways \u00b7 %d links",
                    n_prot, n_pw, nrow(link_df)),
            side = 3, line = 0.4, cex = 0.85, font = 3, col = "grey30")
    } else {
      title(main = display_name, cex.main = 1.4, font.main = 2, line = 0.5)
    }
    circos.clear()

    if (mode == "chordonly") {
      dev.off()
      next
    }

    # --- Legend panel ---
    par(mar = c(2, 0.5, 2, 1))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))

    # NES gradient
    text(0, 0.97, "Pathway NES", adj = 0, cex = 0.8, font = 2)
    n_grad <- 50
    nes_seq <- seq(-3, 3, length.out = n_grad)
    grad_cols <- vapply(nes_seq, function(n) nes_color(n), character(1))
    for (j in seq_len(n_grad)) {
      rect((j - 1) / n_grad * 0.8, 0.90, j / n_grad * 0.8, 0.94,
           col = grad_cols[j], border = NA)
    }
    text(0, 0.88, "-3", adj = 0, cex = 0.6)
    text(0.4, 0.88, "0", adj = 0.5, cex = 0.6)
    text(0.8, 0.88, "+3", adj = 1, cex = 0.6)

    # logFC gradient
    text(0, 0.82, "Protein logFC", adj = 0, cex = 0.8, font = 2)
    lfc_seq <- seq(-2, 2, length.out = n_grad)
    lfc_cols <- vapply(lfc_seq, function(l) logfc_color(l), character(1))
    for (j in seq_len(n_grad)) {
      rect((j - 1) / n_grad * 0.8, 0.74, j / n_grad * 0.8, 0.78,
           col = lfc_cols[j], border = NA)
    }
    text(0, 0.72, "-2", adj = 0, cex = 0.6)
    text(0.4, 0.72, "0", adj = 0.5, cex = 0.6)
    text(0.8, 0.72, "+2", adj = 1, cex = 0.6)

    # Contrast colors
    text(0, 0.64, "Primary Contrast", adj = 0, cex = 0.8, font = 2)
    ctr_show <- CONTRAST_COL
    y_pos <- seq(0.58, by = -0.05, length.out = length(ctr_show))
    for (k in seq_along(ctr_show)) {
      rect(0, y_pos[k] - 0.015, 0.06, y_pos[k] + 0.015,
           col = ctr_show[k], border = NA)
      text(0.09, y_pos[k], names(ctr_show)[k], adj = 0, cex = 0.55)
    }
    # Multi-contrast
    rect(0, y_pos[length(ctr_show)] - 0.065, 0.06,
         y_pos[length(ctr_show)] - 0.035, col = MULTI_COL, border = NA)
    text(0.09, y_pos[length(ctr_show)] - 0.05, "Multi (3+)", adj = 0, cex = 0.55)

    # Stats
    text(0, 0.20, sprintf("%d proteins", n_prot), adj = 0, cex = 0.65)
    text(0, 0.15, sprintf("%d pathways", n_pw), adj = 0, cex = 0.65)
    text(0, 0.10, sprintf("%d links", nrow(link_df)), adj = 0, cex = 0.65)

    dev.off()
  }

  message(sprintf("  Rendered: %s (%d proteins, %d pathways, %d links)",
                  contrast_name, n_prot, n_pw, nrow(link_df)))
}

# --- Render each contrast ---
panels_rendered <- character(0)

for (ctr in CONTRASTS) {
  cd <- chord_data_list[[ctr]]
  if (is.null(cd)) {
    message(sprintf("  %s: skipped (no data)", ctr))
    next
  }

  short <- SHORT_NAMES[ctr]
  fbase <- file.path(RPT_PNG, paste0("SUPP_fgsea_chord_", short))

  render_merged_chord(cd$chord, cd$pw, ctr, fbase, mode = "full")
  render_merged_chord(cd$chord, cd$pw, ctr,
                      paste0(fbase, "_sq"), mode = "chordonly")
  panels_rendered <- c(panels_rendered, ctr)
}

message(sprintf("\nRendered %d merged panels: %s", length(panels_rendered),
                paste(panels_rendered, collapse = ", ")))
message("\nDone: F05 Panel F fGSEA Leading-Edge Chord Diagrams (SUPP)")
