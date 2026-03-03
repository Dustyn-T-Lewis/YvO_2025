# 05_Figures -- Refined Figure Pipeline (F1-F3)

Self-contained rewrites of Figures 1-3 from the YvO proteomics analysis.

## Purpose

Refined versions of 04_Figures/F1-F3 with:
- Consistent styling via shared/style.R
- Trimmed, readable panel scripts
- Self-contained (no runtime dependency on 04_Figures)

## Directory Structure

```
05_Figures/
  shared/
    style.R              # Palettes, sizing constants, FIG_THEME, helpers
    volcano_ring.R       # Volcano ring plot builder (used by F1-F3)
  F1/
    a_script/
      YvO_F1_setup.R             # Setup: packages, data loading, F1 palettes
      panel_A.R                  # CV% violins
      panel_B.R                  # logFC density histograms
      panel_C.R                  # PCA biplot + PERMANOVA
      panel_D.R                  # DEPs per contrast (stacked bar)
      panel_E.R                  # UpSet plot
      panel_F.R                  # fGSEA grouped bar chart
      supp_S1_7.R                # Pi-score distributions (supplementary)
      YvO_figure1_composite.R    # Composite assembly + Excel
    b_reports/                   # PDF/PNG outputs
    c_data/                      # Panel data CSVs + supplementary Excel
  F2/
    a_script/
      YvO_F2_setup.R             # Setup: packages, data loading, F2 palettes
      panel_A.R                  # Volcano ring: Training (Young)
      panel_B.R                  # Volcano ring: Training (Old)
      panel_C.R                  # Volcano ring: Interaction
      panel_D.R                  # Concordance scatter (logFC x logFC)
      panel_E.R                  # RRHO2 + ORA bars
      panel_F.R                  # Heatmap -> Sankey -> pathway bars
      YvO_figure2_composite.R    # Composite assembly + Excel
    b_reports/                   # PDF/PNG outputs
    c_data/                      # Panel data CSVs + supplementary Excel
  F3/
    a_script/
      YvO_F3_setup.R             # Setup: packages, data, reversal tests
      panel_A.R                  # Volcano ring: Aging Effect
      panel_B.R                  # Volcano ring: Reversal Effect
      panel_C.R                  # Reversal scatter (logFC Aging vs Training Old)
      panel_D.R                  # RRHO2 reversal map + ORA bars
      panel_E.R                  # Classification heatmap + Sankey + enrichment bars
      YvO_figure3_composite.R    # Composite assembly + Excel
    b_reports/                   # PDF/PNG outputs
    c_data/                      # Panel data CSVs + supplementary Excel
```

## Running

```bash
cd /path/to/A_YvO_2025
Rscript 05_Figures/F1/a_script/YvO_figure1_composite.R
Rscript 05_Figures/F2/a_script/YvO_figure2_composite.R
Rscript 05_Figures/F3/a_script/YvO_figure3_composite.R
```

Each composite script sources its setup and all panel scripts, then assembles
the final multi-panel figure and supplementary Excel workbook. Individual panel
scripts can also be sourced standalone for development.

## Dependencies

R packages: tidyverse, patchwork, ggplot2, ggrepel, ggforce, ggnewscale,
scales, grid, fgsea, msigdbr, rrvgo, GOSemSim, org.Hs.eg.db, GO.db,
clusterProfiler, ComplexHeatmap, vegan, ggsignif, boot, openxlsx

## Relationship to 04_Figures

04_Figures is the original reference implementation. 05_Figures reads the same
upstream data (01_normalization, 02_Imputation, 03_DEP) and the fGSEA cache
from 04_Figures/F2/c_data/shared/, but does not source any 04_Figures scripts.
