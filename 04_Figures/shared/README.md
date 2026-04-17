# `04_Figures/shared/` — Shared Figure Assets

Cross-figure resources used by F02–F07 panel scripts.

## Files

| File | Role | Active writer | Consumers |
|---|---|---|---|
| `style.R` | Palettes (`DIR_COLORS`, `DB_COLORS`, `PAL_CLASS`), `FIG_THEME`, `AGE_COLORS`, `get_pdf_device()` helper | curated | every panel + stitcher |
| `figure_supplement_helpers.R` | 7-function Excel supplement API: `add_sheet`, `safe_read`, `build_workbook`, `cleanup_after_workbook`, `read_sheet_df`, `read_matrix_sheet`, `matrix_to_df` | curated | every `90_stitch_figure.R` |
| `pathway_utils.R` | `build_pathway_collection()`, `run_fgsea_deduplicated()`, `run_ora_deduplicated()`, `classify_database()`, Jaccard deduplication, `CONSOLIDATED_PATHWAY_ORDER`/`_COLORS` | curated | F04 (panels A/B/D/E + supp chord/ora_chord/pattern_ora/rrho_manual/enrichment_heatmap), F05 (same + prepare_data Reversal append), F06 (`YvO_WGCNA_run.R`, panel_B_module_nes_scatter, supp panel_B_triptych + panel_D_hub) |
| `volcano_ring.R` | `select_ring_terms()`, `build_ring_with_gaps()`, `build_nes_legend_bar()`, `make_volcano_ring()` — volcano-ring rendering primitives | curated | F03 panels A/B/C/D + stitcher; F04/F05 supp `chord_diagrams.R` (shared ring-term selection) |
| `go_slim_categories.R` | GO Slim BP assignment (15 consolidated pathway labels); `bp_slim`, `SLIM_CONSOLIDATED`, `assign_go_slim_consolidated()` | curated | F04 + F05 only (stitchers + `panel_C_pattern_heatmap`); also re-exported through `pathway_utils.R` |
| `fgsea_tstat_all_v2.csv` | **Frozen fGSEA cache** (see below) | — (frozen) | 9 consumers (F02 panel_E; F03 A/B/C/D; F04 panel_B_nes_scatter + supp chord; F05 panel_B_nes_scatter + supp chord) |
| `fgsea_f1_panel_F.csv` | Frozen fGSEA subset for the `panel_F_fgsea_chord` supplementary panel (the "f1" suffix is historical — pre-2026-04-02 figure numbering — not related to the current F01 figure) | — (frozen) | F04 + F05 `supp/panels/panel_F_fgsea_chord.R` |

## `fgsea_tstat_all_v2.csv` — the frozen pathway cache

**Schema (10 columns):** `pathway, pval, padj, log2err, ES, NES, size, leadingEdge, database, contrast`

**Content:** 10,885 data rows covering:
- **Contrasts:** Aging, Training_Young, Training_Old, Interaction, Reversal
- **Databases:** Hallmark, KEGG, Reactome, GO:BP, GO Slim

**Status: frozen in git.** There is no single script that regenerates this file from scratch. It exists here as a cross-figure resource because:

1. Full regeneration requires `msigdbr` v25 + ~10 min of fGSEA across 4 databases × 4 base contrasts, plus a separate Reversal-contrast run.
2. The pre-refactor monolith that originally built it (`04_Figures/F03/a_script/.archive_panel_A_volcanoes.R`) writes to a *different* path (`04_Figures/F03/c_data/shared/`), not this canonical location. It has not been updated since the 2026-04-15 refactor moved the cache to `04_Figures/shared/`.
3. `04_Figures/F05/a_script/panels/prepare_data.R:236-264` is the only script that currently writes to this canonical path, and it only appends the Reversal contrast when missing — it assumes the other four contrasts are already present.

**To restore if accidentally deleted or corrupted:**

```bash
git restore 04_Figures/shared/fgsea_tstat_all_v2.csv
```

**To regenerate fully (known work, not currently wired up):**

1. Run `03_DEP/a_script/01_run_dep.R` to get `03_combined_results.csv` with fresh t-statistics per contrast.
2. For each of the 4 base contrasts (Aging, Training_Young, Training_Old, Interaction): rank by t-statistic, run `run_fgsea_deduplicated()` (from `pathway_utils.R`) against the pathway collection from `build_pathway_collection(min_size = 10, max_size = 500)`, tag with `contrast` and `database`.
3. Run `F05/a_script/panels/prepare_data.R` to append the Reversal contrast and write to this path.

Lifting steps 1–2 into a proper `build_fgsea_cache.R` at this location is the right long-term fix and is tracked in `docs/README.md` open items.

## Pathway-analysis usage across figures

Single reference for "where does pathway analysis show up in the manuscript figures" and which shared asset each consumer uses. Five distinct analysis types appear across F02–F07.

### fGSEA — rank-based enrichment

Input: t-statistics from `03_DEP/c_data/03_combined_results.csv`, ranked per contrast. Output: NES, pval, padj per pathway × database. Cache: `fgsea_tstat_all_v2.csv` (frozen).

| Figure / panel | Role in manuscript |
|---|---|
| F02 panel_E | Stacked-bar summary of N significant pathways per contrast × database (overview) |
| F03 panels A–D | Volcano rings — top-N enriched pathways radiating around each contrast's volcano (Aging, Training_Young, Training_Old, Interaction) |
| F04 panel_B (`panel_B_nes_scatter.R`) | NES scatter Training_Young vs Training_Old — pathway-level concordance |
| F04 supp `chord_diagrams.R` | Pathway-level chord edges (cross-training convergence) |
| F04 supp `panel_F_fgsea_chord.R` | Legacy fGSEA chord (uses `fgsea_f1_panel_F.csv` subset) |
| F05 panel_B (`panel_B_nes_scatter.R`) | NES scatter Aging vs Training_Old — pathway-level reversal |
| F05 supp `chord_diagrams.R` | Pathway-level chord edges (reversal chord) |
| F05 supp `panel_F_fgsea_chord.R` | Legacy fGSEA chord (uses `fgsea_f1_panel_F.csv` subset) |

### ORA — over-representation (hypergeometric) on DEP lists

Input: DEP gene lists (Pi-score < 0.05 or FDR-significant) from each contrast, tested against pathway collections. Uses `run_ora_deduplicated()` from `pathway_utils.R`.

| Figure / panel | Role |
|---|---|
| F04 panel_A (`panel_A_ORA.R`) | Quadrant ORA — 4 concordance quadrants (Up-Up, Up-Down, Down-Up, Down-Down) of Training_Young × Training_Old DEPs |
| F04 panel_D flanking bars | ORA on fry driving-protein sets |
| F04 panel_E (`panel_E_rrho2.R`) | ORA on RRHO2 concordant / discordant quadrants |
| F04 supp `pattern_ora.R`, `ora_chord_supp.R`, `enrichment_heatmap.R` | Pattern-level ORA + chord + multi-method heatmap |
| F05 panel_A | Quadrant ORA — Reversed Up/Down, Exacerbated Up/Down (Aging × Training_Old) |
| F05 panel_D flanking bars | ORA on fry driving proteins |
| F05 panel_E | ORA on reversal quadrants |
| F05 supp `pattern_ora.R`, `ora_chord_supp.R`, `enrichment_heatmap.R` | Mirror F04 supp structure for reversal |

### fry — competitive gene-set rotation test (limma)

Input: full limma fit object (`03_DEP/c_data/01_limma_DAList.rds`). Tests whether a gene set shows a coordinated direction more than expected by chance *relative to* the rest of the measured proteome.

| Figure / panel | Role |
|---|---|
| F04 panel_D (`panel_D_fry.R`) | fry barcode + driving-protein set for Training concordance |
| F05 panel_D (`panel_D_fry.R`) | fry for Aging × Training_Old reversal |

### RRHO2 — rank-rank hypergeometric overlap

Input: two ranked gene lists (by t-stat). Uses a pure-R `phyper()` implementation (the RedRibbon package segfaults on the ~2,100-gene matrix — see MEMORY).

| Figure / panel | Role |
|---|---|
| F04 panel_E (`panel_E_rrho2.R`) | RRHO2 heatmap of Training_Young × Training_Old gene-rank overlap |
| F05 panel_E (`panel_E_rrho2.R`) | RRHO2 heatmap of Aging × Training_Old (reversal geometry) |

### WGCNA module enrichment

Input: per-module gene sets from `04_Figures/F06/c_data/wgcna/`. Multi-database enrichment (Hallmark + KEGG + Reactome + GO:BP; `min_size = 15`, no GO Slim — intentionally broader than the figure-level cache). Uses `pathway_utils.R` inside `YvO_WGCNA_run.R`.

| Figure / panel | Role |
|---|---|
| F06 pipeline runner (`YvO_WGCNA_run.R`) | Produces `wgcna_module_enrichment.csv` (module × pathway NES + padj) |
| F06 panel_A heatmap | Incorporates module-level NES annotations into the module-trait heatmap |
| F06 panel_B_module_nes_scatter | Training concordance + aging reversal at the *module* level (not protein level) |
| F06 supp `panel_B_triptych` | Per-module triptychs with enrichment context |
| F06 supp `panel_D_hub` | Hub-network overlays annotated by pathway enrichment |

### Figures that do NOT use pathway analysis

- **F00 (pipeline QC)** — filter / missingness / PCA / normalization QC only.
- **F01 (phenotype)** — training volume, DXA LBM, VL thickness, deadlift 1RM, Type-II fCSA — phenotype metrics, no protein-level tests.
- **F07 (phenotype prediction)** — uses WGCNA *module eigengenes* from F06 as features for ROC / classifier analysis; no fresh pathway inference.

## Archive

`shared/.archive_orphans_20260416/` contains earlier orphans moved out of the live surface during the 2026-04-16 hygiene pass:

| File | Why archived |
|---|---|
| `build_pipeline_supplement.R` | Generated `S1_pipeline.xlsx` (pipeline-level supplement) but was not sourced by any figure stitcher; standalone utility with no active consumer. |
| `S1_pipeline.xlsx` | Output of the above; not referenced by figure pipeline. |
| `pilot_dedup_comparison.csv` | Exploratory deduplication comparison; not consumed by any live script. |
| `fgsea_gobp_v2.csv` | Earlier GO:BP-only fGSEA cache; superseded by the multi-database `fgsea_tstat_all_v2.csv`. |

Nothing in `.archive_orphans_20260416/` is loaded by any live script. Git history is preserved.

## Ownership

Changes to files in this directory affect every downstream figure. Follow the Operating Rules change types (especially #6 naming-harmonization + #7 validation-enhancement). Do not change pathway-collection assembly, fGSEA parameters, Jaccard cutoff, or colour palettes casually — they cascade across the entire figure set.
