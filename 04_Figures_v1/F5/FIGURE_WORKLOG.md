# Figure 5 — Refinement Worklog

## Figure purpose
**WGCNA Co-Expression Network Architecture.**
Systems-level characterization of the skeletal muscle proteome: 9 co-expression modules mapped to biological function, phenotypic traits, and hub proteins. Demonstrates module preservation across training. Validates all module-based inference in F6.

---

## Current file inventory

### Scripts (`a_script/`)
| File | Role |
|------|------|
| `YvO_F5_setup.R` | Setup: loads WGCNA outputs, MEs, trait correlations, kME, GO enrichment, biological module labels |
| `YvO_WGCNA_run.R` | Upstream WGCNA execution — do NOT modify unless re-running WGCNA |
| `panel_A.R` | Dendrogram + module color bars |
| `panel_B.R` | Module-trait correlation heatmap (now with biological labels) |
| `panel_C.R` | Multi-view module characterization triptych (5 cols x 6 modules) |
| `panel_D.R` | GO enrichment dotplot |
| `panel_E.R` | Hub protein networks (2x3 grid) |
| `panel_F.R` | Module preservation (compacted horizontal bar chart) |
| `panel_supp_overlap.R` | FCM-WGCNA alluvial overlap |
| `YvO_figure5_composite.R` | Assembly script (main + supplementary) |

### Outputs (`b_reports/`)
- `Figure_5.pdf` — main figure (180 x 240 mm)
- `Figure_S5.pdf` — supplementary figure
- `panel_A_dendrogram.pdf`
- `panel_B_heatmap.pdf`
- `panel_C_triptych.pdf`
- `panel_D_go_enrichment.pdf`
- `panel_E_hub_network.pdf`
- `panel_F_preservation.pdf`
- `supp_fcm_wgcna_overlap.pdf`
- `age_gap_closure.pdf`
- `wgcna/01_soft_threshold.pdf`
- `wgcna/02_module_dendrogram.pdf`
- `wgcna/03_module_trait_heatmap.pdf`

### Data (`c_data/`)
- `01_panel_A_dendrogram_data.csv`
- `02_panel_B_*.csv` (heatmap data, correlation CIs)
- `03_panel_C_*.csv` (eigengene, heatmap z-scores, Sankey links, triptych enrichment)
- `04_panel_D_enrichment_data.csv`
- `05_panel_E_*.csv` (hub network, hub CIs, context ring)
- `06_panel_F_preservation.csv`
- `07_module_reversal.csv`
- `08_age_gap_closure.csv`
- `09_fcm_wgcna_overlap.csv`
- `10_cross_method_integration.csv`
- `Figure_5_data.xlsx`
- `wgcna/` — module assignments, correlations, p-values, GO, hub proteins, network.rds

## Panel inventory

| Panel | Status | Main vs Supp | Key message |
|-------|--------|-------------|-------------|
| A | Functional | Main | 9 modules from 2,124 proteins |
| B | Biological labels applied | Main | Blue r=0.72 age, Brown r=0.40 training |
| C | Full triptych in supp; simplified in main | Split done | Per-module functional characterization |
| D | Functional | Supp | GO enrichment for 4 modules |
| E | Dense, illegible at print | Supp | Hub networks for 6 modules |
| F | Compacted to bar chart | Main | All Zsummary > 10 |
| Overlap | Supplementary | Supp | FCM-WGCNA alluvial |
| Age gap | Non-significant (p=0.91) | Supp | No systematic gap closure |

## Biological module label mapping
Centralized in `YvO_F5_setup.R` as `mod_bio_labels` and `mod_display_label()`.

| Color | Biology | Key trait |
|-------|---------|-----------|
| Blue | Catabolic Metabolism | Age r = 0.72 |
| Brown | Cytoskeletal Remodeling | Age r = -0.52, Training r = 0.40 |
| Turquoise | Protein Modification | Age r = -0.46 |
| Green | Muscle Structure | Age r = -0.36, Type II fCSA r = 0.35 |
| Black | Translation Machinery | — |
| Pink | Immune/Complement | — |
| Yellow | Ubiquitin-Proteasome | — |
| Red | (to be characterized) | — |
| Magenta | (to be characterized) | — |

## Known issues

### Resolved
- [x] Composite 600x800 mm — restructured to 180x240 mm main + supplementary
- [x] Panel C: 5 columns x 6 modules too complex for main figure — simplified version in main (eigengene slopes + trait strip only)
- [x] Module names are color defaults — biological labels added via `mod_bio_labels` / `mod_display_label()` in setup
- [x] Panel F visually uninformative — compacted to horizontal bar chart with biological labels
- [x] Assembly script needs refactoring — now produces Figure_5.pdf (main) + Figure_S5.pdf (supp)
- [x] Hub networks illegible in composite — moved to supplementary

### Unresolved
- [ ] Panels C and D show different module subsets without explanation
- [ ] Age gap closure non-significant (p = 0.91) — needs context in legend
- [ ] No module shows significant eigengene reversal after BH correction — noted in data
- [ ] Panel B uses n=62 for all trait correlations (liberal for phenotypes with n=44-50) — documented in subtitle
- [ ] Color scale limits differ between Panel B (-0.8 to 0.8) and Panel C trait strip (-1 to 1)

## Changes made

### Session 2026-03-01
1. **`YvO_F5_setup.R`**: Added `mod_bio_labels` named vector and `mod_display_label()` helper function for biological module labels throughout all panels.
2. **`panel_B.R`**: Module sidebar now displays "Biology (Color)" labels via `mod_display_label()` instead of raw color names.
3. **`panel_F.R`**: Replaced scatter plot with compact horizontal bar chart. Added biological module labels. Reduced save dimensions to 180x120 mm.
4. **`YvO_figure5_composite.R`**: Complete restructure:
   - Main figure (180 x 240 mm): A (dendrogram) / B (heatmap) + simplified C (eigengene slopes + trait strip) / F (bar chart)
   - Supplementary figure S5: Full triptych (C) + GO dotplot (D) + Hub networks (E)
   - Simplified Panel C built inline using `build_simple_row()` function

## Open questions
- Should Red and Magenta modules get biological labels? (need GO enrichment review)
- Should the age gap closure and FCM-WGCNA overlap be included in the supplementary composite?

## Publication-readiness checklist
- [x] Main figure fits journal full-page format (~180 x 240 mm)
- [x] Biological module labels present in B, C (simplified), F
- [ ] All text legible at print scale (needs visual verification after render)
- [x] Heatmap readable
- [x] Preservation panel compacted
- [x] Supplementary content separated
- [ ] Age gap closure result contextualized
- [ ] Module reversal null result noted
- [x] Colors follow shared/palettes.R
- [ ] Panel letters visible (needs verification)

## Next recommended steps
1. Render and visually inspect both main and supplementary figures
2. Characterize Red and Magenta modules for biological labels
3. Add age gap closure context to figure legend
4. Address Panel B/C color scale discrepancy
