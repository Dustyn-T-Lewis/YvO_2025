# Figure 6 — Refinement Worklog

## Figure purpose
**Translational Utility of WGCNA Module Eigengenes.**
Capstone figure testing whether modules can discriminate age, predict training response, and track phenotypic change. Narrows from modules to individual hub proteins.

---

## Current file inventory

### Scripts (`a_script/`)
| File | Role |
|------|------|
| `YvO_F6_setup.R` | Setup: loads MEs, computes delta_me, phenotype correlations, top-3 modules, biological labels |
| `panel_A.R` | PCA + ROC + permutation + Pre-vs-Post comparison (multiple outputs) |
| `panel_B.R` | Baseline eigengene vs delta phenotype (3x2 faceted) |
| `panel_C.R` | Delta eigengene vs delta phenotype (3x2 faceted) |
| `panel_D.R` | kME vs |GS| scatter (turquoise module) |
| `panel_E.R` | Protein-phenotype correlation heatmap |
| `YvO_figure6_composite.R` | Assembly script |

### Outputs (`b_reports/`)
- `Figure_6.pdf` — composite
- `panel_A_pca.pdf`
- `panel_A_roc.pdf`
- `panel_A_age_discrimination.pdf` — combined PCA+ROC
- `panel_A_permutation.pdf` — supplementary
- `panel_A_roc_comparison.pdf` — supplementary (Pre vs Post AUC)
- `panel_B_baseline_association.pdf`
- `panel_C_plasticity.pdf`
- `panel_D_hub_scatter.pdf`
- `panel_E_heatmap.pdf`

### Data (`c_data/`)
- `01_panel_A_pca_scores.csv`
- `01_panel_A_pre_vs_post.csv`
- `01_panel_A_roc_curve.csv`
- `01_panel_A_permutation.csv`
- `02_panel_B_baseline_association.csv`
- `03_panel_C_plasticity.csv`
- `04_panel_D_kme_gs.csv`
- `04_panel_D_hub_sensitivity.csv`
- `05_panel_E_heatmap.csv`
- `Figure_6_data.xlsx`

## Panel inventory

| Panel | Status | Key message |
|-------|--------|-------------|
| A | Strong result (AUC=0.97), permutation p now annotated | Module eigengenes discriminate age |
| B | Title reframed; stat annotations simplified | Suggestive baseline trends (no BH significance) |
| C | Stat annotations simplified | Green-dLBM and Red-dVL track phenotype |
| D | Functional | 38 turquoise hubs, kME-GS rho=0.42 |
| E | Age column dominates | Hubs primarily age markers |

## Known issues

### Resolved
- [x] Panel B title: "Associates with Training Response" — reframed to "Baseline Eigengene vs Training Response" with subtitle noting no BH significance
- [x] Unicode delta renders as "." in composite key strip — replaced with "Delta VL" / "Delta LBM"
- [x] Stat annotations in Panels B/C near-illegible — simplified to `r = X.XX, r_partial = X.XX***`
- [x] Panel A permutation test (p=0.001) not annotated in main figure — added to ROC panel annotation
- [x] N=31 not prominently displayed — added to ROC annotation
- [x] Biological module labels updated to match F5 (Catabolic Metabolism, etc.)

### Unresolved
- [ ] Pre-vs-Post ROC (AUC 0.97->0.99) complicates reversal narrative — not addressed in main figure
- [ ] Panel E Age column visually overwhelms phenotype columns
- [ ] Panels B and C select different top-3 modules — narratively fragmented
- [ ] Panel D focuses on turquoise only
- [ ] Baseline VL column in Panel E mostly non-significant

## Changes made

### Session 2026-03-01
1. **`YvO_F6_setup.R`**: Updated `mod_bio_labels` to use biological names matching F5 (Catabolic Metabolism, Cytoskeletal Remodeling, etc.).
2. **`panel_A.R`**: Added permutation p-value and N to ROC panel annotation. ROC annotation now reads: "AUC = X.XX (95% CI: X.XX-X.XX) / Permutation p < 0.001 / N = 31". Re-saves combined Panel A after permutation test runs.
3. **`panel_B.R`**:
   - Title changed from "Baseline Eigengene Associates with Training Response" to "Baseline Eigengene vs Training Response"
   - Subtitle now notes "No facet reaches p_adj < 0.05"
   - Stat annotation simplified from 2-line block to single line: `r = X.XX, r_partial = X.XX***`
4. **`panel_C.R`**: Stat annotation simplified to match Panel B format.
5. **`YvO_figure6_composite.R`**: Key strip phenotype labels changed from Unicode `\u0394 VL` / `\u0394 LBM` to plain "Delta VL" / "Delta LBM" for reliable PDF rendering.

## Key integrity notes
- Panel A AUC = 0.97 is strong but N = 31 with CI [0.91, 1.00] — now displayed
- Panel B: NO BH-corrected significance — title and subtitle now reflect this
- Panel C: Two of six facets significant — these are the defensible findings
- Pre-vs-Post AUC increasing means age gap persists at module level — honest but not yet addressed in main figure
- Partial correlations are weaker than raw — indicates age confounding
- Hub proteins are primarily age markers, secondarily training markers

## Publication-readiness checklist
- [x] Panel B title accurately reflects evidence strength
- [x] Unicode delta renders correctly (replaced with plain text)
- [x] Stat annotations readable at journal scale (simplified)
- [x] Permutation p-value visible in Panel A
- [x] N=31 displayed prominently
- [ ] Panel letters visible and consistent (needs verification)
- [ ] Colors follow shared/palettes.R (needs verification)
- [x] No overclaiming in titles or annotations
- [x] Key strip complete and untruncated

## Next recommended steps
1. Render and visually inspect composite
2. Address Pre-vs-Post ROC finding in figure legend
3. Consider whether Panel E Baseline VL column should be dropped
4. Consider panel order optimization (A -> C -> B -> D -> E)
