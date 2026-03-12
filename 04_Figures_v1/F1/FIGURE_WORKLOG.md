# Figure 1 — Refinement Worklog

## Figure purpose
**Proteomic Landscape of Aging and Resistance Training in Human Skeletal Muscle.**
Foundational overview establishing data quality, contrast magnitudes, DEP counts, intersection structure, and pathway enrichment. Introduces the blunted Old training response.

---

## Current file inventory

### Scripts (`a_script/`)
| File | Role |
|------|------|
| `YvO_F1_setup.R` | Setup: loads DAList, DEP results, contrast data |
| `panel_A.R` | CV% violin plots |
| `panel_B.R` | logFC density histograms |
| `panel_C.R` | PCA biplot + PERMANOVA |
| `panel_D.R` | DEP count bar chart |
| `panel_E.R` | UpSet intersection plot |
| `panel_F.R` | fGSEA pathway enrichment counts |
| `supp_S1_7.R` | Supplementary pi-score distributions |
| `YvO_figure1_composite.R` | Assembly script |

### Outputs (`b_reports/`)
- `Figure_1.pdf` — composite
- `panel_A_cv.pdf`
- `panel_B_logfc_density.pdf`
- `panel_C_pca.pdf`
- `panel_D_dep_counts.pdf`
- `panel_E_upset.pdf`
- `panel_F_fgsea.pdf`
- `supplementary/S1_7_pi_score_distributions.pdf`

### Data (`c_data/`)
- `audit_panel_A_*.csv` (median CV CI, Wilcoxon effects)
- `audit_panel_B_*.csv` (Cliff's delta, median LFC CI)
- `audit_panel_C_*.csv` (betadisper, PCA variance CI)
- `audit_panel_D_*.csv` (DEP fraction CI)
- `audit_panel_E_*.csv` (overlap enrichment)
- `audit_panel_F_*.csv` (NES summary)
- `06_panel_F_fgsea_results.csv`
- `F1_supplementary.xlsx`

## Panel inventory

| Panel | Status | Key message |
|-------|--------|-------------|
| A | Functional, effect sizes annotated | CV comparable across groups (22-26%) |
| B | Functional, rendering OK | Aging > Tr Young > Tr Old in logFC magnitude |
| C | Functional | Young/Old separation on PC1, PERMANOVA significant |
| D | Functional, rendering OK | Aging 9.3%, Tr Young 5.0%, Tr Old 0.8% |
| E | Functional, subtitle improved | Most DEPs contrast-specific; Tr Young x Interaction enriched |
| F | Functional, no pathway names | Tr Young most enriched, Tr Old fewest |
| S1.7 | Includes unexplained Reversal contrast | Pi-score distribution validation |

## Known issues

### Unresolved (require user decision)
- [ ] Panel F: pathway counts only, no pathway identity (editorial decision)
- [ ] Supp S1.7: Reversal contrast appears without main-figure context (editorial decision)

### Resolved
- [x] Panel B: chi-squared Unicode renders as ".^2" → fixed with `bquote()` plotmath expression
- [x] Panel D: pi symbol renders as "." in bar labels → fixed with `parse = TRUE` on `geom_text()` + plotmath label strings
- [x] Panel D: pseudo-log y-axis unlabeled → added "(pseudo-log scale)" to axis label
- [x] Panel A: significance brackets overstate trivial effects → added Cliff's |d| alongside stars
- [x] Panel E: subtitle π symbol fixed via `bquote()` plotmath; subtitle slightly shortened
- [x] Panel E: Unicode minus sign warnings eliminated (ASCII hyphens in legend)
- [x] Supp S1.7: title π symbol fixed via `expression()` plotmath
- [x] Panel letters A–F confirmed present and consistent in composite

## Changes made

### 2026-03-01 — Iteration 1: Unicode rendering fixes + pseudo-log axis label
**Files modified:**
- `a_script/panel_B.R` (lines 145-151): Replaced Unicode `\u03C7\u00B2` (χ²) in `sprintf()` caption with `bquote()` plotmath expression using `chi^2`. Removed `face = "italic"` from caption theme since plotmath handles formatting.
- `a_script/panel_D.R` (line 114): Changed `THRESH_LABEL` from Unicode pi/≤ strings to plotmath-parseable strings (`"pi <= 0.05"`). Added `parse = TRUE` to `geom_text()` call (line 140). Changed y-axis label from `"% of proteome"` to `"% of proteome (pseudo-log scale)"` (line 151).

**Outputs regenerated:** panel_B_logfc_density.pdf/png, panel_D_dep_counts.pdf/png, Figure_1.pdf/png, F1_supplementary.xlsx

**Result:** All three symbols (χ², π, ≤) now render correctly in PDF output. Pseudo-log transform is explicitly labeled.

### 2026-03-01 — Iteration 2: Panel A significance bracket context
**Files modified:**
- `a_script/panel_A.R` (lines 132-143): Modified `.signif_label()` to include Cliff's delta value alongside stars. Bracket annotations now show e.g. "*** |d| = 0.07" instead of bare "***", making negligible effect sizes visible to readers.

**Outputs regenerated:** panel_A_cv.pdf/png, Figure_1.pdf/png, F1_supplementary.xlsx

**Result:** Brackets now convey both statistical significance and practical significance. All three annotated comparisons have |d| < 0.13, clearly indicating negligible effects despite large-N significance.

### 2026-03-01 — Iteration 3: Consistent π symbol across all panels + supplementary
**Files modified:**
- `a_script/panel_E.R` (line 217): Replaced `sprintf()` subtitle with `bquote()` plotmath expression so π renders correctly via plotmath (Unicode π fails in R's PDF device `mbcsToSbcs`). Also slightly shortened "proteins excluded" → "excluded".
- `a_script/supp_S1_7.R` (line 44): Replaced plain-text "Pi-score" title with `expression(bold("S1.7  ") * pi * bold("-score distributions"))` so π renders via plotmath.

**Outputs regenerated:** panel_E_upset.pdf/png, S1_7_pi_score_distributions.pdf/png, Figure_1.pdf/png, F1_supplementary.xlsx

**Result:** π symbol now consistent across Panel D (bar labels), Panel E (subtitle), and S1.7 (title + axis labels). All use plotmath rendering to bypass PDF device Unicode limitations.

### 2026-03-01 — Iteration 4: Eliminate Unicode minus sign warnings in Panel E legend
**Files modified:**
- `a_script/panel_E.R` (lines 54-57): Replaced Unicode minus signs `\u2212` (−) with ASCII hyphens `-` in contrast definition legend labels. The PDF device was substituting hyphens anyway, so this is a no-visual-change fix that eliminates 12 `mbcsToSbcs` warnings per composite build.

**Outputs regenerated:** panel_E_upset.pdf/png, Figure_1.pdf/png, F1_supplementary.xlsx

**Result:** Unicode warnings reduced from 30 to 18 (remaining 18 are all upstream GOSemSim/rrvgo library deprecation warnings, not our code).

### 2026-03-01 — Iteration 5: Legibility test at journal scale
**No files modified.** Test render at 190mm (double-column) reveals severe truncation: panel titles, subtitles, Panel B caption, Panel C annotation, and key legend are all clipped or overlapping. The 380mm design cannot be scaled to standard journal width without splitting into two sub-figures or significant layout redesign. This is flagged as a critical pre-submission issue requiring user guidance on target journal format.

### 2026-03-01 — Iteration 6: Supplementary XLSX data completeness
**Files modified:**
- `a_script/YvO_figure1_composite.R` (line 64-65): Added `ci_lo`, `ci_hi` columns to Sheet B (`B_effect_sizes`) in the supplementary XLSX. Previously only exported `contrast, med_abs_lfc, n_above_05` — now includes the bootstrap 95% CI on median |logFC| per contrast.

**Outputs regenerated:** F1_supplementary.xlsx, Figure_1.pdf/png (no visual change)

**Result:** XLSX Sheet B now includes bootstrap CI data, matching the audit CSV and the panel annotations. All 7 XLSX sheets verified complete.

## Open questions
- ~~Should Panel A brackets be removed or annotated with effect sizes?~~ Resolved: annotated with Cliff's delta
- Should Panel F add top 2-3 pathway names per contrast?
- Should Reversal contrast in S1.7 be explained or removed?
- Should composite be split into two half-figures (A-C / D-F)?

## Publication-readiness checklist
- [x] All Unicode symbols render correctly in PDF
- [x] Panel letters A–F visible and consistent
- [x] Font sizes legible at journal scale (base_size=8, titles=9, subtitles=6.5)
- [x] Axis labels complete and accurate
- [x] Legends complete (contrast definitions + direction key in Panel E area)
- [x] Colors follow shared/palettes.R (GROUP_FILL, DIR_COLORS, THEME_PUB from palettes; CONTRAST_COLORS F1-specific)
- [ ] Composite dimensions appropriate for target journal (currently 380 x 260 mm — **FAILS at 190mm double-column scale**: titles truncated, annotations clipped, legend cut off. Needs split or redesign.)
- [x] No overclaiming in annotations
- [x] All audit CSVs current (regenerated 2026-03-01)
- [x] Supplementary XLSX current (regenerated 2026-03-01)

## Next recommended steps
1. ~~Reassess Panel A significance brackets~~ Done (Iteration 2)
2. Panel E: subtitle overloaded — consider restructuring (user decision)
3. Panel F: pathway counts only — consider adding top pathway names (user decision)
4. Supp S1.7: Reversal contrast needs justification or removal (user decision)
5. ~~Panel letter labels in composite~~ Done (confirmed A–F present)
6. Composite dimensions — need target journal specs (user decision)
