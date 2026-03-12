# Figure 3 — Refinement Worklog

## Figure purpose
**Training-Mediated Reversal of the Aging Skeletal Muscle Proteome.**
Central claim figure providing convergent evidence (protein, threshold-free, pathway) that resistance training in older adults partially reverses aging-associated proteomic changes. This is the climax of the first half of the paper.

**Status: Currently the strongest figure in the manuscript.** Refinements should preserve strengths and sharpen precision of claims.

---

## Current file inventory

### Scripts (`a_script/`)
| File | Role |
|------|------|
| `YvO_F3_setup.R` | Setup: loads DEP results, computes reversal metrics, runs fGSEA |
| `panel_AB.R` | Volcano rings for Aging (A) and Reversal (B) |
| `panel_C.R` | Reversal scatter (logFC aging vs training old) |
| `panel_D.R` | RRHO heatmap |
| `panel_E.R` | Pathway-level NES reversal scatter |
| `panel_F.R` | Top 30 reversal DEPs classification |
| `YvO_figure3_composite.R` | Assembly script |

### Outputs (`b_reports/`)
- `Figure_3.pdf` — composite
- `panel_A_volcano.pdf`
- `panel_B_volcano.pdf`
- `panel_C_reversal_scatter.pdf`
- `panel_D_rrho2.pdf`
- `panel_E_nes_scatter.pdf`
- `panel_F_classification.pdf`
- `melov_reversal_permutation.pdf` — supplementary

### Data (`c_data/`)
- `fgsea_tstat_all.csv`
- `panel_A/volcano_aging.csv`, `ring_terms.csv`
- `panel_B/volcano_reversal.csv`, `ring_terms.csv`
- `panel_C/reversal_scatter.csv`
- `panel_D/rrho2_matrix.csv`, `rrho2_summary.csv`
- `panel_E/nes_scatter.csv`
- `panel_F/sankey_links.csv`, `classification.csv`, `enrichment_bars.csv`
- `reversal_tests/melov_permutation.csv`, `reversal_contingency.csv`, `signed_reversal_score.csv`
- `F3_supplementary.xlsx`

## Panel inventory

| Panel | Status | Key message |
|-------|--------|-------------|
| A | Publication-ready | 198 Aging DEPs; OXPHOS/catabolism up, translation/myogenesis down |
| B | Publication-ready | 180 Reversal DEPs; enrichment pattern largely inverted vs A |
| C | Publication-ready | r = -0.35 [-0.38, -0.31], ρ = -0.33, 67% reversal [64, 70] |
| D | Publication-ready | Reversal and exacerbation hotspots both present |
| E | Publication-ready | Pathway r = -0.61 [-0.78, -0.36], 79% reversed (n=38) |
| F | Functional but somewhat circular | Top 30: 26 Reversed + 4 Attenuated |
| Melov (supp) | Non-significant (p = 0.14) | 10.1% reversal — tension with main panels |

## Known issues

### Unresolved
- [ ] RRHO exacerbation signal comparable to reversal but visually underemphasized
- [ ] Panel F is circular (filtered to reversal DEPs, then shows they're reversed)
- [ ] Panel E: n=38 pathways, CI [-0.78, -0.36] is wide — inherent limitation
- [ ] Multiple color encodings for "reversal" across panels

### Resolved
- [x] Melov non-significance acknowledged in Panel C → annotated along anti-diagonal: "Magnitude reversal: 10.1%, p = 0.14 (n.s.)" (2026-03-01)
- [x] Unicode characters (ρ, π, −, en-space) rendering as dots in PDF → replaced with ASCII equivalents (2026-03-01)
- [x] Panel F legend text overlaps in composite → shortened labels ("Trn. Old"/"Trn. Young"), tightened layout (2026-03-01)
- [x] Panel F "F" letter not visible in composite → moved title from plot_annotation() to dumbbell labs() (2026-03-01)
- [x] Panel C subtitle truncated in composite → split into 2 lines (2026-03-01)
- [x] Panel E subtitle truncated in composite → split into 2 lines (2026-03-01)
- [x] No panel letter labels (A–F) in composite → added "A  Title" format to all panels (2026-03-01); upgraded to `labs(tag)` system with 16pt bold (2026-03-01)
- [x] Bottom key strip truncated ("Sig Training only" cut off) → added explicit xlim (2026-03-01)
- [x] Typography inconsistent across panels → applied 5-tier system matching F2 (2026-03-01)
- [x] Composite dimensions too dense → expanded 380x350 to 440x400mm, rebalanced row widths/heights (2026-03-01)

## Changes made

### Session 2026-03-01: Publication-readiness fixes (iterations 1–7)

**Iteration 1–3: Subtitle truncation, panel letters, key strip**
- `panel_C.R` — Split subtitle into 2 lines; added "C  " prefix
- `panel_E.R` — Split subtitle into 2 lines; added "E  " prefix
- `panel_AB.R` — Added "A  " / "B  " prefixes to volcano ring titles
- `panel_D.R` — Added "D  " prefix
- `panel_F.R` — Added "F  " prefix via plot_annotation
- `YvO_figure3_composite.R` — Added explicit `xlim` to key strip `coord_cartesian()`

**Iteration 4–5: Panel F title and legend**
- `panel_F.R` — Moved title/subtitle from `plot_annotation()` to dumbbell `labs()`; removed nested `plot_annotation()` (fixes invisible "F" in composite)
- `panel_F.R` — Shortened legend labels: "Training (Old)" → "Trn. Old", "Training (Young)" → "Trn. Young"; moved Response column from x=3.0 to x=2.5; tightened xlim from 7 to 5.5

**Iteration 6–7: Unicode rendering + Melov annotation**
- `panel_C.R` — Replaced `\u03c1` (ρ) with "rho", `\u2002` (en-space) with double-space in quadrant labels
- `panel_E.R` — Replaced `\u2002` (en-space) with double-space in quadrant labels
- `panel_F.R` — Replaced `\u03C0` (π) with "pi" in subtitle
- `panel_AB.R` — Replaced `\u2212` (−) with "-" in Panel B title
- `panel_C.R` — Added Melov magnitude-reversal annotation along anti-diagonal: "Magnitude reversal: 10.1%, p = 0.14 (n.s.) / (Melov permutation, 576 aging-sig. proteins)". Saved Melov variables (`melov_rev_pct`, `melov_p`) before local `reversal_pct` overwrite.

**Outputs regenerated:**
- All standalone panels (A-F PDFs + PNGs)
- `Figure_3.pdf` and `Figure_3.png` (composite)
- `F3_supplementary.xlsx`

**Verification:**
- Panel C: subtitle fully visible, "rho" renders correctly ✓
- Panel E: subtitle fully visible, quadrant labels clean ✓
- Panel F: "F" letter visible, "pi" renders correctly, legend labels legible ✓
- Panel B: "Training - Aging" renders without warnings ✓
- Panel letters A–F: all visible in composite ✓
- Key strip: "Sig Training only" fully visible ✓
- No Unicode conversion warnings in build output ✓

### Session 2026-03-01: Combined layout + typography refinement

**Typography tier system applied (matching F2):**

| Tier | Role | Size | Face |
|------|------|------|------|
| 1 | Panel tags (A-F) | 16pt theme | bold |
| 2 | Axis titles, legend titles | 11pt theme / 10pt geom | bold |
| 3 | Axis tick text, legend text | 10pt theme / 9pt geom | plain |
| 4 | Quadrant labels, key headers, bar counts, Sankey pathway names | 3.5 geom | bold |
| 5 | Gene/pathway repel labels, direction labels, Melov annotation | 2.8 geom | italic (genes) / plain (other) |

**Files modified:**

- `panel_AB.R` — Converted "A  Title" / "B  Title" to `labs(tag = "A/B")` + stripped letter prefix from titles
- `panel_C.R` — Gene labels: 2.2->2.8, bold->italic; max.overlaps 30->40; box.padding 0.5->0.6; quadrant labels 2.5->3.5; Melov annotation 2.0->2.8; legend.text 5.5->9; legend.title 6->10; imp key text 2.0/1.8->2.8; added `labs(tag = "C")`
- `panel_D.R` — Colorbar title 5.5->9 + added bold face; reversed quadrant labels 2.0->3.5; exacerbated quadrant labels 2.5->3.5; axis direction labels 1.8->2.8; added `labs(tag = "D")`
- `panel_E.R` — Pathway labels 2.2->2.8; max.overlaps 30->40; quadrant labels 2.5->3.5; added `labs(tag = "E")`
- `panel_F.R` — Gene names axis.text.y 6.5->9; pattern group labels 2.2->2.8; Sankey pathway labels 2.8->3.5; bar count labels 2.5->3.5; legend titles 2.8->3.5; legend items 2.5->2.8; added `labs(tag = "F")`
- `YvO_figure3_composite.R` — Title 11->14, subtitle 8->10; plot.tag 16pt bold; strip titles/subtitles in composite (kept for standalone); row1 widths c(0.33,0.33,0.34)->c(0.30,0.30,0.40); heights c(0.44,0.44,0.12)->c(0.46,0.46,0.08); dimensions 380x350->440x400mm; axis title bold override on panels C/D/E; pF uses `&` operator to recursively strip titles

**Outputs regenerated:**
- All standalone panels (A-F PDFs + PNGs)
- `Figure_3.pdf` and `Figure_3.png` (composite, 440x400mm)
- `F3_supplementary.xlsx`

**Verification:**
- Panel tags A-F: all visible at 16pt bold in composite
- Quadrant labels: uniformly 3.5pt bold across C/D/E panels
- Gene labels (Panel C): italic at 2.8pt, no excessive overlap with max.overlaps=40
- Pathway labels (Panel E): bold at 2.8pt, max.overlaps=40
- RRHO (Panel D): quadrant labels all 3.5pt (was mixed 2.0/2.5), direction labels 2.8pt
- Panel F: gene names legible at 9pt, Sankey labels 3.5pt, bar counts 3.5pt
- Key strip: compact at 8% height, fully readable (KEY_TITLE/KEY_ITEM from palettes.R)
- Composite title/subtitle: 14/10pt, properly centered
- No regressions in standalone panel PDFs

## Key integrity notes
- **Directional reversal** (r = -0.35 protein, r = -0.61 pathway) is well supported
- **Magnitude reversal** (Melov: 10.1%, p = 0.14) is NOT significant
- This distinction (direction vs magnitude) must be maintained in the figure and caption
- The exacerbation signal in RRHO should be noted, not hidden
- Panel F should ideally show ALL aging DEPs classified, not just top 30 reversal

## Open questions
- Should Melov result be annotated in Panel C or as a separate sub-panel?
- Should Panel C add a loess/regression fit instead of (or alongside) the anti-diagonal?
- Should Panel F show classification of all aging-significant proteins?

## Publication-readiness checklist
- [x] Panel letters A–F visible
- [x] All subtitles visible in composite (Panels C and E fixed)
- [x] Unicode characters render correctly in PDF (ASCII replacements)
- [x] Panel F legend labels legible (shortened to "Trn. Old"/"Trn. Young")
- [x] Melov result acknowledged in Panel C (annotation along anti-diagonal)
- [x] Reversal framed as partial, not complete (subtitle shows r = -0.35, 67%)
- [ ] Exacerbation signal noted
- [x] Volcano ring labels legible
- [x] Colors follow shared/palettes.R
- [x] Key strip fully visible
- [x] Composite dimensions rebalanced (440x400mm with wider Panel C allocation and reduced key strip)
- [x] Typography tiers consistent with F2 (5-tier system applied)

### 2026-03-03 — F2/F3 Panel Harmonization & Standardization

**Objective:** Bring F2 and F3 to consistent publication standard with shared infrastructure, unified key panels, and harmonized themes/typography across all panels.

**Shared infrastructure changes (also documented in F2 worklog):**
- `shared/palettes.R` — added THEME_FIG (base_size=14) and TXT_* typography constants
- `shared/go_slim_categories.R` — **new file**: GO Slim super-category shared utility

**Setup changes (`YvO_F3_setup.R`):**
- Removed local THEME_FIG and TXT_* definitions (now from shared/palettes.R)
- Added `source("04_Figures/shared/go_slim_categories.R")`
- Kept figure-specific colors (SIG_COLORS with Reversal terminology, ORA_QUAD_COLORS with Reversed/Exacerbated framing)

**Panel-level changes:**

| Panel | Changes |
|-------|---------|
| panel_C.R | Removed legend (`legend.position = "none"`), switched to `cairo_pdf` |
| panel_D.R | Removed heatmap legend, ORA text sizes → `rel(0.7)`, added PNG output |
| panel_E.R | Switched THEME_PUB → THEME_FIG, label size → TXT_GENE, quadrant size → TXT_QUADRANT, removed hand-built legend, switched to `cairo_pdf` |
| panel_F.R | **Major methodology change**: replaced ORA+greedy pathway mapping with GO Slim super-categories via shared `assign_go_slim_super()`. Eliminates "Muscle Cell Development absorbs 13/30" problem. Updated pathway colors to SUPER_COLORS, title/subtitle updated, switched to `cairo_pdf` |

**Composite changes (`YvO_figure3_composite.R`):**
- Complete rewrite with new 100×100 virtual grid layout
- Dimensions: 350 × 420 mm (was 480 × 380 mm)
- Added unified key panel at bottom (rows 86-100) with F3-specific sections: Significance (Reversal/Sig Both/Sig Aging/Sig Training/NS), Scatter Quadrants (Reversal/Exacerbation), RRHO Overlap, ORA Quadrants (4-color Reversed/Exacerbated), Dumbbell Contrasts (Aging/Training Old), GO Slim Super-categories (13)
- No title stripping — panels keep their own titles/subtitles
- Panel letter tags via `labs(tag=)` at TXT_TAG size (20pt bold)
- Layout: A|C|E / B|C|E / D(span)|E / Key — 2 volcanos stacked left, scatter mid, classification full right, RRHO+ORA spans left+mid bottom, key full-width bottom
- Supplementary Excel workbook preserved

**Outputs to regenerate:** All F3 panels + Figure_3.pdf/png + F3_supplementary.xlsx

---

### Session: 2026-03-03 — Typography Unification

**Goal:** Unify all annotation text sizes across F2/F3 panels for consistent readability.

**Shared infrastructure (`shared/palettes.R`):**
- Replaced 8 separate TXT_* constants (3.5–4.5 range) with 2 unified constants: `TXT_TAG = 18`, `TXT_LABEL = 3.5`
- Legacy aliases preserve backward compatibility

**F3 panel_D.R:**
- "No enrichment" placeholder: size=3→TXT_LABEL (3.5)

**F3 panel_F.R:**
- Pattern group labels: 4.5→TXT_LABEL, pathway labels: 5.5→TXT_LABEL, count labels: 5.5→TXT_LABEL
- Key section headers: 5.5→TXT_LABEL, key item labels: 4.5–5→TXT_LABEL
- Removed hardcoded axis.title.x size=12 (THEME_FIG provides 14)

**Automatic via aliases (no code changes needed):**
- panel_C.R gene/quadrant/stats labels (TXT_GENE/STATS 4.5→3.5)
- panel_D.R quadrant labels, ORA bar text, direction labels (TXT_ORA_BAR/STATS 4.5→3.5)
- panel_E.R pathway/quadrant labels (TXT_GENE 4.5→3.5)

**Verified:** F3 composite renders successfully with all changes.

### Session: 2026-03-03 — Panel F Text Size Compensation & Layout Widening

**Goal:** Fix Panel F text being too small in composite (root cause: wrap_elements() shrinks 350–450mm standalone to ~200mm, making 3.5mm text → ~1.5mm). Also widen composite to give Panel E (classification) more space.

**Shared infrastructure (`shared/palettes.R`):**
- Added `TXT_PF = 7.0` (annotation text for Panel F, compensates for ~0.48× shrinkage → ~3.4mm effective)
- Added `TXT_PF_GENE = 20` (gene symbol axis text in pt, shrinks to ~10pt in composite)

**F3 panel_F.R:**
- All `TXT_LABEL` (3.5mm) → `TXT_PF` (7.0mm) for: pattern group labels, pathway labels, count labels, all key text
- Gene axis text: `size = 12` → `size = TXT_PF_GENE` (20pt)
- Key relocated from below-bar to overlaying upper-left of bar area (patchwork area overlap)
- Key given white background with grey border for visual separation
- Patchwork design: panels now use full 30 rows (was 26+4 for separate key row)
- Standalone dimensions: 350 → 450mm width, min height 200 → 250mm

**F3 composite (`YvO_figure3_composite.R`):**
- Panel E (classification) allocation: columns 59-100 → 50-100 (51% width, was 42%)
- Other panels adjusted: A,B cols 1-24, C cols 25-49, D cols 1-49
- Composite width: 350 → 420mm (height unchanged at 420mm)

**Outputs regenerated:** All F3 panels + Figure_3.pdf/png (420×420mm)

**Verified:**
- Composite renders without errors
- Panel E gets 51% of 420mm = ~214mm allocation
- Panel F text scales to ~3.4mm effective in composite (matches other panels' TXT_LABEL = 3.5mm)
- Gene symbols scale to ~10pt effective (readable)
- Key is overlaid on bar area, not wasting vertical space

### Session: 2026-03-03 — Panel F Gene Count Bars + Panel G Per-Contrast ORA

**Goal:** (1) Switch Panel F enrichment bars from fold enrichment to gene counts. (2) Rename "super-categories" to "consolidated pathways". (3) Add new Panel G (now Panel F in composite): per-contrast ORA dot plot.

**F3 panel_F.R:**
- Bars now show gene count (not fold enrichment): `bar_total_w <- (pw_row$gene_count / max_count) * S_MAX_LEN`
- Axis label: "Gene count" with integer ticks
- Removed fold=1.0 reference line
- S_MAX_LEN: 2.8 → 3.4 (wider bars)
- Standalone width: 280 → 320mm
- Terminology: "super-category" → "consolidated pathway" in comments, messages, subtitle
- Uses `CONSOLIDATED_COLORS` for pathway bars (via legacy alias `SUPER_COLORS`)

**F3 panel_G.R (new):**
- 2 direction-collapsed gene lists: Aging (Up+Down), Reversal (Rev. Age Up + Rev. Age Down)
- 4 databases: Hallmark (MSigDB H), GO:BP, GO:MF, GO:CC via clusterProfiler
- Top 5 terms per contrast per database, BH correction
- Dot plot: x = -log10(padj), y = term, color = contrast, size = gene count
- Faceted by database
- Exports: `c_data/panel_G/ora_per_contrast.csv`, `b_reports/panel_G_ora.pdf`
- F3 results: 17 significant terms across 5 contrast-database pairs

**F3 composite (`YvO_figure3_composite.R`):**
- Sources panel_G.R → pG, wrapped as Panel F with tag "F"
- Layout: added row 73-90 for Panel F (ORA), key moved to 91-100
- Composite dimensions: 420 × 500mm (was 420 × 420mm)
- Key section: "GO Slim Super-categories" → "Consolidated Pathways"
- Supplementary workbook: added `F_ora_per_contrast` sheet

**Outputs regenerated:** All F3 panels + Panel G + Figure_3.pdf/png + F3_supplementary.xlsx

**Verified:**
- Panel F: 246 proteins, 11 consolidated pathways, gene count bars scale correctly (largest = Metabolism n=56), integer axis ticks, sankey totals match (223 OK)
- Panel G: 17 ORA terms across 5 contrast-database pairs
- Composite renders with all 6 panels (A–F) + key, no errors

## Next recommended steps
1. Assess exacerbation signal framing in RRHO panel
2. Visual review of new 420×500mm composite at print scale
3. Consider whether ORA results warrant further filtering or term deduplication
