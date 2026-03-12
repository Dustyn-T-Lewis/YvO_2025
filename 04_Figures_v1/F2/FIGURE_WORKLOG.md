# Figure 2 — Refinement Worklog

## Figure purpose
**Training-Effect Concordance Across Age Groups.**
Central comparison figure showing the degree of concordance versus divergence in the training response between Young and Old, at protein, threshold-free, and pathway levels. Establishes the primary hypothesis that aging impairs the proteomic training response.

---

## Current file inventory

### Scripts (`a_script/`)
| File | Role |
|------|------|
| `YvO_F2_setup.R` | Setup: loads DEP results, runs fGSEA, prepares contrasts |
| `panel_AB.R` | Volcano rings for Training Young (A) and Training Old (B) |
| `panel_C.R` | logFC concordance scatter |
| `panel_D.R` | RRHO heatmap |
| `panel_E.R` | NES bubble scatter (standalone, removed from composite) |
| `panel_F.R` | Interaction DEPs triptych (standalone, removed from composite) |
| `YvO_figure2_composite.R` | Assembly script |

### Outputs (`b_reports/`)
- `Figure_2.pdf` — composite
- `panel_A_volcano.pdf`
- `panel_B_volcano.pdf`
- `panel_C_concordance.pdf`
- `panel_D_rrho2.pdf`
- `panel_E_nes_bubble.pdf`
- `panel_F_interaction.pdf`

### Data (`c_data/`)
- `panel_A/volcano_young.csv`, `ring_terms.csv`, `ring_terms_old.csv`
- `panel_B/volcano_old.csv`
- `panel_C/concordance.csv`, `concordance_stats.csv`
- `panel_D/rrho2_matrix.csv`, `rrho2_summary.csv`
- `panel_E/nes_scatter.csv`, `nes_scatter_stats.csv`
- `panel_F/interaction_classification.csv`, `interaction_dot.csv`, `interaction_dot_long.csv`, `interaction_patterns.csv`, `pattern_frequency_test.csv`, `sankey_dot.csv`, `sankey_links.csv`, `srplot_input.csv`
- `shared/fgsea_tstat_all_v2.csv`
- `F2_supplementary.xlsx`

## Panel inventory

| Panel | Status | Key message |
|-------|--------|-------------|
| A | Publication-ready | Young: 106 DEPs, broad functional enrichment |
| B | Publication-ready | Old: 17 DEPs, attenuated response |
| C | Publication-ready | r = 0.26, concordance = 65%, symmetric axes |
| D | Publication-ready | Concordant and discordant hotspots of comparable magnitude + ORA |
| E | Removed from composite (standalone still available) | NES r = 0.34 [0.10, 0.54] |
| F | Removed from composite (standalone still available) | 7-category expanded heatmap |

## Known issues

### Unresolved
- [ ] Panels A/B show different pathway selections — direct comparison harder (deliberate design choice, see code comments in volcano_ring.R line 694)
- [ ] Panel D RRHO colorbar legend small in composite (visible and readable; viridis pattern + quadrant labels provide sufficient context)

### Resolved
- [x] Panel C subtitle truncated in composite — split to two lines in panel_C.R (2026-03-01)
- [x] Composite key strip "Sig Old only" truncated — tightened spacing (box_w 0.35->0.30, item_gap 0.6->0.45) and added explicit xlim in composite.R (2026-03-01)
- [x] Panel letters A-F added to composite — wrap_elements for patchwork composites (pA, pB, pF), labs(tag=) for ggplots (pC, pD, pE) (2026-03-01)
- [x] Panel F becomes very small in composite — resolved by 4-over-2 layout (F gets 77% of row 2 width) (2026-03-01)
- [x] Composite was 380 x 350 mm — restructured to 420 x 300 mm landscape with 4-over-2 layout (2026-03-01)
- [x] Volcano rings A/B too small and padded — tightened coord padding, margins, and allocation (2026-03-01)
- [x] Panel D RRHO wasted space with off-panel labels — moved direction cues to axis titles, removed clip="off" (2026-03-01)
- [x] Panel C subtitle too long (2-line with rho and all CIs) — simplified to single line (2026-03-01)
- [x] Unicode rho in Panel C subtitle caused PNG warnings — removed rho from subtitle (detail in concordance_stats.csv) (2026-03-01)

## Changes made

### 2026-03-01 — Loop 1: Panel C subtitle + key strip truncation
**Files modified:** `panel_C.R` (line 236), `YvO_figure2_composite.R` (lines 31/33/96)
**What changed:**
- Panel C subtitle split into 2 lines: line 1 = context (contrast, n), line 2 = statistics (r, rho, concordance with CIs)
- Composite key strip: reduced box_w (0.35->0.30), item_gap (0.6->0.45), added explicit `xlim = c(-0.5, x_cursor + 0.5)` to prevent right-edge truncation
**Outputs regenerated:** all panels + Figure_2.pdf/png + F2_supplementary.xlsx
**Result:** Both subtitle and key strip now fully visible. No analytical changes.

### 2026-03-01 — Loop 2: Panel letters A-F
**Files modified:** `YvO_figure2_composite.R`
**What changed:**
- Added panel letter tags (A-F) to composite assembly
- Used `wrap_elements(full = ...)` for patchwork composites (pA, pB, pF) to reliably tag them
- Used `labs(tag = ...)` for plain ggplot objects (pC, pD, pE)
- All tags positioned top-left, bold 12pt
**Outputs regenerated:** Figure_2.pdf/png
**Result:** All six panel letters visible and consistently positioned.

### 2026-03-01 — Loop 3: Remove Panel C redundant legend in composite
**Files modified:** `YvO_figure2_composite.R`
**What changed:**
- Added `theme(legend.position = "none")` to pC in composite assembly
- Panel C's significance legend was redundant with the unified key strip at the bottom
- Standalone `panel_C_concordance.pdf` still retains its own legend (unaffected)
**Outputs regenerated:** Figure_2.pdf/png
**Result:** Panel C scatter plot now uses full vertical space. More legible gene labels and data points at composite scale. No information lost (significance encoding in unified key strip).

### 2026-03-01 — Loop 4: Composite layout restructured to 4-over-2
**Files modified:** `YvO_figure2_composite.R`
**What changed:**
- Layout changed from 3+3 grid (A|B|C / D|E|F) to 4-over-2 (A|B|C|D / E|F)
- Row 1: A(27%) | B(27%) | C(23%) | D(23%), height 52%
- Row 2: E(25%) | F(75%), height 38%
- Key strip: height 10%
- Volcano ring margins reduced (15mm->5mm) before wrapping for composite (standalone PDFs unaffected)
- Figure dimensions changed from 380x350mm (portrait) to 420x300mm (landscape)
- Panel F now gets 75% of row 2 width (~315mm) vs previous 40% of row 2 (~152mm) -- effective area +53%
- No panel scripts modified; all changes in composite assembly only
**Outputs regenerated:** Figure_2.pdf/png + F2_supplementary.xlsx
**Result:** Volcano rings (A, B) visually larger due to margin recovery. Panel F triptych (dumbbell + Sankey + enrichment) now readable. RRHO (D) moved to top row alongside concordance scatter (C). Panel E compact in bottom-left (appropriate for its simpler content). All panel letters correctly positioned. Key strip untruncated.

### 2026-03-01 — Loop 5: Denser layout — larger volcanoes, tighter RRHO, wider concordance
**Files modified:**
- `../shared/volcano_ring.R` — tighter internal layout
- `panel_D.R` — removed off-panel labels, tighter coord
- `panel_C.R` — simplified subtitle
- `YvO_figure2_composite.R` — new width/height allocations, tighter A/B margins, suppress RRHO legend

**Volcano ring changes (`shared/volcano_ring.R`):**
- Coord padding reduced: `label_r + 3.0` -> `label_r + 1.8` (ring fills ~14% more of canvas)
- Plot margin reduced: `margin(15, 15, 5, 15, "mm")` -> `margin(3, 3, 2, 3, "mm")`
- Docstring corrected: was "shared GO terms from Young" but code does independent per-panel selection

**Panel D changes:**
- Removed 4 off-panel direction annotations (at negative x/y coords) that required `clip = "off"`
- Direction cues moved to axis titles: `"Training (Young) rank  <- Up | Down ->"`
- `coord_cartesian(clip = "off")` -> `coord_fixed(ratio = 1)` with `scale_x/y_continuous(expand = c(0, 0))`
- Added `plot.margin = margin(2, 2, 2, 2, "mm")` to theme
- Standalone PDF: heatmap now fills edge-to-edge

**Panel C changes:**
- Subtitle simplified from 2-line (with rho, all CIs) to single line: `"n = X | r = X [X, X] | concordance = X%"`
- Removes rho (redundant with r for this data) — still in concordance_stats.csv
- Eliminates Unicode rho character that caused PNG renderer warnings

**Composite changes:**
- Top row widths: `c(0.27, 0.27, 0.23, 0.23)` -> `c(0.30, 0.30, 0.24, 0.16)`
- Bottom row widths: `c(0.25, 0.75)` -> `c(0.23, 0.77)`
- Row heights: `c(0.52, 0.38, 0.10)` -> `c(0.57, 0.33, 0.10)`
- A/B composite margins: `margin(5, 5, 2, 5)` -> `margin(1, 1, 1, 1)`
- Panel D legend suppressed in composite (`legend.position = "none"`)

**Outputs regenerated:** all panels + Figure_2.pdf/png + F2_supplementary.xlsx

**Result:**
- A and B are materially larger (30% width + tighter coord/margins = ring fills much more of panel)
- C is wider (24%) with compact single-line subtitle
- D is denser (16% width, edge-to-edge heatmap, no wasted off-panel space)
- F preserved at 77% of bottom row
- E slightly narrower (23%) but functional
- Key strip and panel letters preserved
- Note: volcano_ring.R changes also affect F3 (beneficial — same tighter layout)

### 2026-03-01 — Loop 6: Further tightened volcanoes, taller figure, rebalanced bottom row
**Files modified:**
- `../shared/volcano_ring.R` — further reduced coord padding and margins
- `YvO_figure2_composite.R` — taller figure, row height rebalance, bottom row 1/3 + 2/3

**Volcano ring changes:**
- Coord padding: `label_r + 1.8` -> `label_r + 1.0` (ring fills ~20% more of canvas vs original 3.0)
- Plot margin: `margin(3, 3, 2, 3)` -> `margin(1, 1, 1, 1)` (minimal)
- Labels still fully visible (clip="off" allows extension beyond coord limits)

**Composite changes:**
- Figure height: 300mm -> 360mm (20% taller, more room for row 2)
- Row heights: `c(0.57, 0.33, 0.10)` -> `c(0.48, 0.42, 0.10)` (row 2 now 42% = ~151mm vs 99mm)
- Bottom row widths: `c(0.23, 0.77)` -> `c(1, 2)` (1/3 Panel E, 2/3 Panel F)
- Top row widths unchanged: `c(0.30, 0.30, 0.24, 0.16)`

**Outputs regenerated:** all panels + Figure_2.pdf/png + F2_supplementary.xlsx

**Result:**
- Volcano rings A/B substantially larger — rings dominate their panels with minimal dead space
- Panel F has more vertical room (151mm vs 99mm) — dumbbell gene names, Sankey, and enrichment all more readable
- Panel E at 1/3 bottom row (~140mm) is more spacious than before (was 23% = ~97mm)
- Standalone panel A/B PNGs confirmed: no clipping, all labels readable
- Note: volcano_ring.R changes also affect F3 (beneficial)

### 2026-03-01 — Loop 7: Restructured to 3-row × 2-column landscape layout
**Files modified:** `YvO_figure2_composite.R`
**What changed:**
- Layout changed from 4-over-2 (A|B|C|D / E|F) to 3-row × 2-column (A|C / B|D / E|F)
- Used patchwork `area()` design on a 100-unit virtual grid instead of row stacking
- Row 1 (27%): A (volcano Young) | C (concordance scatter)
- Row 2 (27%): B (volcano Old) | D (RRHO heatmap)
- Row 3 (40%): E (NES bubble, 35% width) | F (interaction triptych, 65% width)
- Footer (6%): unified key strip
- Column split: left 45% (volcanos), right 55% (C, D, F)
- Figure dimensions: 420×360mm → 440×330mm (wider landscape)
- A and B now have identical cell dimensions (27% height × 45% width) guaranteeing equal rendering
- Each volcano gets its own full row height instead of sharing a single row with two other panels
- No panel scripts modified; all changes in composite assembly only

**Outputs regenerated:** all panels + Figure_2.pdf/png + F2_supplementary.xlsx

**Result:**
- A and B stacked vertically in left column, identically sized — visual comparison much easier
- C and D paired in right column (protein-level and threshold-free concordance side by side vertically)
- Bottom row: E (35%) and F (65%) span full width; F triptych gets ~286mm
- Key strip thin footer, not a structural row
- All panel letters visible and correctly positioned
- Volcano rings substantially larger per panel (each gets ~198×89mm cell vs previous ~126×173mm shared row)

### 2026-03-01 — Loop 8: Panel F to right column, three square concordance panels bottom
**Files modified:** `YvO_figure2_composite.R`
**What changed:**
- Layout changed from A|C / B|D / E|F to: A|F(span) / B|F(span) / C|E|D
- Panel F moved to right column spanning rows 1-2 (55% width × 54% height = ~242×205mm)
- Bottom row: C (protein-level) | E (pathway-level) | D (RRHO) — three equal ~square panels
- Bottom panels each 33% width × 38% height → ~145×144mm at 440×380mm ≈ square
- Figure dimensions: 440×330mm → 440×380mm (taller to accommodate square bottom panels)
- Panel addition order changed to match area assignments: pA + pB + pF + pC + pE + pD + pKey
- No panel scripts modified; all changes in composite assembly only

**Outputs regenerated:** all panels + Figure_2.pdf/png + F2_supplementary.xlsx

**Result:**
- F triptych has full right-column height (~205mm) for gene names and Sankey
- C, E, D bottom row nearly square — clean visual hierarchy
- Volcanos A/B stacked left, identically sized
- Bottom row groups all three concordance measures (protein, pathway, threshold-free) together
- Key strip and panel letters preserved

## Open questions
- Should volcano ring pathway selections be harmonized across A and B?
- Should RRHO colorbar be added to the unified key strip?

## Publication-readiness checklist
- [x] Panel C 1:1 concordance line at true 45° (symmetric axes + coord_fixed)
- [x] Key strip complete and untruncated (Direction + Significance + Quadrant)
- [x] Volcano ring labels legible at journal scale
- [x] Panel letters A-D present and consistent (14pt bold)
- [x] Colors follow shared/palettes.R (DIR_COLORS, SIG_COLORS consistent)
- [x] Composite dimensions appropriate (340×300mm landscape, 4 panels)
- [x] Concordance framing honest (r = 0.26 with CI, no overclaiming)
- [x] A and B identically sized (same grid cell dimensions)
- [x] Panel D RRHO fills space (title stripped in composite)
- [x] ORA bars readable in composite (pathway names, -log10 padj values)
- [x] Supplementary Excel updated (E/F worksheets removed, stats sheet added)

### 2026-03-02 — Refinement pass: Panels C, D, F visualization upgrades

**Files modified:**
- `a_script/YvO_F2_setup.R` — added `classify_expanded()`, `EXPANDED_COLORS`, `EXPANDED_ORDER`, `ORA_CONCORDANT_COL`, `ORA_DISCORDANT_COL`, `library(ggnewscale)`
- `a_script/panel_C.R` — replaced single geom_point layer with hexbin (NS) + ggnewscale::new_scale_fill() + geom_point (significant categories)
- `a_script/panel_D.R` — added RRHO hotspot gene extraction, per-quadrant ORA (Hallmark + GO:BP), flanking bar plots, patchwork assembly (RRHO | concordant/discordant bars)
- `a_script/panel_F.R` — complete redesign: expanded 7-category classification, direction-split heatmap with category annotation bar, category-to-pathway aggregate Sankey, stacked composition bars (F4 pattern)
- `a_script/YvO_figure2_composite.R` — Panel D wrapped with wrap_elements (now patchwork), added D_rrho2_ora and F_expanded worksheets to supplementary Excel, updated metadata

**Panel C changes:**
- NS proteins (~1900) rendered as hexbin density background (bins=40, alpha=0.7, grey95→grey50 gradient)
- All 4 significance categories (Interaction, Sig Both, Sig Young only, Sig Old only) rendered as explicit colored points on top of hexbin
- Uses ggnewscale to manage dual fill scales (hexbin gradient + categorical)
- All other elements preserved: quadrant shading, reference lines, gene labels, statistics, CSV exports

**Panel D changes:**
- RRHO heatmap preserved unchanged
- Added hotspot gene extraction: for each quadrant, finds peak cell (max |hmat|), extracts overlap gene set
- Per-quadrant ORA via enricher() on each hotspot gene set (Hallmark + GO:BP, padj < 0.05)
- Flanking bar plots: concordant (warm #FFB74D) and discordant (cool #64B5F6) faceted by specific quadrant
- Patchwork assembly: RRHO | (concordant / discordant bars), widths c(3,1)
- New CSV exports: rrho2_hotspot_genes.csv, rrho2_ora_concordant.csv, rrho2_ora_discordant.csv
- Edge case: empty panel with "No enrichment" note if quadrant has no significant ORA results

**Panel F changes:**
- Expanded from 4 significance categories to 7 direction-split categories via classify_expanded()
- Categories: Interaction, Sig Both Up/Down, Sig Young Up/Down, Sig Old Up/Down
- Heatmap: gene rows grouped by expanded category with small y-gaps between groups, narrow category annotation bar (left), gene×(Young,Old) logFC tiles
- Replaced gene-level Sankey with category-to-pathway aggregate Sankey (ribbons proportional to gene count per category×pathway pair)
- Added stacked composition bars (F4 cumulative slot tracking pattern): one bar per pathway, segments colored by contributing category
- Soft gene cap at 150 (top per category by |logFC|); compress row height if needed
- 2 sub-panels: heatmap (30%) | Sankey+stacked bars (70%)
- New CSV exports: expanded_classification.csv, expanded_sankey_links.csv, stacked_bar_data.csv
- Legacy sankey_links.csv still exported for backward compat

**New dependencies:**
- ggnewscale (Panel C dual fill scales; already available in env from F5)

**Outputs to regenerate:** All panels + Figure_2.pdf/png + F2_supplementary.xlsx

**Remaining issues:**
- Panel E left unchanged in this pass (may need light readability refinements after viewing composite)
- Stacked bar grid interval set to 10 — may need adjustment based on actual gene counts
- Panel D ORA bar pathway names truncated to 30 chars — verify legibility at composite scale

### 2026-03-02 — Pass 2: Layout & analytical refinements

**Files modified:**
- `a_script/panel_D.R` — restructured patchwork from side-by-side to stacked layout
- `a_script/panel_E.R` — added coord_fixed for square aspect ratio
- `a_script/panel_F.R` — removed pathway cap, show all padj < 0.05 pathways
- `a_script/YvO_figure2_composite.R` — rebalanced grid layout, D enlarged to 50% bottom width

**Panel D changes:**
- Layout changed from `pD | (pD_conc / pD_disc)` (side-by-side, widths 3:1) to `pD / (pD_conc | pD_disc)` (stacked, heights 4:1)
- RRHO heatmap now occupies ~80% of panel height (naturally square via coord_fixed)
- ORA bars span full width below RRHO in two columns (concordant left, discordant right)
- UD quadrant confirmed as genuine biology: 0 enrichment at padj < 0.05. Added explicit "No enrichment" note with hotspot gene count when UD is empty
- ORA bar pathway label truncation increased from 30 to 40 chars (wider bars in new layout)
- Custom empty labels with gene counts for concordant/discordant groups
- Standalone PDF dimensions changed from 260×180mm (landscape) to 220×260mm (portrait)

**Panel E changes:**
- Replaced `coord_cartesian(xlim, ylim)` with `coord_fixed(ratio = 1, xlim = c(-nes_lim, nes_lim), ylim = c(-nes_lim, nes_lim))`
- Plot now renders as a square with symmetric axes centered on origin
- Quadrant corner labels updated to use `nes_lim` (data-driven) instead of hardcoded positions
- Repel label limits updated to `±nes_lim * 0.9` for consistent margin

**Panel F changes:**
- Removed `max_total_pws = 15` hard cap on pathway count
- Changed pathway selection from top-5 Hallmark + top-5 GO:BP to ALL pathways with `p.adjust < 0.05`
- Greedy rescue loop increased from 15 to 50 iterations (no cap check)
- Added adaptive text sizing: `pw_label_size` scales down for >15 or >20 pathways
- Adaptive `S_SBAR_H` prevents stacked bar overlap with many pathways
- "Other" category expected to shrink substantially or disappear

**Composite changes:**
- Panel D allocation enlarged: 50% of bottom row width (was 34%)
- Panel E allocation: 25% of bottom row width (was 33%)
- Panel C allocation: 25% of bottom row width (was 33%)
- Figure height increased from 400mm to 480mm
- RRHO at 80% of D's height: ~184mm square (was ~130mm)
- ORA bars at 20% of D's height: ~230mm × 46mm (full pathway names readable)

**Outputs to regenerate:** All panels + Figure_2.pdf/png + F2_supplementary.xlsx

## Next recommended steps
1. ~~Fix Panel C subtitle truncation~~ DONE
2. ~~Fix key strip truncation~~ DONE
3. ~~Add panel letters~~ DONE
4. ~~Make A/B larger, C wider, D denser~~ DONE
5. ~~Further tighten volcanoes, rebalance bottom row~~ DONE
6. ~~Restructure to 3-row × 2-column layout~~ DONE
7. ~~Panel F to right column, square concordance panels bottom~~ DONE
8. ~~Panel C hexbin + category overlays~~ DONE
9. ~~Panel D RRHO + ORA flanks~~ DONE
10. ~~Panel F expanded 7-category heatmap + stacked bars~~ DONE
11. ~~Panel D stacked layout (RRHO over ORA bars)~~ DONE
12. ~~Panel E square aspect ratio~~ DONE
13. ~~Panel F pathway cap removed~~ DONE
14. ~~Composite D enlarged to 50% width~~ DONE
15. ~~Redesign: drop Panels E/F, clean 2×2 layout (A|B / C|D)~~ DONE
16. ~~Panel C: fix asymmetric axes + add coord_fixed~~ DONE
17. ~~Panel D: strip title/subtitle in composite~~ DONE
18. Consider Panel A/B pathway harmonization (deliberate independent selection)
19. ~~Panel CE prototype: concordance scatter + pathway theme hulls~~ REMOVED (user decision)
21. ~~Panel F: replace ORA with GO Slim super-category mapping~~ DONE

### 2026-03-02 — Figure redesign: 4-panel (A|B / C|D) publication-ready layout

**Rationale:**
Figure 2 was redesigned from a 6-panel (A-F) to a focused 4-panel layout. Panels E (NES bubble) and F (expanded heatmap/Sankey) were removed to sharpen the figure's narrative: individual contrasts (A, B) followed by protein-level (C) and threshold-free (D) concordance measures.

**Scientific validation (pre-implementation):**
- Panels A/B: Pi-score threshold documented, independent pathway selection justified, NES fill clamped — no issues
- Panel C: **Fixed** — asymmetric axes (x: [-2.5, 2] vs y: [-1, 2]) without coord_fixed was distorting the 1:1 concordance line. Now uses symmetric data-driven range and coord_fixed(ratio=1)
- Panel D: Standard RRHO, no MTC across grid cells (standard practice, documented), ORA hotspot extraction sound

**Files modified:**
- `a_script/panel_C.R` — replaced hardcoded asymmetric axis ranges with symmetric data-driven `axis_max`; replaced `coord_cartesian` with `coord_fixed(ratio = 1)`
- `a_script/YvO_figure2_composite.R` — complete rewrite: removed Panel E/F sources; simplified to 2×2 grid with patchwork `area()` design; updated key strip (Direction + Significance + Quadrant); stripped Panel D inner titles via `& theme(plot.title = element_blank())`; reduced dimensions from 460×480mm to 340×300mm; simplified supplementary Excel (removed E/F worksheets, added C_concordance_stats)

**Layout (340×300mm):**
- Row 1 (42%): A (volcano Young, 50%) | B (volcano Old, 50%)
- Row 2 (51%): C (concordance, 42%) | D (RRHO + ORA, 58%)
- Footer (6%): Key strip (Direction, Significance, Quadrant)

**Key strip simplified:**
- Direction: Up/Down (volcano rings)
- Significance: Interaction, Sig Both, Sig Young only, Sig Old only (Panel C)
- Quadrant: Concordant/Discordant shading (Panels C and D)

**Panel scripts NOT modified:** panel_AB.R, panel_D.R (standalone PDFs unaffected)
**Panel scripts modified:** panel_C.R (axis symmetry only; standalone PDF benefits)

**Outputs regenerated:** panel_C_concordance.pdf/png, Figure_2.pdf/png, F2_supplementary.xlsx

**Verification:**
- Composite: all 4 panel letters visible (A-D), key strip untruncated
- Panel C: 1:1 line at true 45°, symmetric axes ~±2.2, statistics unchanged (r=0.26 [0.22, 0.30])
- Panel D: RRHO fills panel (title stripped), ORA bars readable, colorbar legend visible
- Panel A/B: unchanged from prior session
- Standalone PDFs: Panel C now has symmetric axes (improvement); all others unchanged
- Unicode warnings in PNG (U+2002 en-space, U+2190/2192 arrows) — cosmetic, PDF unaffected

### 2026-03-02 — Panel CE: Concordance scatter with pathway theme hulls (prototype)

**Rationale:**
Consolidates the information from Panel C (protein-level concordance) and Panel E (pathway-level concordance) into a single panel using soft convex hull annotations. Pathway-level biological themes are overlaid as subtle hull regions on the existing concordance scatter, showing which functional groups the significant proteins belong to without disrupting the core interpretation.

**Design document:** `docs/plans/2026-03-02-panel-CE-concordance-hulls-design.md`

**New file created:** `a_script/panel_CE.R`
- Sources `panel_C.R` to reuse all computed objects (scatter_df, ns_df, sig_df, label_df, statistics)
- Uses `msigdbr` (Hallmark + GO:BP) to map significant proteins to 6 biological themes via `THEME_MAP`
- Each gene assigned to its primary theme (most pathway memberships; ties broken alphabetically)
- Themes with <3 proteins excluded from hull rendering
- `ggforce::geom_mark_hull` for concave hulls; `ggnewscale` for 3 fill scales (hexbin, hulls, significance)
- Toggle: `SHOW_THEME_HULLS <- TRUE` (set FALSE for standard Panel C)

**Two-layer hull approach (key design decision):**
- Layer 1 (before points): per-theme hull fills using a loop with fixed `fill = HULL_COLORS[tn]` — no scale needed, very subtle (alpha=0.10), dashed grey border
- Layer 2 (after gene labels): invisible hulls (`alpha = 0, color = NA`) with `aes(label = theme)` for ggforce's built-in boundary label placement with elbow connectors
- This separation ensures hull fills sit behind points/labels while hull labels are auto-positioned at hull boundaries without overlapping gene labels

**Theme assignment results (75 significant proteins, 4 valid themes):**
| Theme | Proteins | Example genes |
|-------|----------|---------------|
| Proteostasis | 5 | HSPB6, HPRT1, CCT6A, USP9X, TES |
| Metabolism | 21 | CARNS1, NDUFS7, NDUFA3, ALDH1B1, GPT2 |
| Cytoskeleton | 41 | ACTN4, MYL2, ACTN2, COL3A1, FLNC |
| ECM/Remodeling | 8 | COL1A1, COL4A1, LAMC1, FBN1, PREX1 |

Calcium/Signaling and Stress/Inflammation excluded (<3 proteins each).

**Outputs generated:**
- `b_reports/panel_CE_concordance.pdf` (220×220mm)
- `b_reports/panel_CE_concordance.png` (220×220mm, 300 dpi)
- `c_data/panel_CE/theme_assignments.csv` (75 rows: gene, theme, significance, logFC values)
- `c_data/panel_CE/concordance_themed.csv` (full scatter data with optional theme column)

**Hull visual parameters:**
- Colors: Proteostasis=#FFD54F, Metabolism=#A5D6A7, Cytoskeleton=#FFAB91, ECM/Remodeling=#90CAF9
- Fill alpha: 0.10 (very subtle); border: dashed grey40 at alpha 0.4, linewidth 0.4
- Concavity: 2 (moderately tight); expand/radius: 2mm each
- Labels: 8pt bold italic grey25 on white background, elbow connectors

**Iterations during development:**
1. Centroid-based labels → all near origin, overlapping → rejected
2. geom_label_repel at centroids → still clustered → rejected
3. Extreme-point anchoring → GPT2 overlap with gene label → rejected
4. Two-layer geom_mark_hull → ggforce auto-places at hull boundary → **adopted**

**Integration status:** Standalone prototype only. Not yet integrated into composite (Figure_2.pdf still uses original 4-panel A|B / C|D layout). User evaluation needed before deciding whether to replace Panel C or keep as supplementary.

**No existing files modified.** Panel C, composite, and all other scripts unchanged.

### 2026-03-02 — Panel F: Replace ORA with GO Slim super-category mapping

**Rationale:**
Panel F's ORA (Hallmark + GO:BP via `enricher()`) assigned 131 significant proteins to pathways, but 33% (43/131) fell into "Other" because metabolic enzymes, proteasome subunits, ribosomal proteins, etc. received no meaningful label. GO:BP term redundancy also wasted pathway slots on overlapping terms (6+ muscle/cardiac variants).

**Replacement approach:** GO Slim Generic (62 curated, non-redundant GO BP terms) mapped to 13 biologically coherent super-categories via GOBPANCESTOR ancestor traversal.

**Files modified:**
- `a_script/panel_F.R` — Section 2 (lines 69–288) completely replaced:
  - Removed: `msigdbr()` + `enricher()` ORA pipeline, greedy rescue loop, pathway redistribution, pathway cap
  - Added: GO:BP annotation lookup for all universe genes, ancestor-based GO→Slim mapping, specificity-weighted 1:1 gene→super-category assignment, Fisher's exact test per category, singleton merge
  - Stat audit comment updated (Section 0)
  - Subtitle updated from "ORA: Hallmark + GO:BP" to "GO Slim super-categories (Fisher, BH adj.)"
  - Enrichment bar fallback: when all Fisher p.adjust=1 (expected for broad categories), bars scale by gene count instead of -log10(padj)
- `a_script/test_goslim.R` — deleted (prototype no longer needed)

**New dependency:** `GO.db` (for `GOBPANCESTOR` ancestor graph; already installed)

**Super-category results (131 genes, 13 categories):**
| Category | n genes |
|----------|---------|
| Metabolism | 21 |
| Gene Expression | 15 |
| Protein Homeostasis | 13 |
| Signaling | 12 |
| Transport | 12 |
| DNA & Cell Cycle | 10 |
| Cytoskeleton & Motility | 9 |
| ECM & Adhesion | 7 |
| Immune & Inflammation | 7 |
| Mitochondria & Energy | 7 |
| Circulatory System | 5 |
| Muscle & Contractile | 3 |
| Other | 10 |

- "Other" reduced from 43 → 10 genes (7.6% vs 33%)
- "Development" category had only 1 gene → merged into "Other" automatically

**Interface contract preserved:**
- `ora_top`: tibble with `pathway_label`, `pvalue`, `p.adjust`, `ID`, `database`
- `links_1to1`: tibble with `gene`, `pathway`, `pvalue`, `expanded_cat`, `logFC_Training_Young`, `logFC_Training_Old`, `database`
- Sections 3–8 (ordering, plot build, CSV export, save, verification) unchanged

**Outputs regenerated:**
- `b_reports/panel_F_interaction.pdf` / `.png` — 12 super-category pathway labels, gene-count-scaled bars, Sankey ribbons with category colors
- `c_data/panel_F/expanded_sankey_links.csv` — super-category names in pathway column
- `c_data/panel_F/expanded_classification.csv`, `stacked_bar_data.csv`, `sankey_links.csv` — all updated
- `b_reports/Figure_2.pdf` / `.png` — composite renders cleanly (Panel F is standalone, not in 4-panel composite)
- `c_data/F2_supplementary.xlsx` — regenerated

**Verification:**
- Script completes without error
- 131 genes mapped, Sankey source total matches (121 non-Other = 121 OK)
- 43 Sankey ribbons built, 43 enrichment bar segments
- Enrichment bars use gene count scaling (axis label "Gene count") — appropriate for descriptive super-categories
- Composite Figure_2.pdf/png unaffected (Panel F standalone only)

### 2026-03-03 — F2/F3 Panel Harmonization & Standardization

**Objective:** Bring F2 and F3 to consistent publication standard with shared infrastructure, unified key panels, and harmonized themes/typography across all panels.

**Shared infrastructure changes:**
- `shared/palettes.R` — added THEME_FIG (base_size=14) and TXT_* typography constants (TXT_PATHWAY, TXT_GENE, TXT_QUADRANT, TXT_STATS, TXT_TAG, TXT_ORA_BAR, TXT_ORA_AXIS, TXT_ORA_STRIP)
- `shared/go_slim_categories.R` — **new file**: extracted GO Slim super-category system from panel_F.R into shared utility (62 GO Slim IDs, 13 super-categories, GOBPANCESTOR hierarchy, specificity-weighted 1:1 gene assignment)

**Setup changes (`YvO_F2_setup.R`):**
- Removed local THEME_FIG and TXT_* definitions (now from shared/palettes.R)
- Added `source("04_Figures/shared/go_slim_categories.R")`
- Kept figure-specific colors (SIG_COLORS, CONTRAST_COLORS, ORA_QUAD_COLORS, EXPANDED_COLORS)

**Panel-level changes:**

| Panel | Changes |
|-------|---------|
| panel_C.R | Removed legend (`legend.position = "none"`), removed imputed-protein border encoding (all borders now grey75), switched to `cairo_pdf` |
| panel_D.R | Switched THEME_PUB → THEME_FIG, removed heatmap legend, switched to horizontal layout `(pD \| pD_ora)` with widths c(3,2), ORA facet_wrap ncol=2→1, hardcoded sizes → `rel(0.7)`, switched to `cairo_pdf`, added PNG output |
| panel_E.R | Switched THEME_PUB → THEME_FIG, label size → TXT_GENE, quadrant size → TXT_QUADRANT, removed hand-built legend (pE_key/pE_combined), switched to `cairo_pdf` |
| panel_F.R | Replaced ~220 lines of inline GO Slim code with shared `assign_go_slim_super()`, switched theme_minimal → THEME_FIG, removed legend (`legend.position = "none"`), updated title/subtitle, switched to `cairo_pdf` |

**Composite changes (`YvO_figure2_composite.R`):**
- Complete rewrite with new 100×100 virtual grid layout
- Dimensions: 350 × 480 mm (was 440 × 380 mm)
- Added unified key panel at bottom (rows 86-100) with annotate()-based sections: Significance (5-level), Scatter Quadrants, RRHO Overlap, ORA Quadrants, Heatmap logFC, GO Slim Super-categories (13)
- No title stripping — panels keep their own titles/subtitles
- Panel letter tags via `labs(tag=)` at TXT_TAG size (20pt bold)
- Layout: A|D / B|D / C|E / F(span) / Key — 3 volcanos stacked left, scatter and RRHO+ORA mid, heatmap-sankey right, key full-width bottom
- Supplementary Excel workbook preserved

**Outputs to regenerate:** All F2 panels + Figure_2.pdf/png + F2_supplementary.xlsx

---

### Session: 2026-03-03 — Typography Unification

**Goal:** Unify all annotation text sizes across F2/F3 panels for consistent readability.

**Shared infrastructure (`shared/palettes.R`):**
- Replaced 8 separate TXT_* constants (3.5–4.5 range) with 2 unified constants: `TXT_TAG = 18`, `TXT_LABEL = 3.5`
- Legacy aliases (TXT_GENE, TXT_PATHWAY, TXT_QUADRANT, TXT_STATS, TXT_ORA_BAR, TXT_ORA_AXIS, TXT_ORA_STRIP) all point to TXT_LABEL for backward compatibility
- Key changes: TXT_TAG 20→18, TXT_GENE 4.5→3.5, TXT_STATS 4.5→3.5, TXT_ORA_BAR 4.5→3.5, TXT_ORA_AXIS/STRIP 4.0→3.5

**F2 panel_D.R:**
- Colorbar title size: 9→11 (harmonized with F3)

**F2 panel_F.R:**
- Removed custom theme overrides that fought THEME_FIG (plot.title=11→inherited 16, plot.subtitle=8 italic→inherited 11 bold.italic, legend.text=8→inherited 12)
- Replaced all hardcoded annotation sizes with TXT_LABEL: pw_label_size (4–5.5→3.5), count labels (4→3.5), column headers (5→3.5), axis/tick labels (4.5→3.5)
- Colorbar guide title 10→11, label 8→9

**Automatic via aliases (no code changes needed):**
- panel_C.R gene labels (TXT_GENE 4.5→3.5), quadrant labels (TXT_QUADRANT stays 3.5)
- panel_D.R quadrant labels (TXT_QUADRANT stays 3.5), ORA bar text (TXT_ORA_BAR 4.5→3.5)
- panel_E.R pathway labels (TXT_GENE 4.5→3.5), quadrant labels (TXT_QUADRANT stays 3.5)

**Verified:** F2 composite renders successfully with all changes.

### Session: 2026-03-03 — Panel F Text Size Compensation & Layout Widening

**Goal:** Fix Panel F text being too small in composite (root cause: wrap_elements() shrinks 400–500mm standalone to ~220mm, making 3.5mm text → ~2mm). Also widen composite to give Panel F more space.

**Shared infrastructure (`shared/palettes.R`):**
- `TXT_PF`: 5.0 → 7.0 (compensates for ~0.56× shrinkage → ~3.9mm effective in composite)
- `TXT_PF_GENE`: 14 → 20 (gene symbols shrink to ~11pt in composite)

**F2 panel_F.R:**
- All `TXT_LABEL` (3.5mm) → `TXT_PF` (7.0mm) for: pathway labels, count labels, column headers, axis labels, key text
- gene_text_pt doubled: 7/8.5/10 → 14/16/20 (for >100/>60/≤60 genes)
- Key position: moved from right-center (`S_X_BAR + S_MAX_LEN * 0.45`) to upper-left (`S_X_BAR + 0.05`)
- Standalone width: 500 → 400mm

**F2 composite (`YvO_figure2_composite.R`):**
- Panel F allocation: columns 56-100 → 48-100 (53% width, was 45%)
- Other panels adjusted: A,B,C cols 1-23, D,E cols 24-47
- Composite width: 350 → 420mm (height unchanged at 480mm)

**Outputs regenerated:** All F2 panels + Figure_2.pdf/png (420×480mm) + F2_supplementary.xlsx

**Verified:**
- Composite renders without errors
- Panel F gets 53% of 420mm = ~223mm allocation
- Panel F text scales to ~3.9mm effective in composite (close to other panels' TXT_LABEL = 3.5mm)
- Gene symbols scale to ~8pt effective for 131 genes (readable)
- Key relocated to upper-left of enrichment bar area
- Unicode warnings (U+2002, U+03C1, U+2190/2192) are pre-existing cosmetic issues in Panels C and D

### Session: 2026-03-03 — Panel F Gene Count Bars + Panel G Per-Contrast ORA

**Goal:** (1) Switch Panel F enrichment bars from fold enrichment to gene counts. (2) Rename "super-categories" to "consolidated pathways". (3) Add new Panel G: per-contrast ORA dot plot using Hallmark + GO:BP/MF/CC.

**Shared infrastructure (`shared/go_slim_categories.R`):**
- Renamed: `SUPER_CATEGORY_ORDER` → `CONSOLIDATED_PATHWAY_ORDER`, `SUPER_COLORS` → `CONSOLIDATED_COLORS`, `SLIM_SUPER` → `SLIM_CONSOLIDATED`, `assign_go_slim_super()` → `assign_go_slim_consolidated()`
- Return column: `super` → `consolidated`
- Legacy aliases preserved for backward compatibility: `SLIM_SUPER`, `SUPER_COLORS`, `SUPER_CATEGORY_ORDER`, `assign_go_slim_super()` (returns `super` column via rename)

**F2 panel_F.R:**
- Bars now show gene count (not fold enrichment): `bar_total_w <- (pw_row$gene_count / max_count) * S_MAX_LEN`
- Axis label: "Gene count" with integer ticks
- Removed fold=1.0 reference line
- S_MAX_LEN: 2.8 → 3.4 (wider bars)
- Standalone width: 280 → 320mm
- Terminology: "super-category" → "consolidated pathway" in comments, messages, subtitle
- Fold enrichment still computed and exported to CSV for reference

**F2 panel_G.R (new):**
- 4 direction-collapsed gene lists: Interaction, Sig Both, Sig Young, Sig Old
- 4 databases: Hallmark (MSigDB H), GO:BP, GO:MF, GO:CC via clusterProfiler
- Top 5 terms per contrast per database, BH correction (padj < 0.05, qvalue < 0.2)
- Dot plot: x = -log10(padj), y = term, color = contrast, size = gene count
- Faceted by database
- Exports: `c_data/panel_G/ora_per_contrast.csv`, `b_reports/panel_G_ora.pdf`
- F2 results: 30 significant terms across 8 contrast-database pairs

**F2 composite (`YvO_figure2_composite.R`):**
- Sources panel_G.R → pG, wrapped as Panel G with tag "G"
- Layout: added row 73-90 for Panel G, key moved to 91-100
- Composite dimensions: 420 × 560mm (was 420 × 480mm)
- Key section: "GO Slim Super-categories" → "Consolidated Pathways" (uses CONSOLIDATED_* variables)
- Supplementary workbook: added `G_ora_per_contrast` sheet

**Outputs regenerated:** All F2 panels + Panel G + Figure_2.pdf/png + F2_supplementary.xlsx

**Verified:**
- Panel F: 131 proteins, 12 consolidated pathways, gene count bars scale correctly (largest = Metabolism n=22), integer axis ticks, sankey totals match (119 OK)
- Panel G: 30 ORA terms across 8 contrast-database pairs, dot plot rendered at 280 × 300mm
- Composite renders with all 7 panels (A–G) + key, no errors
