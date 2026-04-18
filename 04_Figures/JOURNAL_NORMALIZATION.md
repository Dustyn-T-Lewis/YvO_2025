# Journal of Physiology — Figure Normalization Plan

## Target Venue
**The Journal of Physiology** (Wiley) — ISSN 1469-7793

## Figure Specifications

| Spec | Value |
|---|---|
| Single column width | 85 mm (3.35 in) |
| Double column width | 178 mm (7.01 in) |
| Max page height | ~240 mm (with caption space) |
| Resolution | 300 dpi (submission), 600 dpi (production) |
| Font family | Sans-serif (Arial or Helvetica) |
| Font size | 10–12 pt at reproduction size |
| Panel labels | **Lowercase** (a, b, c), 12 pt sans-serif |
| Color mode | RGB |
| File format | PDF (vector), TIFF |

## Width Assignments

| Figure | Type | Target Width | Current Width | Ratio |
|---|---|---|---|---|
| F01 MAIN | **single-col** (85 mm) | 85 mm | 215 mm | 2.5× |
| F01 SUPP | **single-col** (85 mm) | 85 mm | 170 mm | 2.0× |
| F02 MAIN | double-col (178 mm) | 178 mm | 440 mm | 2.5× |
| F02 SUPP | double-col (178 mm) | 178 mm | 300 mm | 1.7× |
| F03 MAIN | **single-col** (85 mm) | 85 mm | 380 mm | 4.5× |
| F03 SUPP | double-col (178 mm) | 178 mm | 300 mm | 1.7× |
| F04 MAIN | double-col (178 mm) | 178 mm | 720 mm | 4.0× |
| F04 SUPP | double-col (178 mm) | 178 mm | 460 mm | 2.6× |
| F05 MAIN | double-col (178 mm) | 178 mm | 720 mm | 4.0× |
| F05 SUPP | double-col (178 mm) | 178 mm | 460 mm | 2.6× |
| F06 MAIN | double-col (178 mm) | 178 mm | 520 mm | 2.9× |
| F06 SUPP | double-col (178 mm) | 178 mm | varies | — |
| F07 MAIN | double-col (178 mm) | 178 mm | 450 mm | 2.5× |
| F07 SUPP | double-col (178 mm) | 178 mm | varies | — |

## Changes Required (per figure)

### Global (shared/style.R)
- [ ] Switch `FIG_THEME` from `theme_bw()` serif default to sans-serif (Arial/Helvetica) via `base_family = "Helvetica"`
- [ ] Recalibrate `FIG_TITLE_SIZE`, `FIG_SUBTITLE_SIZE`, `FIG_AXIS_TEXT`, etc. for 178 mm canvas (current sizes tuned for 300–720 mm)
- [ ] Update `composite_text_sizes()` scaling function for journal dimensions
- [ ] Panel tag convention: uppercase → **lowercase** (a, b, c)

### Per-figure stitcher changes
For each `90_stitch_main.R`:
- [ ] Set `COMP_W` to 178 mm (or 85 mm for F01/F03)
- [ ] Recalculate `COMP_H` to preserve aspect ratio or fit page height (~240 mm max)
- [ ] Update cowplot `draw_label()` tag positions (x/y normalized coords change when canvas shrinks)
- [ ] Update tag text: `"A"` → `"a"`, `"B"` → `"b"`, etc.
- [ ] Recalibrate `TAG_SZ` for journal size
- [ ] Update panel layout `widths`/`heights` ratios if needed
- [ ] Verify `plot.margin` values don't clip at smaller size

### Per-panel script changes
For each `panel_*.R`:
- [ ] Verify axis text ≥ 8 pt at reproduction size
- [ ] Verify title/subtitle readable at target width
- [ ] Check `geom_text`/`annotate` font sizes scale appropriately
- [ ] Check legend key sizes, spacing
- [ ] Verify `plot.margin` doesn't waste space at smaller canvas

### Single-column figures (F01, F03)
- F01: 3 panels stacked vertically (already has `90_stitch_single_col.R` at 86 mm — close to 85 mm target)
- F03: 4 volcano rings in 2×2 — need aggressive text scaling for 85 mm width

## Approach

Work figure by figure:
1. **F01** (simplest, 3 panels, single-col) — prototype the style changes
2. **F02** (6 panels, double-col) — test the 178 mm width
3. **F03** (4 panels, single-col) — challenging text at 85 mm
4. **F04/F05** (5 panels each, mirror) — complex layouts
5. **F06** (3 main + nested supp) — heatmap-heavy
6. **F07** (2 main + 6 supp) — ROC/scatter grid

For each figure:
1. Update `COMP_W` and `COMP_H`
2. Re-render and visually inspect
3. Adjust text sizes, margins, spacing
4. Check panel labels are lowercase
5. Verify at 300 dpi the output is clean

## Current Text Size Baseline (shared/style.R)

```r
FIG_TITLE_SIZE    <- 12    # panel titles
FIG_SUBTITLE_SIZE <- 9     # panel subtitles
FIG_STRIP_SIZE    <- 10    # facet strips
FIG_AXIS_TEXT     <- 8.5   # axis tick labels
FIG_LEGEND_TITLE  <- 9.5   # legend title
FIG_LEGEND_TEXT   <- 8.5   # legend labels
base_size         <- 10    # theme_bw base
```

These sizes were tuned for 300–720 mm canvases. At 178 mm, they should mostly work but titles/subtitles may need to shrink. At 85 mm, everything needs halving.

## Files to Modify

### Global
- `shared/style.R` — FIG_THEME font family, text sizes, composite_text_sizes()

### Per figure (90_stitch_main.R in each)
- `F01/a_script/main/90_stitch_main.R`
- `F02/a_script/main/90_stitch_main.R`
- `F03/a_script/main/90_stitch_main.R`
- `F04/a_script/main/90_stitch_main.R`
- `F05/a_script/main/90_stitch_main.R`
- `F06/a_script/main/90_stitch_main.R`
- `F07/a_script/main/90_stitch_main.R`

### Per figure (supp stitchers)
- `F01/a_script/supp/90_stitch_supp.R`
- `F02/a_script/supp/90_stitch_supp.R`
- `F03/a_script/supp/90_stitch_supp.R`
- `F04/a_script/supp/90_stitch_supp.R`
- `F05/a_script/supp/90_stitch_supp.R`
- `F06/a_script/supp/` (nested: 01_QC, 02_analysis, 03_module, mega)
- `F07/a_script/supp/90_stitch_supp.R`

### Single-col variants
- `F01/a_script/main/90_stitch_single_col.R` (already ~86 mm — adjust to 85)
- `F03/a_script/main/90_stitch_single_col.R`
