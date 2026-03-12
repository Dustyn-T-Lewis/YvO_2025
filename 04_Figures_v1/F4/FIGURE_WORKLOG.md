# Figure 4 — Refinement Worklog

## Figure purpose
**Proteomic Response Archetypes to Resistance Training.**
Data-driven FCM clustering on per-subject training deltas identifies four distinct response programs with different age-group behaviors. Bridges single-protein analysis (F1–F3) with network analysis (F5).

---

## Current file inventory

### Scripts (`a_script/`)
| File | Role |
|------|------|
| `YvO_F4_setup.R` | Setup: loads imputed data, runs Mfuzz FCM, computes ORA, assigns themes, defines unified text sizes |
| `panel_A.R` | Cluster profile spaghetti + ribbon plots (with biological labels) |
| `panel_B.R` | 4 highlight-on-grey PCA sub-panels (one per cluster) |
| `panel_C.R` | Per-cluster triptych (heatmap + Sankey + enrichment bars) |
| `panel_D.R` | Cluster-to-theme synthesis Sankey + stacked bars |
| `supp_dmin_elbow.R` | Supplementary: Dmin elbow + bootstrap ARI |
| `YvO_figure4_composite.R` | Assembly script |

### Outputs (`b_reports/`)
- `Figure_4.pdf` — composite
- `panel_A_profiles.pdf`
- `panel_B_pca.pdf`
- `panel_C_triptych.pdf`
- `panel_D_synthesis.pdf`
- `supp_dmin_elbow.pdf`

### Data (`c_data/`)
- `01_panel_A_cluster_profiles.csv`
- `02_panel_B_pca_scores.csv`
- `03_panel_C_enrichment.csv`
- `04_panel_C_sankey_links.csv`
- `05_panel_D_cluster_themes.csv`
- `06_mfuzz_assignments.csv` — consumed downstream by F5
- `07_cluster_stability.csv`
- `F4_supplementary.xlsx`

## Panel inventory

| Panel | Status | Key message |
|-------|--------|-------------|
| A | Biological labels added | C1: Glycolysis; C2: OXPHOS; C3: mTORC1 Signaling; C4: Myogenesis |
| B | 4 highlight-on-grey sub-panels | Per user preference (reverted from single combined scatter) |
| C | Compressed scale, improved gaps | shared_x_max=15, heatmap-Sankey gap increased, pathway labels readable |
| D | Two-pass theme redesign | 8 themes, 756 proteins mapped, transparency annotation, Hallmark+GO:BP only |
| Dmin (supp) | Functional | k=4 elbow, ARI = 0.956 |

## Known issues

### Unresolved
- [ ] Panel A spaghetti lines add visual noise (alpha tuning may help)
- [ ] Panel C heatmap strip narrow (4 columns) — inherent to 4-group design
- [ ] Composite dimensions (380x300mm) larger than typical journal full-page (~180x240mm)

### Resolved
- [x] Panel A: Biological cluster labels added (C1: Glycolysis, C2: OXPHOS, C3: mTORC1 Signaling, C4: Myogenesis) — uses `top_hallmark` from setup
- [x] Panel B: Reverted to 4 highlight-on-grey sub-panels per user preference
- [x] Subtitle simplified — technical parameters removed, now shows protein/subject counts only
- [x] Key row removed entirely from composite (per user request)
- [x] COL_WIDTHS rebalanced: `c(0.11, 0.12, 0.44, 0.33)`
- [x] Heatmap: blue-white-red diverging palette (`#2166AC` to `#B2182B`)
- [x] Unified text sizes across all panels (TXT_TITLE=9, TXT_SUBTITLE=7, TXT_AXIS=7, TXT_TICK=6.5, TXT_ANNOT=3.0, TXT_HEADER=10), all bold
- [x] Panel C: Compressed enrichment bar scale (shared_x_max_C: 30→15) with "//" break indicators
- [x] Panel C: Increased heatmap–Sankey gap (heatmap margin r: -5→2, X_PW: 1.55→1.85)
- [x] Panel C–D whitespace reduced (D_X_CL: 1.0→0.4, D_X_TH: 4.0→3.4)
- [x] Panel C pathway labels: all bold, readable at panel scale
- [x] Panel D: "Other" exclusion now reported transparently (annotation on plot + console output)
- [x] Panel D: Coverage maximized — 756 themed proteins (up from 541), 8 literature-grounded themes
- [x] GO:CC dropped from enrichment pipeline (Hallmark + GO:BP only)
- [x] Two-pass rescue assignment added for unmapped proteins

## Changes made

### 2026-03-01 — Session 1: High-priority refinements
**Files modified:**
- `panel_A.R`: Added biological cluster labels from `top_hallmark` object. Title now shows "C1: Glycolysis" etc.
- `panel_B.R`: Temporarily changed to single combined scatter, then reverted to 4 sub-panels per user preference.
- `panel_C.R`: Heatmap color — briefly changed to viridis, reverted to blue-white-red per user request.
- `YvO_F4_setup.R`: Updated `COL_WIDTHS` from `c(0.11, 0.13, 0.43, 0.33)` to `c(0.11, 0.12, 0.44, 0.33)`.
- `YvO_figure4_composite.R`: Removed key row and technical note entirely per user request. Simplified subtitle.

### 2026-03-01 — Session 2: Text unification + Panel C legibility
**Files modified:**
- `YvO_F4_setup.R`: Added unified text size constants (TXT_TITLE, TXT_SUBTITLE, TXT_AXIS, TXT_TICK, TXT_ANNOT, TXT_HEADER).
- `panel_A.R`: All text sizes updated to TXT_* constants, all bold.
- `panel_B.R`: All text sizes updated to TXT_* constants, all bold.
- `panel_C.R`: Compressed enrichment bar scale (shared_x_max_C: 30→15), increased heatmap–Sankey gap (X_PW: 1.55→1.85, coord xlim: 1.65→1.95, heatmap margin r: -5→2). All text bold.
- `panel_D.R`: Shifted coordinates left (D_X_CL: 1.0→0.4, D_X_TH: 4.0→3.4, xlim: 0.5→0.05). All text bold.
- `YvO_figure4_composite.R`: Header text uses TXT_HEADER, y-axis label uses TXT_AXIS.

**Validation:** Full composite regenerated successfully. All panels render. Warnings are expected (geom_line removed rows for spaghetti alpha filtering, Unicode en-dash substitution in PDF).

### 2026-03-02 — Session 3: Enrichment coverage improvement
**Files modified:**
- `YvO_F4_setup.R`: Increased top-N term selection from 3 to 5 per database per cluster (line 645). This roughly doubles the pathway pool and greedy assignment coverage.

**Impact:**
- Total enrichment terms: 48 max → 55 actual (capped by significant terms available)
- Mapped proteins: ~23.5% → ~50% of core proteins (1,039 / 2,064)
- Sankey coverage: 541 proteins across 7 themes (up from previous ~23.5%)
- Per-cluster: C1=254/572, C2=266/537, C3=321/569, C4=198/386

**Validation:** Full composite regenerated successfully. Panel D, Panel C enrichment bars, and supplementary data all updated.

### 2026-03-02 — Session 4: Panel D coverage maximization (two-pass theme redesign)

**Problem:** Panel D coverage was 52% (541/1039 mapped proteins). C3 had 91% "Other" because themes lacked proteostasis/cytoskeletal categories. GO:CC terms were "stealing" proteins from functional GO:BP terms during greedy assignment.

**Design:** Brainstormed and approved a literature-validated redesign:
- Drop GO:CC from enrichment (localization terms don't map to functional themes)
- Increase top-N from 5 to 7 (Hallmark + GO:BP only)
- Add two-pass rescue assignment (Pass 2 searches all significant terms, p.adjust < 0.05)
- Redesign themes from 5 → 8 literature-grounded categories

**Files modified:**
- `YvO_F4_setup.R`:
  - Lines 642–649: Filtered to Hallmark + GO:BP only, N=7
  - Lines 701–742: Inserted rescue pass for unmapped proteins
  - Lines 793–861: Redesigned `assign_theme()` with 8 themes + new `THEME_COLORS`
- `panel_C.R`: Updated stat audit comments (GO:CC removed, top-7 noted)
- `panel_D.R`: Added exclusion count computation + transparency annotation
- `YvO_figure4_composite.R`: Updated supplementary metadata descriptions

**New theme categories (8 total):**
1. Mitochondrial & Energy Metabolism (primary: C2)
2. Muscle Structure & Myogenesis (primary: C4)
3. Proteostasis & Stress Response (primary: C3) — NEW
4. Cytoskeletal & Cell Division (primary: C3) — NEW
5. Immune & Complement (primary: C1)
6. ECM & Tissue Remodeling (primary: C4)
7. Metabolic & Redox Regulation (primary: C1)
8. Intracellular Transport & Signaling (primary: C1)

**Coverage results (before → after):**
- Total themed: 541 → 756 proteins (+40%)
- "Other": 498 → 117 (down 77%)
- Unmapped: 1191 (57.7%) — proteins not in any significant pathway
- Per-cluster themes: C1=Metabolic(97)+Transport(71)+Immune(31), C2=Mitochondrial(229), C3=Proteostasis(127)+Cytoskeletal(67), C4=Muscle(49)+ECM(48)

**Bug found and fixed:** Regex patterns in `assign_theme()` used literal spaces (e.g., `"oxidative phosph"`) but MSigDB pathway names use underscores (e.g., `hallmark_oxidative_phosphorylation`). Changed to `.` (regex any-char) in Themes 1 and 3. This fixed Mitochondrial from 42 → 229 proteins.

**Literature validation:**
- FCM clustering on training deltas: Kumar & Futschik 2007 (Mfuzz)
- ORA for cluster characterization: standard approach for discrete gene sets (clusterProfiler; Wu et al. 2021)
- rrvgo for GO redundancy: programmatic REVIGO successor (Sayols 2023)
- Theme categories driven by enrichment results; consistent with Ubaida-Mohien et al. 2019 (eLife; mitochondrial aging proteomics)

**Validation:** Full composite regenerated twice (before and after regex fix). All panels render without error. Panel D Sankey shows 8 theme ribbons with convergence from multiple clusters. Transparency annotation shows excluded protein count.

### 2026-03-02 — Session 5: Verification audit and corrections

**Audit performed:** Systematic cross-validation of all changes from Session 4.

**Literature citations corrected:**
- Removed Robinson 2017 attribution (paper is about translational machinery, not mTORC1/UPR/chaperones)
- Removed Schild 2015 attribution (paper focused on OXPHOS/TCA, not ECM remodeling)
- Removed MoTrPAC 2024 as FCM precedent (used graphical clustering, not FCM)
- Removed Melov 2007 as cytoskeletal reversal citation (reversal was primarily mitochondrial)
- Theme comments now describe themes as driven by enrichment data, citing only Ubaida-Mohien 2019 (confirmed)

**Regex false negatives fixed (3):**
- Added `coagulat` to Theme 5 (Immune): captures HALLMARK_COAGULATION → C1 Immune 31→35
- Added `peroxide` to Theme 7 (Metabolic & Redox): captures GOBP_HYDROGEN_PEROXIDE_CATABOLIC_PROCESS
- Added `actin` to Theme 4 (Cytoskeletal): captures GOBP_REGULATION_OF_ACTIN_FILAMENT_BASED_PROCESS → C4 Cytoskeletal 0→15

**GO:CC computation removed:**
- Removed GO:CC gene set loading, enricher() call, rrvgo reduction, and DB_COLORS entry
- GO:CC was being enriched and reduced but then filtered out at top-N selection; wasted computation

**Stale comments fixed:**
- panel_C.R line 47: "Top 3 per database" → "Top 7 per database (Hallmark + GO:BP)"
- setup.R line 47: header comment updated to remove GO:CC reference

**Coverage after all fixes:**
- Themed: 775 (up from 756), Other: 98 (down from 117)
- Remaining "Other" sources: HALLMARK_HYPOXIA, HALLMARK_UV_RESPONSE_DN, GOBP_NEURON_PROJECTION_GUIDANCE, GOBP_NEGATIVE_REGULATION_OF_GENE_EXPRESSION, GOBP_ESTABLISHMENT_OF_RNA_LOCALIZATION, GOBP_CELL_CELL_RECOGNITION, GOBP_REGULATION_OF_TELOMERE_MAINTENANCE, GOBP_REGULATION_OF_CELLULAR_COMPONENT_SIZE, GOBP_KETONE_BIOSYNTHETIC_PROCESS — these are biologically generic or don't cleanly map to any theme

**Validation:** Full pipeline re-run, all panels render, no regressions.

## Open questions
- Is the current composite size (380 x 300 mm) appropriate or should it be reduced for journal submission?

## Publication-readiness checklist
- [x] Biological cluster labels present in Panel A
- [x] Panel letters visible and consistent
- [x] Pathway names legible at panel scale (compressed x-axis, bold text)
- [x] Panel B shows distinct cluster geometry (highlight-on-grey)
- [x] "Other" exclusion in Panel D documented (annotation + transparency note)
- [ ] Composite dimensions reviewed for journal requirements
- [x] Colors follow shared/palettes.R
- [x] Subtitle concise
- [x] All text sizes unified and bold across panels
- [x] Key row removed (per user preference)

## Next recommended steps
1. Evaluate composite dimensions — current 380x300mm may be too large for journal submission
2. Consider reducing spaghetti line noise in Panel A (increase alpha threshold or simplify)
3. Review Panel D at composite scale — verify 8-theme Sankey ribbons are distinguishable
4. Remaining "Other" proteins (98) are from biologically generic pathways — acceptable for publication
