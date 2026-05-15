# YvO 2026 — Young vs Old Skeletal Muscle Proteome (Pre / Post Resistance Training)

This repository contains the analysis pipeline and supplementary outputs
for the YvO study: a label-free quantitative proteomic comparison of
skeletal muscle in young and older adults before and after a resistance
training intervention.

The repository is organized as a numbered, sequential pipeline
(`00_input` → `01_normalization` → `02_Imputation` → `03_DEP` →
`04_Figures`). Each stage is self-contained: scripts in `a_script/`,
PDF/HTML reports in `b_reports/`, and tabular/RDS outputs in `c_data/`.

Install dependencies once after cloning: `Rscript setup.R`.

---

## Study design

| Factor    | Levels                            |
|-----------|-----------------------------------|
| Group     | Young, Old                        |
| Timepoint | Pre, Post (resistance training)   |
| Design    | Repeated measures within subject  |

Four planned contrasts are tested at the proteome level:

- `Training_Young = Young_Post − Young_Pre`
- `Training_Old   = Old_Post − Old_Pre`
- `Aging          = Old_Pre − Young_Pre`
- `Interaction    = (Old_Post − Old_Pre) − (Young_Post − Young_Pre)`

Significance threshold: raw p ≤ 0.10 with π-score ≥ 0.05
(Xiao et al. π-score combines magnitude and significance).

---

## Inputs (`00_input/`)

| File                                      | Contents                                          |
|-------------------------------------------|---------------------------------------------------|
| `YvO_raw.xlsx`                            | Label-free intensities + UniProt annotation       |
| `YvO_meta.xlsx`                           | Sample metadata (Col_ID, Group, Timepoint)        |
| `YvO_pheno_calc.xlsx`                     | Phenotypic measures (training volume, DXA, VL)    |
| `HPA_skeletal_muscle_annotations.tsv`     | Human Protein Atlas skeletal-muscle reference     |

---

## Stage 01 — Normalization (`01_normalization/`)

Filtering cascade (top to bottom; counts logged in
`c_data/01_normalization.xlsx`):

1. **HPA tissue filter** — keep only proteins present in the HPA
   skeletal-muscle reference.
2. **Blood-contaminant removal** — Geyer (2016) plasma protein list
   plus HPA-annotated immunoglobulins (HBA/HBB, ALB, TF, complement,
   apolipoproteins, fibrinogens, IGH/IGK/IGL families).
3. **Deduplication by UniProt ID** — when multiple rows share an
   accession, keep the highest mean-intensity row.
4. **Missingness filter** — require detection in **≥ 10 replicates**
   total (`MIN_REPS = 10`).
5. **Outlier detection** — 4-method consensus (sample-level PCA
   distance, Mahalanobis with `MAHAL_P = 0.01`, MAD with `MAD_K = 3`,
   k-NN); samples flagged by ≥ 3 of 4 methods are removed.
6. **Normalization** — `cycloess` (cyclic LOESS) via `proteoDA`.

Outputs:

- `c_data/01_normalization.xlsx` — 8-sheet supplement (filter log,
  blood/Ig removals, missingness summary, outlier consensus, PCA pre/
  post, cycloess diagnostics).
- `c_data/03_DAList_normalized.rds` — `proteoDA` object passed to
  stages 02–04.
- `c_data/02_normalized.csv` — text-serialized matrix used by stages 02
  and 03 to ensure float reproducibility across machines.
- `b_reports/01_qc_filter.pdf`, `02_qc_pre.pdf`, `03_qc_post.pdf`,
  `04_diagnostics.pdf`.

## Stage 02 — Imputation (`02_Imputation/`)

1. **MAR / MNAR classification** — 3-method consensus (proportion
   missing, intensity-rank logistic classifier, `missMethyl`-style
   left-censoring test). Each protein labelled Complete / MAR / MNAR.
2. **Imputation** — `missForest` (random-forest based) on the full
   matrix; per-protein OOB error tracked. Proteins above
   `MISS_UNRELIABLE = 50%` flagged as low-confidence.

Outputs:

- `c_data/01_DAList_imputed.rds` — `proteoDA` object with imputation
  annotations.
- `c_data/02_imputation.xlsx` — multi-sheet supplement (imputed
  matrix, MAR/MNAR classification, imputation mask, MNAR audit,
  summary, and benchmark ranking when the benchmark has been run).
- `c_data/01_imputed.csv`, `c_data/02_mar_mnar_classification.csv` —
  CSV mirrors of the imputed matrix and classification table.
- `b_reports/01_missingness_report.pdf`,
  `02_imputation_report.pdf`.

## Stage 03 — Differential Expression (`03_DEP/`)

`limma` with `duplicateCorrelation` (within-subject blocking).

- Input: cycloess-normalized **non-imputed** matrix; `limma` handles
  per-protein NAs internally.
- Design: `~ 0 + group + (1 | subject)` with the four contrasts above.
- Empirical Bayes shrinkage; π-score computed per protein per contrast.
- Robustness pass: `03_run_robustness.R` repeats the fit with leave-
  one-subject-out and permutation-based p-values (Phipson–Smyth).

Outputs:

- `c_data/01_limma_DAList.rds` — fitted `proteoDA` object.
- `c_data/03_combined_results.csv` — wide-format results consumed by
  Stage 04 figures.
- `c_data/03_DEP_results.xlsx` — multi-sheet workbook (per-contrast
  results, summary, parameter log).
- `b_reports/01_proteoDA/` — `proteoDA` HTML reports + static plots.

---

## Stage 04 — Figures (`04_Figures/`)

Each figure directory follows the same layout
(`a_script/`, `b_reports/`, `c_data/`) and emits both a main composite
PDF and a supplementary workbook (`F0X_supplementary.xlsx`).

**Per-figure run order:**
`a_script/01_main_panels.R` → `02_supp_panels.R` → `90_stitch_F0X.R`.
The stitcher sources the panel scripts in the right order, builds the
composite, writes `F0X_supplementary.xlsx`, and cleans up consumed
panel CSVs.

**Cross-figure dependency:**
F06 (WGCNA) must run before F07 (F07 reads F06 module assignments).
F00 runs last — it reads stage 01/02/03 report intermediates to
summarize the entire pipeline. F01–F05 are independent and can run
in any order. F06's `YvO_WGCNA_run.R` is run separately before
`90_stitch_F06.R` (it's the modeling step, not a panel script).

| Figure | Title                       | Main panels                                                                |
|--------|-----------------------------|----------------------------------------------------------------------------|
| F00    | Pipeline QC                 | Filter cascade, sample missingness, PCA pre/post, MAR/MNAR, imputation     |
| F01    | Phenotype                   | Training volume, DXA lean body mass, VL thickness                          |
| F02    | Proteome overview           | 6-panel grid (DEP counts, p-distributions, Pi-score, contrast overlap)     |
| F03    | Volcano rings               | Four ring volcanoes (Aging, Training_Young, Training_Old, Interaction)     |
| F04    | Training concordance        | Chord diagram + module overlap (Young vs Old training response)            |
| F05    | Aging reversal              | Chord diagram + reversal classification (proteins moving back toward Young)|
| F06    | WGCNA module–trait          | Module–trait heatmap, eigengene trajectories, hub proteins                 |
| F07    | Phenotype prediction        | LOSO classifier AUC, top-feature hero grid, ROC supplements                |

Shared utilities live under `04_Figures/shared/` (style, palettes,
volcano-ring renderer, fgsea cache, GO-slim categories).

---

## Manuscript outputs

### Main figures (Figures 1–7)

Each main composite PDF is at:
`04_Figures/F0X/b_reports/main/pdf/MAIN_F0X_composite.pdf`

| Figure | Path |
|--------|------|
| Fig 1 | `04_Figures/F01/b_reports/main/pdf/MAIN_F01_composite.pdf` |
| Fig 2 | `04_Figures/F02/b_reports/main/pdf/MAIN_F02_composite.pdf` |
| Fig 3 | `04_Figures/F03/b_reports/main/pdf/MAIN_F03_composite.pdf` |
| Fig 4 | `04_Figures/F04/b_reports/main/pdf/MAIN_F04_composite.pdf` |
| Fig 5 | `04_Figures/F05/b_reports/main/pdf/MAIN_F05_composite.pdf` |
| Fig 6 | `04_Figures/F06/b_reports/main/pdf/MAIN_F06_composite.pdf` |
| Fig 7 | `04_Figures/F07/b_reports/main/pdf/MAIN_F07_composite.pdf` |

Individual panels (PNG/PDF) are in the `panels/` subdirectory alongside
each composite.

### Supplementary figures (S1–S9)

Each supplementary composite PDF is at:
`04_Figures/F0X/b_reports/supp/pdf/SUPP_F0X_*.pdf`

| Ref | Source | Path |
|-----|--------|------|
| S1 Fig | F00 | `04_Figures/F00/b_reports/supp/pdf/SUPP_F00_normalization.pdf` |
| S2 Fig | F00 | `04_Figures/F00/b_reports/supp/pdf/SUPP_F00_imputation.pdf` |
| S3 Fig | F01 | `04_Figures/F01/b_reports/supp/pdf/SUPP_F01_composite.pdf` |
| S4 Fig | F02 | `04_Figures/F02/b_reports/supp/pdf/SUPP_F02_composite.pdf` |
| S5 Fig | F03 | `04_Figures/F03/b_reports/supp/pdf/SUPP_F03_composite.pdf` |
| S6 Fig | F04 | `04_Figures/F04/b_reports/supp/pdf/SUPP_F04_diagnostics.pdf` |
| S7 Fig | F05 | `04_Figures/F05/b_reports/supp/pdf/SUPP_F05_diagnostics.pdf` |
| S8 Fig | F06 | `04_Figures/F06/b_reports/supp/pdf/SUPP_F06_composite.pdf` |
| S9 Fig | F07 | `04_Figures/F07/b_reports/supp/pdf/SUPP_F07_composite_main.pdf` |

### Supplementary tables (S1–S11)

Tables S1–S3 are produced by stages 01–03. Tables S4–S11 are produced by
each figure's `90_stitch_F0X.R` stitcher.

| Ref | Source | Pipeline path |
|-----|--------|---------------|
| S1 Table | Stage 01 | `01_normalization/c_data/01_normalization.xlsx` |
| S2 Table | Stage 02 | `02_Imputation/c_data/02_imputation.xlsx` |
| S3 Table | Stage 03 | `03_DEP/c_data/03_DEP_results.xlsx` |
| S4 Table | F00 | `04_Figures/F00/c_data/F00_supplementary.xlsx` |
| S5 Table | F01 | `04_Figures/F01/c_data/F01_supplementary.xlsx` |
| S6 Table | F02 | `04_Figures/F02/c_data/F02_supplementary.xlsx` |
| S7 Table | F03 | `04_Figures/F03/c_data/F03_supplementary.xlsx` |
| S8 Table | F04 | `04_Figures/F04/c_data/F04_supplementary.xlsx` |
| S9 Table | F05 | `04_Figures/F05/c_data/F05_supplementary.xlsx` |
| S10 Table | F06 | `04_Figures/F06/c_data/F06_supplementary.xlsx` |
| S11 Table | F07 | `04_Figures/F07/c_data/F07_supplementary.xlsx` |

---

## Reproducibility notes

- All scripts set `set.seed(42)` before any stochastic step.
- The numeric matrix passed between Stage 01 and Stages 02 / 03 is
  written to CSV (`02_normalized.csv`) deliberately — RDS-binary
  doubles differ from CSV round-trips at machine epsilon, which is
  enough to shift `missForest` tree splits.
- Output regeneration: re-run a stage's `01_*` script first, then the
  reports script, then the figure stitchers. Each figure's
  `90_stitch_F0X.R` rebuilds its supplementary workbook from
  intermediate panel CSVs and cleans them up afterward.

### Box delivery (optional)

The figure stitchers (`90_stitch_F0X.R`) and stage report scripts
copy main composites and supplementary tables/figures to a Box
directory if it exists. The destination is configurable via two
environment variables:

| Variable | Targets | Used by |
|----------|---------|---------|
| `YVO_BOX_DIR`  | Manuscript root (`02_Figures/`, `03_Supplementary/{figures,tables}/`) | All `04_Figures/F0X/a_script/90_stitch_F0X.R` |
| `YVO_BOX_SUPP` | `03_Supplementary/` subdir (tables only) | `01_normalization/a_script/02_generate_reports.R`, `02_Imputation/a_script/02_generate_reports.R`, `03_DEP/a_script/03_run_robustness.R` |

If unset, both fall back to the author's local Box path. If that
path doesn't exist on the running machine, the copy block is
skipped silently and the in-repo outputs are still produced. To
mirror to a different location:

```sh
export YVO_BOX_DIR=/path/to/manuscript_root
export YVO_BOX_SUPP=/path/to/manuscript_root/03_Supplementary
Rscript 04_Figures/F06/a_script/90_stitch_F06.R
```

---

## Repository layout

```
00_input/                   raw data + HPA reference
01_normalization/
  a_script/                 normalize, generate_reports
  b_reports/                QC PDFs
  c_data/                   normalized DAList, supplement xlsx
02_Imputation/
  a_script/                 impute, generate_reports, benchmark/
  b_reports/                missingness + imputation reports
  c_data/                   imputed DAList, supplement xlsx
03_DEP/
  a_script/                 run_dep, generate_reports, run_robustness
  b_reports/                proteoDA HTML + static plots
  c_data/                   limma DAList, combined results, supp xlsx
04_Figures/
  F00 ... F07/              one directory per figure
    a_script/               panel + stitch scripts
    b_reports/main, supp/   composite PDFs + PNG renders
    c_data/                 supplement workbooks
  shared/                   palettes, fgsea cache, helpers
docs/                       audit / methodology notes
```
