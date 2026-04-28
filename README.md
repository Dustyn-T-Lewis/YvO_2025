# YvO 2026 — Young vs Old Skeletal Muscle Proteome (Pre / Post Resistance Training)

This repository contains the analysis pipeline and supplementary outputs
for the YvO study: a label-free quantitative proteomic comparison of
skeletal muscle in young and older adults before and after a resistance
training intervention.

The repository is organized as a numbered, sequential pipeline
(`00_input` → `01_normalization` → `02_Imputation` → `03_DEP` →
`04_Figures`). Each stage is self-contained: scripts in `a_script/`,
PDF/HTML reports in `b_reports/`, and tabular/RDS outputs in `c_data/`.

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
- `c_data/02_imputation.xlsx` — 5-sheet supplement (classification
  table, imputation summary, OOB error, MNAR audit, parameter log).
- `c_data/02_mar_mnar_classification.csv`,
  `03_imputation_mask.csv`, `04_mnar_imputation_audit.csv`.
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
See `04_Figures/RUN_ORDER.md` for the dependency-correct execution
order. F06 (WGCNA) is the bottleneck — F04, F05, and F07 read its
module assignments.

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
  F00 ... F07/              one directory per main figure
    a_script/               panel + stitch scripts
    b_reports/main, supp/   composite PDFs + PNG renders
    c_data/                 supplement workbooks
  shared/                   palettes, fgsea cache, helpers
  RUN_ORDER.md              dependency-correct run order
docs/                       audit / methodology notes
```
