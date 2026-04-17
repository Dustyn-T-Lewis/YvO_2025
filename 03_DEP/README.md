# Stage 03: Differential Expression

Fit the limma model on cycloess-normalized, non-imputed data and emit canonical DEP artifacts consumed by every downstream figure.

## Scripts

| Script | Role |
|--------|------|
| `a_script/01_run_dep.R` | limma + duplicateCorrelation + Pi-score |
| `a_script/02_dep_reports.R` | Per-contrast volcanos, summary bar chart, outlier sensitivity |
| `a_script/03_dep_robustness.R` | Blunting diagnostics, bootstrap CIs, power analysis, imputation sensitivity |

## Outputs

| File | Description |
|------|-------------|
| `c_data/01_limma_DAList.rds` | Fitted DAList with eBayes results |
| `c_data/02_DA_summary.csv` | Per-contrast DEP counts |
| `c_data/03_combined_results.csv` | Wide-format results (all 5 contrasts) — primary input for all figures |
| `c_data/04_per_contrast_results/*.csv` | Per-contrast stats-only CSVs |
| `c_data/05_results.xlsx` | Multi-sheet results workbook |
| `c_data/10_DEP_supplementary.xlsx` | Robustness analyses workbook |
| `b_reports/01_proteoDA/` | proteoDA volcano + p-histogram plots |

## Methodology

- **Input**: Cycloess-normalized, NON-imputed matrix from stage 01 (limma handles per-protein NAs; Karpievitch 2012)
- **Model**: `~ 0 + group` with subject blocking via `duplicateCorrelation()` (Smyth 2005; iterated per Phipson 2016)
- **5 contrasts**: Aging (Old_Pre - Young_Pre), Training_Young, Training_Old, Interaction, Reversal
- **Reversal** = Training_Old - Aging. Non-orthogonal — shares Old_Pre with opposite sign. Computed for ranking only; NOT reported as primary FDR-significant findings.
- **Thresholds**: FDR 0.10 (BH) for exploratory discovery; Pi-score (P^|logFC| < 0.05) as secondary effect-size-weighted filter. Both 0.05 and 0.10 FDR reported.
- **Pi-score**: Xiao et al. 2014, PMID 22321699

## Important: Subject_ID is not globally unique

10 Subject_IDs are shared between Young and Old cohorts (different individuals). Always derive blocking keys from `sub("_(Pre|Post)$", "", Col_ID)`, never from `meta$Subject_ID`.

## Guardrails

- Do not change contrasts, thresholds, or blocking without strong evidence
- Do not switch between imputed and non-imputed input without justification
- Do not promote Reversal to primary-finding status without acknowledging the shared-variance caveat
