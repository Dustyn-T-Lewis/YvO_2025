# Stage 02: Imputation

Characterize missingness, apply missForest imputation, and emit imputation artifacts.

## Scripts

| Script | Role |
|--------|------|
| `a_script/01_apply_missforest.R` | MAR/MNAR classification + missForest imputation |
| `a_script/02_imputation_reports.R` | Missingness + imputation quality report PDFs |
| `a_script/benchmark/` | 23-method imputation benchmark (separate from canonical pipeline) |

## Outputs

| File | Description |
|------|-------------|
| `c_data/01_imputed.csv` | Imputed protein matrix (used by WGCNA in F06, robustness check in stage 03) |
| `c_data/01_DAList_imputed.rds` | DAList with missingness annotations |
| `c_data/02_mar_mnar_classification.csv` | Per-protein MAR/MNAR classification |
| `c_data/03_imputation_mask.csv` | Boolean mask of which values were imputed |
| `c_data/04_mnar_imputation_audit.csv` | MNAR protein shift audit |
| `c_data/05_imputation_summary.txt` | Summary statistics |
| `b_reports/01_missingness_report.pdf` | Missingness characterization |
| `b_reports/02_imputation_report.pdf` | Benchmark + imputation quality |

## Methodology

- **MAR/MNAR classification**: `msImpute::selectFeatures()` EBM method with k-means fallback. Used for flagging and audit only — does NOT gate imputation.
- **Imputation**: `missForest::missForest()` with `set.seed(42)`, treating all missing values as MAR per the algorithm.
- **Benchmark validation**: missForest ranked #1 in 23-method benchmark (composite 0.928, NRMSE-MCAR 0.129, fold-change preservation rho = 0.977, NES preservation rho = 0.991).

Note: The primary DEP analysis (stage 03) runs on **non-imputed** data. The imputed matrix is used for WGCNA module detection (F06) and as an imputation sensitivity check in stage 03's robustness analysis.

## Guardrails

- Do not change imputation method without re-running the benchmark
- Do not collapse benchmark outputs into a single opaque artifact
- The MAR/MNAR classification is diagnostic — do not use it to drive class-specific imputation
