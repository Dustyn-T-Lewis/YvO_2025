# Figure Execution Order

Figures must run in dependency order. Each stage writes to `c_data/` and `b_reports/`.

## Sequence

```
F00  (QC — standalone, reads 01-03 intermediates)
F01  (Phenotype — standalone)
F02  (Proteome overview — reads shared/fgsea_tstat_all_v2.csv)
F03  (Volcano rings — reads shared/fgsea_tstat_all_v2.csv)
F04  (Training concordance — reads shared/fgsea_tstat_all_v2.csv)
F05  (Aging reversal — reads shared/fgsea_tstat_all_v2.csv)
F06  (WGCNA — standalone; generates c_data/wgcna/ + F06_supplementary.xlsx)
F07  (Phenotype prediction — reads F06/c_data/*)
```

## Critical dependencies

- **F02, F03, F04, F05** all read `shared/fgsea_tstat_all_v2.csv` (frozen fGSEA cache)
- **F07** depends on F06: reads WGCNA module data from `F06/c_data/` and sheets from `F06/c_data/F06_supplementary.xlsx`
- **F04, F05** supp chord scripts read module assignments from `F06/c_data/F06_supplementary.xlsx`

## Per-figure execution

Each figure directory contains:
- `a_script/panels/panel_*.R` — individual panel scripts
- `a_script/90_stitch_figure.R` — sources panels, builds composite + supplementary xlsx
- `a_script/supp/` — supplementary panel scripts + supp stitchers

Run each figure by executing `90_stitch_figure.R` (it sources everything else).

## Upstream pipeline (must complete before any figure)

```
01_normalization/a_script/01_run_normalization.R
01_normalization/a_script/02_norm_reports.R
02_Imputation/a_script/01_apply_missforest.R
02_Imputation/a_script/02_imputation_reports.R
03_DEP/a_script/01_run_dep.R
03_DEP/a_script/02_dep_reports.R
03_DEP/a_script/03_dep_robustness.R
```
