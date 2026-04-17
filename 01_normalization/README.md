# Stage 01: Normalization

Ingest raw DIA-MS data, clean, filter, normalize, and emit canonical artifacts for downstream stages.

## Scripts

| Script | Role |
|--------|------|
| `a_script/01_run_normalization.R` | Main pipeline: filter, outlier detection, normalization |
| `a_script/02_norm_reports.R` | Diagnostic plots + supplementary workbook |

## Outputs

| File | Description |
|------|-------------|
| `c_data/02_normalized.csv` | Cycloess-normalized protein matrix (consumed by stages 02 + 03) |
| `c_data/03_DAList_normalized.rds` | proteoDA DAList object with metadata |
| `c_data/01_DAList_prenorm.rds` | Pre-normalization DAList (QC reference) |
| `c_data/00_report_intermediates.rds` | Plot data handoff to 02_norm_reports.R |
| `c_data/05_normalization_supp.xlsx` | Supplementary workbook (3 sheets) |
| `b_reports/01_norm_comparison.pdf` | proteoDA normalization comparison |
| `b_reports/02_qc_pre.pdf` | Pre-normalization QC |
| `b_reports/03_qc_post.pdf` | Post-normalization QC |
| `b_reports/04_diagnostics.pdf` | Custom 4-page diagnostic report |

## Methodology

1. **HPA tissue filter** — retain proteins with Human Protein Atlas skeletal muscle evidence
2. **Blood contaminant removal** — top plasma proteins (Geyer et al. 2016, PMID 27135364) + HPA immunoglobulin genes
3. **UniProt ID deduplication** — keep highest mean intensity per duplicate
4. **Missingness filter** — require >= 10 non-missing samples in at least one Group_Time level
5. **4-method consensus outlier detection** — missingness IQR, PCA Mahalanobis (chi-sq p < 0.01), MAD median intensity (k = 3), inter-sample correlation (k = 3 MAD); remove if >= 3/4 agree
6. **Cyclic loess normalization** via proteoDA (Bolstad 2003; Thurman 2023)

CV is reported on the linear scale where applicable (Brenes 2024).

## Guardrails

- Do not change normalization methodology without strong justification
- Do not rename outputs without updating all downstream references in stages 02, 03, and 04_Figures
- The canonical outputs (`02_normalized.csv`, `03_DAList_normalized.rds`) are consumed by every downstream stage
