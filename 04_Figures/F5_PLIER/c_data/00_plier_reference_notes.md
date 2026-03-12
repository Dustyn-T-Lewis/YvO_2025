# PLIER Implementation Reference Notes

## Function signatures (from wgmao/PLIER source)

### PLIER()
```r
PLIER(data, priorMat, svdres = NULL, k = NULL, L1 = NULL, L2 = NULL, L3 = NULL,
      frac = 0.7, max.iter = 350, trace = FALSE, scale = TRUE, Chat = NULL,
      maxPath = 10, doCrossval = TRUE, penalty.factor = rep(1, ncol(priorMat)),
      glm_alpha = 0.9, minGenes = 10, tol = 1e-06, seed = 123456,
      allGenes = FALSE, rseed = NULL, pathwaySelection = c("complete", "fast"))
```

### num.pc()
```r
num.pc(data, method = "elbow", B = 20, seed = NULL)
# Returns: integer (estimated significant PCs)
```

### Key helpers
```r
commonRows(data1, data2)   # shared rownames
rowNorm(x)                 # z-score each row
combinePaths(...)          # merge binary pathway matrices
nameB(plierRes, top = 1, fdr.cutoff = 0.01)  # name LVs by top pathway
plotU(plierRes, auc.cutoff = 0.6, fdr.cutoff = 0.05, top = 3)
```

## Return object fields
- B: k × p (LV scores, rows=LVs, cols=samples)
- Z: n × k (gene loadings)
- U: m × k (pathway utilization, sparse non-negative)
- C: n × m (filtered prior matrix used)
- Uauc: m × k (cross-validation AUC per pathway-LV)
- Up: m × k (p-values for pathway-LV)
- summary: data.frame (pathway, LV index, AUC, p, FDR)
- withPrior: integer vector (LV indices with significant pathway associations)
- residual: scalar (reconstruction error)
- L1, L2, L3: final lambda values

## Prior matrix best practices
- Reactome + KEGG + Hallmark is a good default for non-immune tissue
- Blood cell markers (IRIS/DMAP, svmMarkers) are irrelevant for skeletal muscle
- PLIER handles redundant/overlapping pathways well
- minGenes=10 default; may need to lower for proteomics
- Build via msigdbr → binary matrix conversion

## k selection
- num.pc() gives estimated significant PCs
- Paper recommends 2x multiplier
- Method is not sensitive to exact k (robust to overestimation)
- For ~2000 genes × 62 samples, expect ~10-20 PCs → k ≈ 20-40

## Proteomics-specific considerations
1. No published PLIER-on-proteomics precedent
2. Gene symbol duplicates must be resolved before fitting
3. Fewer features = fewer surviving pathways
4. Data is already log2-transformed (DIA-MS intensity)
5. scale=TRUE in PLIER will z-score rows (appropriate)

## Pitfalls to avoid
- Wrong matrix orientation (must be genes × samples)
- Duplicate rownames (causes silent errors)
- Not using commonRows() to intersect data and prior
- Running PLIER per group instead of all samples together
- SVD errors from zero-variance genes (filter first)
