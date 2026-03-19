# Interrogating NES Scatter Quadrants: Approaches for Multi-Resolution Pathway Analysis

## Context and Problem Statement

Two NES scatter plots form the pathway-level backbone of the concordance (F04 Panel G) and reversal (F05 Panel E) analyses:

- **Panel G (Concordance)**: Training_Young NES vs Training_Old NES across 69 a priori terms (42 Hallmark + 27 GO Slim). Asks whether the same pathways respond to exercise in both age groups.
- **Panel E (Reversal)**: Aging NES vs Training_Old NES across the same 69 terms. Asks whether training in old adults reverses aging-associated pathway changes.

Both panels plot all 69 terms (NS in grey, significant in color), compute Spearman rho on all terms, and annotate quadrant counts on significant terms only. The a priori collection is deliberately broad -- GO Slim terms are high-level biological process categories (62 terms, 27 passing size filters), and Hallmark gene sets are curated, non-redundant summaries of coherent biological states.

**The question**: Can we trace finer biological dynamics within specific quadrants (e.g., what specific processes drive discordance?) without flooding the plot, violating the a priori principle, or cherry-picking?

### What we already have

The pipeline has already computed fGSEA on a much larger collection per contrast:

| Database | Terms tested per contrast | Significant (padj<0.05), Aging | Sig, Training_Young | Sig, Training_Old |
|----------|--------------------------|-------------------------------|---------------------|-------------------|
| GO:BP | ~998 | 57 | 71 | 30 |
| Reactome | ~432 | 135 | 217 | 41 |
| KEGG | ~35 | 16 | 25 | ~5 |
| WikiPathways | ~191 | 19 | 15 | ~5 |
| GO Slim | ~27 | 5 | 10 | 9 |
| Hallmark | ~42 | 4 | 9 | 4 |

These results are cached in `fgsea_tstat_all_v2.csv` (~2,170 terms per contrast, 10,862 rows total) and `fgsea_gobp_v2.csv` (~998 GO:BP terms per contrast, 4,992 rows). Leading edge genes are stored for every term/contrast combination.

Additionally, protein-level RRHO2 (F04 Panel E, F05 Panel D) already performs per-quadrant ORA on RRHO hotspot genes, identifying specific enriched terms within concordant and discordant protein sets.

---

## Approach 1: Hierarchical Enrichment (Broad Screen then Focused Follow-Up)

### Concept

Run fGSEA on the full GO:BP or Reactome collection (already done), then use the a priori NES scatter quadrants to define "regions of interest." Within each quadrant, filter the full enrichment results to terms that fall in that quadrant (same NES sign pattern) and are significant in at least one contrast. Present these as a focused drill-down.

### Methodological basis

This is conceptually a two-stage analysis: stage 1 is the a priori broad screen (69 terms, main figure), stage 2 is the hypothesis-generating fine-grained screen (restricted to quadrant-consistent terms, supplementary). The two-stage approach avoids the multiple-testing burden of displaying all terms while maintaining the a priori scaffold.

**Key distinction**: The full-collection fGSEA is run on ALL terms with genome-wide BH correction. The quadrant filtering is a post-hoc selection on already-corrected results -- it does NOT re-test or re-correct, so there is no statistical inflation. This is analogous to how Subramanian et al. (2005) describe the leading-edge analysis as a "refinement step" after the primary GSEA.

### Supporting literature

- **Subramanian et al. (2005)** PMID 16199517 -- Original GSEA paper. Section on leading-edge analysis explicitly describes post-hoc refinement of enrichment results as a valid interpretive step, distinct from the primary statistical test.
- **Reimand et al. (2019)** PMID 30664679 -- "Pathway enrichment analysis and visualization of omics data," Nat Protocols. Recommends a workflow of: (1) run enrichment broadly, (2) reduce redundancy, (3) visualize as networks/maps. Explicitly endorses multi-resolution presentation where main figures show consolidated results and supplements show fine-grained results.
- **Savage et al. (2019)** PMID 31142576 -- "Graph Algorithms for Condensing and Consolidating Gene Set Analysis Results," Mol Cell Proteomics. Introduces graph-based methods for navigating between broad and specific enrichment results, directly relevant to the zoom-in concept.
- **Korotkevich et al. (2021)** -- fgsea preprint (bioRxiv 060012). The `collapsePathways()` function in fgsea is itself a hierarchical tool: it identifies "independent" (non-redundant) pathways from a larger set, which is already used in `pathway_utils.R`.

### Feasibility for this dataset

**High**. The data already exists (`fgsea_tstat_all_v2.csv` and `fgsea_gobp_v2.csv`). Implementation would be:
1. Pivot all ~998 GO:BP terms to wide format (same as Panel G/E already does for the 69 a priori terms).
2. Filter to terms where both contrasts have valid NES values.
3. Classify into quadrants by NES sign.
4. Within each quadrant, select terms significant in at least one contrast (padj < 0.05 after genome-wide BH).
5. Apply Jaccard deduplication (already implemented in `deduplicate_enrichment()`).
6. Present as a supplementary table or focused panel.

### Pitfalls

- **Redundancy explosion**: Full GO:BP has massive overlap between terms. Without aggressive deduplication, the drill-down will be uninterpretable. Jaccard dedup at 0.35-0.50 is essential.
- **NES comparability across collections**: NES from GO Slim terms (large, ~50-222 genes) and GO:BP terms (smaller, ~10-50 genes) have different variance. Smaller gene sets produce noisier NES estimates. This is not a statistical error, but a visual interpretation issue.
- **Interpretation trap**: A term significant in Training_Young but not Training_Old with NES signs in the discordant quadrant is NOT the same as a formally tested interaction. The NES scatter provides visual summary, not formal discordance testing.

### Verdict: RECOMMENDED as primary approach

This is the most natural and statistically defensible zoom-in. It preserves the a priori structure in the main figure while providing mechanistic detail in supplements. The data is already computed.

---

## Approach 2: ORA on Quadrant-Specific Proteins

### Concept

Instead of working at the pathway level, go back to the protein level. Take the proteins that drive each quadrant (e.g., proteins significantly up in Training_Young but significantly down in Training_Old) and run ORA on them. This moves from protein-level to pathway, providing a different perspective from the pathway-level NES approach.

### Methodological basis

This is already partially implemented in the RRHO2 panels (F04 Panel E, F05 Panel D), where RRHO hotspot genes are subjected to ORA. The approach here would use the DEP classification (Pi-score < 0.05) rather than RRHO hotspot membership to define the foreground sets.

The key question is: what defines "quadrant-specific proteins"?

**Option A -- Significance-based**: Proteins with Pi_Training_Young < 0.05 AND Pi_Training_Old >= 0.05 AND logFC_Training_Young > 0 AND logFC_Training_Old < 0. This is very strict (few proteins, likely underpowered).

**Option B -- Rank-based**: Proteins in the top/bottom N% of t-statistics for each contrast, classified by the sign combination. This is what RRHO effectively does.

**Option C -- Interaction-filtered**: Proteins with Pi_Interaction < 0.05. These are formally discordant. But with only 34 Interaction DEPs, ORA will be underpowered.

### Supporting literature

- **Huang et al. (2009)** PMID 19033363 -- "Bioinformatics enrichment tools: paths toward the comprehensive functional analysis of large gene lists." Classifies enrichment approaches and discusses the validity of ORA on subsets defined by multiple criteria.
- **Reimand et al. (2019)** PMID 30664679 -- Explicitly warns against ORA on arbitrarily filtered subsets but endorses it when the filtering criteria are biologically motivated and pre-specified.
- **Existing RRHO2 ORA** (F04 Panel E) -- Already demonstrates this approach. The concordant-up quadrant ORA found CCT/TRiC chaperonin complex (padj = 2.2e-7), protein folding, cilium assembly. The discordant ORA found respiratory electron transport (padj = 9.5e-5), striated muscle contraction.

### Feasibility for this dataset

**Moderate**. Feasibility depends on foreground set size:
- Training_Young DEPs (Pi < 0.05): 99 proteins. Reasonable for ORA.
- Training_Old DEPs (Pi < 0.05): 18 proteins. Marginal for ORA.
- Interaction DEPs: 34. Marginal.
- Cross-classification (e.g., sig in Y but not O, with specific sign): likely 10-30 proteins per quadrant. Underpowered for most pathway databases.

The RRHO-based approach (already implemented) may actually be better here, because it uses continuous rank thresholds rather than binary significance cutoffs, giving larger foreground sets in the quadrant tails.

### Pitfalls

- **Foreground set size**: ORA requires a sufficiently large foreground (typically >= 20-50 genes) for adequate power. The strict quadrant definitions may yield too few proteins.
- **Universe definition**: The ORA background must be the full set of quantified proteins (~2,138), not the genome. This is already handled correctly in `run_ora_deduplicated()`.
- **Circularity with RRHO ORA**: If this approach produces the same results as the existing RRHO ORA, it adds no new information. The value is only in using a DIFFERENT quadrant definition (DEP-based vs. RRHO-based).
- **Binary threshold artifacts**: ORA is sensitive to the chosen significance cutoff. A protein with Pi = 0.051 vs 0.049 is biologically equivalent but classified differently.

### Verdict: ALREADY PARTIALLY IMPLEMENTED; limited additional value

The RRHO2 per-quadrant ORA (F04 Panel E, F05 Panel D) already does this with a more powerful continuous-rank approach. Adding a DEP-classification-based ORA would be redundant unless the goal is to compare the two definitions. Better to reference the existing RRHO ORA results when discussing quadrant biology.

---

## Approach 3: Leading Edge Intersection

### Concept

fGSEA computes a "leading edge" for each pathway -- the subset of genes that contribute most to the enrichment signal. For a pathway that is significant in both Training_Young and Training_Old (or Aging and Training_Old), the leading edges may overlap substantially (same genes driving the signal) or diverge (different genes within the same pathway responding differently). Comparing leading edges between contrasts reveals the molecular basis of concordance or discordance.

### Methodological basis

The leading-edge concept originates from Subramanian et al. (2005) PMID 16199517, who describe it as the "core of a gene set that accounts for the enrichment signal." The original GSEA paper explicitly describes a "Leading-Edge Subset Analysis" (Section 3.3) that intersects leading edges across multiple gene sets to find shared driver genes. Extending this to intersect leading edges across contrasts for the SAME gene set is a natural generalization.

### What this would look like

For a pathway like GOSLIM_PROTEIN_FOLDING (significant in both contrasts in Panel G):
1. Extract leading edge for Training_Young: e.g., {CCT2, CCT3, HSP90AB1, HSPA5, ...}
2. Extract leading edge for Training_Old: e.g., {CCT2, CCT3, HSPB1, BAG3, ...}
3. Intersection = shared drivers of concordance
4. Symmetric difference = contrast-specific drivers

For a discordant pathway (significant NES in opposite directions):
1. Leading edge in Training_Young (upregulated): {gene set A}
2. Leading edge in Training_Old (downregulated): {gene set B}
3. These sets are conceptually different -- set A contains genes driving upregulation, set B contains genes driving downregulation. Their overlap reveals genes that the pathway includes but that respond in opposite directions.

### Supporting literature

- **Subramanian et al. (2005)** PMID 16199517 -- Defines the leading-edge analysis. Section 3.3 explicitly performs cross-gene-set leading-edge intersection.
- **The original GSEA paper (Mootha et al. 2003)** PMID 12808457 -- First application of GSEA; leading edge genes from OXPHOS pathway were the key finding (cytochrome c oxidase subunits).
- **Reimand et al. (2019)** PMID 30664679 -- Recommends examining leading-edge genes as part of the interpretation workflow.

### Feasibility for this dataset

**High**. Leading edges are already stored in `fgsea_tstat_all_v2.csv` and `fgsea_gobp_v2.csv` (the `leadingEdge` column contains semicolon-delimited gene lists for every term/contrast). Implementation would be:
1. For each a priori term (or quadrant-interesting term), extract leading edges from both contrasts.
2. Compute overlap (Jaccard index) and present as a heatmap or upset plot.
3. For discordant pathways, highlight genes that appear in opposite-direction leading edges.

### Pitfalls

- **Leading edge is direction-dependent**: For a pathway with NES > 0, the leading edge is the top-ranked genes in the pathway. For NES < 0, it is the bottom-ranked genes. When comparing contrasts with opposite NES signs, the leading edges are drawn from opposite tails of the rank distribution, so low overlap is expected by construction.
- **Set size confound**: Larger pathways have larger leading edges, inflating overlap by chance. Any comparison should use a normalized metric (Jaccard, or overlap coefficient).
- **Not a formal test**: Leading-edge intersection is descriptive, not inferential. It identifies candidates but does not provide p-values for the overlap. A hypergeometric test on the leading-edge overlap could formalize this, but the multiple-testing burden would be substantial.
- **Presentation complexity**: A full leading-edge intersection analysis across all 69 terms and 5 contrasts generates a large matrix. This is inherently supplementary material.

### Verdict: RECOMMENDED as a complementary interpretive tool

This is mechanistically the most informative approach -- it names the actual proteins driving concordance or discordance within a pathway. It should be presented as a supplementary table (for all terms) and selectively in the text (for key discordant pathways). It does not replace the NES scatter but deepens its interpretation.

---

## Approach 4: Pathway-Level RRHO

### Concept

Apply the RRHO2 framework at the pathway level instead of the protein level. Rank pathways by NES in each contrast and compute the hypergeometric overlap matrix.

### Methodological basis

RRHO (Plaisier et al. 2010, PMID 20660011; Cahill et al. 2018, PMID 29942049) was designed for genome-scale ranked lists (thousands of genes). Its statistical power depends on the number of items ranked.

### Feasibility for this dataset

**Very low with a priori collection (69 terms)**. RRHO requires a large ranked list for the hypergeometric test to have meaningful resolution. With 69 terms, the grid would be ~35x35 at most, and each cell would test overlap of very small sets. The resulting heatmap would be dominated by noise.

**Potentially feasible with full GO:BP (~998 terms) or all databases (~2,170 terms)**. This would provide a grid of ~100x100 or ~200x200, approaching the resolution where RRHO can detect meaningful hotspots. However:
- Pathway-level ranks (NES) are noisier than gene-level ranks (t-statistics) because NES is itself a summary statistic.
- Many GO:BP terms are highly overlapping in gene membership, violating the RRHO assumption that ranked items are independent.
- The biological interpretation of "pathway-level RRHO hotspots" is unclear -- what does it mean that a cluster of pathways all rank similarly in two contrasts?

### Supporting literature

No published precedent for pathway-level RRHO was found in the literature search. RRHO has been applied at the gene level (Plaisier et al. 2010), the transcript level (Cahill et al. 2018), and in cross-platform comparisons (Lardone et al. 2016, PMID 26883106), but always on individual molecular features, not on pathway summary statistics.

### Pitfalls

- **Independence violation**: GO:BP terms share genes extensively. Two pathways with 80% gene overlap will have nearly identical NES values, creating artificial concordance in the RRHO.
- **NES is a derived statistic**: RRHO on NES values applies a nonparametric rank test to a statistic that is already a nonparametric enrichment score. This double-ranking may lose information.
- **Interpretability**: The output of pathway-level RRHO is a heatmap of "-log10(p) of overlap between top-N pathways in contrast A and top-M pathways in contrast B." This is harder to interpret biologically than the gene-level version.
- **No software support**: There is no published tool or validated implementation for pathway-level RRHO.

### Verdict: NOT RECOMMENDED

The theoretical foundations are weak (independence violation, derived statistic), there is no literature precedent, and the a priori collection is too small. The existing gene-level RRHO already captures the concordance/discordance structure at higher resolution. If pathway-level cross-contrast comparison is desired, the NES scatter (already implemented) is a more transparent and interpretable visualization.

---

## Approach 5: Reactome Hierarchy Drill-Down

### Concept

Reactome has a natural hierarchy: 27 top-level pathways -> sub-pathways -> reactions. Use the top-level Reactome categories as the "broad view" (analogous to GO Slim), then drill into specific sub-pathways within interesting top-level categories. This preserves the a priori structure at the top level while allowing systematic zoom.

### Reactome hierarchy structure (verified via API)

Example: Metabolism (R-HSA-1430728, top-level) contains:
- Metabolism of carbohydrates and carbohydrate derivatives
- Metabolism of lipids
- Aerobic respiration and respiratory electron transport
- Metabolism of amino acids and derivatives
- Metabolism of vitamins and cofactors
- ...10 total children

And "Aerobic respiration and respiratory electron transport" (R-HSA-1428517) further contains:
- Pyruvate metabolism
- TCA cycle
- Respiratory electron transport
- Formation of ATP by chemiosmotic coupling
- Mitochondrial uncoupling
- ...6 children

### Methodological basis

Using curated pathway hierarchies for multi-level analysis is well-established in the Reactome ecosystem (Rothfels et al. 2023, PMID 37012728 -- "Using the Reactome Database," Current Protocols). The Paley et al. (2017) Omics Dashboard (PMID 28968752) explicitly implements hierarchical browsing of pathway enrichment from high-level categories down to individual reactions. Reactome's own analysis portal provides hierarchical visualization of enrichment results.

### Feasibility for this dataset

**Moderate-to-high**. Reactome terms are already in `fgsea_tstat_all_v2.csv` (~432 Reactome terms tested per contrast, 135-217 significant depending on contrast). The hierarchy can be obtained via the Reactome API or msigdbr metadata. Implementation would be:

1. Map each Reactome term to its top-level parent (27 categories).
2. Run fGSEA at both the top-level parent level (aggregate gene sets) and the sub-pathway level (individual Reactome terms).
3. For top-level categories that fall in interesting NES scatter quadrants, present the sub-pathway breakdown.
4. This creates a natural "zoom" from 27 top-level categories (comparable to the 27 GO Slim terms already used) to hundreds of specific sub-pathways.

### Pitfalls

- **Gene set construction at top level**: Reactome top-level categories are not pre-built as MSigDB gene sets. They would need to be constructed by aggregating genes from all child pathways. This is doable but introduces a gene set that was not in the original fGSEA run.
- **Overlap with GO Slim**: Reactome top-level categories overlap conceptually with GO Slim terms (e.g., "Metabolism" vs "Metabolic Process"). Running both creates redundancy.
- **Reactome bias**: Reactome has deep coverage of metabolism, signaling, and disease pathways but weaker coverage of structural biology, muscle-specific processes, and aging biology. For a skeletal muscle proteomics study, this may miss key biology.
- **Hierarchy depth varies**: Some top-level categories have 2-3 levels of nesting; others have 5-6. The "drill-down" is not uniform.

### Verdict: USEFUL as a supplementary data source, but not as the primary zoom mechanism

Reactome hierarchy provides a natural multi-level structure, but GO Slim + GO:BP already achieves the same conceptual zoom (GO Slim = broad, GO:BP = fine) with better consistency and wider biological coverage. The Reactome hierarchy is most valuable for specific sub-pathway dissection of metabolic findings (e.g., drilling into which specific branch of lipid metabolism is discordant).

---

## Approach 6: Network-Based Approaches (STRING / PPI)

### Concept

Take the proteins within each NES scatter quadrant (or the leading-edge proteins from quadrant-specific pathways) and analyze them as a PPI network using STRING or similar. Identify densely connected subnetworks (modules) that may represent functional complexes driving the quadrant behavior.

### Methodological basis

STRING (Szklarczyk et al. 2023, PMID 36370105) provides both PPI network construction and built-in enrichment analysis. The typical workflow is: (1) input a protein list, (2) identify connected components and clusters, (3) run enrichment on clusters. This is standard in proteomics papers.

### Relationship to existing WGCNA

The pipeline already includes WGCNA (WGCNA_F5, WGCNA_F6), which performs co-expression network analysis on the full proteome. WGCNA modules capture co-expression patterns across all samples, not quadrant-specific patterns. However, there is a natural connection: proteins in a specific NES scatter quadrant may be enriched in specific WGCNA modules, providing convergent evidence.

### Feasibility for this dataset

**Moderate**. The protein lists exist (from leading-edge extraction or DEP classification). STRING analysis is straightforward. But:
- The input sets need to be well-defined (which proteins represent "the discordant quadrant"?).
- With ~2,138 total proteins, the full network is manageable.
- Per-quadrant subnetworks may be small (10-50 proteins), limiting network analysis power.

### Pitfalls

- **Circular enrichment**: If we define proteins by pathway membership (leading edge of a pathway) and then do network enrichment, we will recover the original pathway. The network adds value only if it identifies NEW connections not captured by the pathway.
- **STRING score thresholds**: Results are sensitive to the minimum interaction score. High confidence (>0.7) may fragment the network; medium confidence (>0.4) may introduce false connections.
- **Presentation burden**: PPI networks are visually complex and journal-unfriendly in small panels. They work best as supplementary figures.
- **WGCNA overlap**: If the network modules recapitulate WGCNA modules, the analysis is redundant. If they disagree, interpreting the disagreement is difficult.

### Verdict: LOW PRIORITY; use only if a specific mechanistic question requires it

Network analysis is most valuable when a specific question cannot be answered by pathway enrichment alone (e.g., "do the discordant proteins form a physical complex?"). For the broad question of "what processes drive discordance?", approaches 1 and 3 are more informative and less assumption-laden.

---

## Approach 7: Supplementary Panel / Companion Table

### Concept

Keep the main figure clean (69 a priori terms) but add supplementary material showing finer-grained results. This is a presentation strategy, not an analytical approach -- it can be combined with any of the above.

### Standard practice in high-impact proteomics papers

Examining the presentation practices in recent exercise/aging proteomics studies:

1. **Main figure**: Shows a curated, interpretable visualization (NES scatter, volcano, heatmap) with a manageable number of terms (typically 10-50 labeled pathways).
2. **Supplementary table**: Provides the full enrichment results (all tested terms, all contrasts, NES, padj, leading edge). This is the "complete record" that readers can interrogate.
3. **Supplementary figure**: May include a focused panel on a specific biological theme (e.g., "mitochondrial sub-pathways in the reversal quadrant") or an EnrichmentMap network visualization.

This two-tier approach (curated main + complete supplement) is explicitly recommended by Reimand et al. (2019, PMID 30664679) and is the standard in Cell, Nature, and PNAS-level proteomics papers.

### Recommended implementation

**Tier 1 (Main figure)**: Keep Panels G and E exactly as they are. The 69 a priori terms with significance classification and Spearman rho are the primary result.

**Tier 2 (Supplementary table)**: Export the full GO:BP + Reactome NES scatter data for each comparison. Format:

| pathway | database | set_size | NES_contrast_A | NES_contrast_B | padj_A | padj_B | quadrant | leading_edge_A | leading_edge_B | LE_overlap_jaccard |
|---------|----------|----------|----------------|----------------|--------|--------|----------|----------------|----------------|--------------------|

This allows reviewers and readers to filter to any quadrant and any significance threshold.

**Tier 3 (Optional supplementary panel)**: A focused drill-down panel for the most biologically interesting quadrant. Candidates:
- For Panel G: The discordant quadrant (Training_Young up / Training_Old down) -- which specific metabolic or structural pathways show age-dependent training response?
- For Panel E: The reversed quadrant (Aging up / Training_Old down) -- which specific catabolic or degradative pathways are reversed by training?

Format: A bar plot or dot plot of the top 15-20 deduplicated GO:BP terms within that quadrant, colored by database, sized by gene set size.

### Verdict: ESSENTIAL complement to any analytical approach

Regardless of which analytical approach is chosen, the supplementary table should be provided. This is a journal expectation, not an optional extra.

---

## Integrated Recommendation

Based on the analysis above, the recommended strategy is a **three-tier approach**:

### Tier 1: Main Figure (No Change)
Panels G and E remain as they are: 69 a priori terms (GO Slim + Hallmark), NES scatter with significance classification, Spearman rho, and quadrant counts. This is the defensible, pre-registered pathway-level summary.

### Tier 2: Supplementary Table (Approach 1 + 3 Combined)
A comprehensive supplementary CSV/Excel with:
- All ~2,170 tested pathways in wide format (NES and padj for each contrast pair)
- Quadrant classification based on NES signs
- Leading-edge genes for each term in each contrast
- Leading-edge overlap (Jaccard index) between the two contrasts
- Jaccard deduplication flag (survives J=0.35, J=0.50)

This table is the complete record. It supports Approach 1 (any reader can filter to a quadrant + significance) and Approach 3 (leading-edge intersection is pre-computed).

### Tier 3: Supplementary Panel (Approach 1, Focused)
One supplementary figure panel per NES scatter, showing the GO:BP drill-down for the most biologically interesting quadrant:

**For Panel G (Concordance)**: Focus on the discordant quadrant(s). Take all GO:BP terms with NES_Training_Young and NES_Training_Old in opposite directions AND padj < 0.05 in at least one contrast. Deduplicate (Jaccard 0.35). Present as a horizontal bar/dot plot of NES values (dual bars: Young vs Old), ordered by the magnitude of discordance (|NES_Young - NES_Old|).

**For Panel E (Reversal)**: Focus on the reversed quadrant(s). Take all GO:BP terms with NES_Aging and NES_Training_Old in opposite directions AND padj < 0.05 in at least one contrast. Same visualization.

### What NOT to do
- Do not create a pathway-level RRHO (Approach 4). No precedent, weak foundations.
- Do not add a full PPI network panel (Approach 6) unless a specific mechanistic question demands it.
- Do not run Reactome hierarchy as a parallel a priori collection (Approach 5) -- it overlaps with GO Slim and adds complexity without proportional insight. Use it selectively if a specific metabolic finding warrants sub-pathway dissection.
- Do not run ORA on DEP-classified quadrant proteins (Approach 2) -- the existing RRHO ORA already does this more powerfully with continuous ranks.

---

## Key Literature References

| Reference | PMID | Relevance |
|-----------|------|-----------|
| Subramanian et al. (2005) PNAS | 16199517 | GSEA, leading-edge analysis, post-hoc refinement |
| Mootha et al. (2003) Nat Genet | 12808457 | Original GSEA application, leading-edge interpretation |
| Reimand et al. (2019) Nat Protocols | 30664679 | Best practices for pathway enrichment, multi-tier presentation |
| Savage et al. (2019) MCP | 31142576 | Graph algorithms for condensing enrichment results |
| Supek et al. (2011) PLoS ONE (REVIGO) | 21789182 | GO term semantic similarity and redundancy reduction |
| Merico et al. (2010) PLoS ONE (EnrichmentMap) | 21085593 | Network visualization of enrichment results |
| Korotkevich et al. (2021) bioRxiv (fgsea) | -- | fgsea, collapsePathways, multilevel p-value estimation |
| Hanzelmann et al. (2013) BMC Bioinformatics (GSVA) | 23323831 | Single-sample pathway scores (alternative to NES comparison) |
| Cahill et al. (2018) Sci Rep (RRHO2) | 29942049 | Improved RRHO with stratified quadrant analysis |
| Plaisier et al. (2010) Nucleic Acids Res (RRHO) | 20660011 | Original RRHO method |
| Szklarczyk et al. (2023) NAR (STRING) | 36370105 | PPI network + enrichment analysis |
| Huang et al. (2009) NAR (DAVID) | 19033363 | Enrichment tool survey, ORA methodology |
| Melov et al. (2007) PLoS ONE | 17520024 | Exercise reversal of aging, precedent for reversal analysis |
| Rothfels et al. (2023) Curr Protocols (Reactome) | 37012728 | Reactome hierarchy and pathway analysis |
| Paley et al. (2017) NAR (Omics Dashboard) | 28968752 | Hierarchical pathway browsing for expression data |

---

## Implementation Priority

1. **Supplementary table** (Tier 2): Immediate. Requires only data reshaping of existing fGSEA outputs. No new analyses.
2. **Leading-edge overlap column** (Approach 3): Add to the supplementary table. Requires parsing the existing `leadingEdge` column and computing Jaccard per term.
3. **Focused supplementary panel** (Tier 3): After reviewing the supplementary table to identify the most informative quadrant.
4. **Reactome sub-pathway drill-down** (Approach 5): Only if a specific metabolic finding warrants it, after reviewing Tier 2 results.

---

*Report generated 2026-03-19. No code was written. This is a research and brainstorming document only.*
