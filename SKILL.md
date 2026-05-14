# SKILL.md — Alpha Diversity Pipeline: Data Interpretation Guide

This document is written for an AI analysis agent that will interpret the outputs
of the `alpha-diversity` Nextflow pipeline. Read it carefully before analysing
any pipeline output.

---

## What this pipeline does

The pipeline takes a microbiome feature table (counts of ASVs or OTUs per sample)
and computes alpha diversity — a measure of within-sample biodiversity — for each
sample. It then runs pairwise statistical tests to compare diversity between groups.
All outputs are raw data CSVs. Plotting is left to the user.

---

## Step 1 — Read the manifest before anything else

Every pipeline run writes `pipeline_manifest.json` to the output directory.
Read it first. It records every parameter that was used, the exact command line,
and timestamps. Use this to write the Methods section.

Key fields:
| Field | What it means |
|---|---|
| `params.normalization` | `clr` or `rarefaction` — changes how metrics are defined |
| `params.min_library_size` | Samples below this depth were removed before analysis |
| `params.taxon_rank` | The level at which features were collapsed (e.g. `Feature`, `Genus`) |
| `params.alpha_metrics` | Which diversity metrics were computed |
| `params.test_method` | Statistical test used (`anova`, `kruskal`, `auto`, `none`) |
| `params.p_adjust_method` | Multiple-testing correction method (e.g. `BH`) |
| `params.group` | Metadata column(s) used as grouping variable |

---

## Output files — what each contains

### `Diversity_{taxon_rank}_{label}.csv` — Wide diversity table

One row per sample. The definitive result file.

| Column | Type | Meaning |
|---|---|---|
| `sample` | string | Sample ID |
| `Richness` | float | Observed species count (CLR) or rarefied richness |
| `Shannon` | float | Shannon entropy H = -sum(p·log(p)), 0 to log(S) |
| `GiniSimpson` | float | 1 - sum(p²), range 0–1 (higher = more diverse) |
| `InvSimpson` | float | 1 / sum(p²), range 1 to S |
| `FisherAlpha` | float | Parametric richness estimator; higher = richer |
| `Chao1` | float | Non-parametric richness estimator; higher = richer |
| `PielouEvenness` | float | H / log(observed_richness), range 0–1 (higher = more even) |
| `FaithsPD` | float | Sum of branch lengths in minimum spanning subtree (requires tree) |
| `Groups` | string | Group label for this sample (from `--group` parameter) |
| `normalization` | string | Either `clr` or `rarefaction` — records how metrics were computed |
| `Type` | string | Optional point-style label (only present if `--type` was used) |

Not all metric columns will be present — only those in `params.alpha_metrics`.

---

### `Diversity_long_{taxon_rank}_{label}.csv` — Long diversity table

Same data as the wide table, reshaped to long format. Easier for plotting.

| Column | Type | Meaning |
|---|---|---|
| `sample` | string | Sample ID |
| `value` | float | Diversity value |
| `measure` | string | Metric name (e.g. `Shannon`) |
| `Groups` | string | Group label |
| `Type` | string | Optional point-style label |

---

### `Diversity_pairwise_{taxon_rank}_{label}.csv` — Statistical comparisons

One row per pairwise group comparison per metric.

| Column | Type | Meaning |
|---|---|---|
| `measure` | string | Metric name |
| `group1` | string | First group |
| `group2` | string | Second group |
| `estimate` | float | Mean difference (ANOVA/TukeyHSD) or NA (Kruskal) |
| `conf_low` | float | Lower 95% CI of difference (TukeyHSD) or NA |
| `conf_high` | float | Upper 95% CI of difference (TukeyHSD) or NA |
| `pvalue` | float | Raw p-value |
| `padj` | float | Adjusted p-value — **always use this for reporting** |
| `significance` | string | `***` (≤0.001), `**` (≤0.01), `*` (≤0.05), `""` (ns) |
| `test_method` | string | Requested method (may be `auto`) |
| `actual_test_method` | string | Method actually used (important when `test_method=auto`) |
| `levene_pvalue` | float | Levene's test p-value (only when `test_method=auto`) |

**Important**: For ANOVA mode, TukeyHSD is used. `padj` equals `p adj` from
TukeyHSD directly (FWER already controlled). For Kruskal mode, `padj` is the
BH/bonferroni/holm-adjusted Mann-Whitney U p-value.

---

### `Diversity_n_{taxon_rank}_{label}.csv` — Sample counts per group

| Column | Type | Meaning |
|---|---|---|
| `Groups` | string | Group label |
| `n` | integer | Number of samples in this group |

Always check this file. Groups with n < 5 make statistical comparisons unreliable.

---

### `Diversity_twoway_{taxon_rank}_{label}.csv` — Two-way ANOVA (optional)

Only present when `--type` was used. Reports the ANOVA table for `Groups * Type`
interaction. **This is exploratory only** — no post-hoc for the interaction term
is computed.

| Column | Meaning |
|---|---|
| `measure` | Metric name |
| `term` | ANOVA term: `Groups`, `Type`, `Groups:Type`, `Residuals` |
| `df` | Degrees of freedom |
| `sum_sq`, `mean_sq` | Sums/means of squares |
| `f_value` | F statistic |
| `pvalue` | Raw p-value for this term |

---

### `rarefaction_data_{label}.csv` — Rarefaction curve data

Used to assess whether sequencing depth was sufficient. Plot this to see if
richness curves plateau before `min_library_size`.

| Column | Meaning |
|---|---|
| `sample` | Sample ID |
| `depth` | Sequencing depth at this point on the curve |
| `richness` | Expected species richness at this depth |
| `Groups` | Group label |
| `min_library_size` | Reference depth (the filter threshold) |

A curve that has not plateaued at `min_library_size` suggests insufficient depth.

---

### `Diversity_correlation_{label}.csv` — Inter-metric Spearman correlations

| Column | Meaning |
|---|---|
| `metric1`, `metric2` | Pair of metrics being compared |
| `rho` | Spearman correlation coefficient (−1 to +1) |
| `pvalue` | Raw p-value for the correlation |

High correlation (|rho| > 0.9) between two metrics means they are redundant —
report only one of them.

---

### `Diversity_sensitivity_{taxon_rank}_{label}.csv` — Sensitivity analysis

Only present when `--sensitivity_analysis true --normalization rarefaction`.
Tests whether metric values change substantially at different rarefaction depths.

| Column | Meaning |
|---|---|
| `sample` | Sample ID |
| `measure` | Metric name |
| `value` | Diversity value at this depth |
| `Groups` | Group label |
| `depth_fraction` | Fraction of min_library_size used (e.g. 0.5) |
| `rarefaction_depth` | Actual read count used for rarefaction |

If metric rankings or group differences change substantially across depths,
the results are not robust to the choice of rarefaction depth — note this.

---

### `alpha_diversity_{label}.parquet` — Merged Parquet (optional)

Only present when `--merge_parquet true`. Contains all three CSVs (wide, long,
pairwise) merged with `analysis` and `table` metadata columns. Filter by
`table` column: `diversity_wide`, `diversity_long`, `diversity_pairwise`.

---

### `r_session_info.txt` — R package versions

Lists all R package versions used. Include this in supplementary materials.

---

## Diversity metric interpretation guide

### Richness
- **What it measures**: Number of distinct taxa present in a sample
- **CLR mode**: Observed count of non-zero features
- **Rarefaction mode**: Expected number of species at `min_library_size` reads
- **Range**: 0 to total number of features
- **High value**: Many distinct taxa — high biodiversity
- **Caveat**: Sensitive to sequencing depth (more reads = more rare taxa detected)

### Shannon (H)
- **What it measures**: Entropy — accounts for both richness and evenness
- **Formula**: H = −Σ(pᵢ · ln(pᵢ)) where pᵢ is the relative abundance of taxon i
- **Range**: 0 (one taxon completely dominant) to ln(S) (perfectly even)
- **High value**: Many taxa present at similar abundances
- **Caveat**: More sensitive to rare species than Simpson

### GiniSimpson (= 1 − D)
- **What it measures**: Probability that two randomly chosen individuals belong to different taxa
- **Formula**: 1 − Σ(pᵢ²)
- **Range**: 0 (one taxon dominates) to approaching 1 (many equal taxa)
- **Important**: This is the Gini-Simpson index, NOT the raw Simpson index D = Σ(pᵢ²).
  D and GiniSimpson move in opposite directions — report as GiniSimpson.
- **High value**: High diversity, low dominance

### InvSimpson (= 1/D)
- **What it measures**: Effective number of equally abundant species (Hill number N₂)
- **Formula**: 1/Σ(pᵢ²)
- **Range**: 1 (one taxon) to S (perfectly even)
- **High value**: Equivalent to many equally abundant taxa
- **Caveat**: Less sensitive to rare taxa than Shannon

### FisherAlpha
- **What it measures**: Parametric richness estimator based on log-series distribution
- **Range**: Positive real number, no upper bound
- **High value**: Higher richness
- **Caveat**: Assumes log-series species abundance distribution; may not suit all datasets

### Chao1
- **What it measures**: Non-parametric estimator of true species richness
- **Based on**: Frequency of singletons and doubletons in raw counts
- **High value**: Estimated true richness is high
- **Caveat**: Requires raw integer counts; biased when sample size is very small

### PielouEvenness (J)
- **What it measures**: How evenly abundances are distributed among taxa
- **Formula**: J = Shannon / ln(observed_richness)
- **Range**: 0 (one taxon dominates) to 1 (perfectly even)
- **High value**: Abundances are distributed evenly — no dominant taxon

### FaithsPD
- **What it measures**: Sum of branch lengths connecting all taxa in a sample on the phylogenetic tree
- **Range**: Positive real, in tree branch-length units
- **High value**: Phylogenetically diverse community
- **Requires**: `--tree_file` must have been provided; feature IDs must match tree tip labels
- **Caveat**: Values depend on the tree and its branch length units

---

## Statistical output guide

### Which p-value to report
Always use `padj` (adjusted p-value), not `pvalue`. `padj` accounts for the fact
that you are running multiple comparisons simultaneously.

### Significance symbols
| Symbol | Meaning |
|---|---|
| `***` | padj ≤ 0.001 |
| `**` | padj ≤ 0.01 |
| `*` | padj ≤ 0.05 |
| `""` (empty) | padj > 0.05 — not significant |

### ANOVA vs Kruskal-Wallis
- **ANOVA** (`actual_test_method = "anova"`): parametric, assumes normality; post-hoc
  is Tukey HSD (controls family-wise error rate). `estimate` and confidence intervals
  are meaningful — report them.
- **Kruskal-Wallis** (`actual_test_method = "kruskal"`): non-parametric; post-hoc is
  Mann-Whitney U. `estimate`, `conf_low`, `conf_high` are NA for this test.
- **auto mode**: Check `levene_pvalue`. If < 0.05, heterogeneous variance was detected
  and Kruskal was used. Report both `test_method=auto` and `actual_test_method`.

---

## How to write the Methods section

Use `pipeline_manifest.json` to fill in the blanks:

> Alpha diversity was computed using the alpha-diversity Nextflow pipeline
> (v{params.version}). Samples with fewer than {params.min_library_size} reads
> were excluded. Features were {collapsed to {params.taxon_rank} level / retained
> at feature level}. The following diversity metrics were computed using the
> {params.normalization} normalisation strategy: {params.alpha_metrics}.
> Pairwise group comparisons were performed using {test method} with
> {params.p_adjust_method} correction for multiple testing. R package versions
> are reported in r_session_info.txt.

For CLR normalization, add:
> Diversity metrics were computed on relative abundances without rarefaction,
> consistent with a compositional data analysis approach (Gloor et al. 2017).

For rarefaction normalization, add:
> Each sample was rarefied to {params.min_library_size} reads prior to metric
> computation to equalise sequencing depth across samples.

---

## How to write the Results section

1. Report n per group from `Diversity_n_*.csv` — note any small groups
2. For each metric, report the median (IQR) per group
3. Report significant pairwise comparisons from `Diversity_pairwise_*.csv`
   using `padj` and `significance` columns
4. If sensitivity analysis was run, state that results were consistent across
   rarefaction depths (or note inconsistencies)
5. Note any highly correlated metrics from `Diversity_correlation_*.csv`
   (|rho| > 0.9) and explain why only one was reported

---

## Red flags — always check these

| Flag | What to look for | What to say |
|---|---|---|
| **Small groups** | Any row in `Diversity_n_*.csv` with `n < 5` | "Statistical comparisons for groups with n < 5 should be interpreted with caution" |
| **Null result** | All `padj > 0.05` | A null result is a valid finding — state it clearly |
| **Depth inadequacy** | Rarefaction curves in `rarefaction_data_*.csv` not plateauing | "Sequencing depth may have been insufficient to capture full diversity" |
| **Sensitivity failure** | Large shifts across `depth_fraction` values | "Metric {X} was sensitive to rarefaction depth and should be interpreted cautiously" |
| **Metric redundancy** | `rho > 0.9` between two metrics | Report only the most interpretable metric |
| **Inconsistent direction** | Shannon up but GiniSimpson down between same groups | Check for dominance shifts; explain the biological interpretation |
| **GiniSimpson confusion** | Column is named `GiniSimpson` | Never call it "Simpson's index" — it is 1−D, not D |
