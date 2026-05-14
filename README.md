# alpha-diversity

A Nextflow pipeline for computing alpha diversity metrics from microbiome feature tables, with optional statistical comparisons between groups.

## Introduction

This pipeline accepts a feature table (BIOM, TSV, or GTDB format) and a sample metadata CSV. It collapses features to a chosen taxonomic level, computes a suite of alpha diversity indices per sample, and performs pairwise statistical tests between groups. All outputs are raw data CSV files — plotting is left to the user.

See [SKILL.md](SKILL.md) for a complete guide on interpreting the outputs, aimed at both human analysts and AI data agents.

## Quick start

```bash
nextflow run main.nf \
  --feature_table /path/to/feature_table.biom \
  --meta_table    /path/to/meta_table.csv \
  --group         Treatment \
  --label         my_analysis
```

## Parameters

### Input / Output

| Parameter | Default | Description |
|---|---|---|
| `--feature_table` | *(required)* | Path to feature table (BIOM, TSV, or GTDB format) |
| `--meta_table` | *(required)* | Path to sample metadata CSV (first column = sample IDs) |
| `--taxonomy_table` | `""` | Taxonomy TSV (required for `tsv`/`gtdb` input formats) |
| `--tree_file` | `""` | Newick tree file (required when `FaithsPD` is in `--alpha_metrics`) |
| `--input_format` | `biom` | `biom` \| `tsv` \| `gtdb` |
| `--output_dir` | `results/` | Directory for output files |

### Filtering

| Parameter | Default | Description |
|---|---|---|
| `--min_library_size` | `5000` | Minimum per-sample read depth; samples below this are dropped. Also used as the rarefaction depth under `--normalization rarefaction`. |
| `--exclude_column` | `""` | Metadata column used to identify samples for exclusion |
| `--exclude_values` | `""` | Comma-separated values in `exclude_column` to remove |

### Grouping

| Parameter | Default | Description |
|---|---|---|
| `--group` | `""` | One metadata column, or comma-separated columns pasted together as the group label |
| `--type` | `""` | Optional point-style label (single or comma-separated columns) |

### Alpha diversity

| Parameter | Default | Description |
|---|---|---|
| `--taxon_rank` | `Feature` | Taxonomic level for feature collation (`Feature` \| `Genus` \| `Family` \| `Order` \| `Class` \| `Phylum`) |
| `--label` | `analysis` | Label prepended to output file names |
| `--normalization` | `clr` | `clr` (proportion-based, recommended) \| `rarefaction` (subsample to `min_library_size`) |
| `--alpha_metrics` | `Richness,Shannon,GiniSimpson,FisherAlpha,PielouEvenness` | Comma-separated list; see [Available metrics](#available-metrics) |
| `--test_method` | `anova` | `anova` \| `kruskal` \| `none` \| `auto` |
| `--p_adjust_method` | `BH` | `BH` \| `bonferroni` \| `holm` \| `fdr` \| `none` |

### Output options

| Parameter | Default | Description |
|---|---|---|
| `--merge_parquet` | `false` | Merge all output CSVs into a single Parquet file |
| `--sensitivity_analysis` | `false` | Run sensitivity analysis across rarefaction depths (only active when `--normalization rarefaction`) |
| `--sensitivity_depths` | `0.5,0.75,1.0` | Depth fractions for sensitivity analysis |

## Available metrics

| Metric | Description | Normalization |
|---|---|---|
| `Richness` | Observed species count (CLR) or rarefied richness (rarefaction) | Both |
| `Shannon` | Entropy H = −Σ(p·ln p) | Both |
| `GiniSimpson` | 1 − Σ(p²) — diversity index; **not** the raw Simpson D | Both |
| `InvSimpson` | 1 / Σ(p²) — Hill number N₂ | Both |
| `FisherAlpha` | Parametric richness estimator | Both |
| `Chao1` | Non-parametric richness estimator (requires integer counts) | Both |
| `PielouEvenness` | Shannon / ln(observed richness) | Both |
| `FaithsPD` | Faith's phylogenetic diversity (requires `--tree_file`) | Both |

**Normalization note**: `clr` uses proportions internally for all abundance-weighted
metrics (Shannon, GiniSimpson, InvSimpson, PielouEvenness) and raw counts for
estimators (Chao1, FisherAlpha). No subsampling is performed. This is consistent
with a compositional data analysis framework (Gloor et al. 2017). Use
`--normalization rarefaction` to reproduce the traditional rarefaction approach.

**Statistical testing**: ANOVA mode uses Tukey HSD for pairwise post-hoc tests
(family-wise error rate controlled). Kruskal-Wallis mode uses pairwise
Mann-Whitney U with `--p_adjust_method` correction. The `auto` option runs
Levene's test per metric and selects ANOVA or Kruskal-Wallis accordingly.

## Outputs

All files are written to `--output_dir`.

| File | Always? | Description |
|---|---|---|
| `Diversity_{taxon_rank}_{label}.csv` | Yes | Wide-format diversity table (one row per sample, one column per metric) |
| `Diversity_long_{taxon_rank}_{label}.csv` | Yes | Long-format version with `measure` and `value` columns |
| `Diversity_pairwise_{taxon_rank}_{label}.csv` | Yes | Pairwise comparisons: test statistic, p-value, padj, significance, CI |
| `Diversity_n_{taxon_rank}_{label}.csv` | Yes | Per-group sample counts |
| `Diversity_twoway_{taxon_rank}_{label}.csv` | When `--type` used | Two-way ANOVA table (Groups × Type, exploratory) |
| `rarefaction_data_{label}.csv` | Yes | Rarefaction curve data (sample, depth, richness, Groups) |
| `Diversity_correlation_{label}.csv` | Yes | Spearman correlations between metrics (metric1, metric2, rho, pvalue) |
| `Diversity_sensitivity_{taxon_rank}_{label}.csv` | When `--sensitivity_analysis true` and `--normalization rarefaction` | Metrics at multiple rarefaction depths |
| `alpha_diversity_{label}.parquet` | When `--merge_parquet true` | All three main CSVs merged |
| `pipeline_manifest.json` | Yes | Machine-readable record of all parameters, timestamps, and run metadata |
| `r_session_info.txt` | Yes | R package versions used |
| `pipeline_report.html` | Yes | Nextflow HTML execution report |
| `pipeline_trace.txt` | Yes | Per-process resource usage trace |
| `pipeline_timeline.html` | Yes | Pipeline execution timeline |

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 23.04
- [conda](https://docs.conda.io/) or [mamba](https://mamba.readthedocs.io/) (default — environment built automatically from `environment.yml`)
- **or** Docker with `-profile docker`
- **or** Singularity with `-profile singularity`

Core R dependencies (resolved automatically by conda/Docker): `vegan`, `phyloseq`,
`picante`, `ape`, `car`, `arrow`, `data.table`, `optparse`.

## Docker

```bash
# Build
docker build -t reubenduncan/ecologyflow-alpha:latest .

# Run with Docker profile
nextflow run main.nf -profile docker \
  --feature_table /path/to/table.biom \
  --meta_table    /path/to/meta.csv
```

## Example scripts

The `examples/` directory contains scripts that demonstrate how to visualise
the pipeline's CSV outputs. These are **not called by the pipeline**.

| Script | Description |
|---|---|
| `examples/plot_alpha_diversity.R` | Publication-quality boxplots with significance brackets (CLD, full, or selective) |
| `examples/plot_rarefaction_curves.R` | Rarefaction curve plot from `rarefaction_data_*.csv` |
| `examples/plot_metric_correlation.R` | Spearman correlation heatmap from `Diversity_correlation_*.csv` |
| `examples/standalone_alpha_diversity.R` | Legacy standalone script (reference only) |

## Normalisation: CLR vs rarefaction

By default (`--normalization clr`) the pipeline uses a compositionally sound
approach: no subsampling is performed, abundance-weighted metrics are computed on
relative abundances (scale-invariant), and richness is the observed species count.
This is consistent with the compositional data analysis framework (Gloor et al.
2017) and avoids the well-documented pitfalls of rarefaction (McMurdie & Holmes
2014).

Use `--normalization rarefaction` to reproduce the traditional rarefaction approach.
Note that `--sensitivity_analysis` only functions under rarefaction normalization.

## Workflow of workflows

The data-import and parquet-merge logic is designed to be factored out as shared
Nextflow subworkflows (`subworkflows/local/data_import` and
`subworkflows/local/merge_parquet`) reusable across beta-diversity,
differential-abundance, and other ecology pipelines. This refactor is planned for
the next architectural revision.
