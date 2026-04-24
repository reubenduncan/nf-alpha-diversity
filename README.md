# alpha-diversity

A Nextflow pipeline for computing alpha diversity metrics from microbiome feature tables, with optional statistical comparisons between groups.

## Introduction

This pipeline accepts a feature table (BIOM, TSV, or GTDB format) and a sample metadata CSV, collapses features to a chosen taxonomic level, and computes a suite of alpha diversity indices per sample. It then performs pairwise statistical tests (ANOVA or Kruskal-Wallis) between groups and writes results to CSV. An optional `--merge_parquet` flag consolidates all outputs into a single Parquet file.

## Quick start

```bash
nextflow run main.nf \
  --feature_table /path/to/feature_table.biom \
  --meta_table    /path/to/meta_table.csv \
  --groups_column Treatment \
  --label         my_analysis
```

## Parameters

### Input / Output

| Parameter | Default | Description |
|---|---|---|
| `--feature_table` | *(required)* | Path to feature table (BIOM, TSV, or GTDB format) |
| `--meta_table` | *(required)* | Path to sample metadata CSV (first column = sample IDs) |
| `--taxonomy_table` | `""` | Taxonomy TSV (required for `tsv`/`gtdb` input formats) |
| `--tree_file` | `""` | Newick tree file (accepted for API consistency; not used) |
| `--input_format` | `biom` | `biom` \| `tsv` \| `gtdb` |
| `--output_dir` | `results/` | Directory for output files |
| `--scripts_dir` | `/opt/ecology-scripts` | Path to R scripts (override for local runs) |

### Filtering

| Parameter | Default | Description |
|---|---|---|
| `--min_library_size` | `5000` | Minimum per-sample read depth; samples below this are dropped |
| `--exclude_column` | `""` | Metadata column used to identify samples for exclusion |
| `--exclude_values` | `""` | Comma-separated values in `exclude_column` to remove |

### Grouping

| Parameter | Default | Description |
|---|---|---|
| `--groups_column` | `""` | Metadata column for the primary grouping variable |
| `--groups_paste_columns` | `""` | Comma-separated columns pasted together to form groups |
| `--type_column` | `""` | Secondary metadata column included in output CSVs |
| `--type2_column` | `""` | Tertiary metadata column included in output CSVs |
| `--type2_levels` | `""` | Ordered factor levels for `type2_column` |
| `--connections_column` | `""` | Column identifying paired/repeated-measures connections |

### Alpha diversity

| Parameter | Default | Description |
|---|---|---|
| `--which_level` | `Phylum` | Taxonomic level for feature collation (`Otus` \| `Genus` \| `Family` \| `Order` \| `Class` \| `Phylum`) |
| `--label` | `analysis` | Label prepended to output file names |
| `--alpha_metrics` | `Richness,Shannon,Simpson,FisherAlpha,PielouEvenness` | Comma-separated list of metrics to compute |
| `--test_method` | `anova` | Statistical test: `anova` \| `kruskal` \| `none` |
| `--p_adjust_method` | `BH` | Multiple-testing correction: `BH` \| `bonferroni` \| `holm` \| `fdr` \| `none` |

### Output options

| Parameter | Default | Description |
|---|---|---|
| `--merge_parquet` | `false` | Merge all output CSVs into a single Parquet file |

## Outputs

All files are written to `--output_dir`.

| File | Description |
|---|---|
| `Diversity_{level}_{label}.csv` | Wide-format diversity table (one row per sample, one column per metric) |
| `Diversity_long_{level}_{label}.csv` | Long-format version with `metric` and `value` columns |
| `Diversity_pairwise_{level}_{label}.csv` | Pairwise group comparison results (test statistic, p-value, adjusted p-value) |
| `alpha_diversity_{label}.parquet` | All three CSVs merged with `analysis` and `table` metadata columns (`--merge_parquet` only) |

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 23.04
- Docker (default) **or** a local R installation with: `optparse`, `stringr`, `data.table`, `vegan`, `phyloseq`, `arrow`

## Running without Docker

```bash
nextflow run main.nf \
  -c nextflow.config \
  --scripts_dir "$(pwd)" \
  --feature_table /path/to/table.biom \
  --meta_table    /path/to/meta.csv \
  --groups_column Treatment \
  --label         my_analysis
```

Add `docker.enabled = false` and `process.container = null` to a local override config, or pass them on the command line.
