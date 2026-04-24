#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Alpha Diversity Pipeline
// Standalone workflow — no imports from parent project.
// Outputs: wide CSV, long CSV, pairwise statistics CSV.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Process
// ---------------------------------------------------------------------------
process ALPHA_DIVERSITY {
    tag "${meta.id}"

    publishDir "${params.output_dir}", mode: 'copy'

    input:
    tuple val(meta), path(feature_table), path(meta_table)
    path tree_file   // accepted for API consistency; not used by R script

    output:
    path "Diversity_${params.which_level}_${params.label}.csv",         emit: diversity_wide
    path "Diversity_long_${params.which_level}_${params.label}.csv",    emit: diversity_long
    path "Diversity_pairwise_${params.which_level}_${params.label}.csv", emit: pairwise_csv

    script:
    def tax_arg          = params.taxonomy_table  ? "--taxonomy_table '${params.taxonomy_table}'"  : ""
    def excl_col_arg     = params.exclude_column  ? "--exclude_column '${params.exclude_column}'"  : ""
    def excl_val_arg     = params.exclude_values  ? "--exclude_values '${params.exclude_values}'"  : ""
    def grp_col_arg      = params.groups_column   ? "--groups_column '${params.groups_column}'"    : ""
    def grp_paste_arg    = params.groups_paste_columns ? "--groups_paste_columns '${params.groups_paste_columns}'" : ""
    def type_col_arg     = params.type_column     ? "--type_column '${params.type_column}'"        : ""
    def type2_col_arg    = params.type2_column    ? "--type2_column '${params.type2_column}'"      : ""
    def type2_lvl_arg    = params.type2_levels    ? "--type2_levels '${params.type2_levels}'"      : ""
    def conn_col_arg     = params.connections_column ? "--connections_column '${params.connections_column}'" : ""

    """
    Rscript ${params.scripts_dir}/src/R/alpha_diversity.R \\
        --feature_table   '${feature_table}'             \\
        --meta_table      '${meta_table}'                \\
        --input_format    '${params.input_format}'       \\
        --output_dir      '.'                            \\
        --which_level     '${params.which_level}'        \\
        --label           '${params.label}'              \\
        --min_library_size ${params.min_library_size}    \\
        --alpha_metrics   '${params.alpha_metrics}'      \\
        --test_method     '${params.test_method}'        \\
        --p_adjust_method '${params.p_adjust_method}'    \\
        ${tax_arg}      \\
        ${excl_col_arg} \\
        ${excl_val_arg} \\
        ${grp_col_arg}  \\
        ${grp_paste_arg}\\
        ${type_col_arg} \\
        ${type2_col_arg}\\
        ${type2_lvl_arg}\\
        ${conn_col_arg}
    """
}

// ---------------------------------------------------------------------------
// Process: Merge all CSVs into a single Parquet file
// ---------------------------------------------------------------------------
process MERGE_PARQUET {
    tag "${params.label}"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path csvs

    output:
    path "alpha_diversity_${params.label}.parquet"

    script:
    """
    Rscript ${params.scripts_dir}/src/R/merge_parquet.R \\
        --label      '${params.label}' \\
        --output_dir '.'
    """
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------
workflow {

    // Build a single-element channel carrying feature_table + meta_table
    def meta_map = [id: params.label]

    ch_inputs = Channel.of(
        tuple(meta_map,
              file(params.feature_table, checkIfExists: true),
              file(params.meta_table,    checkIfExists: true))
    )

    // Optional tree file (null path if not provided)
    ch_tree = params.tree_file
        ? Channel.fromPath(params.tree_file, checkIfExists: true)
        : Channel.value(file("NO_FILE"))

    ad_out = ALPHA_DIVERSITY(ch_inputs, ch_tree)

    if (params.merge_parquet) {
        all_csvs = ad_out.diversity_wide
            .mix(ad_out.diversity_long)
            .mix(ad_out.pairwise_csv)
            .collect()
        MERGE_PARQUET(all_csvs)
    }
}
