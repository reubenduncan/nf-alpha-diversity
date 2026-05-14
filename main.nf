#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Alpha Diversity Pipeline
// Outputs: wide CSV, long CSV, pairwise statistics CSV, n CSV,
//          optional two-way ANOVA CSV, rarefaction data CSV,
//          metric correlation CSV, optional sensitivity CSV,
//          optional merged Parquet, pipeline manifest JSON,
//          R session info TXT.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Validate required parameters early to fail fast with a clear message
// ---------------------------------------------------------------------------
if (!params.feature_table)
    error "ERROR: --feature_table is required. Example: --feature_table data/table.biom"
if (!params.meta_table)
    error "ERROR: --meta_table is required. Example: --meta_table data/meta.csv"

// ---------------------------------------------------------------------------
// Process: Dump R package versions for reproducibility
// ---------------------------------------------------------------------------
process DUMP_VERSIONS {
    publishDir "${params.output_dir}", mode: 'copy'

    output:
    path "r_session_info.txt"

    script:
    """
    Rscript -e "writeLines(capture.output(sessionInfo()), 'r_session_info.txt')"
    """
}

// ---------------------------------------------------------------------------
// Process: Core alpha diversity analysis
// ---------------------------------------------------------------------------
process ALPHA_DIVERSITY {
    tag "${meta.id}"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    tuple val(meta), path(feature_table), path(meta_table)
    path taxonomy_table
    path tree_file

    output:
    path "Diversity_${params.taxon_rank}_${params.label}.csv",          emit: diversity_wide
    path "Diversity_long_${params.taxon_rank}_${params.label}.csv",     emit: diversity_long
    path "Diversity_pairwise_${params.taxon_rank}_${params.label}.csv", emit: pairwise_csv
    path "Diversity_n_${params.taxon_rank}_${params.label}.csv",        emit: n_csv
    path "Diversity_twoway_${params.taxon_rank}_${params.label}.csv",   emit: twoway_csv, optional: true

    script:
    def tax_arg      = (taxonomy_table.name != 'NO_TAXONOMY') ? "--taxonomy_table '${taxonomy_table}'" : ""
    def tree_arg     = (tree_file.name      != 'NO_FILE')     ? "--tree_file '${tree_file}'"           : ""
    def excl_col_arg = params.exclude_column ? "--exclude_column '${params.exclude_column}'" : ""
    def excl_val_arg = params.exclude_values ? "--exclude_values '${params.exclude_values}'" : ""
    def grp_arg      = params.group          ? "--group '${params.group}'"                   : ""
    def type_arg     = params.type           ? "--type  '${params.type}'"                    : ""

    """
    Rscript ${projectDir}/src/R/alpha_diversity.R \\
        --feature_table    '${feature_table}'              \\
        --meta_table       '${meta_table}'                 \\
        --input_format     '${params.input_format}'        \\
        --output_dir       '.'                             \\
        --taxon_rank       '${params.taxon_rank}'          \\
        --label            '${params.label}'               \\
        --min_library_size ${params.min_library_size}      \\
        --normalization    '${params.normalization}'       \\
        --alpha_metrics    '${params.alpha_metrics}'       \\
        --test_method      '${params.test_method}'         \\
        --p_adjust_method  '${params.p_adjust_method}'     \\
        ${tax_arg}      \\
        ${tree_arg}     \\
        ${excl_col_arg} \\
        ${excl_val_arg} \\
        ${grp_arg}      \\
        ${type_arg}
    """
}

// ---------------------------------------------------------------------------
// Process: Rarefaction curves (QC data — no plotting)
// ---------------------------------------------------------------------------
process RAREFACTION_CURVES {
    tag "${meta.id}"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    tuple val(meta), path(feature_table), path(meta_table)
    path taxonomy_table

    output:
    path "rarefaction_data_${params.label}.csv"

    script:
    def tax_arg = (taxonomy_table.name != 'NO_TAXONOMY') ? "--taxonomy_table '${taxonomy_table}'" : ""
    def grp_arg = params.group ? "--group '${params.group}'" : ""

    """
    Rscript ${projectDir}/src/R/rarefaction_curves.R \\
        --feature_table    '${feature_table}'          \\
        --meta_table       '${meta_table}'             \\
        --input_format     '${params.input_format}'    \\
        --min_library_size ${params.min_library_size}  \\
        --label            '${params.label}'           \\
        --output_dir       '.'                         \\
        ${tax_arg} \\
        ${grp_arg}
    """
}

// ---------------------------------------------------------------------------
// Process: Sensitivity analysis across rarefaction depths
//          Only runs under --normalization rarefaction
// ---------------------------------------------------------------------------
process SENSITIVITY_RAREFACTION {
    tag "${params.label}"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    tuple val(meta), path(feature_table), path(meta_table)
    path taxonomy_table
    path tree_file

    output:
    path "Diversity_sensitivity_${params.taxon_rank}_${params.label}.csv"

    script:
    def tax_arg  = (taxonomy_table.name != 'NO_TAXONOMY') ? "--taxonomy_table '${taxonomy_table}'" : ""
    def tree_arg = (tree_file.name      != 'NO_FILE')     ? "--tree_file '${tree_file}'"           : ""
    def grp_arg  = params.group ? "--group '${params.group}'" : ""

    """
    Rscript ${projectDir}/src/R/sensitivity_rarefaction.R \\
        --feature_table      '${feature_table}'              \\
        --meta_table         '${meta_table}'                 \\
        --input_format       '${params.input_format}'        \\
        --taxon_rank         '${params.taxon_rank}'          \\
        --alpha_metrics      '${params.alpha_metrics}'       \\
        --min_library_size   ${params.min_library_size}      \\
        --sensitivity_depths '${params.sensitivity_depths}'  \\
        --label              '${params.label}'               \\
        --output_dir         '.'                             \\
        ${tax_arg}  \\
        ${tree_arg} \\
        ${grp_arg}
    """
}

// ---------------------------------------------------------------------------
// Process: Metric correlation (CSV only — no plot)
// ---------------------------------------------------------------------------
process METRIC_CORRELATION {
    tag "${params.label}"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path wide_csv

    output:
    path "Diversity_correlation_${params.label}.csv", emit: correlation_csv

    script:
    """
    Rscript ${projectDir}/src/R/metric_correlation.R \\
        --wide_csv   '${wide_csv}'     \\
        --label      '${params.label}' \\
        --output_dir '.'
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
    Rscript ${projectDir}/src/R/merge_parquet.R \\
        --label      '${params.label}' \\
        --output_dir '.'
    """
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------
workflow {

    def meta_map = [id: params.label]

    ch_inputs = Channel.of(
        tuple(meta_map,
              file(params.feature_table, checkIfExists: true),
              file(params.meta_table,    checkIfExists: true))
    )

    // Optional taxonomy table (staged into work dir for S3 compatibility)
    ch_tax = (params.taxonomy_table && params.taxonomy_table != "")
        ? Channel.fromPath(params.taxonomy_table, checkIfExists: true).first()
        : Channel.value(file("NO_TAXONOMY"))

    // Optional tree file
    ch_tree = (params.tree_file && params.tree_file != "")
        ? Channel.fromPath(params.tree_file, checkIfExists: true).first()
        : Channel.value(file("NO_FILE"))

    DUMP_VERSIONS()

    ad_out = ALPHA_DIVERSITY(ch_inputs, ch_tax, ch_tree)

    RAREFACTION_CURVES(ch_inputs, ch_tax)

    METRIC_CORRELATION(ad_out.diversity_wide)

    // Sensitivity analysis only makes sense with rarefaction normalization
    if (params.sensitivity_analysis) {
        if (params.normalization != "rarefaction") {
            log.warn "sensitivity_analysis=true has no effect under normalization='${params.normalization}' — sensitivity requires normalization='rarefaction'."
        } else {
            SENSITIVITY_RAREFACTION(ch_inputs, ch_tax, ch_tree)
        }
    }

    if (params.merge_parquet) {
        all_csvs = ad_out.diversity_wide
            .mix(ad_out.diversity_long)
            .mix(ad_out.pairwise_csv)
            .collect()
        MERGE_PARQUET(all_csvs)
    }
}

// ---------------------------------------------------------------------------
// Write a machine-readable manifest of all parameters and run metadata.
// This file is consumed by downstream analysis agents (see SKILL.md).
// ---------------------------------------------------------------------------
workflow.onComplete {
    def manifest_data = [
        pipeline        : workflow.manifest.name,
        version         : workflow.manifest.version,
        run_name        : workflow.runName,
        session_id      : workflow.sessionId,
        nextflow_version: nextflow.version.toString(),
        start           : workflow.start.toString(),
        complete        : workflow.complete.toString(),
        duration        : workflow.duration.toString(),
        success         : workflow.success,
        exit_status     : workflow.exitStatus,
        command_line    : workflow.commandLine,
        profile         : workflow.profile,
        params          : params.toMapString()
    ]
    def json = groovy.json.JsonOutput.prettyPrint(
                   groovy.json.JsonOutput.toJson(manifest_data))
    new File("${params.output_dir}/pipeline_manifest.json").text = json
    log.info "Pipeline manifest written to ${params.output_dir}/pipeline_manifest.json"
}
