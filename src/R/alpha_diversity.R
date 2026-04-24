#!/usr/bin/env Rscript
# alpha_diversity.R
# Main analysis script for alpha diversity calculation.
# Supports BIOM, TSV, and GTDB input formats.
# Outputs CSV only — no plots or HTML.

suppressPackageStartupMessages({
  library(optparse)
  library(phyloseq)
  library(stringr)
  library(data.table)
  library(vegan)
})

# ---------------------------------------------------------------------------
# Source ingestion helper (resolved relative to this script's location so it
# works both when called directly and via Nextflow's scripts_dir param).
# ---------------------------------------------------------------------------
script_dir <- tryCatch(
  dirname(normalizePath(sys.frames()[[1]]$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_flag <- grep("^--file=", args, value = TRUE)
    if (length(file_flag)) dirname(normalizePath(sub("^--file=", "", file_flag)))
    else "."
  }
)
source(file.path(script_dir, "load_feature_table.R"))

# ---------------------------------------------------------------------------
# Option parsing
# ---------------------------------------------------------------------------
option_list <- list(
  # Primary input (format-agnostic name)
  make_option("--feature_table",    type = "character", default = NULL,
              help = "Path to feature table (BIOM, TSV, or GTDB) [required]"),
  # Backwards-compatibility alias
  make_option("--biom_file",        type = "character", default = NULL,
              help = "Alias for --feature_table (BIOM input) [backwards compat]"),

  make_option("--meta_table",       type = "character", default = NULL,
              help = "Path to metadata CSV [required]"),
  make_option("--taxonomy_table",   type = "character", default = "",
              help = "Path to taxonomy TSV (required if input_format is tsv or gtdb)"),
  make_option("--input_format",     type = "character", default = "biom",
              help = "Input format: biom | tsv | gtdb [default: biom]"),
  make_option("--output_dir",       type = "character", default = ".",
              help = "Output directory [default: .]"),

  # Filtering
  make_option("--which_level",      type = "character", default = "Phylum",
              help = "Taxonomy level: Otus Genus Family Order Class Phylum [default: Phylum]"),
  make_option("--label",            type = "character", default = "analysis",
              help = "Analysis label used in output filenames [default: analysis]"),
  make_option("--min_library_size", type = "integer",   default = 5000L,
              help = "Minimum reads per sample [default: 5000]"),
  make_option("--exclude_column",   type = "character", default = "",
              help = "Metadata column to filter samples by"),
  make_option("--exclude_values",   type = "character", default = "",
              help = "Comma-separated values to exclude from exclude_column"),

  # Grouping
  make_option("--group", type = "character", default = "",
              help = "One metadata column, or comma-separated columns to paste as the group label"),
  make_option("--type",  type = "character", default = "",
              help = "One metadata column, or comma-separated columns to paste as the point-style label (optional)"),

  # Alpha diversity
  make_option("--alpha_metrics",    type = "character",
              default = "Richness,Shannon,Simpson,FisherAlpha,PielouEvenness",
              help = "Comma-separated metrics: Richness,Shannon,Simpson,FisherAlpha,PielouEvenness,Chao1,InvSimpson"),
  make_option("--test_method",      type = "character", default = "anova",
              help = "Statistical test: anova | kruskal | none [default: anova]"),
  make_option("--p_adjust_method",  type = "character", default = "BH",
              help = "P-value adjustment: BH | bonferroni | holm | fdr | none [default: BH]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---------------------------------------------------------------------------
# Backwards-compat: --biom_file is an alias for --feature_table
# ---------------------------------------------------------------------------
if (is.null(opt$feature_table) && !is.null(opt$biom_file)) {
  opt$feature_table <- opt$biom_file
  message("Note: --biom_file is deprecated; please use --feature_table.")
}

# ---------------------------------------------------------------------------
# Validate required parameters
# ---------------------------------------------------------------------------
if (is.null(opt$feature_table) || opt$feature_table == "")
  stop("--feature_table (or --biom_file) is required.")
if (is.null(opt$meta_table) || opt$meta_table == "")
  stop("--meta_table is required.")
if (!file.exists(opt$feature_table))
  stop("Feature table file not found: ", opt$feature_table)
if (!file.exists(opt$meta_table))
  stop("Metadata file not found: ", opt$meta_table)

valid_formats <- c("biom", "tsv", "gtdb")
if (!tolower(opt$input_format) %in% valid_formats)
  stop("--input_format must be one of: ", paste(valid_formats, collapse = ", "),
       ". Got: ", opt$input_format)

valid_levels <- c("Otus", "Genus", "Family", "Order", "Class", "Phylum")
if (!opt$which_level %in% valid_levels)
  stop("--which_level must be one of: ", paste(valid_levels, collapse = ", "),
       ". Got: ", opt$which_level)

valid_metrics <- c("Richness", "Shannon", "Simpson", "FisherAlpha",
                   "PielouEvenness", "Chao1", "InvSimpson")
requested_metrics <- trimws(strsplit(opt$alpha_metrics, ",")[[1]])
bad_metrics <- setdiff(requested_metrics, valid_metrics)
if (length(bad_metrics) > 0)
  stop("Unknown alpha metrics: ", paste(bad_metrics, collapse = ", "),
       ". Valid options: ", paste(valid_metrics, collapse = ", "))

if (!opt$test_method %in% c("anova", "kruskal", "none"))
  stop("--test_method must be one of: anova, kruskal, none. Got: ", opt$test_method)

valid_adjust <- c("BH", "bonferroni", "holm", "fdr", "none")
if (!opt$p_adjust_method %in% valid_adjust)
  stop("--p_adjust_method must be one of: ", paste(valid_adjust, collapse = ", "),
       ". Got: ", opt$p_adjust_method)

if (!dir.exists(opt$output_dir)) {
  message("Creating output directory: ", opt$output_dir)
  dir.create(opt$output_dir, recursive = TRUE)
}

# ---------------------------------------------------------------------------
# DATA IMPORT via shared helper
# ---------------------------------------------------------------------------
message("Loading feature table (format=", opt$input_format, ") ...")
tax_tbl_arg <- if (opt$taxonomy_table != "") opt$taxonomy_table else NULL
loaded <- load_feature_table(
  feature_table  = opt$feature_table,
  input_format   = opt$input_format,
  taxonomy_table = tax_tbl_arg
)
abund_table  <- loaded$abund_table   # samples x features
OTU_taxonomy <- loaded$OTU_taxonomy  # features x 7

# ---------------------------------------------------------------------------
# METADATA
# ---------------------------------------------------------------------------
message("Loading metadata: ", opt$meta_table)
meta_table <- tryCatch({
  dt <- fread(opt$meta_table, header = TRUE, check.names = FALSE)
  df <- as.data.frame(dt, stringsAsFactors = FALSE)
  rownames(df) <- df[[1]]
  df[, -1, drop = FALSE]
}, error = function(e) stop("Failed to read metadata: ", conditionMessage(e)))

# ---------------------------------------------------------------------------
# LIBRARY SIZE FILTER
# ---------------------------------------------------------------------------
lib_sizes   <- rowSums(abund_table)
keep_samps  <- lib_sizes >= opt$min_library_size
n_dropped   <- sum(!keep_samps)
if (n_dropped > 0)
  message("Dropped ", n_dropped, " sample(s) below min_library_size=", opt$min_library_size)
abund_table <- abund_table[keep_samps, , drop = FALSE]

if (nrow(abund_table) < 2)
  stop("Fewer than 2 samples remain after min_library_size filter (",
       opt$min_library_size, "). Cannot continue.")

# Remove zero-count features
abund_table  <- abund_table[, colSums(abund_table) > 0, drop = FALSE]
OTU_taxonomy <- OTU_taxonomy[colnames(abund_table), , drop = FALSE]

# ---------------------------------------------------------------------------
# ALIGN SAMPLES WITH METADATA
# ---------------------------------------------------------------------------
common_samps <- intersect(rownames(abund_table), rownames(meta_table))
if (length(common_samps) < 2)
  stop("Fewer than 2 samples overlap between feature table and metadata. ",
       "Check that row names of metadata match sample IDs in the feature table.")

abund_table <- abund_table[common_samps, , drop = FALSE]
meta_table  <- meta_table[common_samps, , drop = FALSE]

# Remove any features that are now zero across retained samples
abund_table  <- abund_table[, colSums(abund_table) > 0, drop = FALSE]
OTU_taxonomy <- OTU_taxonomy[colnames(abund_table), , drop = FALSE]

# ---------------------------------------------------------------------------
# HYPOTHESIS SPACE FILTERS
# ---------------------------------------------------------------------------
if (opt$exclude_column != "" && opt$exclude_values != "") {
  if (!opt$exclude_column %in% colnames(meta_table))
    stop("exclude_column '", opt$exclude_column, "' not found in metadata columns.")
  exc_vals   <- trimws(strsplit(opt$exclude_values, ",")[[1]])
  keep_mask  <- !meta_table[[opt$exclude_column]] %in% exc_vals
  n_excl     <- sum(!keep_mask)
  message("Excluding ", n_excl, " sample(s) based on exclude_column='",
          opt$exclude_column, "', exclude_values='", opt$exclude_values, "'")
  meta_table  <- meta_table[keep_mask, , drop = FALSE]
  abund_table <- abund_table[rownames(meta_table), , drop = FALSE]
  abund_table <- abund_table[, colSums(abund_table) > 0, drop = FALSE]
  OTU_taxonomy <- OTU_taxonomy[colnames(abund_table), , drop = FALSE]
  if (nrow(abund_table) < 2)
    stop("Fewer than 2 samples remain after exclusion filter.")
}

# ---------------------------------------------------------------------------
# GROUPING
# ---------------------------------------------------------------------------
resolve_columns <- function(param_val, param_name, meta) {
  cols <- trimws(strsplit(param_val, ",")[[1]])
  missing <- setdiff(cols, colnames(meta))
  if (length(missing) > 0)
    stop(param_name, " references columns not in metadata: ", paste(missing, collapse = ", "))
  if (length(cols) == 1) {
    as.factor(as.character(meta[[cols]]))
  } else {
    as.factor(do.call(paste, c(meta[, cols, drop = FALSE], sep = " ")))
  }
}

if (opt$group != "") {
  meta_table$Groups <- resolve_columns(opt$group, "--group", meta_table)
} else {
  meta_table$Groups <- factor(rep("All", nrow(meta_table)))
  message("No --group specified — all samples assigned to group 'All'.")
}

if (opt$type != "") {
  meta_table$Type <- resolve_columns(opt$type, "--type", meta_table)
} else {
  meta_table$Type <- NULL
}

abund_table  <- abund_table[rownames(meta_table), , drop = FALSE]
abund_table  <- abund_table[, colSums(abund_table) > 0, drop = FALSE]
OTU_taxonomy <- OTU_taxonomy[colnames(abund_table), , drop = FALSE]

# ---------------------------------------------------------------------------
# COLLATE AT TAXONOMIC LEVEL
# ---------------------------------------------------------------------------
which_level <- opt$which_level
message("Collating features at taxonomic level: ", which_level)

if (which_level == "Otus") {
  new_abund_table <- abund_table
} else {
  lvl_list        <- unique(OTU_taxonomy[, which_level])
  new_abund_table <- NULL
  for (i in lvl_list) {
    feat_idx <- rownames(OTU_taxonomy)[OTU_taxonomy[, which_level] == i]
    feat_idx <- intersect(feat_idx, colnames(abund_table))
    if (length(feat_idx) == 0) next
    col_name <- if (is.na(i) || i == "") "__Unknowns__" else i
    tmp <- data.frame(rowSums(abund_table[, feat_idx, drop = FALSE]),
                      check.names = FALSE)
    colnames(tmp) <- col_name
    new_abund_table <- if (is.null(new_abund_table)) tmp else cbind(new_abund_table, tmp)
  }
  if (is.null(new_abund_table))
    stop("No features could be collated at level '", which_level, "'.")
}

abund_table <- as.matrix(new_abund_table)
storage.mode(abund_table) <- "numeric"

# ---------------------------------------------------------------------------
# DIVERSITY CALCULATIONS
# ---------------------------------------------------------------------------
message("Computing alpha diversity metrics: ", paste(requested_metrics, collapse = ", "))

metric_frames <- list()

compute_richness <- function(mat) {
  tryCatch(
    vegan::rarefy(mat, min(rowSums(mat))),
    error = function(e) {
      message("Warning: Richness (rarefy) failed: ", conditionMessage(e))
      setNames(rep(NA_real_, nrow(mat)), rownames(mat))
    }
  )
}

compute_shannon <- function(mat) {
  tryCatch(vegan::diversity(mat, index = "shannon"),
           error = function(e) {
             message("Warning: Shannon failed: ", conditionMessage(e))
             setNames(rep(NA_real_, nrow(mat)), rownames(mat))
           })
}

compute_simpson <- function(mat) {
  tryCatch(vegan::diversity(mat, index = "simpson"),
           error = function(e) {
             message("Warning: Simpson failed: ", conditionMessage(e))
             setNames(rep(NA_real_, nrow(mat)), rownames(mat))
           })
}

compute_invsimpson <- function(mat) {
  tryCatch(vegan::diversity(mat, index = "invsimpson"),
           error = function(e) {
             message("Warning: InvSimpson failed: ", conditionMessage(e))
             setNames(rep(NA_real_, nrow(mat)), rownames(mat))
           })
}

compute_fisher <- function(mat) {
  tryCatch(vegan::fisher.alpha(mat),
           error = function(e) {
             message("Warning: FisherAlpha failed: ", conditionMessage(e))
             setNames(rep(NA_real_, nrow(mat)), rownames(mat))
           })
}

compute_chao1 <- function(mat) {
  tryCatch({
    est <- vegan::estimateR(mat)
    est["S.chao1", ]
  }, error = function(e) {
    message("Warning: Chao1 failed: ", conditionMessage(e))
    setNames(rep(NA_real_, nrow(mat)), rownames(mat))
  })
}

# Pre-compute Shannon for Pielou (shared)
shannon_vals <- NULL

for (metric in requested_metrics) {
  vals <- switch(metric,
    Richness      = compute_richness(abund_table),
    Shannon       = { shannon_vals <<- compute_shannon(abund_table); shannon_vals },
    Simpson       = compute_simpson(abund_table),
    InvSimpson    = compute_invsimpson(abund_table),
    FisherAlpha   = compute_fisher(abund_table),
    Chao1         = compute_chao1(abund_table),
    PielouEvenness = {
      if (is.null(shannon_vals)) shannon_vals <<- compute_shannon(abund_table)
      S <- vegan::specnumber(abund_table)
      tryCatch(shannon_vals / log(S),
               error = function(e) {
                 message("Warning: PielouEvenness failed: ", conditionMessage(e))
                 setNames(rep(NA_real_, nrow(abund_table)), rownames(abund_table))
               })
    }
  )
  # Align names to abund_table rows
  vals_named <- setNames(as.numeric(vals), names(vals))
  metric_frames[[metric]] <- data.frame(
    sample  = rownames(abund_table),
    value   = vals_named[rownames(abund_table)],
    measure = metric,
    stringsAsFactors = FALSE
  )
}

df_long <- do.call(rbind, metric_frames)
rownames(df_long) <- NULL

# ---------------------------------------------------------------------------
# CSV 1: Wide diversity table
# ---------------------------------------------------------------------------
wide_df <- data.frame(row.names = rownames(abund_table))
for (metric in requested_metrics) {
  mf <- metric_frames[[metric]]
  rownames(mf) <- mf$sample
  wide_df[[metric]] <- mf[rownames(wide_df), "value"]
}
wide_df$Groups <- meta_table[rownames(wide_df), "Groups"]
if (!is.null(meta_table$Type))
  wide_df$Type  <- meta_table[rownames(wide_df), "Type"]
# Add sample as explicit column (first col)
wide_df <- cbind(sample = rownames(wide_df), wide_df)

wide_csv <- file.path(opt$output_dir,
  paste0("Diversity_", which_level, "_", opt$label, ".csv"))
write.csv(wide_df, wide_csv, row.names = FALSE)
message("Wrote wide CSV: ", wide_csv)

# ---------------------------------------------------------------------------
# Attach metadata to long-format df
# ---------------------------------------------------------------------------
df_long$Groups <- meta_table[df_long$sample, "Groups"]
if (!is.null(meta_table$Type))
  df_long$Type  <- meta_table[df_long$sample, "Type"]

# ---------------------------------------------------------------------------
# CSV 2: Long diversity table
# ---------------------------------------------------------------------------
long_cols <- c("sample", "value", "measure", "Groups")
if (!is.null(meta_table$Type)) long_cols <- c(long_cols, "Type")

long_csv <- file.path(opt$output_dir,
  paste0("Diversity_long_", which_level, "_", opt$label, ".csv"))
write.csv(df_long[, long_cols, drop = FALSE], long_csv, row.names = FALSE)
message("Wrote long CSV: ", long_csv)

# ---------------------------------------------------------------------------
# CSV 3: Pairwise tests
# ---------------------------------------------------------------------------
grouping_column <- "Groups"
group_levels    <- levels(meta_table$Groups)
n_groups        <- length(group_levels)

pw_csv <- file.path(opt$output_dir,
  paste0("Diversity_pairwise_", which_level, "_", opt$label, ".csv"))

if (opt$test_method == "none" || n_groups < 2) {
  if (opt$test_method == "none") {
    message("test_method='none' — skipping pairwise tests.")
  } else {
    message("Warning: Only ", n_groups, " group(s) found — skipping pairwise tests.")
  }
  empty_pw <- data.frame(
    measure     = character(),
    group1      = character(),
    group2      = character(),
    pvalue      = numeric(),
    padj        = numeric(),
    significance = character(),
    test_method = character(),
    stringsAsFactors = FALSE
  )
  write.csv(empty_pw, pw_csv, row.names = FALSE)
  message("Wrote empty pairwise CSV: ", pw_csv)
} else {
  s       <- combn(group_levels, 2)
  pw_rows <- list()

  for (k in requested_metrics) {
    for (l in seq_len(ncol(s))) {
      g1  <- s[1, l]
      g2  <- s[2, l]
      sub <- df_long[df_long$measure == k &
                     df_long[[grouping_column]] %in% c(g1, g2), ]
      sub <- sub[!is.na(sub$value), ]

      # Check each group has >= 2 samples
      counts <- table(sub[[grouping_column]])
      thin_groups <- names(counts[counts < 2])
      if (length(thin_groups) > 0) {
        message("Warning: Skipping pair (", g1, " vs ", g2, ") for metric '",
                k, "' — group(s) ", paste(thin_groups, collapse = ", "),
                " have fewer than 2 samples.")
        pw_rows[[length(pw_rows) + 1]] <- data.frame(
          measure      = k,
          group1       = g1,
          group2       = g2,
          pvalue       = NA_real_,
          padj         = NA_real_,
          significance = "insufficient_data",
          test_method  = opt$test_method,
          stringsAsFactors = FALSE
        )
        next
      }

      pv <- tryCatch({
        if (opt$test_method == "anova") {
          summary(aov(
            as.formula(paste("value ~", grouping_column)),
            data = sub
          ))[[1]][["Pr(>F)"]][1]
        } else {
          # kruskal
          kt <- kruskal.test(
            as.formula(paste("value ~", grouping_column)),
            data = sub
          )
          kt$p.value
        }
      }, error = function(e) {
        message("Warning: Statistical test failed for metric '", k,
                "' pair (", g1, " vs ", g2, "): ", conditionMessage(e))
        NA_real_
      })

      pw_rows[[length(pw_rows) + 1]] <- data.frame(
        measure      = k,
        group1       = g1,
        group2       = g2,
        pvalue       = pv,
        padj         = NA_real_,   # filled below
        significance = "",
        test_method  = opt$test_method,
        stringsAsFactors = FALSE
      )
    }
  }

  df_pairwise <- do.call(rbind, pw_rows)

  # P-value adjustment (per metric, across all pairs)
  if (opt$p_adjust_method != "none" && nrow(df_pairwise) > 0) {
    for (k in unique(df_pairwise$measure)) {
      idx <- df_pairwise$measure == k
      raw <- df_pairwise$pvalue[idx]
      df_pairwise$padj[idx] <- p.adjust(raw, method = opt$p_adjust_method)
    }
  } else {
    df_pairwise$padj <- df_pairwise$pvalue
  }

  # Significance stars (based on padj)
  sig_label <- function(p) {
    if (is.na(p)) return("")
    if (p <= 0.001) "***"
    else if (p <= 0.01) "**"
    else if (p <= 0.05) "*"
    else ""
  }
  df_pairwise$significance <- sapply(df_pairwise$padj, sig_label)

  write.csv(df_pairwise, pw_csv, row.names = FALSE)
  message("Wrote pairwise CSV: ", pw_csv)
}

message("Alpha diversity analysis complete.")
