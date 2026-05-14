#!/usr/bin/env Rscript
# alpha_diversity.R
# Main analysis script for alpha diversity calculation.
# Supports BIOM, TSV, and GTDB input formats.
# Outputs CSV only — no plots or HTML.

suppressPackageStartupMessages({
  library(optparse)
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
  make_option("--tree_file",        type = "character", default = "",
              help = "Path to Newick tree file (required for FaithsPD metric)"),
  make_option("--input_format",     type = "character", default = "biom",
              help = "Input format: biom | tsv | gtdb [default: biom]"),
  make_option("--output_dir",       type = "character", default = ".",
              help = "Output directory [default: .]"),

  # Filtering
  make_option("--taxon_rank",       type = "character", default = "Feature",
              help = "Taxonomy level: Feature Genus Family Order Class Phylum [default: Feature]"),
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

  # Normalization
  make_option("--normalization",    type = "character", default = "clr",
              help = "clr (default): proportion-based, no rarefaction. rarefaction: subsample to min_library_size."),

  # Alpha diversity
  make_option("--alpha_metrics",    type = "character",
              default = "Richness,Shannon,GiniSimpson,FisherAlpha,PielouEvenness",
              help = paste0("Comma-separated metrics: Richness,Shannon,GiniSimpson,InvSimpson,",
                            "FisherAlpha,PielouEvenness,Chao1,FaithsPD")),
  make_option("--test_method",      type = "character", default = "anova",
              help = paste0("Statistical test: anova | kruskal | none | auto [default: anova]. ",
                            "auto runs Levene's test (Brown-Forsythe) per metric and selects kruskal ",
                            "when variances are heterogeneous (p<=0.05), anova otherwise.")),
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

valid_levels <- c("Feature", "Genus", "Family", "Order", "Class", "Phylum")
if (!opt$taxon_rank %in% valid_levels)
  stop("--taxon_rank must be one of: ", paste(valid_levels, collapse = ", "),
       ". Got: ", opt$taxon_rank)

if (!opt$normalization %in% c("clr", "rarefaction"))
  stop("--normalization must be 'clr' or 'rarefaction'. Got: ", opt$normalization)

valid_metrics <- c("Richness", "Shannon", "GiniSimpson", "InvSimpson",
                   "FisherAlpha", "PielouEvenness", "Chao1", "FaithsPD")
requested_metrics <- trimws(strsplit(opt$alpha_metrics, ",")[[1]])
bad_metrics <- setdiff(requested_metrics, valid_metrics)
if (length(bad_metrics) > 0)
  stop("Unknown alpha metrics: ", paste(bad_metrics, collapse = ", "),
       ". Valid options: ", paste(valid_metrics, collapse = ", "))

if (!opt$test_method %in% c("anova", "kruskal", "none", "auto"))
  stop("--test_method must be one of: anova, kruskal, none, auto. Got: ", opt$test_method)

valid_adjust <- c("BH", "bonferroni", "holm", "fdr", "none")
if (!opt$p_adjust_method %in% valid_adjust)
  stop("--p_adjust_method must be one of: ", paste(valid_adjust, collapse = ", "),
       ". Got: ", opt$p_adjust_method)

if (!dir.exists(opt$output_dir)) {
  message("Creating output directory: ", opt$output_dir)
  dir.create(opt$output_dir, recursive = TRUE)
}

# ---------------------------------------------------------------------------
# Load phylogenetic tree (required only for FaithsPD)
# ---------------------------------------------------------------------------
tree <- NULL
tree_path <- opt$tree_file
if (!is.null(tree_path) && nchar(tree_path) > 0 &&
    tree_path != "NO_FILE" && file.exists(tree_path)) {
  message("Loading phylogenetic tree: ", tree_path)
  tree <- tryCatch(
    ape::read.tree(tree_path),
    error = function(e) stop("Failed to read tree file: ", conditionMessage(e))
  )
}

if ("FaithsPD" %in% requested_metrics && is.null(tree))
  stop("FaithsPD metric requires --tree_file to be provided.")

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
abund_table      <- loaded$abund_table       # samples x features
feature_taxonomy <- loaded$feature_taxonomy  # features x 7

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
abund_table      <- abund_table[, colSums(abund_table) > 0, drop = FALSE]
feature_taxonomy <- feature_taxonomy[colnames(abund_table), , drop = FALSE]

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
abund_table      <- abund_table[, colSums(abund_table) > 0, drop = FALSE]
feature_taxonomy <- feature_taxonomy[colnames(abund_table), , drop = FALSE]

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
  feature_taxonomy <- feature_taxonomy[colnames(abund_table), , drop = FALSE]
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

abund_table      <- abund_table[rownames(meta_table), , drop = FALSE]
abund_table      <- abund_table[, colSums(abund_table) > 0, drop = FALSE]
feature_taxonomy <- feature_taxonomy[colnames(abund_table), , drop = FALSE]

# ---------------------------------------------------------------------------
# COLLATE AT TAXONOMIC LEVEL
# ---------------------------------------------------------------------------
taxon_rank <- opt$taxon_rank
message("Collating features at taxonomic level: ", taxon_rank)

if (taxon_rank == "Feature") {
  new_abund_table <- abund_table
} else {
  lvl_list        <- unique(feature_taxonomy[, taxon_rank])
  new_abund_table <- NULL
  for (i in lvl_list) {
    feat_idx <- rownames(feature_taxonomy)[feature_taxonomy[, taxon_rank] == i]
    feat_idx <- intersect(feat_idx, colnames(abund_table))
    if (length(feat_idx) == 0) next
    col_name <- if (is.na(i) || i == "") "__Unknowns__" else i
    tmp <- data.frame(rowSums(abund_table[, feat_idx, drop = FALSE]),
                      check.names = FALSE)
    colnames(tmp) <- col_name
    new_abund_table <- if (is.null(new_abund_table)) tmp else cbind(new_abund_table, tmp)
  }
  if (is.null(new_abund_table))
    stop("No features could be collated at level '", taxon_rank, "'.")
}

abund_table <- as.matrix(new_abund_table)
storage.mode(abund_table) <- "numeric"

# ---------------------------------------------------------------------------
# DIVERSITY CALCULATIONS
# ---------------------------------------------------------------------------
message("Normalization mode: ", opt$normalization)
message("Computing alpha diversity metrics: ", paste(requested_metrics, collapse = ", "))

# Faith's PD helper — only called when tree is available
compute_faiths_pd <- function(mat, tr) {
  if (!requireNamespace("picante", quietly = TRUE))
    stop("Package 'picante' is required for FaithsPD.")
  common <- intersect(colnames(mat), tr$tip.label)
  if (length(common) < 2)
    stop("Fewer than 2 features are shared between the feature table and the tree. ",
         "Check that feature IDs match tree tip labels.")
  mat_pruned  <- mat[, common, drop = FALSE]
  tree_pruned <- ape::keep.tip(tr, common)
  pd_res      <- picante::pd(mat_pruned, tree_pruned, include.root = TRUE)
  setNames(pd_res$PD, rownames(mat))
}

metric_frames <- list()

# Pre-compute Shannon once (reused by PielouEvenness)
shannon_vals <- NULL

for (metric in requested_metrics) {
  vals <- switch(metric,

    Richness = {
      if (opt$normalization == "clr") {
        # Observed richness: no rarefaction, no data loss
        tryCatch(
          vegan::specnumber(abund_table),
          error = function(e) {
            message("Warning: Richness (specnumber) failed: ", conditionMessage(e))
            setNames(rep(NA_real_, nrow(abund_table)), rownames(abund_table))
          }
        )
      } else {
        tryCatch(
          vegan::rarefy(abund_table, opt$min_library_size),
          error = function(e) {
            message("Warning: Richness (rarefy) failed: ", conditionMessage(e))
            setNames(rep(NA_real_, nrow(abund_table)), rownames(abund_table))
          }
        )
      }
    },

    Shannon = {
      tryCatch({
        # vegan::diversity normalises to proportions internally; scale-invariant
        shannon_vals <<- vegan::diversity(abund_table, index = "shannon")
        shannon_vals
      }, error = function(e) {
        message("Warning: Shannon failed: ", conditionMessage(e))
        setNames(rep(NA_real_, nrow(abund_table)), rownames(abund_table))
      })
    },

    # vegan computes Gini-Simpson = 1 - sum(p^2), NOT the raw Simpson D = sum(p^2)
    GiniSimpson = {
      tryCatch(
        vegan::diversity(abund_table, index = "simpson"),
        error = function(e) {
          message("Warning: GiniSimpson failed: ", conditionMessage(e))
          setNames(rep(NA_real_, nrow(abund_table)), rownames(abund_table))
        }
      )
    },

    InvSimpson = {
      tryCatch(
        vegan::diversity(abund_table, index = "invsimpson"),
        error = function(e) {
          message("Warning: InvSimpson failed: ", conditionMessage(e))
          setNames(rep(NA_real_, nrow(abund_table)), rownames(abund_table))
        }
      )
    },

    FisherAlpha = {
      tryCatch(
        vegan::fisher.alpha(abund_table),
        error = function(e) {
          message("Warning: FisherAlpha failed: ", conditionMessage(e))
          setNames(rep(NA_real_, nrow(abund_table)), rownames(abund_table))
        }
      )
    },

    Chao1 = {
      tryCatch({
        est <- vegan::estimateR(abund_table)
        est["S.chao1", ]
      }, error = function(e) {
        message("Warning: Chao1 failed: ", conditionMessage(e))
        setNames(rep(NA_real_, nrow(abund_table)), rownames(abund_table))
      })
    },

    PielouEvenness = {
      if (is.null(shannon_vals))
        shannon_vals <<- vegan::diversity(abund_table, index = "shannon")
      # Use per-sample observed richness as denominator (correct for any normalization)
      obs_richness <- vegan::specnumber(abund_table)
      tryCatch(
        shannon_vals / log(obs_richness),
        error = function(e) {
          message("Warning: PielouEvenness failed: ", conditionMessage(e))
          setNames(rep(NA_real_, nrow(abund_table)), rownames(abund_table))
        }
      )
    },

    FaithsPD = {
      tryCatch(
        compute_faiths_pd(abund_table, tree),
        error = function(e) {
          message("Warning: FaithsPD failed: ", conditionMessage(e))
          setNames(rep(NA_real_, nrow(abund_table)), rownames(abund_table))
        }
      )
    }
  )

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
wide_df$Groups        <- meta_table[rownames(wide_df), "Groups"]
wide_df$normalization <- opt$normalization
if (!is.null(meta_table$Type))
  wide_df$Type <- meta_table[rownames(wide_df), "Type"]
wide_df <- cbind(sample = rownames(wide_df), wide_df)

wide_csv <- file.path(opt$output_dir,
  paste0("Diversity_", taxon_rank, "_", opt$label, ".csv"))
write.csv(wide_df, wide_csv, row.names = FALSE)
message("Wrote wide CSV: ", wide_csv)

# ---------------------------------------------------------------------------
# Attach metadata to long-format df
# ---------------------------------------------------------------------------
df_long$Groups <- meta_table[df_long$sample, "Groups"]
if (!is.null(meta_table$Type))
  df_long$Type <- meta_table[df_long$sample, "Type"]

# ---------------------------------------------------------------------------
# CSV 2: Long diversity table
# ---------------------------------------------------------------------------
long_cols <- c("sample", "value", "measure", "Groups")
if (!is.null(meta_table$Type)) long_cols <- c(long_cols, "Type")

long_csv <- file.path(opt$output_dir,
  paste0("Diversity_long_", taxon_rank, "_", opt$label, ".csv"))
write.csv(df_long[, long_cols, drop = FALSE], long_csv, row.names = FALSE)
message("Wrote long CSV: ", long_csv)

# ---------------------------------------------------------------------------
# Helper: Brown-Forsythe Levene test (median-based, no extra dependencies)
# Returns p-value from one-way ANOVA on absolute deviations from group medians
# ---------------------------------------------------------------------------
levene_pval <- function(values, groups) {
  groups <- factor(groups)
  grp_med <- tapply(values, groups, median, na.rm = TRUE)
  abs_dev <- abs(values - grp_med[as.character(groups)])
  tryCatch(
    summary(aov(abs_dev ~ groups))[[1]][["Pr(>F)"]][1],
    error = function(e) NA_real_
  )
}

# ---------------------------------------------------------------------------
# CSV 3: Pairwise tests
# ---------------------------------------------------------------------------
grouping_column <- "Groups"
group_levels    <- levels(meta_table$Groups)
n_groups        <- length(group_levels)

pw_csv <- file.path(opt$output_dir,
  paste0("Diversity_pairwise_", taxon_rank, "_", opt$label, ".csv"))

empty_pw <- data.frame(
  measure            = character(),
  group1             = character(),
  group2             = character(),
  estimate           = numeric(),
  conf_low           = numeric(),
  conf_high          = numeric(),
  pvalue             = numeric(),
  padj               = numeric(),
  significance       = character(),
  test_method        = character(),
  actual_test_method = character(),
  levene_pvalue      = numeric(),
  stringsAsFactors   = FALSE
)

if (opt$test_method == "none" || n_groups < 2) {
  if (opt$test_method == "none") {
    message("test_method='none' — skipping pairwise tests.")
  } else {
    message("Warning: Only ", n_groups, " group(s) found — skipping pairwise tests.")
  }
  write.csv(empty_pw, pw_csv, row.names = FALSE)
  message("Wrote empty pairwise CSV: ", pw_csv)
} else {
  pw_rows <- list()

  for (k in requested_metrics) {
    levene_p      <- NA_real_
    actual_method <- opt$test_method

    if (opt$test_method == "auto") {
      df_all   <- df_long[df_long$measure == k & !is.na(df_long$value), ]
      levene_p <- levene_pval(df_all$value, df_all[[grouping_column]])
      actual_method <- if (!is.na(levene_p) && levene_p <= 0.05) "kruskal" else "anova"
      message("  Metric '", k, "': Levene p=", sprintf("%.4g", levene_p),
              " — using ", actual_method)
    }

    if (actual_method == "anova" && n_groups >= 2) {
      # Tukey HSD: standard post-hoc for ANOVA, controls family-wise error rate
      df_met <- df_long[df_long$measure == k & !is.na(df_long$value), ]
      df_met[[grouping_column]] <- factor(df_met[[grouping_column]])

      tryCatch({
        fit     <- aov(as.formula(paste("value ~", grouping_column)), data = df_met)
        tukey   <- TukeyHSD(fit, which = grouping_column)[[grouping_column]]
        for (row_nm in rownames(tukey)) {
          parts <- strsplit(row_nm, "-")[[1]]
          # TukeyHSD names are "b-a" (second minus first); swap to g1 vs g2 convention
          g1 <- trimws(parts[2])
          g2 <- trimws(parts[1])
          pw_rows[[length(pw_rows) + 1]] <- data.frame(
            measure            = k,
            group1             = g1,
            group2             = g2,
            estimate           = tukey[row_nm, "diff"],
            conf_low           = tukey[row_nm, "lwr"],
            conf_high          = tukey[row_nm, "upr"],
            pvalue             = tukey[row_nm, "p adj"],
            padj               = tukey[row_nm, "p adj"],  # already adjusted
            significance       = "",
            test_method        = opt$test_method,
            actual_test_method = actual_method,
            levene_pvalue      = levene_p,
            stringsAsFactors   = FALSE
          )
        }
      }, error = function(e) {
        message("Warning: TukeyHSD failed for metric '", k, "': ", conditionMessage(e))
      })

    } else if (actual_method == "kruskal") {
      # Kruskal-Wallis: pairwise Mann-Whitney U with p-value correction
      s <- combn(group_levels, 2)

      for (l in seq_len(ncol(s))) {
        g1  <- s[1, l]
        g2  <- s[2, l]
        sub <- df_long[df_long$measure == k &
                       df_long[[grouping_column]] %in% c(g1, g2), ]
        sub <- sub[!is.na(sub$value), ]

        counts      <- table(sub[[grouping_column]])
        thin_groups <- names(counts[counts < 2])
        if (length(thin_groups) > 0) {
          message("Warning: Skipping pair (", g1, " vs ", g2, ") for metric '", k,
                  "' — group(s) ", paste(thin_groups, collapse = ", "),
                  " have fewer than 2 samples.")
          pw_rows[[length(pw_rows) + 1]] <- data.frame(
            measure = k, group1 = g1, group2 = g2,
            estimate = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
            pvalue = NA_real_, padj = NA_real_,
            significance = "insufficient_data",
            test_method = opt$test_method, actual_test_method = actual_method,
            levene_pvalue = levene_p, stringsAsFactors = FALSE
          )
          next
        }

        pv <- tryCatch(
          kruskal.test(
            as.formula(paste("value ~", grouping_column)),
            data = sub
          )$p.value,
          error = function(e) {
            message("Warning: Kruskal test failed for metric '", k,
                    "' pair (", g1, " vs ", g2, "): ", conditionMessage(e))
            NA_real_
          }
        )

        pw_rows[[length(pw_rows) + 1]] <- data.frame(
          measure = k, group1 = g1, group2 = g2,
          estimate = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
          pvalue = pv, padj = NA_real_, significance = "",
          test_method = opt$test_method, actual_test_method = actual_method,
          levene_pvalue = levene_p, stringsAsFactors = FALSE
        )
      }

      # Apply multiple-testing correction across all pairs for this metric
      if (opt$p_adjust_method != "none" && length(pw_rows) > 0) {
        k_rows <- sapply(pw_rows, function(r) r$measure == k)
        raw_pv <- sapply(pw_rows[k_rows], function(r) r$pvalue)
        adj_pv <- p.adjust(raw_pv, method = opt$p_adjust_method)
        idx    <- which(k_rows)
        for (j in seq_along(idx)) pw_rows[[idx[j]]]$padj <- adj_pv[j]
      } else {
        for (j in seq_along(pw_rows)) {
          if (pw_rows[[j]]$measure == k)
            pw_rows[[j]]$padj <- pw_rows[[j]]$pvalue
        }
      }
    }
  }

  if (length(pw_rows) == 0) {
    df_pairwise <- empty_pw
  } else {
    df_pairwise <- do.call(rbind, pw_rows)
  }

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

# ---------------------------------------------------------------------------
# CSV 4: Per-group sample counts
# ---------------------------------------------------------------------------
n_csv <- file.path(opt$output_dir,
  paste0("Diversity_n_", taxon_rank, "_", opt$label, ".csv"))
n_df  <- as.data.frame(table(Groups = meta_table$Groups),
                        stringsAsFactors = FALSE)
colnames(n_df) <- c("Groups", "n")
write.csv(n_df, n_csv, row.names = FALSE)
message("Wrote n CSV: ", n_csv)

# ---------------------------------------------------------------------------
# CSV 5: Two-way ANOVA (Groups * Type) — only when --type was provided
# Note: exploratory only; no post-hoc for interaction term is computed.
# ---------------------------------------------------------------------------
if (!is.null(meta_table$Type) && opt$test_method != "none") {
  tw_csv <- file.path(opt$output_dir,
    paste0("Diversity_twoway_", taxon_rank, "_", opt$label, ".csv"))
  tw_rows <- list()

  for (k in requested_metrics) {
    df_sub <- df_long[df_long$measure == k & !is.na(df_long$value), ]
    df_sub[[grouping_column]] <- factor(df_sub[[grouping_column]])
    df_sub$Type               <- factor(df_sub$Type)

    tryCatch({
      fit <- aov(
        as.formula(paste("value ~", grouping_column, "* Type")),
        data = df_sub
      )
      tbl <- summary(fit)[[1]]
      for (term in rownames(tbl)) {
        tw_rows[[length(tw_rows) + 1]] <- data.frame(
          measure  = k,
          term     = trimws(term),
          df       = tbl[term, "Df"],
          sum_sq   = tbl[term, "Sum Sq"],
          mean_sq  = tbl[term, "Mean Sq"],
          f_value  = tbl[term, "F value"],
          pvalue   = tbl[term, "Pr(>F)"],
          stringsAsFactors = FALSE
        )
      }
    }, error = function(e) {
      message("Warning: Two-way ANOVA failed for metric '", k, "': ", conditionMessage(e))
    })
  }

  if (length(tw_rows) > 0) {
    write.csv(do.call(rbind, tw_rows), tw_csv, row.names = FALSE)
    message("Wrote two-way ANOVA CSV: ", tw_csv)
  }
}

message("Alpha diversity analysis complete.")
