#!/usr/bin/env Rscript
# sensitivity_rarefaction.R
# Assess alpha diversity sensitivity to rarefaction depth by computing metrics
# at multiple fractions of min_library_size.
# Only meaningful under --normalization rarefaction; main.nf gates this process.

suppressPackageStartupMessages({
  library(optparse)
  library(vegan)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Resolve script directory (same pattern as alpha_diversity.R)
# ---------------------------------------------------------------------------
script_dir <- tryCatch(
  dirname(normalizePath(sys.frames()[[1]]$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_flag <- grep("^--file=", args, value = TRUE)
    if (length(file_flag))
      dirname(normalizePath(sub("^--file=", "", file_flag)))
    else "."
  }
)
source(file.path(script_dir, "load_feature_table.R"))

# ---------------------------------------------------------------------------
# Option parsing
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--feature_table",      type = "character", default = NULL,
              help = "Path to feature table [required]"),
  make_option("--meta_table",         type = "character", default = NULL,
              help = "Path to sample metadata CSV [required]"),
  make_option("--taxonomy_table",     type = "character", default = "",
              help = "Path to taxonomy TSV (required for tsv/gtdb)"),
  make_option("--tree_file",          type = "character", default = "",
              help = "Newick tree file (required for FaithsPD)"),
  make_option("--input_format",       type = "character", default = "biom",
              help = "Input format: biom | tsv | gtdb [default: biom]"),
  make_option("--taxon_rank",         type = "character", default = "Feature",
              help = "Taxonomic level [default: Feature]"),
  make_option("--group",              type = "character", default = "",
              help = "Metadata grouping column(s) (optional)"),
  make_option("--alpha_metrics",      type = "character",
              default = "Richness,Shannon,PielouEvenness",
              help = "Comma-separated metrics [default: Richness,Shannon,PielouEvenness]"),
  make_option("--min_library_size",   type = "integer",   default = 5000L,
              help = "Reference minimum library size [default: 5000]"),
  make_option("--sensitivity_depths", type = "character", default = "0.5,0.75,1.0",
              help = "Depth fractions to test [default: 0.5,0.75,1.0]"),
  make_option("--label",              type = "character", default = "analysis",
              help = "Label used in output filenames [default: analysis]"),
  make_option("--output_dir",         type = "character", default = ".",
              help = "Output directory [default: .]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------
if (is.null(opt$feature_table) || opt$feature_table == "")
  stop("--feature_table is required.")
if (is.null(opt$meta_table) || opt$meta_table == "")
  stop("--meta_table is required.")
if (!file.exists(opt$feature_table))
  stop("Feature table not found: ", opt$feature_table)
if (!file.exists(opt$meta_table))
  stop("Metadata file not found: ", opt$meta_table)

valid_metrics <- c("Richness", "Shannon", "GiniSimpson", "InvSimpson",
                   "FisherAlpha", "Chao1", "PielouEvenness", "FaithsPD")
requested_metrics <- trimws(strsplit(opt$alpha_metrics, ",")[[1]])
bad_metrics <- setdiff(requested_metrics, valid_metrics)
if (length(bad_metrics) > 0)
  stop("Unknown metrics: ", paste(bad_metrics, collapse = ", "),
       ". Supported: ", paste(valid_metrics, collapse = ", "))

fractions <- suppressWarnings(
  as.numeric(trimws(strsplit(opt$sensitivity_depths, ",")[[1]]))
)
if (any(is.na(fractions)) || any(fractions <= 0))
  stop("--sensitivity_depths must be comma-separated positive numeric fractions")

if (!dir.exists(opt$output_dir)) {
  message("Creating output directory: ", opt$output_dir)
  dir.create(opt$output_dir, recursive = TRUE)
}

# ---------------------------------------------------------------------------
# Load phylogenetic tree (required for FaithsPD only)
# ---------------------------------------------------------------------------
tree <- NULL
tree_path <- opt$tree_file
if (!is.null(tree_path) && nchar(tree_path) > 0 &&
    tree_path != "NO_FILE" && file.exists(tree_path)) {
  message("Loading phylogenetic tree: ", tree_path)
  tree <- tryCatch(
    ape::read.tree(tree_path),
    error = function(e) stop("Failed to read tree: ", conditionMessage(e))
  )
}
if ("FaithsPD" %in% requested_metrics && is.null(tree))
  stop("FaithsPD metric requires --tree_file.")

# ---------------------------------------------------------------------------
# Load feature table
# ---------------------------------------------------------------------------
message("Loading feature table (format=", opt$input_format, ") ...")
tax_tbl_arg <- if (opt$taxonomy_table != "") opt$taxonomy_table else NULL
loaded      <- load_feature_table(
  feature_table  = opt$feature_table,
  input_format   = opt$input_format,
  taxonomy_table = tax_tbl_arg
)
abund_table      <- loaded$abund_table
feature_taxonomy <- loaded$feature_taxonomy

# ---------------------------------------------------------------------------
# Load metadata
# ---------------------------------------------------------------------------
message("Loading metadata: ", opt$meta_table)
meta_table <- tryCatch({
  dt <- fread(opt$meta_table, header = TRUE, check.names = FALSE)
  df <- as.data.frame(dt, stringsAsFactors = FALSE)
  rownames(df) <- df[[1]]
  df[, -1, drop = FALSE]
}, error = function(e) stop("Failed to read metadata: ", conditionMessage(e)))

# ---------------------------------------------------------------------------
# Align samples
# ---------------------------------------------------------------------------
common_samps <- intersect(rownames(abund_table), rownames(meta_table))
if (length(common_samps) < 2)
  stop("Fewer than 2 overlapping samples between feature table and metadata.")

abund_table      <- abund_table[common_samps, , drop = FALSE]
meta_table       <- meta_table[common_samps, , drop = FALSE]
abund_table      <- abund_table[, colSums(abund_table) > 0, drop = FALSE]
feature_taxonomy <- feature_taxonomy[colnames(abund_table), , drop = FALSE]

# ---------------------------------------------------------------------------
# Collate at taxonomic level
# ---------------------------------------------------------------------------
taxon_rank <- opt$taxon_rank
message("Collating features at taxonomic level: ", taxon_rank)

if (taxon_rank == "Feature") {
  collated <- abund_table
} else {
  lvl_list <- unique(feature_taxonomy[, taxon_rank])
  collated  <- NULL
  for (i in lvl_list) {
    feat_idx <- rownames(feature_taxonomy)[feature_taxonomy[, taxon_rank] == i]
    feat_idx <- intersect(feat_idx, colnames(abund_table))
    if (length(feat_idx) == 0) next
    col_name <- if (is.na(i) || i == "") "__Unknowns__" else i
    tmp <- data.frame(
      rowSums(abund_table[, feat_idx, drop = FALSE]),
      check.names = FALSE
    )
    colnames(tmp) <- col_name
    collated <- if (is.null(collated)) tmp else cbind(collated, tmp)
  }
  if (is.null(collated))
    stop("No features could be collated at level '", taxon_rank, "'.")
}
collated <- as.matrix(collated)
storage.mode(collated) <- "numeric"

# ---------------------------------------------------------------------------
# Grouping
# ---------------------------------------------------------------------------
if (opt$group != "") {
  grp_cols <- trimws(strsplit(opt$group, ",")[[1]])
  missing  <- setdiff(grp_cols, colnames(meta_table))
  if (length(missing) > 0)
    stop("--group references columns not in metadata: ",
         paste(missing, collapse = ", "))
  if (length(grp_cols) == 1) {
    group_vec <- as.factor(as.character(meta_table[[grp_cols]]))
  } else {
    group_vec <- as.factor(
      do.call(paste, c(meta_table[, grp_cols, drop = FALSE], sep = " "))
    )
  }
  names(group_vec) <- rownames(meta_table)
} else {
  group_vec <- setNames(
    factor(rep("All", nrow(meta_table))),
    rownames(meta_table)
  )
}

# ---------------------------------------------------------------------------
# Faith's PD helper
# ---------------------------------------------------------------------------
compute_faiths_pd <- function(mat, tr) {
  if (!requireNamespace("picante", quietly = TRUE))
    stop("Package 'picante' is required for FaithsPD.")
  common      <- intersect(colnames(mat), tr$tip.label)
  if (length(common) < 2)
    stop("Fewer than 2 features shared between feature table and tree.")
  mat_pruned  <- mat[, common, drop = FALSE]
  tree_pruned <- ape::keep.tip(tr, common)
  pd_res      <- picante::pd(mat_pruned, tree_pruned, include.root = TRUE)
  setNames(pd_res$PD, rownames(mat))
}

# ---------------------------------------------------------------------------
# Sensitivity loop over depth fractions
# ---------------------------------------------------------------------------
message("Running sensitivity analysis over fractions: ",
        paste(fractions, collapse = ", "))

all_rows <- list()

for (f in fractions) {
  depth_f <- max(1L, floor(f * opt$min_library_size))
  message("  Fraction ", f, " => rarefaction depth ", depth_f)

  lib_sz <- rowSums(collated)
  keep   <- lib_sz >= depth_f
  n_keep <- sum(keep)

  if (n_keep < 2) {
    message("  Warning: Only ", n_keep, " sample(s) have >= ", depth_f,
            " reads at fraction ", f, " — skipping.")
    next
  }

  mat_f <- collated[keep, , drop = FALSE]
  mat_f <- mat_f[, colSums(mat_f) > 0, drop = FALSE]

  # Pre-compute Shannon once if needed by Pielou
  sh_f <- NULL

  for (metric in requested_metrics) {
    vals <- tryCatch({
      switch(metric,
        Richness      = vegan::rarefy(mat_f, depth_f),
        Shannon       = { sh_f <<- vegan::diversity(mat_f, "shannon"); sh_f },
        GiniSimpson   = vegan::diversity(mat_f, "simpson"),
        InvSimpson    = vegan::diversity(mat_f, "invsimpson"),
        FisherAlpha   = vegan::fisher.alpha(mat_f),
        Chao1         = vegan::estimateR(mat_f)["S.chao1", ],
        PielouEvenness = {
          if (is.null(sh_f))
            sh_f <<- vegan::diversity(mat_f, "shannon")
          sh_f / log(vegan::specnumber(mat_f))
        },
        FaithsPD      = compute_faiths_pd(mat_f, tree)
      )
    }, error = function(e) {
      message("  Warning: '", metric, "' failed at depth ", depth_f,
              ": ", conditionMessage(e))
      setNames(rep(NA_real_, nrow(mat_f)), rownames(mat_f))
    })

    samp_names <- rownames(mat_f)
    all_rows[[length(all_rows) + 1]] <- data.frame(
      sample            = samp_names,
      measure           = metric,
      value             = as.numeric(vals[samp_names]),
      Groups            = as.character(group_vec[samp_names]),
      depth_fraction    = f,
      rarefaction_depth = depth_f,
      stringsAsFactors  = FALSE
    )
  }
}

# ---------------------------------------------------------------------------
# Write output CSV
# ---------------------------------------------------------------------------
if (length(all_rows) == 0) {
  message("Warning: No results produced — all fractions were skipped.")
  result_df <- data.frame(
    sample            = character(),
    measure           = character(),
    value             = numeric(),
    Groups            = character(),
    depth_fraction    = numeric(),
    rarefaction_depth = integer(),
    stringsAsFactors  = FALSE
  )
} else {
  result_df <- do.call(rbind, all_rows)
  rownames(result_df) <- NULL
}

out_csv <- file.path(opt$output_dir,
  paste0("Diversity_sensitivity_", taxon_rank, "_", opt$label, ".csv"))
write.csv(result_df, out_csv, row.names = FALSE)
message("Wrote sensitivity CSV: ", out_csv)
message("Sensitivity analysis complete.")
