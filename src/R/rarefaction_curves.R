#!/usr/bin/env Rscript
# rarefaction_curves.R
# Compute rarefaction curves and write the underlying data as a CSV.
# No plots are produced here — see examples/plot_rarefaction_curves.R.

suppressPackageStartupMessages({
  library(optparse)
  library(vegan)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Resolve script directory so load_feature_table.R is found regardless of
# how this script is invoked (directly, via Nextflow, or via Rscript --file=)
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
  make_option("--feature_table",    type = "character", default = NULL,
              help = "Path to feature table (BIOM, TSV, or GTDB) [required]"),
  make_option("--meta_table",       type = "character", default = NULL,
              help = "Path to sample metadata CSV [required]"),
  make_option("--taxonomy_table",   type = "character", default = "",
              help = "Taxonomy TSV (required for tsv/gtdb formats)"),
  make_option("--input_format",     type = "character", default = "biom",
              help = "Input format: biom | tsv | gtdb [default: biom]"),
  make_option("--group",            type = "character", default = "",
              help = "Metadata column(s) to annotate curves (optional)"),
  make_option("--min_library_size", type = "integer",   default = 5000L,
              help = "Reference depth annotated in output CSV [default: 5000]"),
  make_option("--label",            type = "character", default = "analysis",
              help = "Label used in output filename [default: analysis]"),
  make_option("--output_dir",       type = "character", default = ".",
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

if (!dir.exists(opt$output_dir)) {
  message("Creating output directory: ", opt$output_dir)
  dir.create(opt$output_dir, recursive = TRUE)
}

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
abund_table <- loaded$abund_table   # samples x features (numeric matrix)

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
if (length(common_samps) < 1)
  stop("No overlapping samples between feature table and metadata.")

abund_table <- abund_table[common_samps, , drop = FALSE]
meta_table  <- meta_table[common_samps, , drop = FALSE]
abund_table <- abund_table[, colSums(abund_table) > 0, drop = FALSE]

# ---------------------------------------------------------------------------
# Library sizes
# ---------------------------------------------------------------------------
lib_sizes <- rowSums(abund_table)
message("Library sizes range: ", min(lib_sizes), " - ", max(lib_sizes))

# ---------------------------------------------------------------------------
# Rarefaction curves (vegan::rarecurve)
# step=200, tmax = maximum observed library size
# ---------------------------------------------------------------------------
message("Computing rarefaction curves (step=200) ...")
rc <- vegan::rarecurve(
  abund_table,
  step   = 200,
  tmax   = max(lib_sizes),
  sample = opt$min_library_size,
  label  = FALSE
)

# ---------------------------------------------------------------------------
# Convert to long data frame:
#   sample | depth | richness | Groups | min_library_size
# ---------------------------------------------------------------------------
rc_list <- lapply(seq_along(rc), function(i) {
  x      <- rc[[i]]
  depths <- attr(x, "Subsample")
  data.frame(
    sample           = common_samps[i],
    depth            = depths,
    richness         = as.numeric(x),
    stringsAsFactors = FALSE
  )
})
df_rc <- do.call(rbind, rc_list)

# Annotate with group column(s)
if (opt$group != "") {
  grp_cols <- trimws(strsplit(opt$group, ",")[[1]])
  missing  <- setdiff(grp_cols, colnames(meta_table))
  if (length(missing) > 0)
    stop("--group references columns not in metadata: ",
         paste(missing, collapse = ", "))
  if (length(grp_cols) == 1) {
    group_vec <- as.character(meta_table[[grp_cols]])
  } else {
    group_vec <- do.call(
      paste, c(meta_table[, grp_cols, drop = FALSE], sep = " ")
    )
  }
  names(group_vec) <- rownames(meta_table)
  df_rc$Groups <- group_vec[df_rc$sample]
} else {
  df_rc$Groups <- "All"
}

# Record the reference depth so downstream plots can draw the vertical line
df_rc$min_library_size <- opt$min_library_size

# ---------------------------------------------------------------------------
# Write CSV
# ---------------------------------------------------------------------------
out_csv <- file.path(
  opt$output_dir,
  paste0("rarefaction_data_", opt$label, ".csv")
)
write.csv(df_rc, out_csv, row.names = FALSE)
message("Wrote rarefaction data CSV: ", out_csv)
message("Rarefaction curves complete.")
