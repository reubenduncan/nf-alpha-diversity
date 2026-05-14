#!/usr/bin/env Rscript
# metric_correlation.R
# Compute Spearman correlations between alpha diversity metrics.
# Outputs a long-format CSV only — no plots.
# See examples/plot_metric_correlation.R for visualisation.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Option parsing
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--wide_csv",   type = "character", default = NULL,
              help = "Path to the wide-format diversity CSV [required]"),
  make_option("--group_col",  type = "character", default = "Groups",
              help = "Name of the grouping column (excluded from correlation) [default: Groups]"),
  make_option("--label",      type = "character", default = "analysis",
              help = "Label used in output filenames [default: analysis]"),
  make_option("--output_dir", type = "character", default = ".",
              help = "Output directory [default: .]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------
if (is.null(opt$wide_csv) || opt$wide_csv == "")
  stop("--wide_csv is required.")
if (!file.exists(opt$wide_csv))
  stop("Wide CSV not found: ", opt$wide_csv)

if (!dir.exists(opt$output_dir)) {
  message("Creating output directory: ", opt$output_dir)
  dir.create(opt$output_dir, recursive = TRUE)
}

# ---------------------------------------------------------------------------
# Load wide CSV
# ---------------------------------------------------------------------------
message("Loading wide diversity CSV: ", opt$wide_csv)
wide_df <- tryCatch({
  dt <- fread(opt$wide_csv, header = TRUE, check.names = FALSE)
  as.data.frame(dt, stringsAsFactors = FALSE)
}, error = function(e) stop("Failed to read wide CSV: ", conditionMessage(e)))

# ---------------------------------------------------------------------------
# Identify metric columns: numeric, excluding known non-metric columns
# ---------------------------------------------------------------------------
non_metric_cols <- unique(c("sample", "Groups", "Type", "normalization", opt$group_col))
numeric_cols    <- names(wide_df)[sapply(wide_df, is.numeric)]
metric_cols     <- setdiff(numeric_cols, non_metric_cols)

if (length(metric_cols) < 2)
  stop("Need at least 2 numeric metric columns for correlation. Found: ",
       paste(metric_cols, collapse = ", "))

message("Computing Spearman correlations for metrics: ",
        paste(metric_cols, collapse = ", "))

mat <- as.matrix(wide_df[, metric_cols, drop = FALSE])

# ---------------------------------------------------------------------------
# Spearman correlation matrix and p-values
# ---------------------------------------------------------------------------
rho_mat <- cor(mat, method = "spearman", use = "pairwise.complete.obs")

n_met  <- length(metric_cols)
pv_mat <- matrix(NA_real_, nrow = n_met, ncol = n_met,
                 dimnames = list(metric_cols, metric_cols))

for (i in seq_len(n_met)) {
  for (j in seq_len(n_met)) {
    if (i == j) {
      pv_mat[i, j] <- NA_real_
    } else {
      x  <- mat[, i]
      y  <- mat[, j]
      ok <- !is.na(x) & !is.na(y)
      if (sum(ok) >= 3) {
        pv_mat[i, j] <- tryCatch(
          cor.test(x[ok], y[ok], method = "spearman")$p.value,
          error = function(e) NA_real_
        )
      }
    }
  }
}

# ---------------------------------------------------------------------------
# Long-format CSV — lower triangle, excluding diagonal
# ---------------------------------------------------------------------------
corr_rows <- list()
for (i in seq_len(n_met)) {
  for (j in seq_len(i - 1)) {
    corr_rows[[length(corr_rows) + 1]] <- data.frame(
      metric1 = metric_cols[i],
      metric2 = metric_cols[j],
      rho     = rho_mat[i, j],
      pvalue  = pv_mat[i, j],
      stringsAsFactors = FALSE
    )
  }
}

if (length(corr_rows) == 0) {
  corr_df <- data.frame(metric1 = character(), metric2 = character(),
                         rho = numeric(), pvalue = numeric(),
                         stringsAsFactors = FALSE)
} else {
  corr_df <- do.call(rbind, corr_rows)
  rownames(corr_df) <- NULL
}

out_csv <- file.path(opt$output_dir,
  paste0("Diversity_correlation_", opt$label, ".csv"))
write.csv(corr_df, out_csv, row.names = FALSE)
message("Wrote correlation CSV: ", out_csv)
message("Metric correlation complete.")
