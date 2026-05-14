#!/usr/bin/env Rscript
# =============================================================================
# EXAMPLE SCRIPT — not called by the pipeline
# Visualises rarefaction curve data from rarefaction_data_{label}.csv
# produced by the pipeline's RAREFACTION_CURVES process.
#
# Usage:
#   Rscript plot_rarefaction_curves.R \
#     --rarefaction_csv results/rarefaction_data_analysis.csv \
#     --output rarefaction_curves.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(data.table)
})

option_list <- list(
  make_option("--rarefaction_csv", type = "character", default = NULL,
              help = "Path to rarefaction_data_{label}.csv [required]"),
  make_option("--output",          type = "character",
              default = "rarefaction_curves.pdf",
              help = "Output file path [default: rarefaction_curves.pdf]"),
  make_option("--width",  type = "double", default = 10,
              help = "Plot width in inches [default: 10]"),
  make_option("--height", type = "double", default = 6,
              help = "Plot height in inches [default: 6]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$rarefaction_csv) || !file.exists(opt$rarefaction_csv))
  stop("--rarefaction_csv is required and must exist.")

df_rc <- as.data.frame(
  fread(opt$rarefaction_csv, header = TRUE),
  stringsAsFactors = FALSE
)

min_lib <- if ("min_library_size" %in% colnames(df_rc))
  df_rc$min_library_size[1] else NA_real_

has_groups <- "Groups" %in% colnames(df_rc) &&
              length(unique(df_rc$Groups)) > 1

p <- ggplot(df_rc, aes(
    x      = depth,
    y      = richness,
    group  = sample,
    colour = if (has_groups) Groups else NULL
  )) +
  geom_line(alpha = 0.4) +
  labs(
    x       = "Sequencing depth (reads)",
    y       = "Expected species richness",
    colour  = if (has_groups) "Group" else NULL,
    caption = if (!is.na(min_lib))
      paste0("Vertical line = min_library_size (", min_lib, ")")
    else NULL
  ) +
  theme_bw() +
  theme(legend.position = if (has_groups) "right" else "none")

if (!is.na(min_lib)) {
  p <- p + geom_vline(
    xintercept = min_lib,
    linetype   = "dashed",
    colour     = "black",
    linewidth  = 0.8
  )
}

ggsave(opt$output, plot = p, width = opt$width, height = opt$height)
message("Wrote: ", opt$output)
