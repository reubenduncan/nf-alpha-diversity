#!/usr/bin/env Rscript
# =============================================================================
# EXAMPLE SCRIPT — not called by the pipeline
# Visualises the metric correlation CSV produced by METRIC_CORRELATION as a
# Spearman rho heatmap (lower triangle).
#
# Usage:
#   Rscript plot_metric_correlation.R \
#     --correlation_csv results/Diversity_correlation_analysis.csv \
#     --output correlation_heatmap.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(data.table)
})

option_list <- list(
  make_option("--correlation_csv", type = "character", default = NULL,
              help = "Path to Diversity_correlation_{label}.csv [required]"),
  make_option("--label",   type = "character", default = "",
              help = "Label for plot title [default: inferred from filename]"),
  make_option("--output",  type = "character",
              default = "correlation_heatmap.pdf",
              help = "Output file [default: correlation_heatmap.pdf]"),
  make_option("--width",   type = "double", default = 7,
              help = "Width in inches [default: 7]"),
  make_option("--height",  type = "double", default = 6,
              help = "Height in inches [default: 6]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$correlation_csv) || !file.exists(opt$correlation_csv))
  stop("--correlation_csv is required and must exist.")

corr_df <- as.data.frame(
  fread(opt$correlation_csv, header = TRUE),
  stringsAsFactors = FALSE
)

if (!all(c("metric1", "metric2", "rho") %in% colnames(corr_df)))
  stop("CSV must have columns: metric1, metric2, rho")

# Reconstruct symmetric matrix for plotting
metrics <- unique(c(corr_df$metric1, corr_df$metric2))

plot_rows <- list()
for (m1 in metrics) {
  for (m2 in metrics) {
    if (m1 == m2) {
      rho_val   <- NA_real_
      rho_label <- ""
    } else {
      row <- corr_df[corr_df$metric1 == m1 & corr_df$metric2 == m2, ]
      if (nrow(row) == 0)
        row <- corr_df[corr_df$metric1 == m2 & corr_df$metric2 == m1, ]
      # Only show lower triangle
      i1 <- which(metrics == m1)
      i2 <- which(metrics == m2)
      if (i1 > i2 && nrow(row) > 0) {
        rho_val   <- row$rho[1]
        rho_label <- sprintf("%.2f", rho_val)
      } else {
        rho_val   <- NA_real_
        rho_label <- ""
      }
    }
    plot_rows[[length(plot_rows) + 1]] <- data.frame(
      metric1   = m1,
      metric2   = m2,
      rho       = rho_val,
      rho_label = rho_label,
      stringsAsFactors = FALSE
    )
  }
}

plot_df          <- do.call(rbind, plot_rows)
plot_df$metric1  <- factor(plot_df$metric1, levels = metrics)
plot_df$metric2  <- factor(plot_df$metric2, levels = rev(metrics))

lbl <- if (nchar(opt$label) > 0) opt$label else
  sub(".*Diversity_correlation_(.+)\\.csv$", "\\1", basename(opt$correlation_csv))

p <- ggplot(plot_df, aes(x = metric1, y = metric2, fill = rho)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = rho_label), size = 3.5) +
  scale_fill_gradient2(
    low      = "blue",
    mid      = "white",
    high     = "red",
    midpoint = 0,
    limits   = c(-1, 1),
    na.value = "grey95",
    name     = "Spearman rho"
  ) +
  labs(
    title = paste0("Alpha diversity metric correlations — ", lbl),
    x = NULL,
    y = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank()
  )

ggsave(opt$output, plot = p, width = opt$width, height = opt$height)
message("Wrote: ", opt$output)
