#!/usr/bin/env Rscript
# =============================================================================
# EXAMPLE SCRIPT — not called by the pipeline
# Generates publication-quality figures from alpha diversity pipeline outputs.
# =============================================================================
# plot_alpha_diversity.R
# Generates figures from alpha diversity pipeline CSV outputs.
#
# Main figure — boxplots with significance overlay (--sig_style):
#   full       — all significant pairwise brackets
#   cld        — compact letter display (recommended for many groups)
#   selective  — brackets only vs. a single reference group (--reference_group)
#
# Optional supplementary figures (generated when the relevant CSV is supplied):
#   --n_csv           adds per-group sample counts to x-axis labels
#   --wide_csv        generates a Spearman correlation matrix figure
#   --sensitivity_csv generates a sensitivity-vs-rarefaction-depth figure
#
# Note: rarefaction curve plots require the raw feature table and are produced
#       separately by the pipeline (rarefaction_curves.R / RAREFACTION_CURVES).
#
# Panel titles show the global omnibus test result and method, e.g.
# "Richness *** (ANOVA)". This is a single test across all groups — it can be
# significant even when no individual pairwise comparison survives correction.

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(data.table)
  library(patchwork)
  library(ggsignif)
  library(multcompView)
  library(scales)
})

# ---------------------------------------------------------------------------
# Option parsing
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--long_csv",        type = "character", default = NULL,
              help = "Long-format diversity CSV from alpha_diversity.R [required]"),
  make_option("--pairwise_csv",    type = "character", default = NULL,
              help = "Pairwise comparison CSV from alpha_diversity.R [required]"),
  make_option("--n_csv",           type = "character", default = NULL,
              help = "Per-group n CSV from alpha_diversity.R (adds n labels to x-axis)"),
  make_option("--wide_csv",        type = "character", default = NULL,
              help = "Wide-format diversity CSV — used to generate correlation matrix figure"),
  make_option("--sensitivity_csv", type = "character", default = NULL,
              help = "Sensitivity rarefaction CSV — used to generate sensitivity figure"),
  make_option("--metrics",         type = "character", default = "Richness,Shannon",
              help = "Comma-separated metrics to plot [default: Richness,Shannon]"),
  make_option("--sig_style",       type = "character", default = "cld",
              help = "Significance style: full | cld | selective [default: cld]"),
  make_option("--reference_group", type = "character", default = NULL,
              help = "Reference group for selective mode [default: first factor level]"),
  make_option("--tile_matrix",     action = "store_true", default = FALSE,
              help = "Add pairwise significance tile matrix beneath boxplots"),
  make_option("--sig_cutoff",      type = "double",   default = 0.05,
              help = "Adjusted p-value cutoff for significance display [default: 0.05]"),
  make_option("--group_col",       type = "character", default = "Groups",
              help = "Column name for grouping variable [default: Groups]"),
  make_option("--type_col",        type = "character", default = "Type",
              help = "Column name for jitter point shape variable [default: Type]"),
  make_option("--test_method",     type = "character", default = "anova",
              help = "Global omnibus test for panel title: anova | kruskal [default: anova]"),
  make_option("--output_file",     type = "character", default = "alpha_diversity.pdf",
              help = "Output file for main figure (.pdf, .png, .svg, .tiff)"),
  make_option("--width",           type = "double",   default = 12,
              help = "Figure width in inches [default: 12]"),
  make_option("--height",          type = "double",   default = 7,
              help = "Figure height in inches [default: 7]"),
  make_option("--dpi",             type = "integer",  default = 300L,
              help = "Resolution for raster outputs [default: 300]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------
if (is.null(opt$long_csv)     || !file.exists(opt$long_csv))
  stop("--long_csv is required and must exist.")
if (is.null(opt$pairwise_csv) || !file.exists(opt$pairwise_csv))
  stop("--pairwise_csv is required and must exist.")
if (!opt$sig_style %in% c("full", "cld", "selective"))
  stop("--sig_style must be one of: full, cld, selective")
if (!opt$test_method %in% c("anova", "kruskal"))
  stop("--test_method must be one of: anova, kruskal")
for (csv_opt in c("n_csv", "wide_csv", "sensitivity_csv")) {
  val <- opt[[csv_opt]]
  if (!is.null(val) && val != "" && !file.exists(val))
    stop("--", csv_opt, " file not found: ", val)
}

# ---------------------------------------------------------------------------
# Load required data
# ---------------------------------------------------------------------------
df_long     <- as.data.frame(fread(opt$long_csv,     header = TRUE))
df_pairwise <- as.data.frame(fread(opt$pairwise_csv, header = TRUE))

metrics <- trimws(strsplit(opt$metrics, ",")[[1]])
missing <- setdiff(metrics, unique(df_long$measure))
if (length(missing) > 0)
  stop("Metrics not found in long CSV: ", paste(missing, collapse = ", "))

group_col <- opt$group_col
type_col  <- opt$type_col

if (!group_col %in% colnames(df_long))
  stop("Group column '", group_col, "' not found in long CSV.")

has_type     <- type_col %in% colnames(df_long) && any(!is.na(df_long[[type_col]]))
group_levels <- levels(factor(df_long[[group_col]]))
n_groups     <- length(group_levels)

# ---------------------------------------------------------------------------
# Load optional CSVs
# ---------------------------------------------------------------------------
n_df <- if (!is.null(opt$n_csv) && opt$n_csv != "")
  as.data.frame(fread(opt$n_csv, header = TRUE)) else NULL

wide_df <- if (!is.null(opt$wide_csv) && opt$wide_csv != "")
  as.data.frame(fread(opt$wide_csv, header = TRUE)) else NULL

sensitivity_df <- if (!is.null(opt$sensitivity_csv) && opt$sensitivity_csv != "")
  as.data.frame(fread(opt$sensitivity_csv, header = TRUE)) else NULL

# ---------------------------------------------------------------------------
# Reference group (selective mode)
# ---------------------------------------------------------------------------
ref_group <- if (!is.null(opt$reference_group) && opt$reference_group != "")
  opt$reference_group else group_levels[1]
if (!ref_group %in% group_levels)
  stop("Reference group '", ref_group, "' not found. Available: ",
       paste(group_levels, collapse = ", "))

# ---------------------------------------------------------------------------
# Colour palette
# ---------------------------------------------------------------------------
group_colours <- setNames(hue_pal()(n_groups), group_levels)
type_shapes   <- c(16, 15, 17, 3, 4, 8)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Look up the actual test method used for a metric in the pairwise CSV.
# Falls back to opt$test_method if the column is absent (older CSV format).
metric_test_method <- function(df_pw, metric, fallback) {
  if ("actual_test_method" %in% colnames(df_pw)) {
    rows <- df_pw[df_pw$measure == metric &
                  !is.na(df_pw$actual_test_method) &
                  df_pw$actual_test_method != "", ]
    if (nrow(rows) > 0) return(rows$actual_test_method[1])
  }
  fallback
}

# Compute global omnibus test p-value; returns list(stars, label).
# The title note clarifies this is the global test, not pairwise.
overall_pval_result <- function(df_sub, group_col, method) {
  tryCatch({
    p <- if (method == "anova") {
      summary(aov(
        as.formula(paste("value ~", group_col)), data = df_sub
      ))[[1]][["Pr(>F)"]][1]
    } else {
      kruskal.test(
        as.formula(paste("value ~", group_col)), data = df_sub
      )$p.value
    }
    stars <- if (is.na(p)) ""
             else if (p <= 0.001) "***"
             else if (p <= 0.01)  "**"
             else if (p <= 0.05)  "*"
             else                 "ns"
    method_label <- if (method == "anova") "ANOVA" else "Kruskal-Wallis"
    list(stars = stars, method = method_label)
  }, error = function(e) list(stars = "", method = ""))
}

# Stars from adjusted p-value
padj_to_stars <- function(p) {
  if (is.na(p) || p > opt$sig_cutoff) return(NA_character_)
  if (p <= 0.001) "***" else if (p <= 0.01) "**" else "*"
}

# Compact letter display — maps group names to safe IDs to avoid separator issues
build_cld <- function(df_pw_metric, group_levels, sig_cutoff) {
  safe_ids <- paste0("G", seq_along(group_levels))
  id_map   <- setNames(safe_ids, group_levels)

  all_pairs <- combn(group_levels, 2, simplify = FALSE)

  pvec <- vapply(all_pairs, function(pair) {
    g1 <- pair[1]; g2 <- pair[2]
    hit <- df_pw_metric[
      (df_pw_metric$group1 == g1 & df_pw_metric$group2 == g2) |
      (df_pw_metric$group1 == g2 & df_pw_metric$group2 == g1), ]
    if (nrow(hit) == 0 || is.na(hit$padj[1])) 1.0 else hit$padj[1]
  }, numeric(1))

  names(pvec) <- vapply(all_pairs, function(p)
    paste(id_map[p[1]], id_map[p[2]], sep = "-"), character(1))

  tryCatch({
    raw <- multcompLetters(pvec, threshold = sig_cutoff)$Letters
    setNames(raw[safe_ids], group_levels)
  }, error = function(e) {
    message("CLD computation failed: ", conditionMessage(e), "\nFalling back to no letters.")
    setNames(rep("", length(group_levels)), group_levels)
  })
}

# ---------------------------------------------------------------------------
# Build boxplot for one metric
# ---------------------------------------------------------------------------
make_boxplot <- function(metric, df_long, df_pairwise, n_df,
                         group_col, type_col, has_type,
                         group_levels, group_colours, type_shapes,
                         sig_style, ref_group, test_method, sig_cutoff,
                         is_first_panel, is_last_panel) {

  df_sub <- df_long[df_long$measure == metric & !is.na(df_long$value), ]
  df_sub[[group_col]] <- factor(df_sub[[group_col]], levels = group_levels)
  df_pw  <- df_pairwise[df_pairwise$measure == metric, ]

  # Use actual method recorded in pairwise CSV when available
  effective_method <- metric_test_method(df_pw, metric, test_method)
  res <- overall_pval_result(df_sub, group_col, effective_method)

  # Title: "Richness *** (ANOVA)" — stars are from the global omnibus test
  panel_title <- paste0(
    metric,
    if (nchar(res$stars) > 0) paste0(" ", res$stars) else "",
    " (", res$method, ")"
  )

  y_min <- min(df_sub$value, na.rm = TRUE)
  y_max <- max(df_sub$value, na.rm = TRUE)
  y_rng <- y_max - y_min

  sig_pairs <- df_pw[!is.na(df_pw$padj) & df_pw$padj <= sig_cutoff, ]

  n_brax <- if (sig_style == "full") {
    nrow(sig_pairs)
  } else if (sig_style == "selective") {
    nrow(sig_pairs[sig_pairs$group1 == ref_group | sig_pairs$group2 == ref_group, ])
  } else {
    0
  }
  y_top_mult <- if (n_brax == 0) 0.08 else max(0.12, n_brax * 0.07)

  # X-axis labels with n if n_df provided
  if (!is.null(n_df)) {
    n_map  <- setNames(n_df$n, n_df$Groups)
    x_labs <- setNames(
      paste0(group_levels, "\n(n=", n_map[group_levels], ")"),
      group_levels
    )
  } else {
    x_labs <- setNames(group_levels, group_levels)
  }

  p <- ggplot(df_sub,
              aes(x      = .data[[group_col]],
                  y      = value,
                  fill   = .data[[group_col]],
                  colour = .data[[group_col]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.35, width = 0.6,
                 linewidth = 0.45, colour = "grey30") +
    scale_fill_manual(values = group_colours,
                      guide  = if (is_last_panel) "legend" else "none") +
    scale_colour_manual(values = group_colours, guide = "none") +
    scale_x_discrete(labels = x_labs) +
    scale_y_continuous(expand = expansion(mult = c(0.05, y_top_mult))) +
    labs(title = panel_title, x = group_col, y = "Observed Values", fill = "Groups") +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title         = element_text(hjust = 0.5, face = "bold"),
      legend.position    = "bottom"
    )

  # Jitter points coloured by group
  if (has_type) {
    df_sub[[type_col]] <- factor(df_sub[[type_col]])
    type_levels <- levels(df_sub[[type_col]])
    shapes_use  <- type_shapes[seq_len(min(length(type_levels), length(type_shapes)))]
    p <- p +
      geom_jitter(aes(shape = .data[[type_col]]),
                  width = 0.15, size = 1.8, alpha = 0.7) +
      scale_shape_manual(values = shapes_use, name = type_col,
                         guide  = if (is_first_panel) "legend" else "none")
  } else {
    p <- p + geom_jitter(width = 0.15, size = 1.8, alpha = 0.7, shape = 16)
  }

  # -------------------------------------------------------------------------
  # Significance overlay
  # -------------------------------------------------------------------------
  if (sig_style == "cld") {
    letters <- build_cld(df_pw, group_levels, sig_cutoff)
    y_lbl   <- y_max + y_rng * 0.04

    lbl_df        <- data.frame(label = as.character(letters), y = y_lbl,
                                stringsAsFactors = FALSE)
    lbl_df[[group_col]] <- factor(names(letters), levels = group_levels)

    p <- p +
      geom_text(data = lbl_df,
                aes(x = .data[[group_col]], y = y, label = label),
                inherit.aes = FALSE, size = 4, fontface = "bold", vjust = 0) +
      coord_cartesian(clip = "off")

  } else if (sig_style == "selective") {
    ref_sig <- sig_pairs[sig_pairs$group1 == ref_group | sig_pairs$group2 == ref_group, ]

    if (nrow(ref_sig) > 0) {
      g_idx   <- setNames(seq_along(group_levels), group_levels)
      ref_idx <- g_idx[ref_group]
      ref_sig$dist <- abs(g_idx[ref_sig$group2] - ref_idx +
                          g_idx[ref_sig$group1] - ref_idx)
      ref_sig <- ref_sig[order(ref_sig$dist), ]

      cmp_list   <- lapply(seq_len(nrow(ref_sig)),
                           function(i) c(ref_sig$group1[i], ref_sig$group2[i]))
      ann_labels <- vapply(ref_sig$padj, function(pp)
        if (pp <= 0.001) "***" else if (pp <= 0.01) "**" else "*", character(1))

      p <- p + geom_signif(comparisons = cmp_list, annotations = ann_labels,
                           step_increase = 0.08, tip_length = 0.01,
                           textsize = 3.5, vjust = 0.3, color = "black")
    }

  } else {
    # full — all significant pairwise brackets, ordered by span
    if (nrow(sig_pairs) > 0) {
      g_idx <- setNames(seq_along(group_levels), group_levels)
      sig_pairs$span <- abs(g_idx[sig_pairs$group2] - g_idx[sig_pairs$group1])
      sig_pairs      <- sig_pairs[order(sig_pairs$span), ]

      cmp_list   <- lapply(seq_len(nrow(sig_pairs)),
                           function(i) c(sig_pairs$group1[i], sig_pairs$group2[i]))
      ann_labels <- vapply(sig_pairs$padj, function(pp)
        if (pp <= 0.001) "***" else if (pp <= 0.01) "**" else "*", character(1))

      p <- p + geom_signif(comparisons = cmp_list, annotations = ann_labels,
                           step_increase = 0.05, tip_length = 0.01,
                           textsize = 3.0, vjust = 0.3, color = "black")
    }
  }

  p
}

# ---------------------------------------------------------------------------
# Build tile matrix for one metric
# ---------------------------------------------------------------------------
make_tile <- function(metric, df_pairwise, group_levels) {
  df_pw <- df_pairwise[df_pairwise$measure == metric, ]

  half <- df_pw[, c("group1", "group2", "significance")]
  half$significance[is.na(half$significance) | half$significance == ""] <- "ns"
  mirror <- setNames(half[, c("group2", "group1", "significance")],
                     c("group1", "group2", "significance"))
  sym <- rbind(half, mirror)

  sym$group1       <- factor(sym$group1, levels = group_levels)
  sym$group2       <- factor(sym$group2, levels = rev(group_levels))
  sym$significance <- factor(sym$significance, levels = c("***", "**", "*", "ns"))
  sym$txt_col      <- ifelse(sym$significance %in% c("***", "**"), "white", "grey35")

  tile_fills <- c("***" = "#08306b", "**" = "#2171b5", "*" = "#6baed6", "ns" = "#f0f0f0")

  ggplot(sym, aes(x = group1, y = group2, fill = significance)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = significance, colour = txt_col), size = 2.6) +
    scale_fill_manual(values = tile_fills, name = "Significance", drop = FALSE,
                      labels = c("***" = "p≤0.001", "**" = "p≤0.01",
                                 "*"   = "p≤0.05",  "ns" = "ns")) +
    scale_colour_identity() +
    labs(title = paste0(metric, " — pairwise"), x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          panel.grid  = element_blank(),
          plot.title  = element_text(hjust = 0.5, size = 9, face = "bold"),
          legend.position = "right") +
    coord_fixed()
}

# ---------------------------------------------------------------------------
# Sensitivity panel: one metric, placed beneath its boxplot
# ---------------------------------------------------------------------------
make_sensitivity_panel <- function(metric, sensitivity_df, group_levels, group_colours) {
  df_m <- sensitivity_df[sensitivity_df$measure == metric & !is.na(sensitivity_df$value), ]
  df_m$Groups <- factor(df_m$Groups, levels = group_levels)

  sum_rows <- list()
  for (g in levels(df_m$Groups)) {
    for (d in sort(unique(df_m$rarefaction_depth))) {
      sub <- df_m[df_m$Groups == g & df_m$rarefaction_depth == d, ]
      if (nrow(sub) == 0) next
      sum_rows[[length(sum_rows) + 1]] <- data.frame(
        Groups            = g,
        rarefaction_depth = d,
        median_val        = median(sub$value),
        q25               = quantile(sub$value, 0.25),
        q75               = quantile(sub$value, 0.75),
        stringsAsFactors  = FALSE
      )
    }
  }
  if (length(sum_rows) == 0) return(NULL)

  sum_df <- do.call(rbind, sum_rows)
  sum_df$Groups <- factor(sum_df$Groups, levels = group_levels)

  ggplot(sum_df, aes(x = rarefaction_depth, y = median_val,
                     colour = Groups, fill = Groups, group = Groups)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    scale_colour_manual(values = group_colours, guide = "none") +
    scale_fill_manual(values   = group_colours, guide = "none") +
    labs(x = "Rarefaction depth (reads)", y = metric,
         caption = "Median ± IQR") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          plot.caption     = element_text(size = 7, colour = "grey50"))
}

# ---------------------------------------------------------------------------
# Correlation matrix figure (Spearman, lower triangle)
# ---------------------------------------------------------------------------
make_correlation_plot <- function(wide_df, group_col) {
  non_metric <- unique(c("sample", "Groups", "Type", group_col))
  metric_cols <- setdiff(names(wide_df)[sapply(wide_df, is.numeric)], non_metric)

  if (length(metric_cols) < 2) {
    message("Fewer than 2 numeric metric columns found — skipping correlation plot.")
    return(NULL)
  }

  mat     <- as.matrix(wide_df[, metric_cols, drop = FALSE])
  rho_mat <- cor(mat, method = "spearman", use = "pairwise.complete.obs")
  n_met   <- length(metric_cols)

  plot_rows <- list()
  for (i in seq_len(n_met)) {
    for (j in seq_len(n_met)) {
      rho_val   <- if (i >= j) rho_mat[i, j] else NA_real_
      rho_label <- if (i == j) ""
                   else if (!is.na(rho_val)) sprintf("%.2f", rho_val)
                   else ""
      plot_rows[[length(plot_rows) + 1]] <- data.frame(
        metric1 = metric_cols[i], metric2 = metric_cols[j],
        rho = rho_val, rho_label = rho_label,
        stringsAsFactors = FALSE
      )
    }
  }
  plot_df <- do.call(rbind, plot_rows)
  plot_df$metric1 <- factor(plot_df$metric1, levels = metric_cols)
  plot_df$metric2 <- factor(plot_df$metric2, levels = rev(metric_cols))

  ggplot(plot_df, aes(x = metric1, y = metric2, fill = rho)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = rho_label), size = 3.5) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-1, 1),
                         na.value = "grey95", name = "Spearman ρ") +
    labs(x = NULL, y = NULL,
         title = "Alpha diversity metric correlations (Spearman ρ, lower triangle)") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid  = element_blank())
}

# ---------------------------------------------------------------------------
# Build main boxplot panels
# ---------------------------------------------------------------------------
n_metrics <- length(metrics)

box_plots <- lapply(seq_along(metrics), function(i) {
  make_boxplot(
    metric         = metrics[i],
    df_long        = df_long,
    df_pairwise    = df_pairwise,
    n_df           = n_df,
    group_col      = group_col,
    type_col       = type_col,
    has_type       = has_type,
    group_levels   = group_levels,
    group_colours  = group_colours,
    type_shapes    = type_shapes,
    sig_style      = opt$sig_style,
    ref_group      = ref_group,
    test_method    = opt$test_method,
    sig_cutoff     = opt$sig_cutoff,
    is_first_panel = (i == 1),
    is_last_panel  = (i == n_metrics)
  )
})

# ---------------------------------------------------------------------------
# Build per-metric supplementary panels (sensitivity)
# ---------------------------------------------------------------------------
if (!is.null(sensitivity_df)) {
  sensitivity_df <- sensitivity_df[sensitivity_df$measure %in% metrics, ]
  sens_panels <- lapply(metrics, make_sensitivity_panel,
                        sensitivity_df = sensitivity_df,
                        group_levels   = group_levels,
                        group_colours  = group_colours)
} else {
  sens_panels <- NULL
}

if (opt$tile_matrix) {
  tile_plots <- lapply(metrics, make_tile,
                       df_pairwise = df_pairwise, group_levels = group_levels)
}

# ---------------------------------------------------------------------------
# Compose main figure: stack panels per column, then arrange columns
# ---------------------------------------------------------------------------
col_plots <- lapply(seq_along(metrics), function(i) {
  stack   <- list(box_plots[[i]])
  heights <- 3
  if (opt$tile_matrix) {
    stack   <- c(stack, list(tile_plots[[i]]))
    heights <- c(heights, 2)
  }
  if (!is.null(sens_panels) && !is.null(sens_panels[[i]])) {
    stack   <- c(stack, list(sens_panels[[i]]))
    heights <- c(heights, 1)
  }
  if (length(stack) == 1) return(stack[[1]])
  Reduce(`/`, stack) + plot_layout(heights = heights)
})

fig <- wrap_plots(col_plots, nrow = 1)

fig <- fig +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# ---------------------------------------------------------------------------
# Save helper
# ---------------------------------------------------------------------------
save_figure <- function(fig, path, width, height, dpi) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("pdf", "svg")) {
    ggsave(path, fig, width = width, height = height, device = ext)
  } else if (ext %in% c("png", "jpg", "jpeg", "tiff")) {
    ggsave(path, fig, width = width, height = height, dpi = dpi)
  } else {
    fallback <- paste0(tools::file_path_sans_ext(path), ".pdf")
    message("Unrecognised extension '", ext, "' — saving as PDF: ", fallback)
    ggsave(fallback, fig, width = width, height = height, device = "pdf")
  }
  message("Saved: ", path)
}

stem <- tools::file_path_sans_ext(opt$output_file)
ext  <- tools::file_ext(opt$output_file)

save_figure(fig, opt$output_file, opt$width, opt$height, opt$dpi)

# ---------------------------------------------------------------------------
# Supplementary: correlation matrix figure
# ---------------------------------------------------------------------------
if (!is.null(wide_df)) {
  corr_fig  <- make_correlation_plot(wide_df, group_col)
  if (!is.null(corr_fig)) {
    corr_path <- paste0(stem, "_correlation.", ext)
    save_figure(corr_fig, corr_path, width = 7, height = 6, dpi = opt$dpi)
  }
}
