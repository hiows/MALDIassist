# Internal plotting helpers.

#' @name .require_ggplot2
#' @title Require ggplot2 for Plotting Functions
#' @keywords internal
#' @noRd
.require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for plotting. Install it with install.packages('ggplot2').",
      call. = FALSE
    )
  }
  invisible(getNamespace("ggplot2"))
}

#' @name .cluster_order
#' @title Compute Clustered Row or Column Order
#' @keywords internal
#' @noRd
.cluster_order <- function(x) {
  if (nrow(x) < 2L && ncol(x) < 2L) {
    return(list(rows = seq_len(nrow(x)), cols = seq_len(ncol(x))))
  }

  row_ord <- if (nrow(x) >= 2L) {
    stats::hclust(stats::dist(x))$order
  } else {
    seq_len(nrow(x))
  }

  col_ord <- if (ncol(x) >= 2L) {
    stats::hclust(stats::dist(t(x)))$order
  } else {
    seq_len(ncol(x))
  }

  list(rows = row_ord, cols = col_ord)
}

#' @name .with_axis_limits
#' @title Apply Optional xlim and ylim to a ggplot Object
#' @keywords internal
#' @noRd
.with_axis_limits <- function(p, xlim = NULL, ylim = NULL) {
  if (!is.null(xlim) || !is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(
      xlim = xlim,
      ylim = ylim,
      expand = FALSE
    )
  }

  p
}

#' @name .zero_centered_heatmap_scale
#' @title Symmetric breaks and diverging colors centered at zero
#' @keywords internal
#' @noRd
.zero_centered_heatmap_scale <- function(mat, n = 100L) {
  vals <- mat[is.finite(mat)]
  if (length(vals) == 0L) {
    return(NULL)
  }

  max_abs <- max(abs(vals), na.rm = TRUE)
  if (max_abs <= 0) {
    max_abs <- 1
  }

  breaks <- seq(-max_abs, max_abs, length.out = n + 1L)
  colors <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(n)

  list(breaks = breaks, color = colors)
}

#' @name .group_annotation_colors
#' @title Fixed Viridis colors for row group annotation in pheatmap
#' @keywords internal
#' @noRd
.group_annotation_colors <- function(group_values) {
  group_factor <- factor(group_values)
  levels <- levels(group_factor)
  n_levels <- length(levels)

  if (n_levels < 1L) {
    return(NULL)
  }

  if (!requireNamespace("colorspace", quietly = TRUE)) {
    stop(
      "Package 'colorspace' is required for group annotation colors in ",
      "heatmap_spectrum_matrix(). Install it with install.packages('colorspace').",
      call. = FALSE
    )
  }

  list(
    group = stats::setNames(
      colorspace::sequential_hcl(n_levels, palette = "Viridis"),
      levels
    )
  )
}
