#' @name heatmap_matched_matrix
#' @title Heatmap of a Sample-by-Marker Matched-Peak Matrix
#'
#' @description
#' Draws a clustered heatmap of a sample-by-marker matched-peak matrix, such as
#' the `detected_matrix` or `delta_mz_matrix` returned by
#' [build_matched_matrix()], with optional group annotation and a zero-centered
#' diverging palette. Requires the suggested package `pheatmap`.
#'
#' @param matched_matrix A numeric matrix with samples as rows and markers as
#'   columns, such as the `detected_matrix` or `delta_mz_matrix` element of
#'   [build_matched_matrix()].
#' @param row_cluster Logical indicating whether to cluster rows. The default
#'   is `TRUE`.
#' @param col_cluster Logical indicating whether to cluster columns. The
#'   default is `TRUE`.
#' @param groups Optional group labels for sample annotation. When supplied,
#'   group colors use a fixed `colorspace` Viridis palette.
#' @param title Optional plot title.
#' @param center_at_zero Logical. When `TRUE` (default), uses a diverging
#'   blue-white-red palette with limits symmetric around `0`.
#' @param hide_rownames Logical. When `TRUE`, sample (row) labels are hidden.
#'   The default is `FALSE`.
#' @param hide_colnames Logical. When `TRUE`, marker (column) labels are hidden.
#'   The default is `FALSE`.
#'
#' @details
#' Missing values (`NA`) are shown as grey cells via `pheatmap`'s `na_col`.
#' When `center_at_zero = TRUE`, finite values are mapped with equal positive
#' and negative limits so that `0` sits at the neutral center color. This is
#' especially useful for `delta_mz_matrix` matrices.
#' When the matrix contains many `NA` values, consider setting
#' `row_cluster = FALSE` or `col_cluster = FALSE` because clustering may be
#' unstable with sparse data.
#'
#' @return A `pheatmap` object.
#'
#' @seealso [build_matched_matrix()]
#'
#' @examples
#' if (requireNamespace("pheatmap", quietly = TRUE)) {
#'   set.seed(1)
#'   m <- matrix(rnorm(40), nrow = 5)
#'   rownames(m) <- paste0("sample_", 1:5)
#'   colnames(m) <- paste0("mz_", 1:8)
#'
#'   heatmap_matched_matrix(m, title = "Example heatmap")
#' }
#'
#' @export
heatmap_matched_matrix <- function(matched_matrix,
                                   row_cluster = TRUE,
                                   col_cluster = TRUE,
                                   groups = NULL,
                                   title = "Matched peaks heatmap",
                                   center_at_zero = TRUE,
                                   hide_rownames = FALSE,
                                   hide_colnames = FALSE) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop(
      "Package 'pheatmap' is required for heatmap_matched_matrix(). Install it with install.packages('pheatmap').",
      call. = FALSE
    )
  }

  if (!is.matrix(matched_matrix) || !is.numeric(matched_matrix)) {
    stop(
      "'matched_matrix' must be a numeric matrix.",
      call. = FALSE
    )
  }

  mat <- matched_matrix

  if (nrow(mat) < 1L || ncol(mat) < 1L) {
    stop(
      "'matched_matrix' must have at least one row and one column.",
      call. = FALSE
    )
  }

  has_na <- anyNA(mat)

  if (!has_na) {
    if (any(!is.finite(mat))) {
      stop(
        "'matched_matrix' must contain only finite values.",
        call. = FALSE
      )
    }
  } else if (any(!is.na(mat) & !is.finite(mat))) {
    stop(
      "'matched_matrix' must contain only finite or NA values.",
      call. = FALSE
    )
  }

  if (is.null(rownames(mat))) {
    rownames(mat) <- paste0("sample_", seq_len(nrow(mat)))
  }

  if (is.null(colnames(mat))) {
    colnames(mat) <- paste0("feature_", seq_len(ncol(mat)))
  }

  if (!is.logical(row_cluster) || length(row_cluster) != 1L || is.na(row_cluster)) {
    stop(
      "'row_cluster' must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!is.logical(col_cluster) || length(col_cluster) != 1L || is.na(col_cluster)) {
    stop(
      "'col_cluster' must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!is.logical(center_at_zero) ||
      length(center_at_zero) != 1L ||
      is.na(center_at_zero)) {
    stop(
      "'center_at_zero' must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!is.logical(hide_rownames) ||
      length(hide_rownames) != 1L ||
      is.na(hide_rownames)) {
    stop(
      "'hide_rownames' must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!is.logical(hide_colnames) ||
      length(hide_colnames) != 1L ||
      is.na(hide_colnames)) {
    stop(
      "'hide_colnames' must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (nrow(mat) < 2L) {
    row_cluster <- FALSE
  }

  if (ncol(mat) < 2L) {
    col_cluster <- FALSE
  }

  annotation_row <- NULL

  if (!is.null(groups)) {
    if (length(groups) != nrow(mat)) {
      stop(
        "'groups' must have the same length as the number of rows in 'matched_matrix'.",
        call. = FALSE
      )
    }

    group_values <- if (!is.null(names(groups))) {
      if (!all(rownames(mat) %in% names(groups))) {
        stop(
          "When 'groups' is named, its names must include all row names in 'matched_matrix'.",
          call. = FALSE
        )
      }
      groups[rownames(mat)]
    } else {
      groups
    }

    annotation_row <- data.frame(
      group = factor(group_values),
      stringsAsFactors = FALSE
    )
    rownames(annotation_row) <- rownames(mat)
  }

  annotation_colors <- NULL
  if (!is.null(annotation_row)) {
    annotation_colors <- .group_annotation_colors(annotation_row$group)
  }

  pheatmap_args <- list(
    mat = mat,
    cluster_rows = row_cluster,
    cluster_cols = col_cluster,
    annotation_row = annotation_row,
    annotation_colors = annotation_colors,
    annotation_names_row = FALSE,
    main = title,
    show_rownames = !hide_rownames,
    show_colnames = !hide_colnames,
    silent = FALSE
  )

  if (has_na) {
    pheatmap_args$na_col <- "grey90"
  }

  if (isTRUE(center_at_zero)) {
    zero_scale <- .zero_centered_heatmap_scale(mat)
    if (!is.null(zero_scale)) {
      pheatmap_args$breaks <- zero_scale$breaks
      pheatmap_args$color <- zero_scale$color
    }
  }

  do.call(pheatmap::pheatmap, pheatmap_args)
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
      "heatmap_matched_matrix(). Install it with install.packages('colorspace').",
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
