#' @name plot_spectrum
#' @title Plot a Single Mass Spectrum with Optional Peak Segments
#'
#' @description
#' Draws a single mass spectrum as a line plot using `ggplot2`. When a peak
#' table is supplied, vertical segments are drawn from zero to each peak
#' intensity. Requires the suggested package `ggplot2`.
#'
#' @param spectrum A numeric matrix or data frame with m/z in the first column
#'   and intensity in the second column.
#' @param peaks Optional peak matrix or data frame. When supplied, vertical peak
#'   segments are drawn from zero to each peak intensity.
#' @param xlim Optional numeric vector of length two limiting the x-axis range.
#' @param ylim Optional numeric vector of length two limiting the y-axis range.
#' @param title Optional plot title.
#'
#' @return A `ggplot2` object.
#'
#' @seealso [plot_spectra()], [find_peaks()]
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   x <- seq(1000, 2000, length.out = 1000)
#'   y <- dnorm(x, mean = 1500, sd = 25) * 100
#'   spectrum <- data.frame(mz = x, intensity = y)
#'
#'   plot_spectrum(spectrum, title = "Example spectrum")
#' }
#'
#' @export
plot_spectrum <- function(spectrum,
                          peaks = NULL,
                          xlim = NULL,
                          ylim = NULL,
                          title = NULL) {
  .require_ggplot2()

  if (!is.matrix(spectrum) && !is.data.frame(spectrum)) {
    stop(
      "'spectrum' must be a numeric matrix or data frame.",
      call. = FALSE
    )
  }

  if (ncol(spectrum) < 2L) {
    stop(
      "'spectrum' must have at least two columns: m/z and intensity.",
      call. = FALSE
    )
  }

  df <- data.frame(
    mz = spectrum[, 1],
    intensity = spectrum[, 2],
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$mz, y = .data$intensity)) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::labs(
      title = title,
      x = "m/z",
      y = "Intensity"
    ) +
    ggplot2::theme_minimal()

  if (!is.null(peaks)) {
    if (!is.matrix(peaks) && !is.data.frame(peaks)) {
      stop(
        "'peaks' must be a numeric matrix or data frame, or NULL.",
        call. = FALSE
      )
    }

    if (ncol(peaks) >= 2L && nrow(peaks) > 0L) {
      pk <- data.frame(
        mz = peaks[, 1],
        intensity = peaks[, 2],
        stringsAsFactors = FALSE
      )

      p <- p + ggplot2::geom_segment(
        data = pk,
        ggplot2::aes(
          x = .data$mz,
          xend = .data$mz,
          y = 0,
          yend = .data$intensity
        ),
        inherit.aes = FALSE,
        linewidth = 0.6,
        color = "firebrick"
      )
    }
  }

  .with_axis_limits(p, xlim = xlim, ylim = ylim)
}

#' @name plot_spectra
#' @title Plot One or More Mass Spectra
#'
#' @description
#' Draws one or more mass spectra as overlaid line plots, colored by sample.
#' Accepts a single spectrum or a named list of spectra. Requires the suggested
#' package `ggplot2`.
#'
#' @param spectra A spectrum matrix/data frame or named list of spectra.
#' @param xlim Optional numeric vector of length two limiting the x-axis range.
#' @param ylim Optional numeric vector of length two limiting the y-axis range.
#' @param title Optional plot title.
#'
#' @return A `ggplot2` object.
#'
#' @seealso [plot_spectrum()]
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   x <- seq(1000, 2000, length.out = 1000)
#'   spectra <- list(
#'     sample_1 = data.frame(mz = x, intensity = dnorm(x, 1400, 25) * 100),
#'     sample_2 = data.frame(mz = x, intensity = dnorm(x, 1600, 25) * 80)
#'   )
#'
#'   plot_spectra(spectra, title = "Example spectra")
#' }
#'
#' @export
plot_spectra <- function(spectra,
                         xlim = NULL,
                         ylim = NULL,
                         title = NULL) {
  .require_ggplot2()

  if (is.matrix(spectra) || is.data.frame(spectra)) {
    spectra <- list(spectrum = spectra)
  }

  if (!is.list(spectra) || length(spectra) < 1L) {
    stop(
      "'spectra' must be a spectrum object or a non-empty list of spectra.",
      call. = FALSE
    )
  }

  if (is.null(names(spectra))) {
    names(spectra) <- as.character(seq_along(spectra))
  }

  df <- do.call(
    rbind,
    lapply(names(spectra), function(nm) {
      sp <- spectra[[nm]]
      data.frame(
        sample = nm,
        mz = sp[, 1],
        intensity = sp[, 2],
        stringsAsFactors = FALSE
      )
    })
  )

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$mz, y = .data$intensity, color = .data$sample)
  ) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::labs(
      title = title,
      x = "m/z",
      y = "Intensity",
      color = "Sample"
    ) +
    ggplot2::theme_minimal()

  .with_axis_limits(p, xlim = xlim, ylim = ylim)
}

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
