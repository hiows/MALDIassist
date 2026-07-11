#' @name visualize_spectrum
#' @title Visualize a Single Mass Spectrum with Optional Peak Segments
#'
#' @description
#' Draws a single mass spectrum as a base-R line plot. When a peak table is
#' supplied, vertical segments are drawn from zero to each peak intensity, with
#' an option to annotate the most intense peaks.
#'
#' @param spectrum A numeric matrix or data frame with m/z in the first column
#'   and intensity in the second column.
#' @param peaks Optional peak matrix or data frame. When supplied, vertical peak
#'   segments are drawn from zero to each peak intensity.
#' @param interest_range Optional numeric vector giving an m/z range to restrict
#'   the spectrum (and peaks) to before plotting.
#' @param annotate_topN Logical. When `TRUE`, the most intense peaks are labelled
#'   with their m/z values. The default is `FALSE`.
#' @param topN Number of peaks to annotate when `annotate_topN = TRUE`.
#' @param xlim Optional numeric vector of length two limiting the x-axis range.
#' @param ylim Optional numeric vector of length two limiting the y-axis range.
#' @param main Optional plot title.
#' @param cex.axis Axis annotation scaling factor.
#' @param cex.lab Axis label scaling factor.
#' @param cex.main Title scaling factor.
#' @param lwd Line width of the spectrum trace.
#' @param col Line color of the spectrum trace.
#' @param peaks_lwd Line width of the peak segments.
#' @param peaks_col Color of the peak segments.
#' @param peaks_lty Line type of the peak segments.
#'
#' @return Invisibly `NULL`; called for the side effect of drawing a plot.
#'
#' @seealso [visualize_spectra()], [find_peaks()]
#'
#' @examples
#' x <- seq(1000, 2000, length.out = 1000)
#' y <- dnorm(x, mean = 1500, sd = 25) * 100
#' spectrum <- data.frame(mz = x, intensity = y)
#'
#' visualize_spectrum(spectrum, main = "Example spectrum")
#'
#' @export
visualize_spectrum <- function(spectrum,
                               peaks = NULL,
                               interest_range = NULL,
                               annotate_topN = FALSE,
                               topN = 10,
                               xlim = NULL,
                               ylim = NULL,
                               main = NULL,
                               cex.axis = 1,
                               cex.lab = 1,
                               cex.main = 1.3,
                               lwd = 1,
                               col = "black",
                               peaks_lwd = 2,
                               peaks_col = "red",
                               peaks_lty = 3) {

  colnames(spectrum)[1:2] <- c("mz", "intensity")
  if (!is.null(interest_range)) {
    keep <- spectrum$mz >= min(interest_range) & spectrum$mz <= max(interest_range)
    spectrum <- spectrum[keep, , drop = FALSE]
  }

  x <- spectrum$mz
  y <- spectrum$intensity

  if (is.null(ylim)) {
    ylim <- c(min(0, min(y)), max(y))
    ylim <- sort(ylim)
    offset <- diff(ylim) * 0.1
    ylim <- c(min(ylim), max(ylim) + offset)
  }

  plot(
    x, y,
    type = "l",
    xlim = xlim,
    ylim = ylim,
    main = main,
    xlab = "m/z",
    ylab = "Intensity",
    lwd = lwd,
    cex.axis = cex.axis,
    cex.lab = cex.lab,
    cex.main = cex.main,
    col = col
  )

  if (!is.null(peaks)) {
    colnames(peaks)[1:2] <- c("mz", "intensity")
    if (!is.null(interest_range)) {
      keep <- peaks$mz >= min(interest_range) & peaks$mz <= max(interest_range)
      peaks <- peaks[keep, , drop = FALSE]
    }

    x_peaks <- peaks$mz
    y_peaks <- peaks$intensity
    segments(
      x_peaks, 0, x_peaks, y_peaks,
      lwd = peaks_lwd,
      lty = peaks_lty,
      col = peaks_col
    )

    if (annotate_topN) {
      peaks_topN <- peaks[order(peaks$intensity, decreasing = TRUE), ]
      peaks_topN <- head(peaks_topN, topN)
      text(
        peaks_topN$mz,
        peaks_topN$intensity,
        labels = round(peaks_topN$mz, 2),
        col = "blue",
        font = 2,
        adj = c(0.5, -0.3)
      )
    }
  }

  invisible(NULL)
}


#' @name visualize_spectra
#' @title Visualize One or More Overlaid Mass Spectra
#'
#' @description
#' Draws one or more mass spectra as overlaid base-R line plots, colored by a
#' sequential palette. No legend is drawn.
#'
#' @param spectra A list of spectrum matrices/data frames, each with m/z in the
#'   first column and intensity in the second.
#' @param interest_range Optional numeric vector giving an m/z range to restrict
#'   the spectra to before plotting.
#' @param xlim Optional numeric vector of length two limiting the x-axis range.
#' @param ylim Optional numeric vector of length two limiting the y-axis range.
#' @param main Optional plot title.
#' @param cex.axis Axis annotation scaling factor.
#' @param cex.lab Axis label scaling factor.
#' @param cex.main Title scaling factor.
#' @param lwd Line width of the spectra traces.
#' @param col_platte Sequential `colorspace` palette used to color the traces.
#'
#' @return Invisibly `NULL`; called for the side effect of drawing a plot.
#'
#' @seealso [visualize_spectrum()]
#'
#' @examples
#' x <- seq(1000, 2000, length.out = 1000)
#' spectra <- list(
#'   sample_1 = data.frame(mz = x, intensity = dnorm(x, 1400, 25) * 100),
#'   sample_2 = data.frame(mz = x, intensity = dnorm(x, 1600, 25) * 80)
#' )
#'
#' visualize_spectra(spectra, main = "Example spectra")
#'
#' @export
visualize_spectra <- function(spectra,
                              interest_range = NULL,
                              xlim = NULL,
                              ylim = NULL,
                              main = NULL,
                              cex.axis = 1,
                              cex.lab = 1,
                              cex.main = 1.3,
                              lwd = 1.5,
                              col_platte = c("Viridis", "YlOrRd")) {

  if (!requireNamespace("colorspace", quietly = TRUE)) {
    stop(
      "Package 'colorspace' is required for visualize_spectra(). Install it with install.packages('colorspace').",
      call. = FALSE
    )
  }

  col_platte <- match.arg(col_platte)

  spectra <- lapply(spectra, function(spectrum) {
    colnames(spectrum)[1:2] <- c("mz", "intensity")
    if (!is.null(interest_range)) {
      keep <- spectrum$mz >= min(interest_range) & spectrum$mz <= max(interest_range)
      spectrum <- spectrum[keep, , drop = FALSE]
    }
    spectrum
  })

  pooled_y <- unlist(lapply(spectra, "[[", 2))

  if (is.null(ylim)) {
    ylim <- c(min(0, min(pooled_y)), max(pooled_y))
    ylim <- sort(ylim)
    offset <- diff(ylim) * 0.1
    ylim <- c(min(ylim), max(ylim) + offset)
  }

  cols <- rev(colorspace::sequential_hcl(n = length(spectra) + 1, col_platte))
  cols <- rev(cols[-1])

  x <- spectra[[1]]$mz
  y <- spectra[[1]]$intensity

  plot(
    x, y,
    type = "l",
    xlim = xlim,
    ylim = ylim,
    main = main,
    xlab = "m/z",
    ylab = "Intensity",
    lwd = lwd,
    cex.axis = cex.axis,
    cex.lab = cex.lab,
    cex.main = cex.main,
    col = cols[1]
  )
  for (i in seq_along(spectra)) {
    x <- spectra[[i]]$mz
    y <- spectra[[i]]$intensity

    lines(
      x, y,
      col = cols[i],
      lwd = lwd
    )
  }

  invisible(NULL)
}
