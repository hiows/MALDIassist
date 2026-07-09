#' @name build_kde_spectrum
#' @title Build a Gaussian KDE-Smoothed Spectrum
#'
#' @description
#' Builds a Gaussian kernel density estimate (KDE) representation of a mass
#' spectrum. The returned object contains the KDE-smoothed spectrum evaluated on
#' the observed m/z grid, together with the reusable KDE closure and the
#' bandwidth used. This prepares data that matches the Gaussian KDE-based peak
#' detection performed by [find_peaks()] and [find_peaks_spectra()].
#'
#' @param spectrum A numeric matrix or data frame with at least two columns.
#'   The first column must contain strictly increasing x-axis values, typically
#'   m/z values, and the second column must contain the corresponding
#'   intensities.
#'
#' @param bw A positive numeric scalar specifying the bandwidth of the Gaussian
#'   KDE. The default is the median interval between adjacent x-axis values.
#'
#' @details
#' The KDE closure is created with [get_gauss_kde()] using derivative order
#' `0`. The smoothed spectrum is obtained by evaluating this closure on the
#' observed x-axis grid, so the output spectrum has the same x-axis values as
#' the input. The column names of the input spectrum are preserved in the
#' output spectrum. When the input has no usable column names, `x` and `y` are
#' used as a fallback.
#'
#' @return
#' A list with three elements:
#'
#' - `spectrum`: A two-column `data.frame` containing the observed x-axis values
#'   and the KDE-smoothed intensities. The column names are inherited from
#'   `spectrum`.
#' - `gauss_kde`: The Gaussian KDE closure returned by [get_gauss_kde()], which
#'   can be re-evaluated at arbitrary x-axis values.
#' - `bw`: The bandwidth used to build the KDE.
#'
#' @seealso [get_gauss_kde()], [find_peaks()], [find_peaks_spectra()]
#'
#' @examples
#' x <- seq(1000, 2000, length.out = 2000)
#' y <- dnorm(x, mean = 1400, sd = 20) * 100 +
#'      dnorm(x, mean = 1700, sd = 25) * 80
#' spectrum <- data.frame(mz = x, intensity = y)
#'
#' kde_spectrum <- build_kde_spectrum(spectrum, bw = 1)
#'
#' head(kde_spectrum$spectrum)
#' kde_spectrum$gauss_kde(c(1400, 1700))
#'
#' @export
build_kde_spectrum <- function(spectrum,
                               bw = stats::median(diff(spectrum[, 1]))) {

  # -------------------------------------------------------------------------
  # 1. Validate input spectrum
  # -------------------------------------------------------------------------
  if (!is.data.frame(spectrum) && !is.matrix(spectrum)) {
    stop(
      "'spectrum' must be a data frame or matrix.",
      call. = FALSE
    )
  }

  if (ncol(spectrum) < 2L) {
    stop(
      "'spectrum' must have at least two columns: x-axis values and intensities.",
      call. = FALSE
    )
  }

  if (nrow(spectrum) < 3L) {
    stop(
      "'spectrum' must contain at least three rows.",
      call. = FALSE
    )
  }

  x <- spectrum[, 1]
  y <- spectrum[, 2]

  if (!is.numeric(x) || !is.numeric(y)) {
    stop(
      "The first and second columns of 'spectrum' must be numeric.",
      call. = FALSE
    )
  }

  if (anyNA(x) || anyNA(y)) {
    stop(
      "The first and second columns of 'spectrum' must not contain NA values.",
      call. = FALSE
    )
  }

  if (any(!is.finite(x)) || any(!is.finite(y))) {
    stop(
      "The first and second columns of 'spectrum' must contain only finite values.",
      call. = FALSE
    )
  }

  if (is.unsorted(x, strictly = TRUE)) {
    stop(
      "The x-axis values in the first column of 'spectrum' must be strictly increasing.",
      call. = FALSE
    )
  }

  # -------------------------------------------------------------------------
  # 2. Validate bandwidth
  # -------------------------------------------------------------------------
  if (!is.numeric(bw) ||
      length(bw) != 1L ||
      is.na(bw) ||
      !is.finite(bw) ||
      bw <= 0) {
    stop(
      "'bw' must be a positive finite numeric scalar.",
      call. = FALSE
    )
  }

  # -------------------------------------------------------------------------
  # 3. Determine output-column names (preserve input names when available)
  # -------------------------------------------------------------------------
  input_colnames <- colnames(spectrum)

  if (is.null(input_colnames) || length(input_colnames) < 2L) {
    output_colnames <- c("x", "y")
  } else {
    output_colnames <- input_colnames[1:2]

    if (anyNA(output_colnames) || any(output_colnames == "")) {
      output_colnames <- c("x", "y")
    }
  }

  # -------------------------------------------------------------------------
  # 4. Build the Gaussian KDE closure and evaluate on the observed grid
  # -------------------------------------------------------------------------
  f <- get_gauss_kde(data = spectrum, bw = bw, d = 0)

  x_eval <- as.double(x)
  y_kde <- f(x_eval)

  kde_spectrum <- data.frame(
    x = x_eval,
    y = y_kde,
    stringsAsFactors = FALSE
  )

  colnames(kde_spectrum)[1:2] <- output_colnames

  # -------------------------------------------------------------------------
  # 5. Assemble result
  # -------------------------------------------------------------------------
  result <- list(
    spectrum = kde_spectrum,
    gauss_kde = f,
    bw = bw
  )

  return(result)
}


#' @rdname build_kde_spectrum
#'
#' @param spectra A two-column `data.frame` or numeric `matrix` representing a
#'   single spectrum, or a `list` of such objects. The first column must contain
#'   strictly increasing x-axis values and the second column must contain
#'   intensity values.
#' @param n_cores A positive integer specifying the number of worker processes
#'   used when `spectra` is a list. The default is `1L`, which preserves
#'   sequential behaviour.
#'
#' @details
#' `build_kde_spectra()` applies [build_kde_spectrum()] to a single spectrum or
#' to each spectrum in a list. When `bw` is `NULL`, the per-spectrum default
#' bandwidth of [build_kde_spectrum()] (the median interval between adjacent
#' x-axis values) is used for every spectrum. Supplying a non-`NULL` `bw`
#' applies the same bandwidth to all spectra.
#'
#' @return
#' `build_kde_spectrum()` returns a single KDE-spectrum list.
#' `build_kde_spectra()` returns a single KDE-spectrum list when `spectra` is a
#' single spectrum, or a named list of KDE-spectrum lists with the same list
#' structure and names as `spectra`.
#'
#' @examples
#' \dontrun{
#' spectra <- list(sample_1 = spectrum_1, sample_2 = spectrum_2)
#'
#' kde_spectra <- build_kde_spectra(
#'   spectra = spectra,
#'   bw = 1,
#'   n_cores = 2
#' )
#' }
#'
#' @export
build_kde_spectra <- function(spectra,
                              bw = NULL,
                              n_cores = 1L) {

  .assert_integer_scalar(
    x = n_cores,
    arg_name = "n_cores",
    lower_bound = 1L
  )

  build_args <- list()

  if (!is.null(bw)) {
    build_args$bw <- bw
  }

  is_single_spectrum <- is.data.frame(spectra) || is.matrix(spectra)

  if (is_single_spectrum) {
    return(
      do.call(
        build_kde_spectrum,
        c(list(spectrum = spectra), build_args)
      )
    )
  }

  if (!is.list(spectra)) {
    stop(
      "'spectra' must be a two-column data.frame, a two-column numeric ",
      "matrix, or a list containing such objects.",
      call. = FALSE
    )
  }

  if (length(spectra) == 0L) {
    return(spectra)
  }

  kde_list <- .parallel_lapply(
    X = seq_along(spectra),
    FUN = function(i) {
      do.call(
        build_kde_spectrum,
        c(list(spectrum = spectra[[i]]), build_args)
      )
    },
    n_cores = n_cores
  )

  names(kde_list) <- names(spectra)

  kde_list
}
