#' @name find_peaks_fast
#' @title Fast Local Peak Detection from a Mass Spectrum
#'
#' @description
#' Detects local peaks from a two-column mass spectrum using a C++ backend.
#' The function first identifies local maxima, including optional plateau
#' handling, and then keeps only the highest local maximum within each
#' \code{hws_peaks}-based window.
#'
#' @param data A data frame or matrix containing the mass spectrum. The first
#'   column is assumed to contain x values, usually m/z values, and the second
#'   column is assumed to contain y values, usually intensity values.
#' @param hws_peaks Numeric. Half-window size used to filter neighboring local
#'   maxima. For each local maximum, only the highest local maximum within
#'   \code{x_peak - hws_peaks} and \code{x_peak + hws_peaks} is retained.
#'   Default is \code{10}.
#' @param plateau Character. Method used to represent flat local maxima.
#'   One of \code{"middle"}, \code{"first"}, or \code{"last"}.
#'   \itemize{
#'     \item \code{"middle"}: use the middle point of a plateau.
#'     \item \code{"first"}: use the first point of a plateau.
#'     \item \code{"last"}: use the last point of a plateau.
#'   }
#' @param na_rm Logical. If \code{TRUE}, rows with missing values in the first
#'   two columns of \code{data} are removed before peak detection. Default is
#'   \code{TRUE}.
#'
#' @return A data frame with two columns. The column names are inherited from
#'   the first two columns of \code{data}. Each row represents one detected
#'   peak.
#'
#' @details
#' This function is an R wrapper around the C++ backend
#' \code{cpp_find_peaks_fast()}.
#'
#' The input \code{data} is converted to a two-column numeric matrix before
#' calling the C++ backend. The first column must be strictly increasing and
#' must not contain duplicated values.
#'
#' Local maxima are first detected from the intensity profile. Then, if multiple
#' local maxima exist within the same \code{hws_peaks}-based x-window, only the
#' local maximum with the highest intensity is retained.
#'
#' This function is intended for fast peak candidate detection from mass
#' spectrum-like one-dimensional signals.
#'
#' @examples
#' x <- seq(1000, 2000, length.out = 5000)
#' y <- dnorm(x, mean = 1400, sd = 20) * 100 +
#'      dnorm(x, mean = 1450, sd = 10) * 30 +
#'      dnorm(x, mean = 1700, sd = 25) * 80
#'
#' data <- data.frame(mz = x, intensity = y)
#'
#' peaks <- find_peaks_fast(
#'   data = data,
#'   hws_peaks = 10,
#'   plateau = "middle"
#' )
#'
#' head(peaks)
#'
#' @export
find_peaks_fast <- function(data,
                            hws_peaks = 10,
                            plateau = c("middle", "first", "last"),
                            na_rm = TRUE) {
  
  plateau <- match.arg(plateau)
  
  # ----------------------------------------
  # Validate data
  # ----------------------------------------
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("'data' must be a data frame or matrix.", call. = FALSE)
  }
  
  if (ncol(data) < 2L) {
    stop("'data' must have at least two columns: x and y.", call. = FALSE)
  }
  
  x <- data[, 1]
  y <- data[, 2]
  
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("The first and second columns of 'data' must be numeric.",
         call. = FALSE)
  }
  
  if (length(x) != length(y)) {
    stop("The first and second columns of 'data' must have the same length.",
         call. = FALSE)
  }
  
  if (length(x) < 3L) {
    stop("'data' must contain at least 3 points.", call. = FALSE)
  }
  
  # ----------------------------------------
  # Validate hws_peaks
  # ----------------------------------------
  if (!is.numeric(hws_peaks) ||
      length(hws_peaks) != 1L ||
      is.na(hws_peaks) ||
      !is.finite(hws_peaks)) {
    stop("'hws_peaks' must be a single finite numeric value.",
         call. = FALSE)
  }
  
  if (hws_peaks <= 0) {
    stop("'hws_peaks' must be greater than 0.", call. = FALSE)
  }
  
  # ----------------------------------------
  # Validate na_rm
  # ----------------------------------------
  if (!is.logical(na_rm) || length(na_rm) != 1L || is.na(na_rm)) {
    stop("'na_rm' must be a single non-missing logical value.",
         call. = FALSE)
  }
  
  # ----------------------------------------
  # Handle missing values
  # ----------------------------------------
  if (na_rm) {
    keep <- complete.cases(x, y)
    x <- x[keep]
    y <- y[keep]
  } else {
    if (anyNA(x) || anyNA(y)) {
      stop("'data' contains missing values. Use 'na_rm = TRUE' to remove them.",
           call. = FALSE)
    }
  }
  
  if (length(x) < 3L) {
    stop("After removing NA values, at least 3 points are required.",
         call. = FALSE)
  }
  
  if (any(!is.finite(x)) || any(!is.finite(y))) {
    stop("The first two columns of 'data' must contain finite numeric values.",
         call. = FALSE)
  }
  
  if (anyDuplicated(x)) {
    stop("'x' values in 'data' must not be duplicated.", call. = FALSE)
  }
  
  if (is.unsorted(x, strictly = TRUE)) {
    stop("'x' values in 'data' must be strictly increasing.", call. = FALSE)
  }
  
  # ----------------------------------------
  # Convert plateau option to C++ integer code
  # ----------------------------------------
  plateau_id <- switch(
    plateau,
    middle = 0L,
    first = 1L,
    last = 2L
  )
  
  # ----------------------------------------
  # Prepare input matrix
  # ----------------------------------------
  data_mat <- cbind(x, y)
  
  # ----------------------------------------
  # Detect peaks using C++
  # ----------------------------------------
  peaks <- cpp_find_peaks_fast(
    data = data_mat,
    hws_peaks = hws_peaks,
    plateau = plateau_id
  )
  
  peaks <- as.data.frame(peaks)
  colnames(peaks) <- colnames(data)[1:2]
  rownames(peaks) <- NULL
  
  peaks$type = "peak"
  
  return(peaks)
}

#' @rdname find_peaks_fast
#'
#' @param spectra A two-column `data.frame` or numeric `matrix` representing a
#'   single spectrum, or a `list` of such objects. The first column is assumed
#'   to contain x values, usually m/z values, and the second column is assumed
#'   to contain y values, usually intensity values.
#' @param n_cores A positive integer specifying the number of worker processes
#'   used when `spectra` is a list. The default is `1L`, which preserves
#'   sequential behaviour.
#'
#' @details
#' `find_peaks_spectra_fast()` applies [find_peaks_fast()] to a single spectrum
#' or to each spectrum in a list.
#'
#' @return
#' `find_peaks_fast()` returns a single peak table.
#' `find_peaks_spectra_fast()` returns a single peak table when `spectra` is a
#' single spectrum, or a named list of peak tables with the same list structure
#' and names as `spectra`.
#'
#' @seealso [find_peaks_spectra()], [preprocess_maldi_spectra()],
#'   [filter_peaks_spectra()]
#'
#' @examples
#' make_spectrum <- function() {
#'   x <- seq(1000, 2000, length.out = 2000)
#'   y <- dnorm(x, mean = 1400, sd = 20) * 100 +
#'        dnorm(x, mean = 1700, sd = 25) * 80
#'   data.frame(mz = x, intensity = y)
#' }
#'
#' spectra <- list(sample_1 = make_spectrum(), sample_2 = make_spectrum())
#'
#' peaks_list <- find_peaks_spectra_fast(spectra, hws_peaks = 10)
#'
#' lapply(peaks_list, head)
#'
#' @export
find_peaks_spectra_fast <- function(spectra,
                                    hws_peaks = 10,
                                    plateau = c("middle", "first", "last"),
                                    na_rm = TRUE,
                                    n_cores = 1L) {
  
  plateau <- match.arg(plateau)
  
  .assert_integer_scalar(
    x = n_cores,
    arg_name = "n_cores",
    lower_bound = 1L
  )
  
  find_args <- list(
    hws_peaks = hws_peaks,
    plateau = plateau,
    na_rm = na_rm
  )
  
  is_single_spectrum <- is.data.frame(spectra) || is.matrix(spectra)
  
  if (is_single_spectrum) {
    return(
      do.call(
        find_peaks_fast,
        c(list(data = spectra), find_args)
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
  
  peaks_list <- .parallel_lapply(
    X = seq_along(spectra),
    FUN = function(i) {
      do.call(
        find_peaks_fast,
        c(list(data = spectra[[i]]), find_args)
      )
    },
    n_cores = n_cores
  )
  
  names(peaks_list) <- names(spectra)
  
  peaks_list
}