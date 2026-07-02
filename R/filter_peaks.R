#' @name filter_peaks
#' @title Filter Detected Peaks Using Intensity, Prominence, and Peak Strength
#' 
#' @description
#' Filters candidate peaks detected from a MALDI-TOF mass spectrum using
#' minimum thresholds for peak intensity, peak prominence, and peak strength.
#' Peak prominence is calculated using [calculate_peak_prominence()], and peak
#' strength is estimated using [estimate_peak_strength()].
#'
#' @param data A numeric matrix or data frame containing the mass spectrum.
#'   The first column must contain strictly increasing x-axis values, typically
#'   m/z values, and the second column must contain the corresponding
#'   intensities.
#'
#' @param peaks A numeric matrix or data frame containing candidate peaks.
#'   The first column must contain peak positions and the second column must
#'   contain the corresponding peak intensities. Additional columns are
#'   preserved in the returned object.
#'
#' @param cutoff_peak_intensity A non-negative finite numeric scalar specifying
#'   the minimum peak intensity. Peaks with intensities less than or equal to
#'   this threshold are removed. If `NULL`, the median absolute deviation
#'   (`MAD`) of the spectrum intensities is used.
#'
#' @param cutoff_peak_prominence A non-negative finite numeric scalar specifying
#'   the minimum peak prominence. Peaks with prominences less than or equal to
#'   this threshold are removed. If `NULL`, the median absolute deviation
#'   (`MAD`) of the spectrum intensities is used.
#'
#' @param cutoff_peak_strength A non-negative finite numeric scalar specifying
#'   the minimum peak-strength score. Peaks with strength values less than or
#'   equal to this threshold are removed. The default is `0.2`.
#'
#' @param k A positive finite numeric scalar passed to
#'   [estimate_peak_strength()]. This parameter controls the contribution of
#'   the noise-related term when calculating peak strength. The default is
#'   `1`.
#'
#' @param normalization_type Character string specifying the intensity
#'   transformation passed to [estimate_peak_strength()]. Available options
#'   are `"raw"`, `"sqrt"`, and `"log10"`.
#'
#' @details
#' A candidate peak is retained only when all three criteria are satisfied:
#'
#' \preformatted{
#' peak intensity  > cutoff_peak_intensity
#' peak prominence > cutoff_peak_prominence
#' peak strength   > cutoff_peak_strength
#' }
#'
#' Peaks with non-finite intensity, prominence, or strength values are removed.
#' Such values may occur, for example, when a boundary peak does not have a
#' well-defined neighboring valley.
#'
#' The default intensity and prominence thresholds are calculated from the
#' spectrum intensity vector rather than from the detected peak intensities.
#'
#' @return An object of the same class as `peaks`, containing only peaks that
#'   satisfy all filtering criteria. All original columns in `peaks` are
#'   preserved. If no peaks remain, an empty object with the same columns is
#'   returned.
#'
#' @seealso [filter_peaks_spectra()]
#'
#' @examples
#' x <- seq(1000, 2000, length.out = 2000)
#' y <- dnorm(x, mean = 1400, sd = 20) * 100 +
#'      dnorm(x, mean = 1700, sd = 25) * 80
#' processed_spectrum <- data.frame(mz = x, intensity = y)
#'
#' detected_peaks <- find_peaks_fast(processed_spectrum, hws_peaks = 10)
#'
#' filtered_peaks <- filter_peaks(
#'   data = processed_spectrum,
#'   peaks = detected_peaks,
#'   cutoff_peak_intensity = NULL,
#'   cutoff_peak_prominence = NULL,
#'   cutoff_peak_strength = 0.2,
#'   k = 1,
#'   normalization_type = "raw"
#' )
#'
#' filtered_peaks
#'
#' @export
filter_peaks <- function(
    data,
    peaks,
    cutoff_peak_intensity = NULL,
    cutoff_peak_prominence = NULL,
    cutoff_peak_strength = 0.2,
    k = 1,
    normalization_type = c("raw", "sqrt", "log10")) {
  
  # ----------------------------------------
  # Match arguments
  # ----------------------------------------
  normalization_type <- match.arg(normalization_type)
  
  # ----------------------------------------
  # Validate spectrum data
  # ----------------------------------------
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop(
      "'data' must be a numeric matrix or data frame.",
      call. = FALSE
    )
  }
  
  if (ncol(data) < 2L) {
    stop(
      "'data' must have at least two columns: x-axis values and intensities.",
      call. = FALSE
    )
  }
  
  if (nrow(data) < 3L) {
    stop(
      "'data' must contain at least three rows.",
      call. = FALSE
    )
  }
  
  x <- data[, 1]
  y <- data[, 2]
  
  if (!is.numeric(x) || !is.numeric(y)) {
    stop(
      "The first and second columns of 'data' must be numeric.",
      call. = FALSE
    )
  }
  
  if (anyNA(x) || anyNA(y) || any(!is.finite(x)) || any(!is.finite(y))) {
    stop(
      "The first and second columns of 'data' must contain only finite values.",
      call. = FALSE
    )
  }
  
  if (anyDuplicated(x)) {
    stop(
      "The x-axis values in 'data' must not contain duplicates.",
      call. = FALSE
    )
  }
  
  if (is.unsorted(x, strictly = TRUE)) {
    stop(
      "The x-axis values in 'data' must be strictly increasing.",
      call. = FALSE
    )
  }
  
  # ----------------------------------------
  # Validate peaks
  # ----------------------------------------
  if (!is.matrix(peaks) && !is.data.frame(peaks)) {
    stop(
      "'peaks' must be a numeric matrix or data frame.",
      call. = FALSE
    )
  }
  
  if (ncol(peaks) < 2L) {
    stop(
      "'peaks' must have at least two columns: peak positions and intensities.",
      call. = FALSE
    )
  }
  
  # Return early when no candidate peaks exist.
  if (nrow(peaks) == 0L) {
    return(peaks)
  }
  
  x_peaks <- peaks[, 1]
  intensity <- peaks[, 2]
  
  if (!is.numeric(x_peaks) || !is.numeric(intensity)) {
    stop(
      "The first and second columns of 'peaks' must be numeric.",
      call. = FALSE
    )
  }
  
  if (anyNA(x_peaks) || any(!is.finite(x_peaks))) {
    stop(
      "Peak positions in the first column of 'peaks' must be finite values.",
      call. = FALSE
    )
  }
  
  if (any(x_peaks < min(x) | x_peaks > max(x))) {
    stop(
      "All peak positions must lie within the x-axis range of 'data'.",
      call. = FALSE
    )
  }
  
  # ----------------------------------------
  # Set and validate thresholds
  # ----------------------------------------
  if (is.null(cutoff_peak_intensity)) {
    cutoff_peak_intensity <- stats::mad(y)
  }
  
  if (is.null(cutoff_peak_prominence)) {
    cutoff_peak_prominence <- stats::mad(y)
  }
  
  validate_nonnegative_scalar <- function(value, name) {
    if (!is.numeric(value) ||
        length(value) != 1L ||
        is.na(value) ||
        !is.finite(value) ||
        value < 0) {
      stop(
        sprintf("'%s' must be a non-negative finite numeric scalar.", name),
        call. = FALSE
      )
    }
  }
  
  validate_nonnegative_scalar(
    cutoff_peak_intensity,
    "cutoff_peak_intensity"
  )
  
  validate_nonnegative_scalar(
    cutoff_peak_prominence,
    "cutoff_peak_prominence"
  )
  
  validate_nonnegative_scalar(
    cutoff_peak_strength,
    "cutoff_peak_strength"
  )
  
  if (!is.numeric(k) ||
      length(k) != 1L ||
      is.na(k) ||
      !is.finite(k) ||
      k <= 0) {
    stop(
      "'k' must be a positive finite numeric scalar.",
      call. = FALSE
    )
  }
  
  # ----------------------------------------
  # Calculate peak metrics
  # ----------------------------------------
  prominence <- calculate_peak_prominence(
    data = data,
    peaks = peaks,
    valley_type = "lower"
  )
  
  strength <- estimate_peak_strength(
    data = data,
    peaks = peaks,
    k = k,
    normalization_type = normalization_type
  )
  
  # ----------------------------------------
  # Validate calculated metrics
  # ----------------------------------------
  if (!is.numeric(prominence) || length(prominence) != nrow(peaks)) {
    stop(
      paste0(
        "'calculate_peak_prominence()' must return a numeric vector ",
        "with the same length as the number of peaks."
      ),
      call. = FALSE
    )
  }
  
  if (!is.numeric(strength) || length(strength) != nrow(peaks)) {
    stop(
      paste0(
        "'estimate_peak_strength()' must return a numeric vector ",
        "with the same length as the number of peaks."
      ),
      call. = FALSE
    )
  }
  
  # ----------------------------------------
  # Filter peaks
  # ----------------------------------------
  keep <- is.finite(intensity) &
    is.finite(prominence) &
    is.finite(strength) &
    intensity > cutoff_peak_intensity &
    prominence > cutoff_peak_prominence &
    strength > cutoff_peak_strength
  
  filtered_peaks <- peaks[keep, , drop = FALSE]
  
  rownames(filtered_peaks) <- NULL
  
  return(filtered_peaks)
}


#' @rdname filter_peaks
#'
#' @param spectra A two-column `data.frame` or numeric `matrix` representing a
#'   single spectrum, or a `list` of such objects aligned with `peaks_list`.
#' @param peaks_list A peak table (`data.frame` or `matrix`) when `spectra` is a
#'   single spectrum, or a `list` of peak tables aligned with `spectra`.
#' @param n_cores A positive integer specifying the number of worker processes
#'   used when `spectra` is a list. The default is `1L`, which preserves
#'   sequential behaviour.
#'
#' @details
#' `filter_peaks_spectra()` applies [filter_peaks()] to a single spectrum and
#' its peak table, or to each spectrum-peak pair in a list. When `spectra` and
#' `peaks_list` are lists, they are matched by name: the intersection of
#' `names(spectra)` and `names(peaks_list)` determines which samples are
#' processed. Unnamed lists are matched positionally by assigning sequential
#' names.
#'
#' @return
#' `filter_peaks()` returns a single filtered peak table.
#' `filter_peaks_spectra()` returns a single filtered peak table when `spectra`
#' is a single spectrum, or a named list of filtered peak tables for the shared
#' sample names when `spectra` and `peaks_list` are lists.
#'
#' @seealso [find_peaks_spectra()], [preprocess_maldi_spectra()]
#'
#' @examples
#' \dontrun{
#' pp_spectra <- preprocess_maldi_spectra(spectra, n_cores = 2)
#' peaks_list <- find_peaks_spectra(pp_spectra, n_cores = 2)
#'
#' filtered <- filter_peaks_spectra(
#'   spectra = pp_spectra,
#'   peaks_list = peaks_list,
#'   cutoff_peak_strength = 0.2,
#'   n_cores = 2
#' )
#' }
#'
#' @export
filter_peaks_spectra <- function(spectra,
                                 peaks_list,
                                 cutoff_peak_intensity = NULL,
                                 cutoff_peak_prominence = NULL,
                                 cutoff_peak_strength = 0.2,
                                 k = 1,
                                 normalization_type = c("raw", "sqrt", "log10"),
                                 n_cores = 1L) {
  
  normalization_type <- match.arg(normalization_type)
  
  .assert_integer_scalar(
    x = n_cores,
    arg_name = "n_cores",
    lower_bound = 1L
  )
  
  filter_args <- list(
    cutoff_peak_intensity = cutoff_peak_intensity,
    cutoff_peak_prominence = cutoff_peak_prominence,
    cutoff_peak_strength = cutoff_peak_strength,
    k = k,
    normalization_type = normalization_type
  )
  
  is_single_spectrum <- is.data.frame(spectra) || is.matrix(spectra)
  
  if (is_single_spectrum) {
    if (!is.data.frame(peaks_list) && !is.matrix(peaks_list)) {
      stop(
        "When 'spectra' is a single spectrum, 'peaks_list' must be a single ",
        "peak table (data.frame or matrix).",
        call. = FALSE
      )
    }
    
    return(
      do.call(
        filter_peaks,
        c(list(data = spectra, peaks = peaks_list), filter_args)
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
  
  if (!is.list(peaks_list)) {
    stop(
      "'peaks_list' must be a peak table (data.frame or matrix) or a list of ",
      "such objects.",
      call. = FALSE
    )
  }
  
  if (length(spectra) == 0L) {
    return(spectra)
  }
  
  if (is.null(names(spectra))) {
    names(spectra) <- as.character(seq_along(spectra))
  }
  
  if (is.null(names(peaks_list))) {
    names(peaks_list) <- as.character(seq_along(peaks_list))
  }
  
  common <- intersect(names(spectra), names(peaks_list))
  
  if (length(common) < 1L) {
    stop(
      "'spectra' and 'peaks_list' must share at least one sample name.",
      call. = FALSE
    )
  }
  
  filtered <- .parallel_lapply(
    X = common,
    FUN = function(nm) {
      do.call(
        filter_peaks,
        c(
          list(
            data = spectra[[nm]],
            peaks = peaks_list[[nm]]
          ),
          filter_args
        )
      )
    },
    n_cores = n_cores
  )
  
  names(filtered) <- common
  
  filtered
}