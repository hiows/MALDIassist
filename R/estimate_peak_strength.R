#' @name estimate_peak_strength
#' @title Estimate Peak Strength
#'
#' @description
#' Estimate a bounded peak strength score for each detected peak.
#'
#' This function validates input spectrum data and detected peaks, then calls
#' the Rcpp backend `cpp_estimate_peak_strength()`. The peak strength score is
#' calculated from the local intensity scale around each peak relative to the
#' global intensity scale of the whole spectrum.
#'
#' The intensity can be transformed before calculation using one of three
#' normalization types:
#'
#' \itemize{
#'   \item `"raw"`: use raw intensity `y`
#'   \item `"sqrt"`: use `sqrt(y)`
#'   \item `"log10"`: use `log10(1 + y)`
#' }
#'
#' @param data A data frame or matrix with at least two columns. The first
#'   column must be x values, such as m/z, and the second column must be
#'   non-negative intensity values.
#' @param peaks A data frame or matrix with at least two columns. The first
#'   column must be peak x positions and the second column must be peak
#'   intensities. Only the first two columns are used.
#' @param k A single non-negative finite numeric value. MAD multiplier used
#'   when estimating local and global strength scales. Default is 1.
#' @param normalization_type Character string specifying the intensity
#'   transformation method. One of `"raw"`, `"sqrt"`, or `"log10"`.
#'   Default is `"raw"`.
#'
#' @return A numeric vector of peak strength scores. The returned vector has
#'   the same length as the number of input peaks. Peaks for which local
#'   valleys cannot be defined may return `NA_real_`.
#'
#' @examples
#' x <- seq(1000, 2000, length.out = 2000)
#' y <- dnorm(x, mean = 1400, sd = 20) * 100 +
#'      dnorm(x, mean = 1700, sd = 25) * 80
#' spectrum_data <- data.frame(mz = x, intensity = y)
#'
#' detected_peaks <- find_peaks_fast(spectrum_data, hws_peaks = 10)
#'
#' peak_strength <- estimate_peak_strength(
#'   data = spectrum_data,
#'   peaks = detected_peaks,
#'   k = 1,
#'   normalization_type = "log10"
#' )
#'
#' peak_strength
#'
#' @export
estimate_peak_strength <- function(data,
                                   peaks,
                                   k = 1,
                                   normalization_type = c("raw", "sqrt", "log10")) {
  
  # ----------------------------------------
  # 0. Validate backend function
  # ----------------------------------------
  if (!exists("cpp_estimate_peak_strength", mode = "function")) {
    stop(
      "'cpp_estimate_peak_strength()' was not found. ",
      "Please compile or load the Rcpp backend first.",
      call. = FALSE
    )
  }
  
  # ----------------------------------------
  # 1. Validate normalization_type
  # ----------------------------------------
  normalization_type <- match.arg(normalization_type)
  
  # ----------------------------------------
  # 2. Validate data
  # ----------------------------------------
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("'data' must be a data frame or matrix.", call. = FALSE)
  }
  
  if (ncol(data) < 2L) {
    stop("'data' must have at least two columns: x and y.", call. = FALSE)
  }
  
  if (nrow(data) < 3L) {
    stop("'data' must contain at least 3 rows.", call. = FALSE)
  }
  
  x_raw <- suppressWarnings(as.numeric(data[, 1L]))
  y_raw <- suppressWarnings(as.numeric(data[, 2L]))
  
  if (length(x_raw) != length(y_raw)) {
    stop("The first and second columns of 'data' must have the same length.", call. = FALSE)
  }
  
  if (anyNA(x_raw) || anyNA(y_raw)) {
    stop("'data' contains NA values in the first or second column.", call. = FALSE)
  }
  
  if (any(!is.finite(x_raw)) || any(!is.finite(y_raw))) {
    stop("'data' contains non-finite values in the first or second column.", call. = FALSE)
  }
  
  if (anyDuplicated(x_raw)) {
    stop("'x' values in 'data' must not be duplicated.", call. = FALSE)
  }
  
  if (is.unsorted(x_raw, strictly = TRUE)) {
    stop("'x' values in 'data' must be strictly increasing.", call. = FALSE)
  }
  
  if (any(y_raw < 0)) {
    stop(
      "'y' values in 'data' must be non-negative when using ",
      "'raw', 'sqrt', or 'log10' normalization.",
      call. = FALSE
    )
  }
  
  data <- cbind(
    x = x_raw,
    y = y_raw
  )
  
  # ----------------------------------------
  # 3. Validate peaks
  # ----------------------------------------
  if (!is.data.frame(peaks) && !is.matrix(peaks)) {
    stop("'peaks' must be a data frame or matrix.", call. = FALSE)
  }
  
  if (ncol(peaks) < 2L) {
    stop("'peaks' must have at least two columns: peak x and peak y.", call. = FALSE)
  }
  
  if (nrow(peaks) == 0L) {
    return(numeric(0))
  }
  
  peak_x <- suppressWarnings(as.numeric(peaks[, 1L]))
  peak_y <- suppressWarnings(as.numeric(peaks[, 2L]))
  
  if (length(peak_x) != length(peak_y)) {
    stop("The first and second columns of 'peaks' must have the same length.", call. = FALSE)
  }
  
  if (anyNA(peak_x) || anyNA(peak_y)) {
    stop("'peaks' contains NA values in the first or second column.", call. = FALSE)
  }
  
  if (any(!is.finite(peak_x)) || any(!is.finite(peak_y))) {
    stop("'peaks' contains non-finite values in the first or second column.", call. = FALSE)
  }
  
  if (anyDuplicated(peak_x)) {
    stop("'peak x' values in 'peaks' must not be duplicated.", call. = FALSE)
  }
  
  if (any(peak_y < 0)) {
    stop("'peak y' values in 'peaks' must be non-negative.", call. = FALSE)
  }
  
  if (any(peak_x < min(x_raw) | peak_x > max(x_raw))) {
    stop("'peak x' values in 'peaks' must be within the x range of 'data'.", call. = FALSE)
  }
  
  if (any(peak_x == min(x_raw) | peak_x == max(x_raw))) {
    warning(
      "Some peak x values are located at the boundary of 'data'. ",
      "Their local valleys may not be defined and may return NA.",
      call. = FALSE
    )
  }
  
  peaks <- cbind(
    x = peak_x,
    y = peak_y
  )
  
  peaks <- peaks[order(peaks[, 1L]), , drop = FALSE]
  
  # ----------------------------------------
  # 4. Validate parameters
  # ----------------------------------------
  if (!is.numeric(k) ||
      length(k) != 1L ||
      is.na(k) ||
      !is.finite(k)) {
    stop("'k' must be a single finite numeric value.", call. = FALSE)
  }
  
  if (k < 0) {
    stop("'k' must be greater than or equal to 0.", call. = FALSE)
  }
  
  # ----------------------------------------
  # 5. Estimate peak strength
  # ----------------------------------------
  peak_strength <- cpp_estimate_peak_strength(
    data = data,
    peaks = peaks,
    k = k,
    normalization = normalization_type
  )
  
  # ----------------------------------------
  # 6. Validate output
  # ----------------------------------------
  if (!is.numeric(peak_strength)) {
    stop("'cpp_estimate_peak_strength()' must return a numeric vector.", call. = FALSE)
  }
  
  if (length(peak_strength) != nrow(peaks)) {
    stop("'cpp_estimate_peak_strength()' must return a vector with the same ", 
      "length as the number of peaks.",
      call. = FALSE
    )
  }
  
  return(peak_strength)
}