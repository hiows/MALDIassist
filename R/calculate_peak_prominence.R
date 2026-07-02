#' @name calculate_peak_prominence
#' @title Calculate Peak Prominence
#'
#' @description
#' Calculate peak prominence for detected peaks.
#'
#' This function validates spectrum data and detected peaks, then calls the
#' Rcpp backend `cpp_calculate_peak_prominence()`. Peak prominence is calculated
#' using the left and right valleys around each peak. The reference valley can
#' be selected by `valley_type`.
#'
#' \itemize{
#'   \item `"higher"`: use the higher of the left and right valleys.
#'   This is the more conservative and standard-like prominence definition.
#'   \item `"lower"`: use the lower of the left and right valleys.
#'   This gives a larger prominence value.
#' }
#'
#' @param data A data frame or matrix with at least two columns. The first
#'   column must be x values, such as m/z, and the second column must be
#'   intensity values.
#' @param peaks A data frame or matrix with at least two columns. The first
#'   column must be peak x positions and the second column must be peak
#'   intensities. Only the first two columns are used.
#' @param valley_type Character string specifying which reference valley should
#'   be used. One of `"higher"` or `"lower"`. Default is `"higher"`.
#' @param zero_tol A single non-negative finite numeric value. Prominence values
#'   with absolute value smaller than or equal to `zero_tol` are treated as zero.
#'   Default is `sqrt(.Machine$double.eps)`.
#'
#' @return A numeric vector of peak prominence values. The returned vector has
#'   the same length as the number of input peaks.
#'
#' @examples
#' x <- seq(1000, 2000, length.out = 2000)
#' y <- dnorm(x, mean = 1400, sd = 20) * 100 +
#'      dnorm(x, mean = 1700, sd = 25) * 80
#' spectrum_data <- data.frame(mz = x, intensity = y)
#'
#' detected_peaks <- find_peaks_fast(spectrum_data, hws_peaks = 10)
#'
#' peak_prominence <- calculate_peak_prominence(
#'   data = spectrum_data,
#'   peaks = detected_peaks,
#'   valley_type = "higher"
#' )
#'
#' peak_prominence
#'
#' @export
calculate_peak_prominence <- function(data,
                                      peaks,
                                      valley_type = c("higher", "lower"),
                                      zero_tol = sqrt(.Machine$double.eps)) {
  
  # ----------------------------------------
  # 0. Validate backend function
  # ----------------------------------------
  if (!exists("cpp_calculate_peak_prominence", mode = "function")) {
    stop(
      "'cpp_calculate_peak_prominence()' was not found. ",
      "Please compile or load the Rcpp backend first.",
      call. = FALSE
    )
  }
  
  # ----------------------------------------
  # 1. Validate valley_type
  # ----------------------------------------
  valley_type <- match.arg(valley_type)
  
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
    stop(
      "The first and second columns of 'data' must have the same length.",
      call. = FALSE
    )
  }
  
  if (anyNA(x_raw) || anyNA(y_raw)) {
    stop(
      "'data' contains NA values in the first or second column.",
      call. = FALSE
    )
  }
  
  if (any(!is.finite(x_raw)) || any(!is.finite(y_raw))) {
    stop(
      "'data' contains non-finite values in the first or second column.",
      call. = FALSE
    )
  }
  
  if (anyDuplicated(x_raw)) {
    stop("'x' values in 'data' must not be duplicated.", call. = FALSE)
  }
  
  # C++ function assumes sorted data.
  # Unlike estimate_peak_strength(), this wrapper can sort data internally.
  data <- cbind(
    x = x_raw,
    y = y_raw
  )
  
  data <- data[order(data[, 1L]), , drop = FALSE]
  
  if (is.unsorted(data[, 1L], strictly = TRUE)) {
    stop("'x' values in 'data' must be strictly increasing.", call. = FALSE)
  }
  
  # ----------------------------------------
  # 3. Validate peaks
  # ----------------------------------------
  if (!is.data.frame(peaks) && !is.matrix(peaks)) {
    stop("'peaks' must be a data frame or matrix.", call. = FALSE)
  }
  
  if (ncol(peaks) < 2L) {
    stop(
      "'peaks' must have at least two columns: peak x and peak y.",
      call. = FALSE
    )
  }
  
  if (nrow(peaks) == 0L) {
    return(numeric(0))
  }
  
  peak_x <- suppressWarnings(as.numeric(peaks[, 1L]))
  peak_y <- suppressWarnings(as.numeric(peaks[, 2L]))
  
  if (length(peak_x) != length(peak_y)) {
    stop(
      "The first and second columns of 'peaks' must have the same length.",
      call. = FALSE
    )
  }
  
  if (anyNA(peak_x) || anyNA(peak_y)) {
    stop(
      "'peaks' contains NA values in the first or second column.",
      call. = FALSE
    )
  }
  
  if (any(!is.finite(peak_x)) || any(!is.finite(peak_y))) {
    stop(
      "'peaks' contains non-finite values in the first or second column.",
      call. = FALSE
    )
  }
  
  if (anyDuplicated(peak_x)) {
    stop("'peak x' values in 'peaks' must not be duplicated.", call. = FALSE)
  }
  
  if (any(peak_x < min(data[, 1L]) | peak_x > max(data[, 1L]))) {
    stop(
      "'peak x' values in 'peaks' must be within the x range of 'data'.",
      call. = FALSE
    )
  }
  
  if (any(peak_x == min(data[, 1L]) | peak_x == max(data[, 1L]))) {
    warning(
      "Some peak x values are located at the boundary of 'data'. ",
      "Their prominence may be underestimated because one side valley may be poorly defined.",
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
  if (!is.numeric(zero_tol) ||
      length(zero_tol) != 1L ||
      is.na(zero_tol) ||
      !is.finite(zero_tol)) {
    stop(
      "'zero_tol' must be a single finite numeric value.",
      call. = FALSE
    )
  }
  
  if (zero_tol < 0) {
    stop(
      "'zero_tol' must be greater than or equal to 0.",
      call. = FALSE
    )
  }
  
  # ----------------------------------------
  # 5. Calculate peak prominence
  # ----------------------------------------
  peak_prominence <- cpp_calculate_peak_prominence(
    data = data,
    peaks = peaks,
    valley_type = valley_type,
    zero_tol = zero_tol
  )
  
  # ----------------------------------------
  # 6. Validate output
  # ----------------------------------------
  if (!is.numeric(peak_prominence)) {
    stop(
      "'cpp_calculate_peak_prominence()' must return a numeric vector.",
      call. = FALSE
    )
  }
  
  if (length(peak_prominence) != nrow(peaks)) {
    stop(
      "'cpp_calculate_peak_prominence()' must return a vector with the same ",
      "length as the number of peaks.",
      call. = FALSE
    )
  }
  
  # names(peak_prominence) <- as.character(peaks[, 1L])
  
  return(peak_prominence)
}