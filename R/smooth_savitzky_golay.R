#' @name smooth_savitzky_golay
#' @title Smooth a Spectrum Using the Savitzky-Golay Filter
#'
#' @description
#' Applies Savitzky-Golay smoothing to a two-column spectrum using a C++ backend.
#'
#' The first column of \code{data} is treated as the x-axis, usually m/z values,
#' and the second column is treated as the signal intensity. Only the intensity
#' values are smoothed; the x-axis values are returned unchanged.
#'
#' @param data A matrix or data frame containing a two-column spectrum. The
#'   first column is assumed to contain x values, usually m/z values, and the
#'   second column is assumed to contain y values, usually intensity values.
#' @param hws Integer-like numeric value. Half-window size for the
#'   Savitzky-Golay filter. The full window size is \code{2 * hws + 1}.
#'   Default is \code{10L}.
#' @param pno Integer-like numeric value. Polynomial order used for local
#'   polynomial fitting. Default is \code{3L}.
#'
#' @return A data frame with two columns. The first column contains the original
#'   x values, and the second column contains the smoothed y values. Column names
#'   are inherited from the first two columns of \code{data}.
#'
#' @details
#' This function is an R wrapper around the C++ backend
#' \code{cpp_savitzkyGolay_filter()}.
#'
#' Boundary regions are handled using boundary-specific Savitzky-Golay
#' coefficients. The center region is smoothed using the central filter
#' coefficients.
#'
#' The full window size, \code{2 * hws + 1}, must be larger than the polynomial
#' order \code{pno}. The number of rows in \code{data} must be greater than or
#' equal to the full window size.
#'
#' @examples
#' x <- seq(1000, 2000, length.out = 5000)
#' y <- dnorm(x, mean = 1500, sd = 30) * 100 +
#'      rnorm(length(x), sd = 2)
#'
#' data <- data.frame(mz = x, intensity = y)
#'
#' smoothed <- smooth_savitzky_golay(
#'   data = data,
#'   hws = 10,
#'   pno = 3
#' )
#'
#' head(smoothed)
#'
#' @export
smooth_savitzky_golay <- function(data,
                                  hws = 10L,
                                  pno = 3L) {
  
  # ------------------------------
  # 1. Input validation
  # ------------------------------
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("'data' must be a matrix or data frame.", call. = FALSE)
  }
  
  if (ncol(data) < 2L) {
    stop("'data' must have at least two columns: x and y.",
         call. = FALSE)
  }
  
  if (nrow(data) < 1L) {
    stop("'data' must contain at least one row.", call. = FALSE)
  }
  
  if (!is.numeric(hws) ||
      length(hws) != 1L ||
      is.na(hws) ||
      !is.finite(hws)) {
    stop("'hws' must be a single finite numeric value.", call. = FALSE)
  }
  
  if (!is.numeric(pno) ||
      length(pno) != 1L ||
      is.na(pno) ||
      !is.finite(pno)) {
    stop("'pno' must be a single finite numeric value.", call. = FALSE)
  }
  
  if (hws < 1L) {
    stop("'hws' must be greater than or equal to 1.", call. = FALSE)
  }
  
  if (pno < 0L) {
    stop("'pno' must be greater than or equal to 0.", call. = FALSE)
  }
  
  if (hws != as.integer(hws)) {
    stop("'hws' must be an integer-like value.", call. = FALSE)
  }
  
  if (pno != as.integer(pno)) {
    stop("'pno' must be an integer-like value.", call. = FALSE)
  }
  
  hws <- as.integer(hws)
  pno <- as.integer(pno)
  
  x <- as.double(data[, 1L])
  y <- as.double(data[, 2L])
  
  if (anyNA(x) || anyNA(y)) {
    stop("'data' contains NA values in the first or second column.",
         call. = FALSE)
  }
  
  if (any(!is.finite(x)) || any(!is.finite(y))) {
    stop("'data' contains non-finite values in the first or second column.",
         call. = FALSE)
  }
  
  n <- length(y)
  window_size <- 2L * hws + 1L
  
  if (window_size <= pno) {
    stop("The full window size, '2 * hws + 1', must be larger than 'pno'.",
         call. = FALSE)
  }
  
  if (n < window_size) {
    stop("The number of rows in 'data' must be greater than or equal to the full window size.",
         call. = FALSE)
  }
  
  # ------------------------------
  # 2. Smooth intensity using C++
  # ------------------------------
  y_smooth <- cpp_savitzkyGolay_filter(
    y = y,
    hws = hws,
    pno = pno
  )
  
  if (!is.numeric(y_smooth) || length(y_smooth) != n) {
    stop("'cpp_savitzkyGolay_filter()' must return a numeric vector with the same length as 'y'.",
         call. = FALSE)
  }
  
  if (anyNA(y_smooth) || any(!is.finite(y_smooth))) {
    stop("'cpp_savitzkyGolay_filter()' returned non-finite values.",
         call. = FALSE)
  }
  
  # ------------------------------
  # 3. Return smoothed spectrum
  # ------------------------------
  df <- data.frame(
    x = x,
    y = y_smooth,
    stringsAsFactors = FALSE
  )
  
  colnames(df) <- colnames(data)[1:2]
  rownames(df) <- NULL
  
  return(df)
}