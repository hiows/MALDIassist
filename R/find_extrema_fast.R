#' @name find_extrema_fast
#' @title Fast Local Extrema Detection by Sign Changes
#'
#' @description
#' Detect local extrema from discrete signal data using sign changes in the
#' first difference.
#'
#' This function identifies local maxima and minima directly from discrete
#' \code{x} and \code{y} values. A local maximum is detected when the sign of
#' the first difference changes from positive to negative, and a local minimum
#' is detected when the sign changes from negative to positive.
#'
#' Flat-top and flat-bottom plateau regions are also handled. A flat-top
#' plateau, represented by a positive slope followed by one or more zero
#' differences and then a negative slope, is classified as a local maximum.
#' A flat-bottom plateau, represented by a negative slope followed by one or
#' more zero differences and then a positive slope, is classified as a local
#' minimum.
#'
#' Compared with derivative-based extrema detection, this function does not
#' require continuous derivative functions. It is intended as a fast candidate
#' detection method for discrete signals such as mass spectrometry spectra.
#'
#' @param x A numeric vector representing the x-axis values, such as m/z values.
#' @param y A numeric vector representing the signal intensity values.
#' @param plateau A character string specifying which point should represent a
#'   plateau-type extremum. One of \code{"middle"}, \code{"first"}, or
#'   \code{"last"}. \code{"first"} uses the first point of the plateau,
#'   \code{"last"} uses the last point of the plateau, and \code{"middle"}
#'   uses the middle point of the plateau. Default is \code{"middle"}.
#' @param na_rm Logical. If \code{TRUE}, pairs with \code{NA} in either
#'   \code{x} or \code{y} are removed before extrema detection. Default is
#'   \code{TRUE}.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{\code{local_min}}{A numeric vector containing the x-coordinates of
#'   detected local minima. Flat-bottom plateau regions are included as local
#'   minima. If no local minimum is detected, \code{NA_real_} is returned.}
#'   \item{\code{local_max}}{A numeric vector containing the x-coordinates of
#'   detected local maxima. Flat-top plateau regions are included as local
#'   maxima. If no local maximum is detected, \code{NA_real_} is returned.}
#' }
#'
#' @details
#' The detection rule is based on the sign of \code{diff(y)}:
#'
#' \itemize{
#'   \item positive to negative sign change: local maximum
#'   \item negative to positive sign change: local minimum
#'   \item positive to zero to negative sign pattern: flat-top local maximum
#'   \item negative to zero to positive sign pattern: flat-bottom local minimum
#' }
#'
#' Boundary plateaus are not classified because their left or right slope cannot
#' be determined safely.
#'
#' This function is fast because it only uses first differences and sign changes.
#' However, it can be sensitive to local noise. For noisy signals, smoothing or
#' additional filtering using intensity, signal-to-noise ratio, prominence, peak
#' width, or statistical significance is recommended after candidate detection.
#'
#' The output format is intentionally matched to \code{find_extrema()} as a named
#' list. Local minima and local maxima are returned separately as
#' \code{local_min} and \code{local_max}. Elements with no detected values are
#' returned as \code{NA_real_}.
#'
#' @examples
#' x <- seq(0, 2 * pi, length.out = 101)
#' y <- sin(x)
#'
#' find_extrema_fast(
#'   x = x,
#'   y = y
#' )
#'
#' # Flat-top plateau example
#' x2 <- 1:6
#' y2 <- c(1, 3, 5, 5, 5, 2)
#'
#' find_extrema_fast(
#'   x = x2,
#'   y = y2,
#'   plateau = "middle"
#' )
#'
#' # With noisy signal
#' set.seed(1)
#' y_noise <- sin(x) + rnorm(length(x), sd = 0.05)
#'
#' find_extrema_fast(
#'   x = x,
#'   y = y_noise,
#'   plateau = "middle"
#' )
#'
#' @export
find_extrema_fast <- function(x, y,
                              plateau = c("middle", "first", "last"),
                              na_rm = TRUE) {
  
  plateau <- match.arg(plateau)
  
  # ----------------------------------------
  # Input validation
  # ----------------------------------------
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("'x' and 'y' must be numeric vectors.", call. = FALSE)
  }
  
  if (length(x) != length(y)) {
    stop("'x' and 'y' must have the same length.", call. = FALSE)
  }
  
  if (length(y) < 3L) {
    stop("'x' and 'y' must have at least 3 points.", call. = FALSE)
  }
  
  if (!is.logical(na_rm) || length(na_rm) != 1L || is.na(na_rm)) {
    stop("'na_rm' must be a single TRUE or FALSE value.", call. = FALSE)
  }
  
  if (na_rm) {
    keep <- complete.cases(x, y)
    x <- x[keep]
    y <- y[keep]
  }
  
  if (length(y) < 3L) {
    stop("After removing NA values, at least 3 points are required.", call. = FALSE)
  }
  
  if (anyNA(x) || anyNA(y)) {
    stop("'x' and 'y' must not contain NA values when 'na_rm = FALSE'.", call. = FALSE)
  }
  
  if (any(!is.finite(x)) || any(!is.finite(y))) {
    stop("'x' and 'y' must contain only finite numeric values.", call. = FALSE)
  }
  
  # ----------------------------------------
  # First difference and sign
  # ----------------------------------------
  dy <- diff(y)
  s <- sign(dy)
  
  local_max_idx <- integer(0)
  local_min_idx <- integer(0)
  
  # ----------------------------------------
  # Case 1: sharp extrema without plateau
  # ----------------------------------------
  # Maxima: + to -
  sharp_max_idx <- which(s[-length(s)] > 0 & s[-1L] < 0) + 1L
  
  # Minima: - to +
  sharp_min_idx <- which(s[-length(s)] < 0 & s[-1L] > 0) + 1L
  
  local_max_idx <- c(local_max_idx, sharp_max_idx)
  local_min_idx <- c(local_min_idx, sharp_min_idx)
  
  # ----------------------------------------
  # Case 2: plateau extrema
  # ----------------------------------------
  zero_idx <- which(s == 0)
  
  if (length(zero_idx) > 0L) {
    
    r <- rle(s == 0)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1L
    
    plateau_runs <- which(r$values)
    
    for (k in plateau_runs) {
      
      start_s <- starts[k]
      end_s <- ends[k]
      
      left_s_idx <- start_s - 1L
      right_s_idx <- end_s + 1L
      
      # Cannot classify boundary plateaus safely
      if (left_s_idx < 1L || right_s_idx > length(s)) {
        next
      }
      
      left_sign <- s[left_s_idx]
      right_sign <- s[right_s_idx]
      
      # dy[start_s:end_s] == 0 means
      # y[start_s] to y[end_s + 1] are flat
      start_y <- start_s
      end_y <- end_s + 1L
      
      idx_one <- switch(
        plateau,
        first = start_y,
        last = end_y,
        middle = round((start_y + end_y) / 2)
      )
      
      # + 0...0 - : flat-top local maximum
      if (left_sign > 0 && right_sign < 0) {
        local_max_idx <- c(local_max_idx, idx_one)
      }
      
      # - 0...0 + : flat-bottom local minimum
      if (left_sign < 0 && right_sign > 0) {
        local_min_idx <- c(local_min_idx, idx_one)
      }
    }
  }
  
  # ----------------------------------------
  # Clean duplicated indices
  # ----------------------------------------
  local_max_idx <- sort(unique(local_max_idx))
  local_min_idx <- sort(unique(local_min_idx))
  
  local_max <- x[local_max_idx]
  local_min <- x[local_min_idx]
  
  if (length(local_max) == 0L) local_max <- NA_real_
  if (length(local_min) == 0L) local_min <- NA_real_
  
  return(list(
    local_min = local_min,
    local_max = local_max
  ))
}