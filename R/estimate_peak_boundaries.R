#' @name estimate_peak_boundaries
#' @title Estimate Adaptive Boundary Points Between Adjacent Peaks
#'
#' @description
#' Estimates adaptive boundary points between adjacent peaks using a
#' user-supplied continuous spectrum function. Each initial boundary is
#' calculated from an ensemble of intensity-weighted partition points based on
#' raw, square-root-transformed, and log-transformed peak intensities. When a
#' sufficiently deep valley exists between adjacent peaks, the weighted
#' boundary is replaced with the observed valley position.
#'
#' @param fun A function representing a continuous mass-spectrum curve, such as
#'   a Gaussian kernel regression function. The function must accept a numeric
#'   vector of x-axis values and return a numeric vector of the same length.
#'
#' @param x A numeric vector containing strictly increasing x-axis values used
#'   as the evaluation grid for valley detection. For MALDI-TOF MS spectra,
#'   these values typically represent m/z values.
#'
#' @param peaks A numeric matrix or data frame containing detected peaks.
#'   The first column must contain peak positions and the second column must
#'   contain the corresponding peak intensities. Additional columns are
#'   preserved in the returned object.
#'
#' @details
#' For each pair of adjacent peaks, three intensity-weighted partition points
#' are calculated using:
#'
#' \preformatted{
#' raw intensity
#' square-root-transformed intensity
#' log10(1 + intensity)-transformed intensity
#' }
#'
#' Let `x_left` and `x_right` denote the positions of two adjacent peaks and
#' let `w_left` and `w_right` denote the corresponding transformed
#' intensities. Each weighted partition point is calculated as:
#'
#' \deqn{
#' x_{bound} =
#' x_{left} +
#' (x_{right} - x_{left})
#' \frac{w_{left}}{w_{left} + w_{right}}
#' }
#'
#' The three partition points are averaged to reduce excessive dependence on a
#' single intensity scale.
#'
#' The supplied spectrum function is then evaluated over the x-axis grid
#' between the adjacent peaks. If the minimum function value within this
#' interval is less than or equal to one-half of the function value at the
#' weighted boundary, the observed valley position is used as the final
#' boundary. Otherwise, the weighted boundary is retained.
#'
#' The returned boundary height is calculated as:
#'
#' \deqn{
#' y_{bound} = \frac{fun(x_{bound})}{2}
#' }
#'
#' @return A list containing:
#'
#' \describe{
#'   \item{peaks}{
#'     The input peak table sorted by peak position. All original columns are
#'     preserved.
#'   }
#'   \item{boundary_points}{
#'     A data frame containing one row for each pair of adjacent peaks:
#'     \describe{
#'       \item{x_bound}{
#'         Estimated x-axis position of the boundary point.
#'       }
#'       \item{y_bound}{
#'         One-half of the spectrum-function value evaluated at `x_bound`.
#'       }
#'     }
#'   }
#' }
#'
#' If fewer than two peaks are supplied, `boundary_points` is returned as an
#' empty data frame.
#'
#' @examples
#' \dontrun{
#' gauss_kde <- get_gauss_kde(
#'   data = spectrum,
#'   bw = 1,
#'   d = 0
#' )
#'
#' boundary_points <- estimate_peak_boundaries(
#'   fun = gauss_kde,
#'   x = spectrum[, 1],
#'   peaks = peaks
#' )
#' }
#'
#' @export
estimate_peak_boundaries <- function(fun, x, peaks) {
  
  # ----------------------------------------
  # Validate spectrum function
  # ----------------------------------------
  if (!is.function(fun)) {
    stop(
      "'fun' must be a function.",
      call. = FALSE
    )
  }
  
  # ----------------------------------------
  # Validate x-axis grid
  # ----------------------------------------
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector.", call. = FALSE)
  }
  
  if (length(x) < 2L) {
    stop("'x' must contain at least two values.", call. = FALSE)
  }
  
  if (anyNA(x) || any(!is.finite(x))) {
    stop("'x' must contain only finite values.", call. = FALSE)
  }
  
  if (anyDuplicated(x)) {
    stop("'x' must not contain duplicates.", call. = FALSE)
  }
  
  if (is.unsorted(x, strictly = TRUE)) {
    stop("'x' must be strictly increasing.", call. = FALSE)
  }
  
  # ----------------------------------------
  # Validate peaks
  # ----------------------------------------
  if (!is.matrix(peaks) && !is.data.frame(peaks)) {
    stop("'peaks' must be a numeric matrix or data frame.", call. = FALSE)
  }
  
  if (ncol(peaks) < 2L) {
    stop("'peaks' must have at least two columns: peak positions and intensities.", call. = FALSE)
  }
  
  # Sort peaks while preserving all original columns.
  if (nrow(peaks) > 0L) {
    peaks <- peaks[order(peaks[, 1]), , drop = FALSE]
    rownames(peaks) <- NULL
  }
  
  # Return early when fewer than two peaks are available.
  if (nrow(peaks) < 2L) {
    boundary_points <- data.frame(
      x_bound = numeric(0),
      y_bound = numeric(0)
    )
    
    return(
      list(
        peaks = peaks,
        boundary_points = boundary_points
      )
    )
  }
  
  x_peaks <- peaks[, 1]
  y_peaks <- peaks[, 2]
  
  if (!is.numeric(x_peaks) || !is.numeric(y_peaks)) {
    stop("The first and second columns of 'peaks' must be numeric.", call. = FALSE)
  }
  
  if (anyNA(x_peaks) ||
      anyNA(y_peaks) ||
      any(!is.finite(x_peaks)) ||
      any(!is.finite(y_peaks))) {
    stop("Peak positions and peak intensities must contain only finite values.", call. = FALSE)
  }
  
  if (any(y_peaks < 0)) {
    stop("Peak intensities must be non-negative.", call. = FALSE)
  }
  
  if (anyDuplicated(x_peaks)) {
    stop("Peak positions must not contain duplicates.", call. = FALSE)
  }
  
  if (any(x_peaks < min(x) | x_peaks > max(x))) {
    stop("All peak positions must lie within the range of 'x'.", call. = FALSE)
  }
  
  # ----------------------------------------
  # Evaluate supplied spectrum function
  # ----------------------------------------
  y_fun <- fun(x)
  
  if (!is.numeric(y_fun) ||
      length(y_fun) != length(x) ||
      anyNA(y_fun) ||
      any(!is.finite(y_fun))) {
    stop(
      paste0(
        "'fun' must return a finite numeric vector with the same length ",
        "as its input."
      ),
      call. = FALSE
    )
  }
  
  # ----------------------------------------
  # Prepare adjacent peak pairs
  # ----------------------------------------
  n_peaks <- length(x_peaks)
  
  x_left <- x_peaks[-n_peaks]
  y_left <- y_peaks[-n_peaks]
  
  x_right <- x_peaks[-1L]
  y_right <- y_peaks[-1L]
  
  # ----------------------------------------
  # Internal weighted-boundary helper
  # ----------------------------------------
  calculate_weighted_boundary <- function(
    x_left_i,
    x_right_i,
    weight_left_i,
    weight_right_i) {
    
    denominator <- weight_left_i + weight_right_i
    
    if (!is.finite(denominator) || denominator <= 0) {
      return((x_left_i + x_right_i) / 2)
    }
    
    x_left_i +
      (x_right_i - x_left_i) *
      weight_left_i /
      denominator
  }
  
  # ----------------------------------------
  # Estimate initial weighted boundaries
  # ----------------------------------------
  x_bound_weighted <- vapply(
    seq_along(x_left),
    function(i) {
      
      x_left_i <- x_left[i]
      x_right_i <- x_right[i]
      
      y_left_i <- max(0, y_left[i])
      y_right_i <- max(0, y_right[i])
      
      sqrt_y_left_i <- sqrt(y_left_i)
      sqrt_y_right_i <- sqrt(y_right_i)
      
      log10_y_left_i <- log10(1 + y_left_i)
      log10_y_right_i <- log10(1 + y_right_i)
      
      x_bound_raw <- calculate_weighted_boundary(
        x_left_i = x_left_i,
        x_right_i = x_right_i,
        weight_left_i = y_left_i,
        weight_right_i = y_right_i
      )
      
      x_bound_sqrt <- calculate_weighted_boundary(
        x_left_i = x_left_i,
        x_right_i = x_right_i,
        weight_left_i = sqrt_y_left_i,
        weight_right_i = sqrt_y_right_i
      )
      
      x_bound_log10 <- calculate_weighted_boundary(
        x_left_i = x_left_i,
        x_right_i = x_right_i,
        weight_left_i = log10_y_left_i,
        weight_right_i = log10_y_right_i
      )
      
      mean(
        c(
          x_bound_raw,
          x_bound_sqrt,
          x_bound_log10
        )
      )
    },
    numeric(1)
  )
  
  # ----------------------------------------
  # Evaluate initial boundary heights
  # ----------------------------------------
  y_bound_weighted <- fun(x_bound_weighted)
  
  if (!is.numeric(y_bound_weighted) ||
      length(y_bound_weighted) != length(x_bound_weighted) ||
      anyNA(y_bound_weighted) ||
      any(!is.finite(y_bound_weighted))) {
    stop(
      paste0(
        "'fun' must return finite numeric values when evaluated at ",
        "weighted boundary points."
      ),
      call. = FALSE
    )
  }
  
  # ----------------------------------------
  # Replace weighted boundaries with valleys
  # when sufficiently deep valleys exist
  # ----------------------------------------
  valley_depth_ratio <- 0.5
  
  x_bound <- vapply(
    seq_along(x_bound_weighted),
    function(i) {
      
      x_left_i <- x_left[i]
      x_right_i <- x_right[i]
      
      keep <- x >= x_left_i & x <= x_right_i
      
      if (!any(keep)) {
        return(x_bound_weighted[i])
      }
      
      x_range <- x[keep]
      y_range <- y_fun[keep]
      
      if (length(y_range) == 0L ||
          anyNA(y_range) ||
          any(!is.finite(y_range))) {
        return(x_bound_weighted[i])
      }
      
      valley_idx <- which.min(y_range)
      
      x_valley <- x_range[valley_idx]
      y_valley <- y_range[valley_idx]
      
      if (y_valley <= y_bound_weighted[i] * valley_depth_ratio) {
        return(x_valley)
      }
      
      x_bound_weighted[i]
    },
    numeric(1)
  )
  
  # ----------------------------------------
  # Calculate final boundary heights
  # ----------------------------------------
  y_bound <- fun(x_bound) / 2
  
  if (!is.numeric(y_bound) ||
      length(y_bound) != length(x_bound) ||
      anyNA(y_bound) ||
      any(!is.finite(y_bound))) {
    stop(
      paste0(
        "'fun' must return finite numeric values when evaluated at ",
        "final boundary points."
      ),
      call. = FALSE
    )
  }
  
  boundary_points <- data.frame(
    x_bound = x_bound,
    y_bound = y_bound
  )
  colnames(boundary_points)[1:2] = colnames(peaks)[1:2]
  
  rownames(boundary_points) <- NULL
  
  # ----------------------------------------
  # Return results
  # ----------------------------------------
  return(boundary_points)
}