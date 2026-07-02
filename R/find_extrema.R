#' @name find_extrema
#' @title Find Local Extrema from First and Second Derivative Functions
#'
#' @description
#' Finds local extrema by locating the roots of a first derivative function and
#' classifying the stationary points using the corresponding second derivative
#' values.
#'
#' This function first identifies points where the first derivative is zero
#' using the internal C++ root finder. The second derivative is then evaluated
#' at each root position to classify the stationary points as local minima,
#' local maxima, or plateau-like points.
#'
#' @param first_deriv A function representing the first derivative of the
#'   original curve. It must accept a numeric vector and return a numeric vector
#'   of the same length.
#' @param second_deriv A function representing the second derivative of the
#'   original curve. It must accept a numeric vector and return a numeric vector
#'   of the same length.
#' @param x A strictly increasing numeric vector defining the search range for
#'   root detection.
#' @param tol A positive numeric scalar used as the numerical tolerance for root
#'   finding and second-derivative-based classification. Default is `1e-5`.
#' @param max_iter A positive numeric scalar specifying the maximum number of
#'   iterations used by the root-finding procedure. Default is `100L`.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{local_min}{Numeric vector of x positions classified as local minima.}
#'   \item{local_max}{Numeric vector of x positions classified as local maxima.}
#'   \item{plateau}{Numeric vector of x positions classified as plateau-like
#'   stationary points.}
#' }
#' If no root is found, or if no point belongs to a given category, `NA_real_`
#' is returned for that element.
#'
#' @details
#' Stationary points are classified according to the second derivative value at
#' each root of the first derivative:
#' \itemize{
#'   \item `second_deriv(x_root) > tol`: local minimum
#'   \item `second_deriv(x_root) < -tol`: local maximum
#'   \item `abs(second_deriv(x_root)) <= tol`: plateau-like stationary point
#' }
#'
#' Root finding is performed by the package's internal C++ backend, which
#' returns numeric root positions within the range of `x`.
#'
#' @examples
#' f_d1 <- function(x) 2 * x
#' f_d2 <- function(x) rep(2, length(x))
#'
#' x <- seq(-5, 5, length.out = 101)
#'
#' find_extrema(
#'   first_deriv = f_d1,
#'   second_deriv = f_d2,
#'   x = x
#' )
#'
#' @export
find_extrema <- function(first_deriv,
                         second_deriv,
                         x,
                         tol = 1e-5,
                         max_iter = 100L) {
  
  # ----------------------------------------
  # Input validation
  # ----------------------------------------
  if (!is.function(first_deriv)) {
    stop("'first_deriv' must be a function.", call. = FALSE)
  }
  
  if (!is.function(second_deriv)) {
    stop("'second_deriv' must be a function.", call. = FALSE)
  }
  
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector.", call. = FALSE)
  }
  
  if (length(x) < 2L) {
    stop("'x' must contain at least two values.", call. = FALSE)
  }
  
  if (anyNA(x)) {
    stop("'x' must not contain NA values.", call. = FALSE)
  }
  
  if (any(!is.finite(x))) {
    stop("'x' must contain only finite numeric values.", call. = FALSE)
  }
  
  if (anyDuplicated(x)) {
    stop("'x' must not contain duplicated values.", call. = FALSE)
  }
  
  if (is.unsorted(x, strictly = TRUE)) {
    stop("'x' must be strictly increasing.", call. = FALSE)
  }
  
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || !is.finite(tol)) {
    stop("'tol' must be a single finite numeric value.", call. = FALSE)
  }
  
  if (tol <= 0) {
    stop("'tol' must be greater than 0.", call. = FALSE)
  }
  
  if (!is.numeric(max_iter) || length(max_iter) != 1L ||
      is.na(max_iter) || !is.finite(max_iter)) {
    stop("'max_iter' must be a single finite numeric value.", call. = FALSE)
  }
  
  if (max_iter < 1L) {
    stop("'max_iter' must be greater than or equal to 1.", call. = FALSE)
  }
  
  max_iter <- as.integer(max_iter)
  
  
  # ----------------------------------------
  # Find roots of the first derivative
  # ----------------------------------------
  x_roots <- tryCatch(
    cpp_find_roots(
      fun = first_deriv,
      x = x,
      yi = 0,
      tol = tol,
      max_iter = max_iter
    ),
    error = function(e) {
      stop("'cpp_find_roots' failed while finding roots of 'first_deriv'.\n",
           "Original error: ", conditionMessage(e),
           call. = FALSE)
    }
  )
  
  if (is.null(x_roots) || length(x_roots) == 0L) {
    return(list(
      local_min = NA_real_,
      local_max = NA_real_,
      plateau = NA_real_
    ))
  }
  
  if (!is.numeric(x_roots)) {
    stop("'cpp_find_roots' must return a numeric vector.", call. = FALSE)
  }
  
  if (anyNA(x_roots) || any(!is.finite(x_roots))) {
    stop("'cpp_find_roots' returned NA, NaN, Inf, or -Inf values.", call. = FALSE)
  }
  
  
  # ----------------------------------------
  # Evaluate second derivative at roots
  # ----------------------------------------
  d2y_roots <- tryCatch(
    second_deriv(x_roots),
    error = function(e) {
      stop("'second_deriv' failed when evaluated at root positions.\n",
           "Original error: ", conditionMessage(e),
           call. = FALSE)
    }
  )
  
  if (!is.numeric(d2y_roots)) {
    stop("'second_deriv' must return a numeric vector.", call. = FALSE)
  }
  
  if (length(d2y_roots) != length(x_roots)) {
    stop("'second_deriv' must return a vector with the same length as 'x_roots'.",
         call. = FALSE)
  }
  
  if (anyNA(d2y_roots) || any(!is.finite(d2y_roots))) {
    stop("'second_deriv' returned NA, NaN, Inf, or -Inf at root positions.",
         call. = FALSE)
  }
  
  
  # ----------------------------------------
  # Classify stationary points
  # ----------------------------------------
  idx_localmax <- d2y_roots < -tol
  idx_localmin <- d2y_roots >  tol
  idx_plateau   <- abs(d2y_roots) <= tol
  
  local_max <- x_roots[idx_localmax]
  local_min <- x_roots[idx_localmin]
  plateau    <- x_roots[idx_plateau]
  
  if (length(local_max) == 0L) local_max <- NA_real_
  if (length(local_min) == 0L) local_min <- NA_real_
  if (length(plateau)    == 0L) plateau    <- NA_real_
  
  return(list(
    local_min = local_min,
    local_max = local_max,
    plateau = plateau
  ))
}