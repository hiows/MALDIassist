#' @name get_curvature_fun
#' @title Create a Curvature Function from First and Second Derivative Functions
#' 
#' @description
#' Creates a closure function that calculates the curvature of a one-dimensional
#' curve \eqn{y = f(x)} using the first and second derivative functions.
#'
#' The returned function evaluates curvature at arbitrary x values. By default,
#' it calculates the absolute curvature:
#'
#' \deqn{
#' \kappa(x) = \frac{|f''(x)|}{\left(1 + [f'(x)]^2\right)^{3/2}}
#' }
#'
#' If \code{absolute = FALSE}, the sign of the second derivative is retained:
#'
#' \deqn{
#' \kappa_s(x) = \frac{f''(x)}{\left(1 + [f'(x)]^2\right)^{3/2}}
#' }
#'
#' @param first_deriv A function representing the first derivative \eqn{f'(x)}.
#'   It must take a numeric vector as input and return a numeric vector of
#'   the same length.
#' @param second_deriv A function representing the second derivative \eqn{f''(x)}.
#'   It must take a numeric vector as input and return a numeric vector of
#'   the same length.
#' @param absolute Logical. If \code{TRUE}, the absolute value of the second
#'   derivative is used in the numerator, returning non-negative curvature
#'   values. If \code{FALSE}, the sign of the second derivative is retained,
#'   returning a signed curvature-like value. Default is \code{TRUE}.
#'
#' @return
#' A function that takes a numeric vector \code{x} and returns curvature values
#' at \code{x}. If \code{absolute = TRUE}, the returned values are
#' non-negative. If \code{absolute = FALSE}, the returned values retain the sign
#' of the second derivative.
#'
#' @details
#' This function returns a closure. The derivative functions
#' \code{first_deriv} and \code{second_deriv}, as well as the \code{absolute}
#' option, are stored inside the returned function environment. Therefore, the
#' returned curvature function can be evaluated repeatedly without passing the
#' derivative functions again.
#'
#' The actual curvature calculation is performed by the internal C++ function
#' \code{cpp_curvature()}.
#'
#' In MALDI-TOF spectrum analysis, curvature can be used to quantify local
#' peak sharpness or shape changes in a smoothed continuous representation
#' of the spectrum. Absolute curvature is useful when only the magnitude of
#' local bending is required, whereas signed curvature-like values can be useful
#' when distinguishing concave-up and concave-down regions.
#'
#' @examples
#' first_deriv <- function(x) 2 * x
#' second_deriv <- function(x) rep(2, length(x))
#'
#' curvature_fun <- get_curvature_fun(
#'   first_deriv = first_deriv,
#'   second_deriv = second_deriv
#' )
#'
#' curvature_fun(seq(-2, 2, length.out = 5))
#'
#' signed_curvature_fun <- get_curvature_fun(
#'   first_deriv = first_deriv,
#'   second_deriv = second_deriv,
#'   absolute = FALSE
#' )
#'
#' signed_curvature_fun(seq(-2, 2, length.out = 5))
#'
#' @export
get_curvature_fun <- function(first_deriv, second_deriv, absolute = TRUE) {
  fun <- function(x) {
    cpp_curvature(
      x = x,
      first_deriv = first_deriv,
      second_deriv = second_deriv,
      absolute = absolute
    )
  }
  return(fun)
}