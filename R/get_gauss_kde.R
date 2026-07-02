#' @name get_gauss_kde
#' @title Create a Gaussian KDE Function for a Mass Spectrum
#'
#' @description
#' Creates a closure function for Gaussian kernel regression-based smoothing
#' of a two-column mass spectrum. The returned function can be evaluated at
#' arbitrary m/z values and can also return the first, second, or third
#' derivative of the smoothed function.
#'
#' The input spectrum is assumed to contain m/z values in the first column and
#' intensity values in the second column. The observed m/z values, observed
#' intensities, bandwidth, and derivative order are stored inside the returned
#' function environment. Therefore, the returned function can be saved with
#' \code{save()} and reused after \code{load()}, as long as the MALDIassist
#' package is loaded.
#'
#' @param data A matrix or data frame with at least two columns. The first
#'   column should contain m/z values and the second column should contain
#'   intensity values.
#' @param bw A positive numeric scalar specifying the bandwidth of the Gaussian
#'   kernel. Default is \code{1}.
#' @param d An integer-like scalar specifying the derivative order to return.
#'   \code{0} returns the smoothed Gaussian KDE function, \code{1} returns the
#'   first derivative, \code{2} returns the second derivative, and \code{3}
#'   returns the third derivative. Default is \code{0}.
#'
#' @return
#' A function that takes a numeric vector \code{x} and returns the estimated
#' Gaussian KDE-smoothed intensity or its derivative at \code{x}.
#'
#' @details
#' This function returns a closure. The returned function stores
#' \code{x_obs}, \code{y_obs}, \code{bw}, and \code{d} internally, so the user
#' does not need to provide the original spectrum again when evaluating the
#' function.
#'
#' The actual numerical evaluation is performed by internal Rcpp functions:
#' \code{cpp_gaussKDE()}, \code{cpp_gaussKDE_1st_deriv()},
#' \code{cpp_gaussKDE_2nd_deriv()}, or \code{cpp_gaussKDE_3rd_deriv()}.
#'
#' @examples
#' mz <- seq(100, 200, by = 1)
#' intensity <- dnorm(mz, mean = 150, sd = 10) * 100
#' spectrum <- data.frame(mz = mz, intensity = intensity)
#'
#' kde_fun <- get_gauss_kde(data = spectrum, bw = 1, d = 0)
#' kde_fun(c(140, 150, 160))
#'
#' d1_fun <- get_gauss_kde(data = spectrum, bw = 1, d = 1)
#' d1_fun(c(140, 150, 160))
#'
#' @export
get_gauss_kde <- function(data,
                          bw = 1,
                          d = 0) {
  
  # ------------------------------
  # 1. Input validation
  # ------------------------------
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("'data' must be a matrix or data frame.", call. = FALSE)
  }
  
  if (ncol(data) < 2L) {
    stop("'data' must have at least two columns: m/z and intensity.",
         call. = FALSE)
  }
  
  if (nrow(data) < 1L) {
    stop("'data' must contain at least one row.", call. = FALSE)
  }
  
  if (!is.numeric(bw) || length(bw) != 1L ||
      is.na(bw) || !is.finite(bw)) {
    stop("'bw' must be a single finite numeric value.", call. = FALSE)
  }
  
  if (bw <= 0) {
    stop("'bw' must be greater than 0.", call. = FALSE)
  }
  
  if (!is.numeric(d) || length(d) != 1L ||
      is.na(d) || !is.finite(d)) {
    stop("'d' must be a single finite numeric value.", call. = FALSE)
  }
  
  if (d != as.integer(d)) {
    stop("'d' must be an integer-like value.", call. = FALSE)
  }
  
  if (!d %in% 0:3) {
    stop("'d' must be one of 0, 1, 2, or 3.", call. = FALSE)
  }
  
  bw <- as.double(bw)
  d <- as.integer(d)
  
  x_obs <- as.double(data[, 1L])
  y_obs <- as.double(data[, 2L])
  
  if (anyNA(x_obs) || anyNA(y_obs)) {
    stop("'data' contains NA values in m/z or intensity.", call. = FALSE)
  }
  
  if (!all(is.finite(x_obs)) || !all(is.finite(y_obs))) {
    stop("'data' contains non-finite values in m/z or intensity.",
         call. = FALSE)
  }
  
  # ------------------------------
  # 2. Internal helper for evaluating x
  # ------------------------------
  validate_x <- function(x) {
    
    if (!is.numeric(x)) {
      stop("'x' must be a numeric vector.", call. = FALSE)
    }
    
    if (anyNA(x) || any(!is.finite(x))) {
      stop("'x' must contain only finite numeric values.", call. = FALSE)
    }
    
    as.double(x)
  }
  
  # ------------------------------
  # 3. Create closure function
  # ------------------------------
  fun <- switch(
    as.character(d),
    
    "0" = function(x) {
      x <- validate_x(x)
      
      cpp_gaussKDE(
        x = x,
        x_obs = x_obs,
        y_obs = y_obs,
        bw = bw
      )
    },
    
    "1" = function(x) {
      x <- validate_x(x)
      
      cpp_gaussKDE_1st_deriv(
        x = x,
        x_obs = x_obs,
        y_obs = y_obs,
        bw = bw
      )
    },
    
    "2" = function(x) {
      x <- validate_x(x)
      
      cpp_gaussKDE_2nd_deriv(
        x = x,
        x_obs = x_obs,
        y_obs = y_obs,
        bw = bw
      )
    },
    
    "3" = function(x) {
      x <- validate_x(x)
      
      cpp_gaussKDE_3rd_deriv(
        x = x,
        x_obs = x_obs,
        y_obs = y_obs,
        bw = bw
      )
    }
  )
  
  return(fun)
}