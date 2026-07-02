#' @name subtract_baseline
#' @title Subtract Baseline from a Mass Spectrum
#'
#' @description
#' Estimates and subtracts the baseline from a two-column mass spectrum using
#' either the SNIP or TopHat baseline correction method.
#'
#' The first column of \code{data} is treated as the x-axis, usually m/z values,
#' and the second column is treated as the signal intensity. The original column
#' names of the first two columns in \code{data} are preserved in the returned
#' raw and baseline-subtracted spectra.
#'
#' @param data A matrix or data frame containing a two-column spectrum. The
#'   first column is assumed to contain x values, usually m/z values, and the
#'   second column is assumed to contain y values, usually intensity values.
#' @param baseline_type Character. Baseline correction method. One of
#'   \code{"snip"} or \code{"tophat"}. Default is \code{"snip"}.
#' @param iter_snip Integer-like numeric value. Number of iterations used for
#'   SNIP baseline estimation. Used only when \code{baseline_type = "snip"}.
#'   Default is \code{50L}.
#' @param hws_tophat Integer-like numeric value. Half-window size used for
#'   TopHat baseline estimation. Used only when
#'   \code{baseline_type = "tophat"}. Default is \code{50L}.
#' @param nonnegative_baseline Logical. If \code{TRUE}, the estimated baseline
#'   is constrained to be non-negative inside the C++ baseline estimation
#'   backend. Default is \code{TRUE}.
#'
#' @return A list with four elements:
#'   \itemize{
#'     \item \code{raw_data}: A data frame containing the original spectrum.
#'     Column names are inherited from the first two columns of \code{data}.
#'     \item \code{subtracted_data}: A data frame containing the
#'     baseline-subtracted spectrum. Column names are inherited from the first
#'     two columns of \code{data}.
#'     \item \code{baseline}: A numeric vector containing the estimated baseline.
#'     \item \code{param}: A list containing the parameters used for baseline
#'     subtraction.
#'   }
#'
#' @details
#' This function is an R wrapper around the C++ baseline estimation functions
#' \code{cpp_SNIP()} and \code{cpp_TopHat()}.
#'
#' If \code{baseline_type = "snip"}, the baseline is estimated using
#' \code{cpp_SNIP()}. If \code{baseline_type = "tophat"}, the baseline is
#' estimated using \code{cpp_TopHat()}.
#'
#' The argument \code{nonnegative_baseline} controls whether the estimated
#' baseline itself is constrained to be non-negative. This is different from
#' the final baseline-subtracted signal, which is clipped to be non-negative by
#' applying:
#'
#' \deqn{
#' y_{corrected} = max(y - baseline, 0)
#' }
#'
#' Therefore, the returned baseline-subtracted intensity is always
#' non-negative.
#'
#' This function assumes that \code{cpp_SNIP()} and \code{cpp_TopHat()} are
#' available.
#'
#' @examples
#' x <- seq(1000, 2000, length.out = 5000)
#' y <- dnorm(x, mean = 1500, sd = 30) * 100 +
#'      seq(20, 50, length.out = length(x))
#'
#' data <- data.frame(mz = x, intensity = y)
#'
#' result <- subtract_baseline(
#'   data = data,
#'   baseline_type = "snip",
#'   iter_snip = 50,
#'   nonnegative_baseline = TRUE
#' )
#'
#' head(result$subtracted_data)
#'
#' @export
subtract_baseline <- function(data,
                              baseline_type = c("snip", "tophat"),
                              iter_snip = 50L,
                              hws_tophat = 50L,
                              nonnegative_baseline = TRUE) {
  
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
  
  baseline_type <- match.arg(baseline_type)
  
  if (!is.numeric(iter_snip) ||
      length(iter_snip) != 1L ||
      is.na(iter_snip) ||
      !is.finite(iter_snip)) {
    stop("'iter_snip' must be a single finite numeric value.",
         call. = FALSE)
  }
  
  if (!is.numeric(hws_tophat) ||
      length(hws_tophat) != 1L ||
      is.na(hws_tophat) ||
      !is.finite(hws_tophat)) {
    stop("'hws_tophat' must be a single finite numeric value.",
         call. = FALSE)
  }
  
  if (iter_snip < 1L) {
    stop("'iter_snip' must be greater than or equal to 1.",
         call. = FALSE)
  }
  
  if (hws_tophat < 1L) {
    stop("'hws_tophat' must be greater than or equal to 1.",
         call. = FALSE)
  }
  
  if (iter_snip != as.integer(iter_snip)) {
    stop("'iter_snip' must be an integer-like value.", call. = FALSE)
  }
  
  if (hws_tophat != as.integer(hws_tophat)) {
    stop("'hws_tophat' must be an integer-like value.", call. = FALSE)
  }
  
  if (!is.logical(nonnegative_baseline) ||
      length(nonnegative_baseline) != 1L ||
      is.na(nonnegative_baseline)) {
    stop("'nonnegative_baseline' must be a single non-missing logical value.",
         call. = FALSE)
  }
  
  iter_snip <- as.integer(iter_snip)
  hws_tophat <- as.integer(hws_tophat)
  
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
  
  input_colnames <- colnames(data)[1:2]
  
  if (is.null(input_colnames) ||
      anyNA(input_colnames) ||
      any(input_colnames == "")) {
    input_colnames <- c("x", "y")
  }
  
  raw_data <- data.frame(
    x = x,
    y = y,
    stringsAsFactors = FALSE
  )
  
  colnames(raw_data) <- input_colnames
  rownames(raw_data) <- NULL
  
  # ------------------------------
  # 2. Estimate baseline
  # ------------------------------
  baseline <- switch(
    baseline_type,
    snip = {
      tryCatch(
        cpp_SNIP(
          intensity = y,
          iterations = iter_snip,
          decreasing = TRUE,
          nonnegative = nonnegative_baseline
        ),
        error = function(e) {
          stop("'cpp_SNIP()' failed.\n",
               "Original error: ", conditionMessage(e),
               call. = FALSE)
        }
      )
    },
    tophat = {
      tryCatch(
        cpp_TopHat(
          intensity = y,
          half_window = hws_tophat,
          nonnegative = nonnegative_baseline
        ),
        error = function(e) {
          stop("'cpp_TopHat()' failed.\n",
               "Original error: ", conditionMessage(e),
               call. = FALSE)
        }
      )
    }
  )
  
  # ------------------------------
  # 3. Validate baseline
  # ------------------------------
  if (!is.numeric(baseline)) {
    stop("The estimated baseline must be a numeric vector.", call. = FALSE)
  }
  
  if (length(baseline) != length(y)) {
    stop("The estimated baseline must have the same length as the intensity vector.",
         call. = FALSE)
  }
  
  if (anyNA(baseline) || any(!is.finite(baseline))) {
    stop("The estimated baseline contains NA, NaN, Inf, or -Inf values.",
         call. = FALSE)
  }
  
  # ------------------------------
  # 4. Subtract baseline
  # ------------------------------
  corrected_y <- pmax(y - baseline, 0)
  
  subtracted_data <- data.frame(
    x = x,
    y = corrected_y,
    stringsAsFactors = FALSE
  )
  
  colnames(subtracted_data) <- input_colnames
  rownames(subtracted_data) <- NULL
  
  # ------------------------------
  # 5. Return result
  # ------------------------------
  result <- list(
    raw_data = raw_data,
    subtracted_data = subtracted_data,
    baseline = baseline,
    param = list(
      baseline_type = baseline_type,
      iter_snip = iter_snip,
      hws_tophat = hws_tophat,
      nonnegative_baseline = nonnegative_baseline
    )
  )
  
  return(result)
}