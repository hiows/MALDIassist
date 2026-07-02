# Internal preprocessing helpers for MALDI-TOF spectra.

#' Preprocess a single MALDI-TOF mass spectrum
#'
#' @param spectrum A two-column `data.frame` or numeric `matrix`.
#' @param spectrum_label A character string used to identify the spectrum in
#'   error messages.
#' @param hws_sg A positive integer specifying the Savitzky-Golay half-window
#'   size.
#' @param pno_sg A non-negative integer specifying the Savitzky-Golay
#'   polynomial order.
#' @param baseline_type A baseline subtraction method.
#' @param iter_snip A positive integer specifying the number of SNIP iterations.
#' @param hws_tophat A positive integer specifying the TopHat half-window size.
#'
#' @return A baseline-subtracted spectrum.
#'
#' @keywords internal
#' @noRd
.preprocess_single_maldi_spectrum <- function(spectrum,
                                              spectrum_label,
                                              hws_sg,
                                              pno_sg,
                                              baseline_type,
                                              iter_snip,
                                              hws_tophat) {
  
  spectrum <- .validate_maldi_spectrum(
    spectrum = spectrum,
    spectrum_label = spectrum_label,
    min_nrow = 2L * hws_sg + 1L
  )
  
  smoothed_spectrum <- tryCatch(
    smooth_savitzky_golay(
      data = spectrum,
      hws = hws_sg,
      pno = pno_sg
    ),
    error = function(e) {
      stop(
        "Savitzky-Golay smoothing failed for ", spectrum_label, ": ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  subtraction_result <- tryCatch(
    subtract_baseline(
      data = smoothed_spectrum,
      baseline_type = baseline_type,
      iter_snip = iter_snip,
      hws_tophat = hws_tophat,
      nonnegative_baseline = TRUE
    ),
    error = function(e) {
      stop(
        "Baseline subtraction failed for ", spectrum_label, ": ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  if (!is.list(subtraction_result) ||
      is.null(subtraction_result$subtracted_data)) {
    stop(
      "'subtract_baseline()' must return a list containing ",
      "'subtracted_data'. The error occurred while processing ",
      spectrum_label, ".",
      call. = FALSE
    )
  }
  
  subtraction_result$subtracted_data
}


#' Validate a MALDI-TOF mass spectrum
#'
#' @param spectrum A spectrum object.
#' @param spectrum_label A character string used in error messages.
#' @param min_nrow The minimum required number of observations.
#'
#' @return A validated two-column `data.frame`.
#'
#' @keywords internal
#' @noRd
.validate_maldi_spectrum <- function(spectrum,
                                     spectrum_label,
                                     min_nrow) {
  
  if (!is.data.frame(spectrum) && !is.matrix(spectrum)) {
    stop(
      spectrum_label,
      " must be a two-column data.frame or numeric matrix.",
      call. = FALSE
    )
  }
  
  if (ncol(spectrum) != 2L) {
    stop(
      spectrum_label,
      " must contain exactly two columns: m/z and intensity.",
      call. = FALSE
    )
  }
  
  if (nrow(spectrum) < min_nrow) {
    stop(
      spectrum_label,
      " must contain at least ", min_nrow,
      " observations when 'hws_sg' is used.",
      call. = FALSE
    )
  }
  
  spectrum <- as.data.frame(spectrum)
  
  mz <- spectrum[[1L]]
  intensity <- spectrum[[2L]]
  
  if (!is.numeric(mz) || !is.numeric(intensity)) {
    stop(
      "The first and second columns of ",
      spectrum_label,
      " must both be numeric.",
      call. = FALSE
    )
  }
  
  if (anyNA(mz) || anyNA(intensity)) {
    stop(
      spectrum_label,
      " must not contain missing values.",
      call. = FALSE
    )
  }
  
  if (any(!is.finite(mz)) || any(!is.finite(intensity))) {
    stop(
      spectrum_label,
      " must not contain infinite or non-finite values.",
      call. = FALSE
    )
  }
  
  if (any(diff(mz) <= 0)) {
    stop(
      "The m/z values in ",
      spectrum_label,
      " must be strictly increasing without duplicated values.",
      call. = FALSE
    )
  }
  
  spectrum
}


#' Validate an integer-like scalar argument
#'
#' @param x An object to validate.
#' @param arg_name The argument name used in the error message.
#' @param lower_bound The minimum permitted value.
#'
#' @return `invisible(TRUE)` when the input is valid.
#'
#' @keywords internal
#' @noRd
.assert_integer_scalar <- function(x,
                                   arg_name,
                                   lower_bound) {
  
  is_valid <- is.numeric(x) &&
    length(x) == 1L &&
    !is.na(x) &&
    is.finite(x) &&
    x == as.integer(x) &&
    x >= lower_bound
  
  if (!is_valid) {
    stop(
      "'", arg_name, "' must be a single integer greater than or equal to ",
      lower_bound, ".",
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}
