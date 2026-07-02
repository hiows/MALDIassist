#' @name preprocess_maldi_spectra
#' @title Preprocess MALDI-TOF Mass Spectra
#'
#' @description
#' Applies Savitzky-Golay smoothing and baseline subtraction to one or more
#' MALDI-TOF mass spectra.
#'
#' @param spectra A two-column `data.frame` or numeric `matrix` representing
#'   a single spectrum, or a `list` of such objects. The first column must
#'   contain strictly increasing m/z values and the second column must contain
#'   intensity values.
#' @param hws_sg A positive integer specifying the half-window size used for
#'   Savitzky-Golay smoothing. The full smoothing window size is
#'   `2 * hws_sg + 1`.
#' @param pno_sg A non-negative integer specifying the polynomial order used
#'   for Savitzky-Golay smoothing. It must be smaller than the full smoothing
#'   window size.
#' @param baseline_type A character string specifying the baseline subtraction
#'   method. One of `"snip"` or `"tophat"`.
#' @param iter_snip A positive integer specifying the number of iterations used
#'   for SNIP baseline estimation.
#' @param hws_tophat A positive integer specifying the half-window size used
#'   for TopHat baseline estimation.
#' @param n_cores Number of worker processes used when `spectra` is a list.
#'   The default is `1L`, which preserves sequential behaviour.
#'
#' @return A preprocessed spectrum or a list of preprocessed spectra. The
#'   returned object has the same list structure and names as `spectra`.
#'
#' @details
#' Each spectrum is smoothed using [smooth_savitzky_golay()] and subsequently
#' baseline-corrected using [subtract_baseline()]. Negative values produced
#' during baseline subtraction are truncated according to
#' `nonnegative_baseline = TRUE`.
#'
#' @examples
#' set.seed(1)
#' spectrum <- data.frame(
#'   mz = seq(2000, 8000, by = 1),
#'   intensity = runif(6001)
#' )
#'
#' pp_spectrum <- preprocess_maldi_spectra(spectrum)
#'
#' spectra <- list(sample_1 = spectrum, sample_2 = spectrum)
#'
#' pp_spectra <- preprocess_maldi_spectra(
#'   spectra = spectra,
#'   hws_sg = 10L,
#'   pno_sg = 3L,
#'   baseline_type = "snip",
#'   iter_snip = 50L
#' )
#'
#' @export
preprocess_maldi_spectra <- function(spectra,
                                     hws_sg = 10L,
                                     pno_sg = 3L,
                                     baseline_type = c("snip", "tophat"),
                                     iter_snip = 50L,
                                     hws_tophat = 50L,
                                     n_cores = 1L) {
  
  baseline_type <- match.arg(baseline_type)
  
  .assert_integer_scalar(
    x = hws_sg,
    arg_name = "hws_sg",
    lower_bound = 1L
  )
  
  .assert_integer_scalar(
    x = pno_sg,
    arg_name = "pno_sg",
    lower_bound = 0L
  )
  
  .assert_integer_scalar(
    x = iter_snip,
    arg_name = "iter_snip",
    lower_bound = 1L
  )
  
  .assert_integer_scalar(
    x = hws_tophat,
    arg_name = "hws_tophat",
    lower_bound = 1L
  )
  
  .assert_integer_scalar(
    x = n_cores,
    arg_name = "n_cores",
    lower_bound = 1L
  )
  
  window_size <- 2L * hws_sg + 1L
  
  if (pno_sg >= window_size) {
    stop(
      "'pno_sg' must be smaller than the full Savitzky-Golay window size ",
      "('2 * hws_sg + 1' = ", window_size, ").",
      call. = FALSE
    )
  }
  
  is_single_spectrum <- is.data.frame(spectra) || is.matrix(spectra)
  
  if (is_single_spectrum) {
    return(
      .preprocess_single_maldi_spectrum(
        spectrum = spectra,
        spectrum_label = "spectra",
        hws_sg = hws_sg,
        pno_sg = pno_sg,
        baseline_type = baseline_type,
        iter_snip = iter_snip,
        hws_tophat = hws_tophat
      )
    )
  }
  
  if (!is.list(spectra)) {
    stop(
      "'spectra' must be a two-column data.frame, a two-column numeric ",
      "matrix, or a list containing such objects.",
      call. = FALSE
    )
  }
  
  if (length(spectra) == 0L) {
    return(spectra)
  }
  
  pp_spectra <- .parallel_lapply(
    X = seq_along(spectra),
    FUN = function(i) {
      
      spectrum_name <- names(spectra)[i]
      
      if (is.null(spectrum_name) ||
          is.na(spectrum_name) ||
          !nzchar(spectrum_name)) {
        spectrum_label <- paste0("spectra[[", i, "]]")
      } else {
        spectrum_label <- paste0("spectra[['", spectrum_name, "']]")
      }
      
      .preprocess_single_maldi_spectrum(
        spectrum = spectra[[i]],
        spectrum_label = spectrum_label,
        hws_sg = hws_sg,
        pno_sg = pno_sg,
        baseline_type = baseline_type,
        iter_snip = iter_snip,
        hws_tophat = hws_tophat
      )
    },
    n_cores = n_cores
  )
  
  names(pp_spectra) <- names(spectra)
  
  pp_spectra
}
