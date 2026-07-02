#' MALDIassist: Mathematical utilities for MALDI-TOF mass spectrometry
#'
#' Tools for loading Bruker MALDI-TOF spectra, preprocessing, Gaussian KDE-based
#' peak detection, peak-quality metrics, and cohort feature analysis workflows.
#'
#' @section Main functions:
#' * [load_maldi_spectra()]: load Bruker raw spectra
#' * [preprocess_maldi_spectra()]: smooth and baseline-correct spectra
#' * [find_peaks()] / [find_peaks_spectra()]: detect ordinary and shoulder peaks
#'   in one spectrum or a list of spectra
#' * [find_peaks_fast()] / [find_peaks_spectra_fast()]: fast local peak detection
#'   in one spectrum or a list of spectra
#' * [filter_peaks()] / [filter_peaks_spectra()]: filter detected peaks by
#'   intensity, prominence, and strength
#' * [find_frequent_mz()]: find frequent m/z values across a cohort
#' * [build_matched_peaks_matrix()]: assemble a cohort peak intensity matrix
#'
#' @section Suggested packages:
#' Some functions require optional packages that are only loaded when used:
#' * [plot_spectrum()], [plot_spectra()]: `ggplot2`
#' * [heatmap_spectrum_matrix()]: `pheatmap`
#'
#' @useDynLib MALDIassist, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom rlang .data
#' @importFrom stats hclust dist complete.cases
"_PACKAGE"
