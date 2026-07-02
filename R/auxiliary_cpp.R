# Internal roxygen documentation for C++ backend functions (see RcppExports.R).

#' SNIP baseline estimation (C++ backend)
#'
#' @param intensity Numeric intensity vector.
#' @param iterations Number of SNIP iterations.
#' @param decreasing Logical. Use decreasing clipping windows when `TRUE`.
#' @param nonnegative Logical. Clip negative baseline values to zero.
#' @return Numeric vector of estimated baseline values.
#' @keywords internal
#' @name cpp_SNIP
#' @aliases cpp_SNIP
#' @noRd
NULL

#' SNIP baseline estimation for a two-column spectrum (C++ backend)
#'
#' @param spectrum Numeric matrix with x and intensity columns.
#' @param iterations Number of SNIP iterations.
#' @param decreasing Logical. Use decreasing clipping windows when `TRUE`.
#' @param nonnegative Logical. Clip negative baseline values to zero.
#' @param corrected_nonnegative Logical. Truncate subtracted intensities at zero.
#' @return A data frame with baseline-corrected intensities.
#' @keywords internal
#' @name cpp_SNIP_spectrum
#' @aliases cpp_SNIP_spectrum
#' @noRd
NULL

#' TopHat baseline estimation (C++ backend)
#'
#' @param intensity Numeric intensity vector.
#' @param half_window Half-window size for moving minimum/maximum filters.
#' @param nonnegative Logical. Clip negative baseline values to zero.
#' @return Numeric vector of estimated baseline values.
#' @keywords internal
#' @name cpp_TopHat
#' @aliases cpp_TopHat
#' @noRd
NULL

#' TopHat baseline estimation for a two-column spectrum (C++ backend)
#'
#' @param spectrum Numeric matrix with x and intensity columns.
#' @param half_window Half-window size for moving minimum/maximum filters.
#' @param nonnegative Logical. Clip negative baseline values to zero.
#' @param corrected_nonnegative Logical. Truncate subtracted intensities at zero.
#' @return A data frame with baseline-corrected intensities.
#' @keywords internal
#' @name cpp_TopHat_spectrum
#' @aliases cpp_TopHat_spectrum
#' @noRd
NULL

#' Calculate peak prominence (C++ backend)
#'
#' @param data Numeric matrix with x and intensity columns.
#' @param peaks Numeric matrix with peak positions and intensities.
#' @param valley_type Valley selection method: `"higher"` or `"lower"`.
#' @param zero_tol Tolerance below which prominence is treated as zero.
#' @return Numeric vector of peak prominences.
#' @keywords internal
#' @name cpp_calculate_peak_prominence
#' @aliases cpp_calculate_peak_prominence
#' @noRd
NULL

#' Estimate peak strength (C++ backend)
#'
#' @param data Numeric matrix or data frame with x and intensity columns.
#' @param peaks Numeric matrix or data frame with peak x and height columns.
#' @param k Non-negative MAD multiplier.
#' @param normalization Intensity transform: `"raw"`, `"sqrt"`, or `"log10"`.
#' @return Numeric vector of peak strength scores in `[0, 1]`, `NA` when local
#'   valleys cannot be defined.
#' @keywords internal
#' @name cpp_estimate_peak_strength
#' @aliases cpp_estimate_peak_strength
#' @noRd
NULL

#' Fast local peak detection (C++ backend)
#'
#' @param data Numeric matrix with x and intensity columns.
#' @param hws_peaks Half-window size in x-axis units.
#' @param plateau Plateau handling code: `0` middle, `1` first, `2` last.
#' @return Numeric matrix of detected peaks.
#' @keywords internal
#' @name cpp_find_peaks_fast
#' @aliases cpp_find_peaks_fast
#' @noRd
NULL

#' Savitzky-Golay smoothing filter (C++ backend)
#'
#' @param y Numeric intensity vector.
#' @param hws Half-window size.
#' @param pno Polynomial order.
#' @return Smoothed numeric vector.
#' @keywords internal
#' @name cpp_savitzkyGolay_filter
#' @aliases cpp_savitzkyGolay_filter
#' @noRd
NULL

#' Gaussian kernel density estimate (C++ backend)
#'
#' @param x Evaluation points.
#' @param x_obs Observed x values.
#' @param y_obs Observed intensities.
#' @param bw Gaussian kernel bandwidth.
#' @return Estimated KDE values at `x`.
#' @keywords internal
#' @name cpp_gaussKDE
#' @aliases cpp_gaussKDE
#' @noRd
NULL

#' First derivative of Gaussian KDE (C++ backend)
#'
#' @inheritParams cpp_gaussKDE
#' @return First derivative values at `x`.
#' @keywords internal
#' @name cpp_gaussKDE_1st_deriv
#' @aliases cpp_gaussKDE_1st_deriv
#' @noRd
NULL

#' Second derivative of Gaussian KDE (C++ backend)
#'
#' @inheritParams cpp_gaussKDE
#' @return Second derivative values at `x`.
#' @keywords internal
#' @name cpp_gaussKDE_2nd_deriv
#' @aliases cpp_gaussKDE_2nd_deriv
#' @noRd
NULL

#' Third derivative of Gaussian KDE (C++ backend)
#'
#' @inheritParams cpp_gaussKDE
#' @return Third derivative values at `x`.
#' @keywords internal
#' @name cpp_gaussKDE_3rd_deriv
#' @aliases cpp_gaussKDE_3rd_deriv
#' @noRd
NULL

#' Gaussian KDE and derivatives (C++ backend)
#'
#' @inheritParams cpp_gaussKDE
#' @return Numeric matrix with columns `kde`, `d1`, `d2`, and `d3`.
#' @keywords internal
#' @name cpp_gaussKDE_all
#' @aliases cpp_gaussKDE_all
#' @noRd
NULL

#' Bisection root finder (C++ backend)
#'
#' @param fun Function to evaluate.
#' @param x_left Left interval endpoint.
#' @param x_right Right interval endpoint.
#' @param yi Target y-value.
#' @param tol Numerical tolerance.
#' @param max_iter Maximum number of bisection iterations.
#' @return Estimated root position.
#' @keywords internal
#' @name cpp_bisection_root
#' @aliases cpp_bisection_root
#' @noRd
NULL

#' Root finder on a function grid (C++ backend)
#'
#' @inheritParams cpp_bisection_root
#' @param x Strictly monotonic numeric grid.
#' @return Numeric vector of estimated roots.
#' @keywords internal
#' @name cpp_find_roots
#' @aliases cpp_find_roots
#' @noRd
NULL

#' Root finder on a precomputed derivative grid (C++ backend)
#'
#' @param x Strictly monotonic numeric grid.
#' @param y_grid Precomputed function values on `x`.
#' @param x_obs Observed x values used for KDE refinement.
#' @param y_obs Observed intensities used for KDE refinement.
#' @param bw Gaussian kernel bandwidth.
#' @param yi Target y-value.
#' @param tol Numerical tolerance.
#' @param max_iter Maximum number of bisection iterations.
#' @return Numeric vector of estimated roots.
#' @keywords internal
#' @name cpp_find_roots_on_grid
#' @aliases cpp_find_roots_on_grid
#' @noRd
NULL

#' Extrema finder from a precomputed first-derivative grid (C++ backend)
#'
#' @param x Strictly monotonic numeric grid.
#' @param d1_grid Precomputed first derivative on `x`.
#' @param x_obs Observed x values used for KDE refinement.
#' @param y_obs Observed intensities used for KDE refinement.
#' @param bw Gaussian kernel bandwidth.
#' @param tol Numerical tolerance.
#' @param max_iter Maximum number of bisection iterations.
#' @return List with `local_min`, `local_max`, and `plateau`.
#' @keywords internal
#' @name cpp_find_extrema_from_grid
#' @aliases cpp_find_extrema_from_grid
#' @noRd
NULL

#' Curvature from derivative functions (C++ backend)
#'
#' @param x Evaluation points.
#' @param first_deriv First derivative function.
#' @param second_deriv Second derivative function.
#' @param absolute Logical. Use absolute second derivative when `TRUE`.
#' @return Numeric vector of curvature values.
#' @keywords internal
#' @name cpp_curvature
#' @aliases cpp_curvature
#' @noRd
NULL

#' Curvature from precomputed derivatives (C++ backend)
#'
#' @param d1 Precomputed first derivative values.
#' @param d2 Precomputed second derivative values.
#' @param absolute Logical. Use absolute second derivative when `TRUE`.
#' @return Numeric vector of curvature values.
#' @keywords internal
#' @name cpp_curvature_from_derivs
#' @aliases cpp_curvature_from_derivs
#' @noRd
NULL
