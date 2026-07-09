#' @name find_peaks
#' @title Detect Ordinary Peaks and Shoulder Peaks Using Weighted Curvature
#' 
#' @description
#' Detects ordinary peaks and shoulder peaks from a smoothed spectrum.
#' Ordinary peaks are identified from the local maxima of a Gaussian kernel
#' density estimate (KDE). Shoulder peaks are inferred from strong local maxima
#' of the weighted reverse-signed curvature that are not located near ordinary
#' KDE peaks.
#'
#' @param data A numeric data frame or matrix with at least two columns.
#'   The first column must contain strictly increasing x-axis values, such as
#'   m/z values. The second column must contain intensity values.
#'
#' @param bw A positive numeric scalar specifying the bandwidth of the Gaussian
#'   KDE. The default is the median interval between adjacent x-axis values.
#'
#' @param hws_peaks A positive numeric scalar specifying the half-window size
#'   used by [find_peaks_fast()] to filter neighboring local maxima from the
#'   curvature profile. The value is interpreted as an x-axis distance. For
#'   MALDI-TOF MS spectra, the unit is typically m/z. For each local maximum,
#'   only the highest local maximum within
#'   `[x_peak - hws_peaks, x_peak + hws_peaks]` is retained.
#'
#' @param merge_tol A positive numeric scalar specifying the maximum x-axis
#'   distance within which two detected candidates are treated as nearby
#'   peaks. For MALDI-TOF MS spectra, the unit is typically m/z. If `NULL`,
#'   `hws_peaks` is used.
#'
#' @param tol A positive numeric scalar specifying the numerical tolerance used
#'   by [find_extrema()].
#'
#' @param max_iter A positive integer specifying the maximum number of
#'   iterations used by [find_extrema()].
#'
#' @param weight_type A character string specifying the intensity-weighting
#'   method applied to the reverse-signed curvature. One of:
#'
#'   - `"none"`: do not apply intensity weighting.
#'   - `"raw"`: multiply the curvature by the non-negative KDE intensity.
#'   - `"sqrt"`: multiply the curvature by the square root of the non-negative
#'     KDE intensity.
#'   - `"log10"`: multiply the curvature by
#'     `log10(non-negative KDE intensity + 1)`.
#'
#' @param hws_grid_kappa_smooth A numeric vector of positive integers
#'   specifying the half-window sizes of the Savitzky-Golay filters applied to
#'   the weighted curvature profile. Each value must be at least `2`, because a
#'   third-order polynomial is used.
#'
#' @param cutoff_kappa_peak_strength A non-negative numeric scalar specifying
#'   the upper bound of the curvature-peak-strength threshold. A curvature peak
#'   is retained as a shoulder candidate only when its strength is greater than
#'   the final cutoff.
#'
#' @param peak_retention_fraction A numeric scalar in the interval `(0, 1]`.
#'   The curvature-peak-strength cutoff is calculated as the smaller value
#'   between `cutoff_kappa_peak_strength` and the
#'   `(1 - peak_retention_fraction)` quantile of the curvature peak strengths.
#'
#' @details
#' The function performs the following steps:
#'
#' 1. Construct a Gaussian KDE and its first and second derivatives.
#' 2. Calculate the signed curvature of the KDE.
#' 3. Retain only the negative signed-curvature component and optionally apply
#'    intensity weighting.
#' 4. Smooth the weighted curvature using multiple Savitzky-Golay filters and
#'    average the resulting profiles.
#' 5. Detect strong local maxima from the averaged curvature profile as
#'    shoulder-peak candidates.
#' 6. Detect ordinary KDE peaks from derivative-based extrema.
#' 7. Remove curvature candidates that are located near ordinary KDE peaks.
#' 8. Remove weaker candidates located near stronger candidates.
#'
#' The reverse-signed curvature is defined as:
#'
#' `abs(pmin(0, curvature))`
#'
#' Therefore, only the concave-down component of the signed curvature
#' contributes to shoulder-peak detection.
#'
#' `hws_peaks` and `merge_tol` are both expressed as x-axis distances but serve
#' different purposes. `hws_peaks` controls the local-maximum filtering window
#' used by [find_peaks_fast()], whereas `merge_tol` controls the distance used
#' to merge nearby ordinary and shoulder-peak candidates.
#'
#' @return A data frame with three columns:
#'
#' - The first column contains the x-axis positions of the detected peaks. Its
#'   name is inherited from the first column of `data`.
#' - The second column contains KDE-estimated intensities. Its name is
#'   inherited from the second column of `data`.
#' - `type` indicates whether each detected feature is an ordinary KDE peak
#'   (`"peak"`) or a curvature-derived shoulder peak (`"shoulder"`).
#'
#' The rows are sorted in ascending order of the x-axis values.
#'
#' @seealso [find_peaks_fast()], [find_extrema()],
#'   [estimate_peak_strength()]
#'
#' @examples
#' x <- seq(1000, 2000, length.out = 2000)
#' y <- dnorm(x, mean = 1400, sd = 20) * 100 +
#'      dnorm(x, mean = 1450, sd = 10) * 30 +
#'      dnorm(x, mean = 1700, sd = 25) * 80
#' spectrum <- data.frame(mz = x, intensity = y)
#'
#' peaks <- find_peaks(
#'   data = spectrum,
#'   hws_peaks = 10,
#'   weight_type = "raw",
#'   hws_grid_kappa_smooth = 3:20,
#'   cutoff_kappa_peak_strength = 0.5,
#'   peak_retention_fraction = 0.25
#' )
#'
#' head(peaks)
#'
#' @export
find_peaks <- function(data,
                       bw = stats::median(diff(data[, 1])),
                       hws_peaks = 10,
                       merge_tol = NULL,
                       tol = 1e-5,
                       max_iter = 100L,
                       weight_type = c("raw", "sqrt", "log10", "none"),
                       hws_grid_kappa_smooth = 3:20,
                       cutoff_kappa_peak_strength = 0.5,
                       peak_retention_fraction = 0.25) {
  
  # -------------------------------------------------------------------------
  # 0. Validate input data
  # -------------------------------------------------------------------------
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop(
      "'data' must be a data frame or matrix.",
      call. = FALSE
    )
  }
  
  if (ncol(data) < 2L) {
    stop(
      "'data' must have at least two columns: x-axis values and intensities.",
      call. = FALSE
    )
  }
  
  if (nrow(data) < 3L) {
    stop(
      "'data' must contain at least three rows.",
      call. = FALSE
    )
  }
  
  x <- data[, 1]
  y <- data[, 2]
  
  if (!is.numeric(x) || !is.numeric(y)) {
    stop(
      "The first and second columns of 'data' must be numeric.",
      call. = FALSE
    )
  }
  
  if (anyNA(x) || anyNA(y)) {
    stop(
      "The first and second columns of 'data' must not contain NA values.",
      call. = FALSE
    )
  }
  
  if (any(!is.finite(x)) || any(!is.finite(y))) {
    stop(
      "The first and second columns of 'data' must contain only finite values.",
      call. = FALSE
    )
  }
  
  if (is.unsorted(x, strictly = TRUE)) {
    stop(
      "The x-axis values in the first column of 'data' must be strictly increasing.",
      call. = FALSE
    )
  }
  
  # -------------------------------------------------------------------------
  # 1. Validate scalar parameters
  # -------------------------------------------------------------------------
  if (!is.numeric(bw) ||
      length(bw) != 1L ||
      is.na(bw) ||
      !is.finite(bw) ||
      bw <= 0) {
    stop(
      "'bw' must be a positive finite numeric scalar.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(hws_peaks) ||
      length(hws_peaks) != 1L ||
      is.na(hws_peaks) ||
      !is.finite(hws_peaks) ||
      hws_peaks <= 0) {
    stop(
      "'hws_peaks' must be a positive finite numeric scalar.",
      call. = FALSE
    )
  }
  
  if (is.null(merge_tol)) {
    merge_tol <- hws_peaks
  }
  
  if (!is.numeric(merge_tol) ||
      length(merge_tol) != 1L ||
      is.na(merge_tol) ||
      !is.finite(merge_tol) ||
      merge_tol <= 0) {
    stop(
      "'merge_tol' must be a positive finite numeric scalar or NULL.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(tol) ||
      length(tol) != 1L ||
      is.na(tol) ||
      !is.finite(tol) ||
      tol <= 0) {
    stop(
      "'tol' must be a positive finite numeric scalar.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(max_iter) ||
      length(max_iter) != 1L ||
      is.na(max_iter) ||
      !is.finite(max_iter) ||
      max_iter < 1 ||
      max_iter != as.integer(max_iter)) {
    stop(
      "'max_iter' must be a positive integer.",
      call. = FALSE
    )
  }
  
  max_iter <- as.integer(max_iter)
  
  weight_type <- match.arg(weight_type)
  
  if (!is.numeric(cutoff_kappa_peak_strength) ||
      length(cutoff_kappa_peak_strength) != 1L ||
      is.na(cutoff_kappa_peak_strength) ||
      !is.finite(cutoff_kappa_peak_strength) ||
      cutoff_kappa_peak_strength < 0) {
    stop(
      "'cutoff_kappa_peak_strength' must be a non-negative finite numeric scalar.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(peak_retention_fraction) ||
      length(peak_retention_fraction) != 1L ||
      is.na(peak_retention_fraction) ||
      !is.finite(peak_retention_fraction) ||
      peak_retention_fraction <= 0 ||
      peak_retention_fraction > 1) {
    stop(
      "'peak_retention_fraction' must be in the interval (0, 1].",
      call. = FALSE
    )
  }
  
  # -------------------------------------------------------------------------
  # 2. Validate smoothing-window grid
  # -------------------------------------------------------------------------
  if (!is.numeric(hws_grid_kappa_smooth) ||
      length(hws_grid_kappa_smooth) == 0L ||
      anyNA(hws_grid_kappa_smooth) ||
      any(!is.finite(hws_grid_kappa_smooth)) ||
      any(hws_grid_kappa_smooth != as.integer(hws_grid_kappa_smooth)) ||
      any(hws_grid_kappa_smooth < 2)) {
    stop(
      paste0(
        "'hws_grid_kappa_smooth' must contain integers greater than or ",
        "equal to 2. A minimum half-window size of 2 is required because ",
        "a third-order Savitzky-Golay polynomial is used."
      ),
      call. = FALSE
    )
  }
  
  hws_grid_kappa_smooth <- sort(unique(as.integer(hws_grid_kappa_smooth)))
  
  sg_window_sizes <- 2L * hws_grid_kappa_smooth + 1L
  
  if (any(sg_window_sizes > nrow(data))) {
    stop(
      paste0(
        "Every Savitzky-Golay window size, calculated as ",
        "'2 * hws_grid_kappa_smooth + 1', must be less than or equal to ",
        "the number of rows in 'data'."
      ),
      call. = FALSE
    )
  }
  
  # -------------------------------------------------------------------------
  # 3. Define output-column names and an empty-result constructor
  # -------------------------------------------------------------------------
  input_colnames <- colnames(data)
  
  if (is.null(input_colnames) || length(input_colnames) < 2L) {
    output_colnames <- c("x", "y")
  } else {
    output_colnames <- input_colnames[1:2]
    
    if (anyNA(output_colnames) || any(output_colnames == "")) {
      output_colnames <- c("x", "y")
    }
  }
  
  create_empty_result <- function() {
    out <- data.frame(
      x = numeric(0),
      y = numeric(0),
      type = character(0),
      stringsAsFactors = FALSE
    )
    
    colnames(out)[1:2] <- output_colnames
    
    return(out)
  }
  
  # -------------------------------------------------------------------------
  # 4. Evaluate Gaussian KDE and derivatives on the full x-grid
  # -------------------------------------------------------------------------
  x_obs <- as.double(x)
  y_obs <- as.double(y)
  
  kde_all <- cpp_gaussKDE_all(
    x = x_obs,
    x_obs = x_obs,
    y_obs = y_obs,
    bw = bw
  )
  
  y_est <- kde_all[, "kde"]
  d1_grid <- kde_all[, "d1"]
  d2_grid <- kde_all[, "d2"]
  
  # -------------------------------------------------------------------------
  # 5. Calculate the weighted reverse-signed curvature
  # -------------------------------------------------------------------------
  
  if (!is.numeric(y_est) ||
      length(y_est) != length(x) ||
      anyNA(y_est) ||
      any(!is.finite(y_est))) {
    stop(
      "'cpp_gaussKDE_all()' must return one finite KDE value per x-axis value.",
      call. = FALSE
    )
  }
  
  weights <- switch(
    weight_type,
    none = rep(1, length(y_est)),
    raw = pmax(y_est, 0),
    sqrt = sqrt(pmax(y_est, 0)),
    log10 = log10(pmax(y_est, 0) + 1)
  )
  
  curvature <- cpp_curvature_from_derivs(
    d1 = d1_grid,
    d2 = d2_grid,
    absolute = FALSE
  )
  
  if (!is.numeric(curvature) ||
      length(curvature) != length(x) ||
      anyNA(curvature) ||
      any(!is.finite(curvature))) {
    stop(
      "'cpp_curvature_from_derivs()' must return one finite numeric value per x-axis value.",
      call. = FALSE
    )
  }
  
  weighted_kappa <- abs(pmin(0, curvature)) * weights
  
  # -------------------------------------------------------------------------
  # 6. Smooth and average the weighted-curvature profiles
  # -------------------------------------------------------------------------
  smoothed_kappa <- lapply(
    hws_grid_kappa_smooth,
    function(hws_i) {
      current_kappa <- cpp_savitzkyGolay_filter(
        y = weighted_kappa,
        hws = hws_i,
        pno = 3
      )
      
      if (!is.numeric(current_kappa) ||
          length(current_kappa) != length(weighted_kappa) ||
          anyNA(current_kappa) ||
          any(!is.finite(current_kappa))) {
        stop(
          paste0(
            "'cpp_savitzkyGolay_filter()' returned an invalid result ",
            "for hws = ",
            hws_i,
            "."
          ),
          call. = FALSE
        )
      }
      
      return(pmax(0, current_kappa))
    }
  )
  
  avg_kappa <- colMeans(
    do.call(
      rbind,
      c(list(weighted_kappa), smoothed_kappa)
    )
  )
  
  data_kappa <- data.frame(
    x = x,
    kappa = avg_kappa
  )
  
  # -------------------------------------------------------------------------
  # 7. Detect shoulder-peak candidates from the curvature profile
  # -------------------------------------------------------------------------
  peaks_kappa <- find_peaks_fast(
    data = data_kappa,
    hws_peaks = hws_peaks,
    plateau = "middle",
    na_rm = TRUE
  )
  
  if (!is.data.frame(peaks_kappa) && !is.matrix(peaks_kappa)) {
    stop(
      "'find_peaks_fast()' must return a data frame or matrix.",
      call. = FALSE
    )
  }
  
  if (ncol(peaks_kappa) < 2L) {
    stop(
      "'find_peaks_fast()' must return at least two columns.",
      call. = FALSE
    )
  }
  
  if (nrow(peaks_kappa) == 0L) {
    
    x_shoulder_candidates <- numeric(0)
    
  } else {
    
    kappa_peaks_strength <- estimate_peak_strength(
      data = data_kappa,
      peaks = peaks_kappa,
      k = 1,
      normalization_type = "raw"
    )
    
    if (!is.numeric(kappa_peaks_strength) ||
        length(kappa_peaks_strength) != nrow(peaks_kappa)) {
      stop(
        paste0(
          "'estimate_peak_strength()' must return one numeric value per ",
          "curvature peak."
        ),
        call. = FALSE
      )
    }
    
    valid_strength <- is.finite(kappa_peaks_strength)
    
    if (!any(valid_strength)) {
      x_shoulder_candidates <- numeric(0)
    } else {
      peaks_kappa <- peaks_kappa[valid_strength, , drop = FALSE]
      kappa_peaks_strength <- kappa_peaks_strength[valid_strength]
      
      cutoff <- min(
        cutoff_kappa_peak_strength,
        stats::quantile(
          kappa_peaks_strength,
          probs = 1 - peak_retention_fraction,
          names = FALSE
        )
      )
      
      x_shoulder_candidates <- peaks_kappa[, 1][
        kappa_peaks_strength > cutoff
      ]
    }
  }
  
  # -------------------------------------------------------------------------
  # 8. Detect ordinary KDE peaks
  # -------------------------------------------------------------------------
  extrema <- cpp_find_extrema_from_grid(
    x = x_obs,
    d1_grid = d1_grid,
    x_obs = x_obs,
    y_obs = y_obs,
    bw = bw,
    tol = tol,
    max_iter = max_iter
  )
  
  if (!is.list(extrema) || is.null(extrema$local_max)) {
    stop(
      "'find_extrema()' must return a list containing 'local_max'.",
      call. = FALSE
    )
  }
  
  x_localmax <- extrema$local_max
  
  if (!is.numeric(x_localmax)) {
    stop(
      "'extrema$local_max' must be a numeric vector.",
      call. = FALSE
    )
  }
  
  # find_extrema() returns NA_real_ when no ordinary local maximum exists.
  # Convert that sentinel value into an empty vector so that a flat spectrum
  # or a spectrum segment without ordinary KDE peaks can be handled safely.
  if (length(x_localmax) == 1L && is.na(x_localmax)) {
    x_localmax <- numeric(0)
  }
  
  if (anyNA(x_localmax) || any(!is.finite(x_localmax))) {
    stop(
      "'extrema$local_max' must contain only finite values or a single NA sentinel.",
      call. = FALSE
    )
  }
  
  # -------------------------------------------------------------------------
  # 9. Remove shoulder candidates located near ordinary KDE peaks
  # -------------------------------------------------------------------------
  overlaps_with_localmax <- vapply(
    x_shoulder_candidates,
    function(xi) {
      any(abs(xi - x_localmax) <= merge_tol)
    },
    logical(1)
  )
  
  x_shoulder <- x_shoulder_candidates[!overlaps_with_localmax]
  
  # -------------------------------------------------------------------------
  # 10. Combine ordinary peaks and shoulder peaks
  # -------------------------------------------------------------------------
  kde_at_peaks <- function(x_peaks) {
    if (length(x_peaks) == 0L) {
      return(numeric(0))
    }
    
    cpp_gaussKDE(
      x = as.double(x_peaks),
      x_obs = x_obs,
      y_obs = y_obs,
      bw = bw
    )
  }
  
  peaks <- data.frame(
    x = x_localmax,
    y = kde_at_peaks(x_localmax),
    type = rep("peak", length(x_localmax)),
    stringsAsFactors = FALSE
  )
  
  shoulder <- data.frame(
    x = x_shoulder,
    y = kde_at_peaks(x_shoulder),
    type = rep("shoulder", length(x_shoulder)),
    stringsAsFactors = FALSE
  )
  
  peaks <- rbind(peaks, shoulder)
  
  if (nrow(peaks) == 0L) {
    return(create_empty_result())
  }
  
  peaks <- peaks[
    order(peaks$x),
    ,
    drop = FALSE
  ]
  
  rownames(peaks) <- NULL
  
  # -------------------------------------------------------------------------
  # 11. Remove weaker nearby candidates
  # -------------------------------------------------------------------------
  n_peaks <- nrow(peaks)
  x_peaks <- peaks$x
  y_peaks <- peaks$y
  
  type_priority <- ifelse(
    peaks$type == "peak",
    1L,
    2L
  )
  
  keep <- vapply(
    seq_len(n_peaks),
    function(i) {
      nearby_indices <- which(
        x_peaks >= x_peaks[i] - merge_tol &
          x_peaks <= x_peaks[i] + merge_tol
      )
      
      nearby_y <- y_peaks[nearby_indices]
      max_y <- max(nearby_y)
      
      strongest_indices <- nearby_indices[
        nearby_y == max_y
      ]
      
      if (length(strongest_indices) == 1L) {
        return(i == strongest_indices)
      }
      
      # In the rare case of an exact intensity tie, prioritize ordinary
      # KDE peaks over shoulder peaks. If the tie remains, retain the
      # leftmost candidate to guarantee deterministic output.
      preferred_priority <- min(type_priority[strongest_indices])
      
      preferred_indices <- strongest_indices[
        type_priority[strongest_indices] == preferred_priority
      ]
      
      selected_index <- preferred_indices[
        which.min(x_peaks[preferred_indices])
      ]
      
      return(i == selected_index)
    },
    logical(1)
  )
  
  peaks <- peaks[
    keep,
    ,
    drop = FALSE
  ]
  
  peaks <- peaks[
    order(peaks$x),
    ,
    drop = FALSE
  ]
  
  rownames(peaks) <- NULL
  
  # -------------------------------------------------------------------------
  # 12. Restore the original input-column names
  # -------------------------------------------------------------------------
  colnames(peaks)[1:2] <- output_colnames
  
  return(peaks)
}

#' @rdname find_peaks
#'
#' @param spectra A two-column `data.frame` or numeric `matrix` representing a
#'   single spectrum, or a `list` of such objects. The first column must contain
#'   strictly increasing x-axis values and the second column must contain
#'   intensity values.
#' @param n_cores A positive integer specifying the number of worker processes
#'   used when `spectra` is a list. The default is `1L`, which preserves
#'   sequential behaviour.
#'
#' @details
#' `find_peaks_spectra()` applies [find_peaks()] to a single spectrum or to each
#' spectrum in a list. When `bw` is `NULL`, the per-spectrum default bandwidth
#' of [find_peaks()] (the median interval between adjacent x-axis values) is
#' used for every spectrum. Supplying a non-`NULL` `bw` applies the same
#' bandwidth to all spectra.
#'
#' @return
#' `find_peaks()` returns a single peak table. `find_peaks_spectra()` returns a
#' single peak table when `spectra` is a single spectrum, or a named list of
#' peak tables with the same list structure and names as `spectra`.
#'
#' @seealso [find_peaks_spectra_fast()], [preprocess_maldi_spectra()],
#'   [filter_peaks_spectra()]
#'
#' @examples
#' \dontrun{
#' spectra <- list(sample_1 = spectrum_1, sample_2 = spectrum_2)
#'
#' peaks_list <- find_peaks_spectra(
#'   spectra = spectra,
#'   hws_peaks = 10,
#'   weight_type = "raw",
#'   n_cores = 2
#' )
#' }
#'
#' @export
find_peaks_spectra <- function(spectra,
                               bw = NULL,
                               hws_peaks = 10,
                               merge_tol = NULL,
                               tol = 1e-5,
                               max_iter = 100L,
                               weight_type = c("raw", "sqrt", "log10", "none"),
                               hws_grid_kappa_smooth = 3:20,
                               cutoff_kappa_peak_strength = 0.5,
                               peak_retention_fraction = 0.25,
                               n_cores = 1L) {
  
  weight_type <- match.arg(weight_type)
  
  .assert_integer_scalar(
    x = n_cores,
    arg_name = "n_cores",
    lower_bound = 1L
  )
  
  find_args <- list(
    hws_peaks = hws_peaks,
    merge_tol = merge_tol,
    tol = tol,
    max_iter = max_iter,
    weight_type = weight_type,
    hws_grid_kappa_smooth = hws_grid_kappa_smooth,
    cutoff_kappa_peak_strength = cutoff_kappa_peak_strength,
    peak_retention_fraction = peak_retention_fraction
  )
  
  if (!is.null(bw)) {
    find_args$bw <- bw
  }
  
  is_single_spectrum <- is.data.frame(spectra) || is.matrix(spectra)
  
  if (is_single_spectrum) {
    return(
      do.call(
        find_peaks,
        c(list(data = spectra), find_args)
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
  
  peaks_list <- .parallel_lapply(
    X = seq_along(spectra),
    FUN = function(i) {
      do.call(
        find_peaks,
        c(list(data = spectra[[i]]), find_args)
      )
    },
    n_cores = n_cores
  )
  
  names(peaks_list) <- names(spectra)
  
  peaks_list
}