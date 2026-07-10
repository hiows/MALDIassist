#' @name align_spectra
#' @title Alignment of MALDI-TOF Spectra to Internal Standards
#'
#' @description
#' Aligns a cohort of MALDI-TOF spectra along the m/z axis using internal
#' standard m/z values selected from frequently observed, high-intensity peaks
#' via [find_frequent_mz()]. Two alignment strategies are provided:
#'
#'   - `"linear"`: a two-point linear alignment. A low-mass and a high-mass
#'     anchor are selected, matched in each spectrum with [match_peaks()], and a
#'     linear map (`aligned_mz = slope * mz + intercept`) that sends the detected
#'     anchors onto the reference standard m/z is applied.
#'   - `"lowess"`: a multi-point non-linear alignment. All frequent standard m/z
#'     values are matched in each spectrum, a locally weighted regression
#'     (`lowess`) of the m/z shift against the detected m/z is fitted per
#'     spectrum, and the interpolated shift is added to the m/z axis.
#'
#' In both cases the alignment is applied to the first column (m/z) of the
#' spectrum and its peak table.
#'
#' @param spectra A list of spectra. Each element must be a matrix or data frame
#'   with m/z in the first column and intensity in the second column. When named,
#'   elements are matched to `peaks_list` by name; otherwise they are matched by
#'   position.
#'
#' @param peaks_list A list of peak tables aligned with `spectra`. Each element
#'   must be a matrix or data frame with m/z in the first column and intensity
#'   in the second column.
#'
#' @param bin_width A positive finite numeric scalar passed to
#'   [find_frequent_mz()] when locating candidate standard m/z values. The
#'   default is `20`.
#'
#' @param alignment_mode Character string specifying the alignment strategy,
#'   either `"linear"` (default) or `"lowess"`. See Details.
#'
#' @param lowess_span A finite numeric scalar in `(0, 1]` giving the smoother
#'   span (`f`) passed to [stats::lowess()] when `alignment_mode = "lowess"`.
#'   Larger values give smoother fits. The default is `2/3`.
#'
#' @param freq_ratio_cutoff A finite numeric scalar in `[0, 1]`. Only frequent
#'   m/z values with `freq_ratio > freq_ratio_cutoff` are considered as standard
#'   candidates. The default is `0.9`.
#'
#' @param hws_alignment A positive finite numeric scalar specifying the
#'   half-window size, in m/z units, passed to [match_peaks()] when matching
#'   each spectrum's peaks to the standard m/z values. The default is `50`.
#'
#' @details
#' Standard candidates are obtained once for the whole cohort from
#' [find_frequent_mz()] and filtered by `freq_ratio_cutoff`. At least two
#' candidates must remain for either mode.
#'
#' \strong{Linear mode}
#'
#' 1. Two standard m/z values are chosen. The retained m/z range is split into a
#'    low-mass and a high-mass half, and within each half the m/z with the
#'    largest `median_intensity` is used as an anchor. These are named `std_1`
#'    (low mass) and `std_2` (high mass).
#' 2. For each spectrum, the two anchors are matched to detected peaks with
#'    [match_peaks()] using `peak_selection_mode = "maximum_intensity"`, so that
#'    when several peaks fall within the matching window the most intense one is
#'    used.
#' 3. A two-point linear map is fitted from the two detected anchors to the two
#'    reference standard m/z values and applied to the first column (m/z) of the
#'    spectrum and its peak table.
#'
#' A spectrum is aligned only when both anchors are matched and the two detected
#' anchor m/z values differ (so that the slope is finite). Otherwise the
#' spectrum and its peaks are returned unchanged and flagged with
#' `is_aligned = FALSE`.
#'
#' \strong{Lowess mode}
#'
#' 1. All retained frequent m/z values are used as standards (`std_1`, ...,
#'    `std_N`).
#' 2. For each spectrum, standards are matched to detected peaks with
#'    [match_peaks()] using `peak_selection_mode = "nearest_mz"`. For every
#'    matched standard the shift `standard_mz - detected_mz` is computed.
#' 3. A `lowess` smoother of the shift against the detected m/z is fitted per
#'    spectrum, and the shift interpolated with [stats::approx()]
#'    (`rule = 2`) is added, independently, to the m/z axis of the spectrum and
#'    of the peak table.
#'
#' A spectrum is aligned only when at least two distinct matched standards are
#' available for the fit. Otherwise the spectrum and its peaks are returned
#' unchanged and flagged with `is_aligned = FALSE`.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{alignment_results}{A named list with one element per sample, each a
#'       list with components `spectrum` and `peaks` containing the (possibly)
#'       m/z-aligned spectrum and peak table.}
#'     \item{standard_mz}{A named numeric vector with the reference standard m/z
#'       values (`std_1`, ...). Length two in `"linear"` mode and one per
#'       retained frequent m/z in `"lowess"` mode.}
#'     \item{matched_mz}{A data frame with one row per sample and columns
#'       `spectrum_name`, one `std_*` column per standard, and `is_aligned`. The
#'       `std_*` columns hold the detected m/z for that sample, or `NA` when the
#'       standard was not matched. `is_aligned` indicates whether the alignment
#'       was applied.}
#'     \item{alignment_mode}{The alignment strategy used.}
#'   }
#'
#' @seealso [find_frequent_mz()], [match_peaks()], [build_matched_matrix()]
#'
#' @examples
#' \dontrun{
#' pp_spectra <- preprocess_maldi_spectra(spectra, n_cores = 2)
#' peaks_list <- find_peaks_spectra(pp_spectra, n_cores = 2)
#'
#' aligned <- align_spectra(
#'   spectra = pp_spectra,
#'   peaks_list = peaks_list,
#'   bin_width = 20,
#'   alignment_mode = "linear",
#'   freq_ratio_cutoff = 0.9,
#'   hws_alignment = 50
#' )
#'
#' aligned$standard_mz
#' aligned$matched_mz
#'
#' aligned_lowess <- align_spectra(
#'   spectra = pp_spectra,
#'   peaks_list = peaks_list,
#'   alignment_mode = "lowess",
#'   lowess_span = 2/3
#' )
#' }
#'
#' @importFrom stats lowess approx setNames
#' @export
align_spectra <- function(spectra,
                          peaks_list,
                          bin_width = 20,
                          alignment_mode = c("linear", "lowess"),
                          lowess_span = 2/3,
                          freq_ratio_cutoff = 0.9,
                          hws_alignment = 50) {

  # ----------------------------------------
  # Input validation
  # ----------------------------------------
  if (missing(spectra) || !is.list(spectra) || is.data.frame(spectra)) {
    stop(
      "'spectra' must be a list of spectra (matrix or data frame each).",
      call. = FALSE
    )
  }

  if (missing(peaks_list) || !is.list(peaks_list) || is.data.frame(peaks_list)) {
    stop(
      "'peaks_list' must be a list of peak tables (matrix or data frame each).",
      call. = FALSE
    )
  }

  if (length(spectra) == 0L || length(peaks_list) == 0L) {
    stop(
      "'spectra' and 'peaks_list' must be non-empty lists.",
      call. = FALSE
    )
  }

  if (!is.numeric(bin_width) ||
      length(bin_width) != 1L ||
      is.na(bin_width) ||
      !is.finite(bin_width) ||
      bin_width <= 0) {
    stop(
      "'bin_width' must be a positive finite numeric scalar.",
      call. = FALSE
    )
  }

  alignment_mode <- match.arg(alignment_mode)

  if (!is.numeric(lowess_span) ||
      length(lowess_span) != 1L ||
      is.na(lowess_span) ||
      !is.finite(lowess_span) ||
      lowess_span <= 0 ||
      lowess_span > 1) {
    stop(
      "'lowess_span' must be a finite numeric scalar in (0, 1].",
      call. = FALSE
    )
  }

  if (!is.numeric(freq_ratio_cutoff) ||
      length(freq_ratio_cutoff) != 1L ||
      is.na(freq_ratio_cutoff) ||
      !is.finite(freq_ratio_cutoff) ||
      freq_ratio_cutoff < 0 ||
      freq_ratio_cutoff > 1) {
    stop(
      "'freq_ratio_cutoff' must be a finite numeric scalar in [0, 1].",
      call. = FALSE
    )
  }

  if (!is.numeric(hws_alignment) ||
      length(hws_alignment) != 1L ||
      is.na(hws_alignment) ||
      !is.finite(hws_alignment) ||
      hws_alignment <= 0) {
    stop(
      "'hws_alignment' must be a positive finite numeric scalar.",
      call. = FALSE
    )
  }

  # ----------------------------------------
  # Match spectra and peak tables by name (or position)
  # ----------------------------------------
  if (is.null(names(spectra))) {
    names(spectra) <- as.character(seq_along(spectra))
  }

  if (is.null(names(peaks_list))) {
    names(peaks_list) <- as.character(seq_along(peaks_list))
  }

  common <- intersect(names(spectra), names(peaks_list))

  if (length(common) < 1L) {
    stop(
      "'spectra' and 'peaks_list' must share at least one sample name.",
      call. = FALSE
    )
  }

  spectra <- spectra[common]
  peaks_list <- peaks_list[common]

  for (i in seq_along(peaks_list)) {
    peaks_i <- peaks_list[[i]]
    if ((!is.matrix(peaks_i) && !is.data.frame(peaks_i)) || ncol(peaks_i) < 2L) {
      stop(
        sprintf(
          paste0(
            "Element '%s' of 'peaks_list' must be a matrix or data frame ",
            "with at least two columns (m/z and intensity)."
          ),
          common[i]
        ),
        call. = FALSE
      )
    }
  }

  for (i in seq_along(spectra)) {
    spectrum_i <- spectra[[i]]
    if ((!is.matrix(spectrum_i) && !is.data.frame(spectrum_i)) ||
        ncol(spectrum_i) < 2L) {
      stop(
        sprintf(
          paste0(
            "Element '%s' of 'spectra' must be a matrix or data frame ",
            "with at least two columns (m/z and intensity)."
          ),
          common[i]
        ),
        call. = FALSE
      )
    }
  }

  # ----------------------------------------
  # Select frequent standard m/z candidates (shared by both modes)
  # ----------------------------------------
  freq_df <- find_frequent_mz(
    peaks_list = peaks_list,
    bin_width = bin_width,
    exclude_mz = NULL
  )

  if (nrow(freq_df) == 0L) {
    stop(
      "No frequent m/z values found; cannot select standard m/z values.",
      call. = FALSE
    )
  }

  freq_df <- freq_df[freq_df$freq_ratio > freq_ratio_cutoff, , drop = FALSE]

  if (nrow(freq_df) < 2L) {
    stop(
      paste0(
        "Fewer than two frequent m/z values pass 'freq_ratio_cutoff'; ",
        "cannot define standard anchors. Try lowering 'freq_ratio_cutoff'."
      ),
      call. = FALSE
    )
  }

  # ----------------------------------------
  # Helper: reorder a match_peaks() result to the standard order
  # ----------------------------------------
  detected_in_order <- function(matching, std_names) {
    ord <- match(std_names, matching$reference_name)
    stats::setNames(matching$detected_mz[ord], std_names)
  }

  # ----------------------------------------
  # Helper: assemble the matched_mz table (unmatched standards shown as NA)
  # ----------------------------------------
  build_matched_mz <- function(detected_by_sample, std_names, is_aligned) {
    detected_matrix <- do.call(rbind, detected_by_sample)
    out <- data.frame(
      spectrum_name = common,
      stringsAsFactors = FALSE
    )
    for (j in seq_along(std_names)) {
      out[[std_names[j]]] <- detected_matrix[, j]
    }
    out$is_aligned <- is_aligned
    rownames(out) <- NULL
    out
  }

  # ----------------------------------------
  # Alignment
  # ----------------------------------------
  results <- switch(
    alignment_mode,

    # ======================================
    # Linear (two-point) alignment
    # ======================================
    linear = {
      mz <- freq_df$mz
      intensity <- freq_df$median_intensity

      min_mz <- min(mz)
      max_mz <- max(mz)
      sel_parts <- seq(min_mz, max_mz, length.out = 3L)

      std_mz <- vapply(seq_len(2L), function(i) {
        idx <- which(mz >= sel_parts[i] & mz <= sel_parts[i + 1L])
        if (length(idx) == 0L) {
          return(NA_real_)
        }
        mz_i <- mz[idx]
        intensity_i <- intensity[idx]
        mz_i[which.max(intensity_i)]
      }, numeric(1))

      std_names <- paste0("std_", seq_len(2L))
      names(std_mz) <- std_names

      if (anyNA(std_mz) || std_mz[1L] == std_mz[2L]) {
        stop(
          paste0(
            "Could not select two distinct standard anchors from the frequent ",
            "m/z values; cannot perform linear alignment."
          ),
          call. = FALSE
        )
      }

      std_matching <- lapply(
        peaks_list,
        match_peaks,
        reference_mz = std_mz,
        reference_names = std_names,
        hws_match = hws_alignment,
        peak_selection_mode = "maximum_intensity"
      )

      detected_by_sample <- lapply(std_matching, detected_in_order, std_names = std_names)

      linear_alignment_fun <- function(std_mz, matched_mz) {
        slope <- (std_mz[2L] - std_mz[1L]) / (matched_mz[2L] - matched_mz[1L])
        intercept <- std_mz[2L] - (slope * matched_mz[2L])
        function(x) slope * x + intercept
      }

      is_aligned <- logical(length(common))

      alignment_results <- lapply(seq_along(common), function(i) {
        matched_mz <- detected_by_sample[[i]]
        spectrum <- spectra[[i]]
        peaks <- peaks_list[[i]]

        can_align <- !anyNA(matched_mz) && (matched_mz[2L] != matched_mz[1L])

        if (can_align) {
          alignment_fun <- linear_alignment_fun(std_mz = std_mz, matched_mz = matched_mz)
          spectrum[, 1L] <- alignment_fun(spectrum[, 1L])
          peaks[, 1L] <- alignment_fun(peaks[, 1L])
          is_aligned[i] <<- TRUE
        } else {
          is_aligned[i] <<- FALSE
        }

        list(spectrum = spectrum, peaks = peaks)
      })
      names(alignment_results) <- common

      matched_mz <- build_matched_mz(detected_by_sample, std_names, is_aligned)

      list(
        alignment_results = alignment_results,
        standard_mz = std_mz,
        matched_mz = matched_mz,
        alignment_mode = alignment_mode
      )
    },

    # ======================================
    # Lowess (multi-point) alignment
    # ======================================
    lowess = {
      std_mz <- freq_df$mz
      std_names <- paste0("std_", seq_along(std_mz))
      names(std_mz) <- std_names

      std_matching <- lapply(
        peaks_list,
        match_peaks,
        reference_mz = std_mz,
        reference_names = std_names,
        hws_match = hws_alignment,
        peak_selection_mode = "nearest_mz"
      )

      detected_by_sample <- lapply(std_matching, detected_in_order, std_names = std_names)

      is_aligned <- logical(length(common))

      alignment_results <- lapply(seq_along(common), function(i) {
        matched_mz <- detected_by_sample[[i]]
        spectrum <- spectra[[i]]
        peaks <- peaks_list[[i]]

        # One (detected_mz, shift) point per matched standard.
        df <- data.frame(
          detected_mz = matched_mz,
          shift = std_mz - matched_mz,
          stringsAsFactors = FALSE
        )
        df <- df[!is.na(df$detected_mz), , drop = FALSE]
        df <- df[!duplicated(df$detected_mz), , drop = FALSE]

        can_align <- nrow(df) >= 2L

        if (can_align) {
          fit <- stats::lowess(x = df$detected_mz, y = df$shift, f = lowess_span)

          spectrum_shift <- stats::approx(
            x = fit$x,
            y = fit$y,
            xout = spectrum[, 1L],
            rule = 2
          )$y
          spectrum[, 1L] <- spectrum[, 1L] + spectrum_shift

          peaks_shift <- stats::approx(
            x = fit$x,
            y = fit$y,
            xout = peaks[, 1L],
            rule = 2
          )$y
          peaks[, 1L] <- peaks[, 1L] + peaks_shift

          is_aligned[i] <<- TRUE
        } else {
          is_aligned[i] <<- FALSE
        }

        list(spectrum = spectrum, peaks = peaks)
      })
      names(alignment_results) <- common

      matched_mz <- build_matched_mz(detected_by_sample, std_names, is_aligned)

      list(
        alignment_results = alignment_results,
        standard_mz = std_mz,
        matched_mz = matched_mz,
        alignment_mode = alignment_mode
      )
    }
  )

  return(results)
}
