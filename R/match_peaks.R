#' @name match_peaks
#' @title Match Detected Peaks to Reference m/z Values
#'
#' @description
#' For each reference m/z value, searches detected peaks within a symmetric
#' half-window and selects one peak per reference. Returns one row per reference
#' m/z with match status and detected peak properties.
#'
#' @param peaks A numeric matrix or data frame containing detected peaks.
#'   The first column must contain peak positions, typically m/z values, and
#'   the second column must contain the corresponding peak intensities.
#'
#' @param reference_mz A numeric vector of reference m/z values to match
#'   against.
#'
#' @param reference_names A character vector of names corresponding to
#'   `reference_mz`. If `NULL`, `NA` is used for all reference names. The
#'   default is `NULL`.
#'
#' @param hws_match A positive finite numeric scalar specifying the half-window
#'   size used for matching. For each reference m/z, only peaks within
#'   `[reference_mz - hws_match, reference_mz + hws_match]` are considered.
#'   The default is `10`.
#'
#' @param peak_selection_mode Character string specifying how a peak is chosen
#'   when multiple peaks fall within the matching window. Either
#'   `"nearest_mz"` (default) or `"maximum_intensity"`.
#'
#' @details
#' When multiple peaks fall within the matching window, selection depends on
#' `peak_selection_mode`:
#'
#'   - `"nearest_mz"`: the peak whose m/z is closest to `reference_mz`.
#'   - `"maximum_intensity"`: the peak with the highest intensity in the
#'     second column of `peaks`.
#'
#' If `peaks` is not sorted by m/z, it is sorted internally before matching.
#'
#' @return A data frame with one row per reference m/z, sorted in ascending
#'   order of `reference_mz`. Columns are:
#'
#'   - `reference_name`: reference label from `reference_names`, or `NA`.
#'   - `reference_mz`: reference m/z value.
#'   - `is_matched`: logical flag indicating whether a peak was found in the
#'     matching window.
#'   - `detected_mz`: m/z of the matched peak, or `NA` when unmatched.
#'   - `detected_intensity`: intensity of the matched peak, or `NA` when
#'     unmatched.
#'   - `delta_mz`: signed m/z difference, `detected_mz - reference_mz`, or
#'     `NA` when unmatched.
#'
#' @examples
#' peaks <- data.frame(
#'   mz = c(1000, 2000.5, 3000),
#'   intensity = c(10, 50, 20)
#' )
#'
#' match_peaks(
#'   peaks = peaks,
#'   reference_mz = c(999, 2005, 4000),
#'   reference_names = c("A", "B", "C"),
#'   hws_match = 10
#' )
#'
#' @export
match_peaks <- function(peaks,
                        reference_mz,
                        reference_names = NULL,
                        hws_match = 10,
                        peak_selection_mode = c("nearest_mz", "maximum_intensity")) {
  if (!is.matrix(peaks) && !is.data.frame(peaks)) {
    stop(
      "'peaks' must be a numeric matrix or data frame.",
      call. = FALSE
    )
  }

  if (ncol(peaks) < 2L) {
    stop(
      "'peaks' must have at least two columns: peak positions and intensities.",
      call. = FALSE
    )
  }

  if (!is.numeric(reference_mz) || length(reference_mz) < 1L) {
    stop(
      "'reference_mz' must be a numeric vector with at least one value.",
      call. = FALSE
    )
  }

  if (anyNA(reference_mz) || any(!is.finite(reference_mz))) {
    stop(
      "'reference_mz' must contain only finite values.",
      call. = FALSE
    )
  }

  if (is.null(reference_names)) {
    reference_names <- rep(NA_character_, length(reference_mz))
  } else {
    if (!is.character(reference_names)) {
      stop(
        "'reference_names' must be a character vector or NULL.",
        call. = FALSE
      )
    }

    if (length(reference_names) != length(reference_mz)) {
      stop(
        "'reference_names' must have the same length as 'reference_mz'.",
        call. = FALSE
      )
    }
  }

  if (!is.numeric(hws_match) ||
      length(hws_match) != 1L ||
      is.na(hws_match) ||
      !is.finite(hws_match) ||
      hws_match <= 0) {
    stop(
      "'hws_match' must be a positive finite numeric scalar.",
      call. = FALSE
    )
  }

  peak_selection_mode <- match.arg(peak_selection_mode)

  n_ref <- length(reference_mz)

  if (nrow(peaks) == 0L) {
    out <- data.frame(
      reference_name = reference_names,
      reference_mz = reference_mz,
      is_matched = rep(FALSE, n_ref),
      detected_mz = rep(NA_real_, n_ref),
      detected_intensity = rep(NA_real_, n_ref),
      delta_mz = rep(NA_real_, n_ref),
      stringsAsFactors = FALSE
    )

    out <- out[order(out$reference_mz), , drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }

  mz <- peaks[, 1]
  intensity <- peaks[, 2]

  if (!is.numeric(mz) || !is.numeric(intensity)) {
    stop(
      "The first and second columns of 'peaks' must be numeric.",
      call. = FALSE
    )
  }

  if (anyNA(mz) || any(!is.finite(mz))) {
    stop(
      "Peak positions in the first column of 'peaks' must be finite values.",
      call. = FALSE
    )
  }

  if (anyNA(intensity) || any(!is.finite(intensity))) {
    stop(
      "Peak intensities in the second column of 'peaks' must be finite values.",
      call. = FALSE
    )
  }

  if (is.unsorted(mz)) {
    ord <- order(mz)
    mz <- mz[ord]
    intensity <- intensity[ord]
  }

  left <- findInterval(reference_mz - hws_match, mz, left.open = TRUE) + 1L
  right <- findInterval(reference_mz + hws_match, mz)

  is_matched <- left <= right
  detected_mz <- rep(NA_real_, n_ref)
  detected_intensity <- rep(NA_real_, n_ref)

  matched_idx <- which(is_matched)
  for (i in matched_idx) {
    idx <- left[i]:right[i]
    j <- switch(
      peak_selection_mode,
      nearest_mz = idx[which.min(abs(mz[idx] - reference_mz[i]))],
      maximum_intensity = idx[which.max(intensity[idx])]
    )
    detected_mz[i] <- mz[j]
    detected_intensity[i] <- intensity[j]
  }

  delta_mz <- detected_mz - reference_mz

  out <- data.frame(
    reference_name = reference_names,
    reference_mz = reference_mz,
    is_matched = is_matched,
    detected_mz = detected_mz,
    detected_intensity = detected_intensity,
    delta_mz = delta_mz,
    stringsAsFactors = FALSE
  )

  out <- out[order(out$reference_mz), , drop = FALSE]
  rownames(out) <- NULL

  return(out)
}
