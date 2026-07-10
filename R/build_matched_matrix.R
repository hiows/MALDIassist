#' @name build_matched_matrix
#' @title Build a Sample-by-Marker Peak Intensity Matrix
#'
#' @description
#' Matches detected peaks from multiple samples to a shared set of reference
#' m/z values and assembles a sample-by-marker feature matrix.
#'
#' @param peaks_list A named list of peak tables, typically from [find_peaks()].
#' @param reference_mz A numeric vector of reference m/z values, for example
#'   from [find_frequent_mz()].
#' @param reference_names Optional character vector of marker names.
#' @param hws_match Half-window size passed to [match_peaks()]. The default is
#'   `10`.
#' @param peak_selection_mode Peak selection rule passed to [match_peaks()].
#'   The default is `"nearest_mz"`.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{detected_matrix}{Integer matrix with samples as rows and markers
#'       as columns, holding `1` where a peak was matched and `0` otherwise.
#'       Column names are `reference_names` when supplied, otherwise
#'       `paste0("mz_", round(reference_mz, 3))`. Columns with no matched peak
#'       in any sample are removed.}
#'     \item{delta_mz_matrix}{Numeric matrix with the same dimensions as
#'       `detected_matrix`, holding the signed m/z difference
#'       (`detected_mz - reference_mz`). Unmatched cells remain `NA`.}
#'     \item{reference_mz}{Reference m/z values retained in the matrices.}
#'     \item{reference_names}{Marker names retained in the matrices, if
#'       supplied.}
#'     \item{sample_names}{Sample names from `peaks_list`.}
#'     \item{matches}{Named list of per-sample [match_peaks()] results.}
#'   }
#'
#' @seealso [match_peaks()], [find_frequent_mz()], [heatmap_matched_matrix()]
#'
#' @examples
#' peaks_list <- list(
#'   sample_1 = data.frame(mz = c(1000, 1500, 2000), intensity = c(10, 30, 20)),
#'   sample_2 = data.frame(mz = c(1002, 2001), intensity = c(12, 25))
#' )
#'
#' result <- build_matched_matrix(
#'   peaks_list = peaks_list,
#'   reference_mz = c(1000, 1500, 2000),
#'   reference_names = c("marker_a", "marker_b", "marker_c"),
#'   hws_match = 10
#' )
#'
#' result$detected_matrix
#'
#' @export
build_matched_matrix <- function(peaks_list,
                                 reference_mz,
                                 reference_names = NULL,
                                 hws_match = 10,
                                 peak_selection_mode = c("nearest_mz", "maximum_intensity")) {
  if (!is.list(peaks_list) || length(peaks_list) < 1L) {
    stop(
      "'peaks_list' must be a non-empty list.",
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

  peak_selection_mode <- match.arg(peak_selection_mode)

  if (!is.null(reference_names)) {
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

  marker_colnames <- if (is.null(reference_names)) {
    paste0("mz_", round(reference_mz, 3))
  } else {
    reference_names
  }

  if (is.null(names(peaks_list))) {
    names(peaks_list) <- as.character(seq_along(peaks_list))
  }

  sample_names <- names(peaks_list)
  n_samples <- length(sample_names)
  n_markers <- length(reference_mz)

  matches <- lapply(
    peaks_list,
    match_peaks,
    reference_mz = reference_mz,
    reference_names = reference_names,
    hws_match = hws_match,
    peak_selection_mode = peak_selection_mode
  )

  detected_mat <- matrix(
    0L,
    nrow = n_samples,
    ncol = n_markers,
    dimnames = list(sample_names, marker_colnames)
  )
  delta_mat <- matrix(
    NA_real_,
    nrow = n_samples,
    ncol = n_markers,
    dimnames = list(sample_names, marker_colnames)
  )

  for (i in seq_len(n_samples)) {
    m <- matches[[i]]
    idx <- match(reference_mz, m$reference_mz)
    if (anyNA(idx)) {
      stop(
        "Per-sample match results do not cover all 'reference_mz' values.",
        call. = FALSE
      )
    }
    detected_mat[i, ] <- as.integer(!is.na(m$detected_mz[idx]))
    delta_mat[i, ] <- m$delta_mz[idx]
  }

  keep <- colSums(detected_mat) > 0L
  if (any(!keep)) {
    detected_mat <- detected_mat[, keep, drop = FALSE]
    delta_mat <- delta_mat[, keep, drop = FALSE]
    reference_mz <- reference_mz[keep]
    if (!is.null(reference_names)) {
      reference_names <- reference_names[keep]
    }
    matches <- lapply(
      matches,
      function(m) {
        m[match(reference_mz, m$reference_mz), , drop = FALSE]
      }
    )
  }

  list(
    detected_matrix = detected_mat,
    delta_mz_matrix = delta_mat,
    reference_mz = reference_mz,
    reference_names = reference_names,
    sample_names = sample_names,
    matches = matches
  )
}
