#' @name find_frequent_mz
#' @title Find Frequent m/z Values from a List of Peak Tables
#'
#' @description
#' Identifies frequently observed m/z values across a cohort of samples. Each
#' sample is represented by a peak table with m/z in the first column and
#' intensity in the second column. The pooled m/z values are scanned with
#' overlapping bins, and each bin location is refined to the strongest local
#' maximum of a Gaussian kernel density estimate (KDE) using [get_gauss_kde()]
#' and [find_extrema()].
#'
#' @param peaks_list A list of peak tables. Each element must be a matrix
#'   or data frame with at least two columns: m/z in the first column and
#'   intensity in the second column.
#'
#' @param bin_width A positive finite numeric scalar specifying the bin width
#'   used for the initial scan. The default is `20`.
#'
#' @param exclude_mz A numeric vector of m/z values to remove from the result,
#'   for example known contaminant or calibrant peaks. For each value, refined
#'   m/z values within `[exclude_mz - hws_exclude, exclude_mz + hws_exclude]`
#'   are removed. If `NULL`, no values are excluded. The default is `NULL`.
#'
#' @param hws_exclude A non-negative finite numeric scalar specifying the
#'   half-window size, in m/z units, used to exclude values around each
#'   `exclude_mz`. A value of `0` removes only exact matches. The default is
#'   `bin_width / 2`.
#'
#' @details
#' The function performs the following steps:
#'
#' 1. Pool the m/z and intensity values across all peak tables and drop any
#'    non-finite observations.
#' 2. Scan the pooled m/z range with three overlapping sets of bins of width
#'    `bin_width` to reduce edge effects.
#' 3. For each bin, build a histogram, refine the m/z location as the strongest
#'    local maximum of a Gaussian KDE, and re-count pooled observations within
#'    `[refined_mz - bin_width, refined_mz + bin_width]`.
#' 4. Drop duplicate refined m/z values produced by the overlapping bins.
#' 5. Remove refined m/z values within
#'    `[exclude_mz - hws_exclude, exclude_mz + hws_exclude]` of any value listed
#'    in `exclude_mz`.
#'
#' The returned table is sorted in ascending order of `mz`.
#'
#' @return A data frame with columns:
#'
#'   - `mz`: density-refined frequent m/z value.
#'   - `median_intensity`: geometric median intensity (via `log10`) of pooled
#'     observations within `[mz - bin_width, mz + bin_width]`.
#'   - `count`: number of pooled observations within
#'     `[mz - bin_width, mz + bin_width]`.
#'   - `freq_ratio`: `count / length(peaks_list)`.
#'
#'   If no frequent m/z values remain, an empty data frame with the same columns
#'   is returned.
#'
#' @examples
#' set.seed(1)
#' peaks_list <- lapply(seq_len(20), function(i) {
#'   mz <- c(3000 + rnorm(1, sd = 1), 5000 + rnorm(1, sd = 1))
#'   data.frame(mz = mz, intensity = runif(2, 100, 1000))
#' })
#'
#' find_frequent_mz(
#'   peaks_list = peaks_list,
#'   bin_width = 20,
#'   exclude_mz = NULL
#' )
#'
#' @seealso [find_peaks()], [get_gauss_kde()], [find_extrema()]
#'
#' @export
find_frequent_mz <- function(peaks_list,
                             bin_width = 20,
                             exclude_mz = NULL,
                             hws_exclude = bin_width / 2) {

  # ----------------------------------------
  # Empty-result helper
  # ----------------------------------------
  empty_out <- function() {
    out <- data.frame(
      mz = numeric(0),
      median_intensity = numeric(0),
      count = numeric(0),
      freq_ratio = numeric(0),
      stringsAsFactors = FALSE
    )
    rownames(out) <- NULL
    return(out)
  }

  # ----------------------------------------
  # Input validation
  # ----------------------------------------
  if (missing(peaks_list) || !is.list(peaks_list)) {
    stop(
      "'peaks_list' must be a list of peak tables (matrix or data frame).",
      call. = FALSE
    )
  }

  n_samples <- length(peaks_list)
  if (n_samples == 0L) {
    return(empty_out())
  }

  for (i in seq_len(n_samples)) {
    peaks_i <- peaks_list[[i]]

    if (!is.matrix(peaks_i) && !is.data.frame(peaks_i)) {
      stop(
        sprintf(
          "Element %d of 'peaks_list' must be a matrix or data frame.",
          i
        ),
        call. = FALSE
      )
    }

    if (ncol(peaks_i) < 2L) {
      stop(
        sprintf(
          paste0(
            "Element %d of 'peaks_list' must have at least two ",
            "columns: m/z and intensity."
          ),
          i
        ),
        call. = FALSE
      )
    }
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

  if (!is.numeric(hws_exclude) ||
      length(hws_exclude) != 1L ||
      is.na(hws_exclude) ||
      !is.finite(hws_exclude) ||
      hws_exclude < 0) {
    stop(
      "'hws_exclude' must be a non-negative finite numeric scalar.",
      call. = FALSE
    )
  }

  if (!is.null(exclude_mz)) {
    if (!is.numeric(exclude_mz) || length(exclude_mz) < 1L) {
      stop(
        "'exclude_mz' must be a numeric vector or NULL.",
        call. = FALSE
      )
    }

    if (anyNA(exclude_mz) || any(!is.finite(exclude_mz))) {
      stop(
        "'exclude_mz' must contain only finite values.",
        call. = FALSE
      )
    }
  }

  # ----------------------------------------
  # Pool m/z and intensity across samples
  # ----------------------------------------
  pooled_mz <- unlist(
    lapply(peaks_list, function(peaks_i) as.numeric(peaks_i[, 1L])),
    use.names = FALSE
  )

  pooled_intensity <- unlist(
    lapply(peaks_list, function(peaks_i) as.numeric(peaks_i[, 2L])),
    use.names = FALSE
  )

  keep <- is.finite(pooled_mz) & is.finite(pooled_intensity)
  pooled_mz <- pooled_mz[keep]
  pooled_intensity <- pooled_intensity[keep]

  if (length(pooled_mz) == 0L) {
    return(empty_out())
  }

  # Sort once so window ranges can be extracted by binary search
  # (findInterval) instead of scanning the full vector per bin. The results
  # are identical: hist() bins values regardless of order, while count and
  # median() are order-independent, and the pairing between m/z and intensity
  # is preserved by reordering both with the same permutation.
  ord <- order(pooled_mz)
  pooled_mz <- pooled_mz[ord]
  pooled_intensity <- pooled_intensity[ord]

  # ----------------------------------------
  # Build overlapping bin breaks
  # ----------------------------------------
  mz_min <- min(pooled_mz)
  mz_max <- max(pooled_mz)
  hws <- bin_width / 2

  brk_0     <- seq(mz_min - bin_width,       mz_max + bin_width, by = bin_width)
  brk_front <- seq(mz_min - bin_width - hws, mz_max + bin_width, by = bin_width)
  brk_back  <- seq(mz_min - bin_width + hws, mz_max + bin_width, by = bin_width)

  n_bins <- length(brk_0) - 1L
  if (n_bins < 1L) {
    return(empty_out())
  }

  # Precompute each window's index bounds in the sorted m/z vector.
  # `nlt_all` counts values strictly below the lower edge and `nle_all`
  # counts values at or below the upper edge, so the window contains the
  # sorted indices `(nlt + 1):nle`. Windows with NA edges (from unequal
  # `brk_*` lengths) yield NA here and are skipped below, matching the
  # original behavior where such bins produced an all-NA `x`.
  bin_idx <- seq_len(n_bins)
  win_lo_all <- pmin(brk_0[bin_idx],      brk_front[bin_idx],      brk_back[bin_idx])
  win_hi_all <- pmax(brk_0[bin_idx + 1L], brk_front[bin_idx + 1L], brk_back[bin_idx + 1L])
  nlt_all <- findInterval(win_lo_all, pooled_mz, left.open = TRUE)
  nle_all <- findInterval(win_hi_all, pooled_mz)

  # ----------------------------------------
  # Refine each bin to the strongest KDE local maximum
  # ----------------------------------------
  freq_list <- lapply(bin_idx, function(i) {

    nlt <- nlt_all[i]
    nle <- nle_all[i]
    if (is.na(nlt) || is.na(nle) || nle <= nlt) {
      return(NULL)
    }

    x <- pooled_mz[(nlt + 1L):nle]

    # A histogram needs a non-zero range; skip degenerate windows.
    if (length(unique(x)) < 2L) {
      return(NULL)
    }

    h <- tryCatch(
      graphics::hist(x, breaks = bin_width * 2, plot = FALSE),
      error = function(e) NULL
    )
    if (is.null(h)) {
      return(NULL)
    }

    df <- data.frame(mid = h$mids, dens = h$density, stringsAsFactors = FALSE)

    if (nrow(df) < 2L) {
      return(NULL)
    }

    m <- df$mid
    bw <- min(diff(df$mid))
    if (!is.finite(bw) || bw <= 0) {
      return(NULL)
    }

    # KDE, first and second derivatives; guard against numeric failure.
    fitted <- tryCatch({
      f  <- get_gauss_kde(data = df, bw = bw)
      d1 <- get_gauss_kde(data = df, bw = bw, d = 1)
      d2 <- get_gauss_kde(data = df, bw = bw, d = 2)
      list(
        f = f,
        ext = find_extrema(first_deriv = d1, second_deriv = d2, x = m)
      )
    }, error = function(e) NULL)

    if (is.null(fitted)) {
      return(NULL)
    }

    x_localmax <- fitted$ext$local_max
    if (all(is.na(x_localmax))) {
      return(NULL)
    }

    x_localmax <- x_localmax[is.finite(x_localmax)]
    if (length(x_localmax) == 0L) {
      return(NULL)
    }

    y_localmax <- fitted$f(x_localmax)
    idx <- which.max(y_localmax)
    if (length(idx) == 0L) {
      return(NULL)
    }

    xi <- x_localmax[idx]

    # Count observations within [xi - bin_width, xi + bin_width] via the same
    # sorted binary search (inclusive on both ends).
    c_lt <- findInterval(xi - bin_width, pooled_mz, left.open = TRUE)
    c_le <- findInterval(xi + bin_width, pooled_mz)
    count <- c_le - c_lt

    sub_intensity <- if (count > 0L) {
      pooled_intensity[(c_lt + 1L):c_le]
    } else {
      numeric(0)
    }
    sub_intensity <- sub_intensity[is.finite(sub_intensity) & sub_intensity > 0]
    median_intensity <- if (length(sub_intensity) == 0L) {
      NA_real_
    } else {
      10^stats::median(log10(sub_intensity))
    }

    data.frame(
      mz = xi,
      median_intensity = median_intensity,
      count = count,
      freq_ratio = count / n_samples,
      stringsAsFactors = FALSE
    )
  })

  freq_list <- freq_list[!vapply(freq_list, is.null, logical(1))]
  if (length(freq_list) == 0L) {
    return(empty_out())
  }

  freq_df <- do.call(rbind, freq_list)
  freq_df <- as.data.frame(freq_df, stringsAsFactors = FALSE)

  # ----------------------------------------
  # Remove duplicate m/z from overlapping bins
  # ----------------------------------------
  freq_df <- freq_df[!duplicated(freq_df$mz), , drop = FALSE]

  # ----------------------------------------
  # Exclude requested m/z regions
  # ----------------------------------------
  if (!is.null(exclude_mz)) {
    drop <- logical(nrow(freq_df))
    for (e in exclude_mz) {
      drop <- drop | abs(freq_df$mz - e) <= hws_exclude
    }
    freq_df <- freq_df[!drop, , drop = FALSE]
  }

  if (nrow(freq_df) == 0L) {
    return(empty_out())
  }

  freq_df <- freq_df[order(freq_df$mz), , drop = FALSE]
  rownames(freq_df) <- NULL

  return(freq_df)
}
