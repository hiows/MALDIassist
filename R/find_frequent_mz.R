#' @name find_frequent_mz
#' @title Find Frequent m/z Values in a Pooled Peak Set
#'
#' @description
#' Identifies frequently observed m/z values across a cohort by binning a pooled
#' m/z vector, retaining high-frequency bins, and refining each bin location with
#' Gaussian kernel density estimation.
#'
#' @param pooled_mz A numeric vector of pooled peak m/z values collected across
#'   multiple samples.
#'
#' @param exclude_mz A numeric vector of m/z values to remove from `pooled_mz`
#'   before analysis, for example known contaminant or calibrant peaks. If
#'   `NULL`, no values are excluded. The default is `NULL`.
#'
#' @param bin_width A positive finite numeric scalar specifying the bin width
#'   used for the initial frequency count. The default is `20`.
#'
#' @param min_freq A numeric scalar specifying the minimum bin count required
#'   to retain a candidate m/z region. Values less than zero are treated as
#'   zero. The default is `30`.
#'
#' @details
#' The function performs the following steps:
#'
#' 1. Remove values listed in `exclude_mz` from `pooled_mz`.
#' 2. Bin the pooled m/z values with width `bin_width` and retain bins with
#'    count greater than `min_freq`.
#' 3. For each retained bin, refine the m/z location as the mode of a Gaussian
#'    kernel density estimate computed on m/z values within
#'    `[bin_midpoint - bin_width, bin_midpoint + bin_width]`.
#' 4. Re-count observations within `[refined_mz - hws, refined_mz + hws]`,
#'    where `hws = bin_width / 2`.
#'
#' The returned table is sorted in descending order of `count`.
#'
#' @return A data frame with columns:
#'
#'   - `mz`: density-refined frequent m/z value.
#'   - `count`: number of pooled observations within `[mz - hws, mz + hws]`.
#'
#'   If no frequent m/z values remain, an empty data frame with the same
#'   columns is returned.
#'
#' @examples
#' set.seed(1)
#' pooled_mz <- c(
#'   3000 + rnorm(40, sd = 1),
#'   5000 + rnorm(35, sd = 1)
#' )
#'
#' find_frequent_mz(
#'   pooled_mz = pooled_mz,
#'   exclude_mz = NULL,
#'   bin_width = 20,
#'   min_freq = 30
#' )
#'
#' @seealso [find_peaks()]
#'
#' @export
find_frequent_mz <- function(pooled_mz,
                             exclude_mz = NULL,
                             bin_width = 20,
                             min_freq = 30) {
  empty_out <- function() {
    out <- data.frame(
      mz = numeric(0),
      count = numeric(0),
      stringsAsFactors = FALSE
    )
    rownames(out) <- NULL
    return(out)
  }

  if (!is.numeric(pooled_mz)) {
    stop(
      "'pooled_mz' must be a numeric vector.",
      call. = FALSE
    )
  }

  pooled_mz <- as.numeric(pooled_mz)

  if (length(pooled_mz) == 0L) {
    return(empty_out())
  }

  if (anyNA(pooled_mz) || any(!is.finite(pooled_mz))) {
    stop(
      "'pooled_mz' must contain only finite values.",
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

  if (!is.numeric(min_freq) ||
      length(min_freq) != 1L ||
      is.na(min_freq) ||
      !is.finite(min_freq)) {
    stop(
      "'min_freq' must be a finite numeric scalar.",
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

    pooled_mz <- pooled_mz[!pooled_mz %in% exclude_mz]

    if (length(pooled_mz) == 0L) {
      return(empty_out())
    }
  }

  mz_pool <- pooled_mz
  hws <- bin_width / 2
  min_freq <- max(0, min_freq)

  bins <- seq(
    min(mz_pool) - bin_width,
    max(mz_pool) + bin_width,
    by = bin_width
  )
  n_bins <- length(bins)

  idx <- findInterval(mz_pool, bins)
  mz_counts <- tabulate(idx, nbins = n_bins - 1L)
  mz_mids <- (bins[-n_bins] + bins[-1L]) / 2

  keep_bins <- mz_counts > min_freq
  if (!any(keep_bins)) {
    return(empty_out())
  }

  df_counts <- data.frame(
    mz = mz_mids[keep_bins],
    count = mz_counts[keep_bins],
    stringsAsFactors = FALSE
  )

  df_counts_opt <- lapply(seq_len(nrow(df_counts)), function(i) {
    mz_i <- df_counts$mz[i]

    x <- mz_pool[
      mz_pool >= mz_i - bin_width &
        mz_pool <= mz_i + bin_width
    ]

    d <- stats::density(
      x,
      bw = hws,
      kernel = "gaussian",
      n = 4096L
    )
    xi <- d$x[which.max(d$y)]

    data.frame(
      mz = xi,
      count = sum(
        mz_pool >= xi - hws &
          mz_pool <= xi + hws
      ),
      stringsAsFactors = FALSE
    )
  })

  df_counts_opt <- do.call(rbind, df_counts_opt)
  rownames(df_counts_opt) <- NULL

  if (nrow(df_counts_opt) == 0L) {
    return(empty_out())
  }

  df <- df_counts_opt[order(df_counts_opt$count, decreasing = TRUE), , drop = FALSE]
  rownames(df) <- NULL

  return(df)
}
