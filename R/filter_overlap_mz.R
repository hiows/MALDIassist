#' @name filter_overlap_mz
#' @title Filter Overlapping Frequent m/z Features
#'
#' @description
#' Removes redundant frequent m/z values that fall within a local half-window
#' of each other. For every candidate, the function retains the row only when
#' its `freq_ratio` equals the maximum `freq_ratio` among all candidates in
#' `[mz - hws_selection, mz + hws_selection]`. This is a post-processing step
#' for tables returned by [find_frequent_mz()], which can otherwise report
#' nearby m/z values (for example, differing by only about 1 Da) as separate
#' features.
#'
#' @param freq_data A data frame or matrix containing at least the columns
#'   `mz` and `freq_ratio`, typically the output of [find_frequent_mz()].
#'   Additional columns are preserved in the returned object.
#'
#' @param hws_selection A non-negative finite numeric scalar specifying the
#'   half-window size, in m/z units, used to define the local neighborhood for
#'   each candidate. The default is `5`.
#'
#' @details
#' Ties are retained: if two or more candidates share the same maximum
#' `freq_ratio` inside a neighborhood, all of them are kept.
#'
#' The returned table keeps the original row order of retained candidates.
#' Row names are reset.
#'
#' @return An object of the same class as `freq_data`, containing only the
#'   non-overlapping frequent m/z rows. If no rows remain, an empty object with
#'   the same columns is returned.
#'
#' @seealso [find_frequent_mz()], [build_matched_matrix()]
#'
#' @examples
#' freq_data <- data.frame(
#'   mz = c(3000, 3001, 5000, 5010),
#'   median_intensity = c(2.0, 1.8, 2.5, 2.1),
#'   count = c(18, 10, 20, 12),
#'   freq_ratio = c(0.9, 0.5, 1.0, 0.6)
#' )
#'
#' filter_overlap_mz(freq_data, hws_selection = 5)
#'
#' @export
filter_overlap_mz <- function(freq_data, hws_selection = 5) {

  # ----------------------------------------
  # Input validation
  # ----------------------------------------
  if (missing(freq_data) || (!is.data.frame(freq_data) && !is.matrix(freq_data))) {
    stop(
      "'freq_data' must be a data frame or matrix.",
      call. = FALSE
    )
  }

  required_cols <- c("mz", "freq_ratio")
  missing_cols <- setdiff(required_cols, colnames(freq_data))
  if (length(missing_cols) > 0L) {
    stop(
      paste0(
        "'freq_data' must contain columns: ",
        paste(required_cols, collapse = ", "),
        ". Missing: ",
        paste(missing_cols, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  if (!is.numeric(hws_selection) ||
      length(hws_selection) != 1L ||
      is.na(hws_selection) ||
      !is.finite(hws_selection) ||
      hws_selection < 0) {
    stop(
      "'hws_selection' must be a non-negative finite numeric scalar.",
      call. = FALSE
    )
  }

  n <- nrow(freq_data)
  if (n == 0L) {
    rownames(freq_data) <- NULL
    return(freq_data)
  }

  mz <- as.numeric(freq_data[, "mz"])
  freq_ratio <- as.numeric(freq_data[, "freq_ratio"])

  if (anyNA(mz) || any(!is.finite(mz))) {
    stop(
      "'freq_data$mz' must contain only finite values.",
      call. = FALSE
    )
  }

  if (anyNA(freq_ratio) || any(!is.finite(freq_ratio))) {
    stop(
      "'freq_data$freq_ratio' must contain only finite values.",
      call. = FALSE
    )
  }

  # ----------------------------------------
  # Keep local freq_ratio maxima within +/- hws_selection
  # ----------------------------------------
  keep <- vapply(seq_len(n), function(i) {
    xi <- mz[i]
    yi <- freq_ratio[i]
    in_window <- mz >= (xi - hws_selection) & mz <= (xi + hws_selection)
    yi == max(freq_ratio[in_window])
  }, logical(1))

  out <- freq_data[keep, , drop = FALSE]
  rownames(out) <- NULL
  return(out)
}
