#' @name estimate_significance
#' @title Estimate Significant m/z Features Between Two Groups
#'
#' @description
#' Performs a per-feature two-group comparison on a sample-by-feature matrix,
#' such as the intensity or detection matrix produced by
#' [build_matched_matrix()]. For each feature (column), a two-sided test is run
#' between the two sample groups and the resulting p-values are optionally
#' adjusted for multiple comparisons.
#'
#' @param matched_matrix A numeric matrix or data frame with samples in rows and
#'   features (m/z markers) in columns. A data frame is coerced with
#'   `as.matrix()` and must contain only numeric columns.
#'
#' @param group A vector of group labels with one entry per row of
#'   `matched_matrix`. It must contain exactly two distinct (non-`NA`) levels.
#'   The first level in factor order is treated as the reference group.
#'
#' @param feat_names Optional character vector of feature names, one per column
#'   of `matched_matrix`. If `NULL`, the column names of `matched_matrix` are
#'   used when available, otherwise `paste0("feat_", seq_len(ncol))`.
#'
#' @param stat_method The statistical test to apply per feature, either
#'   `"t.test"` (Welch two-sample t-test) or `"wilcox"` (Wilcoxon rank-sum test
#'   with `exact = FALSE`). The default is `"t.test"`.
#'
#' @param adj_method The multiple-comparison adjustment passed to
#'   [stats::p.adjust()], one of `"none"`, `"BH"`, or `"bonferroni"`. The
#'   default is `"none"`.
#'
#' @details
#' The two group levels are determined from `as.factor(group)`. Rows are split
#' into the reference group (the first factor level) and the other group, and a
#' two-sided test is applied to each feature column:
#'
#' - For `"t.test"`, if both groups are constant for a feature (zero standard
#'   deviation), the p-value is set to `NaN`, matching the behavior of a
#'   degenerate comparison.
#' - Tests that fail (for example, on constant data) return `NaN` for that
#'   feature instead of raising an error, so a single problematic feature does
#'   not abort the whole scan.
#'
#' @return A data frame with one row per feature and columns:
#'
#'   - `feat_names`: feature name.
#'   - `pvalue`: raw two-sided p-value.
#'   - `adj_pvalue`: p-value after `adj_method` adjustment.
#'
#' @examples
#' set.seed(1)
#' n_per_group <- 10
#' matched_matrix <- rbind(
#'   matrix(rnorm(n_per_group * 3, mean = 0), nrow = n_per_group),
#'   matrix(rnorm(n_per_group * 3, mean = 1), nrow = n_per_group)
#' )
#' colnames(matched_matrix) <- paste0("mz_", c(1000, 2000, 3000))
#' group <- rep(c("control", "case"), each = n_per_group)
#'
#' estimate_significance(
#'   matched_matrix = matched_matrix,
#'   group = group,
#'   stat_method = "t.test",
#'   adj_method = "BH"
#' )
#'
#' @seealso [build_matched_matrix()], [stats::p.adjust()]
#'
#' @export
estimate_significance <- function(matched_matrix,
                                  group,
                                  feat_names = NULL,
                                  stat_method = c("t.test", "wilcox"),
                                  adj_method = c("none", "BH", "bonferroni")) {

  stat_method <- match.arg(stat_method)
  adj_method <- match.arg(adj_method)

  # ----------------------------------------
  # Validate matched_matrix
  # ----------------------------------------
  if (missing(matched_matrix) ||
      (!is.matrix(matched_matrix) && !is.data.frame(matched_matrix))) {
    stop(
      "'matched_matrix' must be a numeric matrix or data frame.",
      call. = FALSE
    )
  }

  if (is.data.frame(matched_matrix)) {
    matched_matrix <- as.matrix(matched_matrix)
  }

  if (!is.numeric(matched_matrix)) {
    stop(
      "'matched_matrix' must contain only numeric values.",
      call. = FALSE
    )
  }

  n_samples <- nrow(matched_matrix)
  n_feats <- ncol(matched_matrix)

  if (is.null(n_samples) || is.null(n_feats) ||
      n_samples < 1L || n_feats < 1L) {
    stop(
      "'matched_matrix' must have at least one row and one column.",
      call. = FALSE
    )
  }

  # ----------------------------------------
  # Validate group
  # ----------------------------------------
  if (missing(group) || length(group) != n_samples) {
    stop(
      "'group' must have one entry per row of 'matched_matrix'.",
      call. = FALSE
    )
  }

  if (!is.factor(group)) {
    group <- as.factor(group)
  } else {
    group <- droplevels(group)
  }

  group_levels <- levels(group)
  if (length(group_levels) != 2L) {
    stop(
      sprintf(
        "'group' must have exactly two distinct levels, but %d were found.",
        length(group_levels)
      ),
      call. = FALSE
    )
  }

  is_ref <- group == group_levels[1L]

  # ----------------------------------------
  # Validate feat_names
  # ----------------------------------------
  if (is.null(feat_names)) {
    feat_names <- colnames(matched_matrix)
    if (is.null(feat_names)) {
      feat_names <- paste0("feat_", seq_len(n_feats))
    }
  } else {
    if (!is.character(feat_names) || length(feat_names) != n_feats) {
      stop(
        "'feat_names' must be a character vector with one name per column.",
        call. = FALSE
      )
    }
  }

  # ----------------------------------------
  # Per-feature two-group test
  # ----------------------------------------
  p <- switch(
    stat_method,
    t.test = vapply(seq_len(n_feats), function(i) {
      x <- matched_matrix[is_ref, i]
      y <- matched_matrix[!is_ref, i]
      if (stats::sd(x) == 0 && stats::sd(y) == 0) {
        return(NaN)
      }
      tryCatch(
        stats::t.test(x, y, alternative = "two.sided")$p.value,
        error = function(e) NaN
      )
    }, numeric(1)),
    wilcox = vapply(seq_len(n_feats), function(i) {
      x <- matched_matrix[is_ref, i]
      y <- matched_matrix[!is_ref, i]
      tryCatch(
        stats::wilcox.test(
          x, y,
          alternative = "two.sided",
          exact = FALSE
        )$p.value,
        error = function(e) NaN
      )
    }, numeric(1))
  )

  adj_p <- stats::p.adjust(p = p, method = adj_method)

  out <- data.frame(
    feat_names = feat_names,
    pvalue = p,
    adj_pvalue = adj_p,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL

  return(out)
}
