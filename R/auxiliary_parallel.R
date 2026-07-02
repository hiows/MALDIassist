# Internal parallel execution helpers.

#' @name .parallel_lapply
#' @title Apply a function to list elements with optional parallel execution
#'
#' @param X A list or vector to iterate over.
#' @param FUN Function applied to each element.
#' @param n_cores Number of worker processes. Values greater than 1 enable
#'   parallel execution.
#'
#' @return A list of results in the same order as `X`.
#'
#' @keywords internal
#' @noRd
.parallel_lapply <- function(X, FUN, n_cores = 1L) {
  n_cores <- as.integer(n_cores)
  
  if (length(X) == 0L || n_cores <= 1L) {
    return(lapply(X, FUN))
  }
  
  n_workers <- min(n_cores, length(X))
  
  if (.Platform$OS.type == "unix") {
    return(
      parallel::mclapply(
        X = X,
        FUN = FUN,
        mc.cores = n_workers
      )
    )
  }
  
  cl <- parallel::makeCluster(n_workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::parLapply(
    cl = cl,
    X = X,
    fun = FUN
  )
}
