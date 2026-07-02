#' @name load_maldi_spectra
#' @title Load Bruker MALDI-TOF MS Spectra
#'
#' @description
#' Loads Bruker MALDI-TOF MS spectra from a directory containing Bruker flex
#' data files.
#'
#' This function searches for spectrum directories containing Bruker \code{fid}
#' files, loads each spectrum using the internal \code{.load_spectrum()}
#' function, assigns sample names based on the directory structure, and returns
#' a named list of raw spectra.
#'
#' The loading procedure is designed for Bruker MALDI-TOF MS directory
#' structures and follows the general concept of \code{readBrukerFlexData}.
#'
#' @param spectra_dir Character. Path to the directory containing Bruker
#'   MALDI-TOF MS spectra.
#' @param return_dir Character or \code{NULL}. Directory to return to after
#'   loading spectra. If \code{NULL}, the current working directory is used.
#'   Default is \code{NULL}.
#'
#' @return A named list of raw spectra. Each element corresponds to one loaded
#'   MALDI-TOF MS spectrum. The list names are generated from the sample
#'   directory structure and sorted using \code{gtools::mixedsort()}.
#'
#' @details
#' The function first normalizes \code{spectra_dir} and \code{return_dir}, then
#' searches recursively for Bruker \code{fid} directories using the internal
#' \code{.search_fid_dir()} function.
#'
#' Each spectrum is loaded after temporarily setting the working directory to
#' the corresponding \code{fid} directory. After loading is completed, the
#' working directory is changed back to \code{return_dir}.
#'
#' Sample names are inferred from directories named \code{"1SLin"} in the
#' Bruker directory hierarchy. If duplicate sample names are detected, parent
#' directory names are appended to make the names more distinguishable.
#'
#' This function assumes that the internal helper functions
#' \code{.search_fid_dir()} and \code{.load_spectrum()} are available.
#'
#' @seealso \code{gtools::mixedsort}
#'
#' @examples
#' \dontrun{
#' raw_spectra <- load_maldi_spectra(
#'   spectra_dir = "path/to/bruker/spectra"
#' )
#'
#' names(raw_spectra)
#' }
#'
#' @importFrom gtools mixedsort
#' @export
load_maldi_spectra = function(spectra_dir, return_dir = NULL){
  
  # STEP 1 ----
  if (is.null(return_dir)) {
    return_dir <- getwd()
  }
  
  return_dir <- normalizePath(return_dir, mustWork = TRUE)
  spectra_dir <- normalizePath(spectra_dir, mustWork = TRUE)
  
  
  # STEP 2 ----
  fid_dirs <- .search_fid_dir(dir = spectra_dir)
  fid_dirs <- normalizePath(fid_dirs, mustWork = TRUE)
  
  raw_spectra = lapply(fid_dirs, function(dir){
    setwd(dir)
    .load_spectrum()
  })
  
  
  # STEP 3 ----
  dir_1SLin <- list.files(
    path = spectra_dir,
    recursive = TRUE,
    full.names = TRUE,
    include.dirs = TRUE
  )
  dir_1SLin <- dir_1SLin[basename(dir_1SLin) == "1SLin"]
  list_names = basename(dirname(dirname(dirname(dir_1SLin))))
  
  i = 2
  while (any(table(list_names) > 1) & i > 0) {
    tmp = basename(Reduce(function(z, i) dirname(z), x = seq_len(i), dir_1SLin))
    replaced = names(which(table(list_names) > 1))
    list_names = ifelse(list_names %in% replaced, paste(replaced, tmp, sep = "/"), list_names)
    i = i - 1
  }
  
  names(raw_spectra) = list_names
  raw_spectra = raw_spectra[mixedsort(names(raw_spectra))]
  
  # STEP 4 ----
  setwd(return_dir)
  
  return(raw_spectra)
}
