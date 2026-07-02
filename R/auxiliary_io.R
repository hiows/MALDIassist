# Internal I/O helpers for Bruker MALDI-TOF spectrum loading.

#' @name .search_fid_dir
#' @title Search Bruker fid directories
#'
#' @description
#' Recursively searches for Bruker `1SLin` marker directories beneath `dir`.
#' The search strategy follows the general concept of `readBrukerFlexData`.
#'
#' @param dir Character. Root directory containing Bruker MALDI-TOF data.
#'
#' @return A character vector of paths to spectrum directories.
#'
#' @keywords internal
#' @noRd
.search_fid_dir <- function(dir){
  path <- list.files(dir, "1SLin", T,T,T,T,T)
  return(path)
}


#' @name .open_fid_file
#' @title Read Bruker acqu parameters and compute m/z values
#'
#' @description
#' Reads the local `acqu` file in the current working directory and converts
#' time-of-flight values to m/z using Bruker calibration parameters.
#'
#' @return A list containing `mass`, `number`, `TimeDelay`, `TimeDelta`, `tof`,
#'   `c1`, `c2`, and `c3`.
#'
#' @keywords internal
#' @noRd
.open_fid_file = function(){
  
  # Step 1. read fid file
  acquFile = sub(pattern = "fid$", x = "fid", replacement = "acqu")
  con = file(acquFile, "rt")
  acquLines = readLines(con, n=-1)
  close(con)
  
  
  # Step 2. extract data
  ## 2.1) nIntensity
  srcStr = acquLines
  patternStr_number = "##\\$TD="
  tmpLine_number = grep(pattern = patternStr_number, x = srcStr, value = TRUE)
  # Locate the acquisition metadata line containing the total point count.
  
  res_number = gsub(pattern="(^.*= *<?)|(>? *$)", replacement = "", x = tmpLine_number)
  res_number = as.double(gsub(",", replacement = ".", res_number))
  
  ## 2.2) TimeDelay
  patternStr_TimeDelay = "##\\$DELAY="
  tmpLine_TimeDelay = grep(pattern = patternStr_TimeDelay, x = srcStr, value = TRUE)
  res_TimeDelay = gsub(pattern = "(^.*= *<?)|(>? *$)", replacement = "", x = tmpLine_TimeDelay)
  res_TimeDelay = as.double(gsub(",", replacement = ".", res_TimeDelay))
  
  ## 2.3) TimeDelta
  patternStr_TimeDelta = "##\\$DW="
  tmpLine_TimeDelta = grep(pattern=patternStr_TimeDelta, x = srcStr, value = TRUE)
  res_TimeDelta = gsub(pattern = "(^.*= *<?)|(>? *$)", replacement = "", x = tmpLine_TimeDelta)
  res_TimeDelta = as.double(gsub(",", replacement = ".", res_TimeDelta))
  
  ## 2.4) tof
  tof = as.double(res_TimeDelay + ((0:(res_number-1)) * res_TimeDelta))
  
  ## 2.5) c1, c2, c3
  temp_c1 = grep(pattern = "##\\$ML1=", x = srcStr, value = TRUE)
  c1 = gsub(pattern = "##\\$ML1=", replacement = "", x = temp_c1) |> as.double(gsub(",", replacement = ".", temp_c1))
  temp_c2 = grep(pattern = "##\\$ML2=", x = srcStr, value = TRUE)
  c2 = gsub(pattern = "##\\$ML2=", replacement = "", x = temp_c2) |> as.double(gsub(",", replacement = ".", temp_c2))
  temp_c3 = grep(pattern = "##\\$ML3=", x = srcStr, value = TRUE)
  c3 = gsub(pattern = "##\\$ML3=", replacement = "", x = temp_c3) |> as.double(gsub(",", replacement = ".", temp_c3))
  
  
  # Step 3. calcuate m/z
  A <- c3
  B <- sqrt(1e+12/c1)
  C <- c2 - tof
  
  if (A == 0) {
    # linear: 0 = B * sqrt(m/z) + C(times)
    mass = (C * C)/(B * B)
  } else {
    # quadratic: 0 = A * (sqrt(m/z))^2 + B * sqrt(m/z) + C(times)
    mass = ((-B + sqrt((B * B) - (4 * A * C)))/(2 * A))^2
  }
  
  
  # Step 4. output
  res = list(mass = mass, 
             number = res_number, TimeDelay = res_TimeDelay, TimeDelta = res_TimeDelta,
             tof = tof, c1 = c1, c2 = c2, c3 = c3)
  
  return(res)
}


#' @name .load_spectrum
#' @title Load one Bruker raw spectrum from the current fid directory
#'
#' @description
#' Reads intensity values from the local `fid` file and pairs them with m/z
#' values obtained from `.open_fid_file()`. The working directory must already
#' point to a Bruker `fid` directory.
#'
#' @return A two-column `data.frame` with `mz` and `intensity`.
#'
#' @keywords internal
#' @noRd
.load_spectrum = function(){
  
  # Step 1. open fid file
  con_fid = file("fid", "rb")
  on.exit(close(con_fid))
  t_data = .open_fid_file()
  
  
  # Step 2. extract data
  ## 2.1) m/z
  mz = t_data$mass
  
  ## 2.2) intensity
  nIntensity = t_data$number
  # intensity = as.double(readBin(file("fid", "rb"), integer(), n=nIntensity))
  intensity = as.double(readBin(con_fid, what = integer(), n = nIntensity, size = 4, endian = "little"))
  intensity[intensity < 0] = 0
  
  
  # Step 3. output
  raw_spectrum = data.frame(mz = mz, intensity = intensity)
  
  return(raw_spectrum)
}
