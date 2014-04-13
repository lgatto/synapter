## these functions are only wrapper around some MSnbase functions and
## should never be part of MSnbase

#' get spectra from MSnExp objects
#' same functionality like "get" but doesn't throw an error if the key is
#' not found
#' @param key character, key
#' @param msnexp MSnExp object
#' @return Spectrum2
.getSpectrum <- function(key, msnexp) {
  if (exists(key, envir=assayData(msnexp))) {
    return(msnexp[[key]])
  } else {
    return(.createEmptyMsnbaseSpectrum2(key=key))
  }
}

.getSpectra <- function(keys, spectralist) {
  mapply(.getSpectrum, keys, spectralist)
}


