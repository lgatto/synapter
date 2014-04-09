## maybe these functions should become part of MSnbase

#' common peaks
#' @param x spectrum1 (MSnbase::Spectrum2)
#' @param y spectrum2 (MSnbase::Spectrum2)
#' @param tolerance double, allowed deviation
#' @return logical, common peaks in y
.commonPeaks <- function(x, y, tolerance=25e-6) {

  if (peaksCount(x) == 0) {
    return(logical(peaksCount(y)))
  } else if (peaksCount(y) == 0) {
    return(logical(peaksCount(x)))
  }

  mx <- mz(x)
  my <- mz(y)

  if (length(mx) == 1) {
    return(abs(mx-my)/my < tolerance)
  }

  ## adopted from MALDIquant:::.which.closest
  ## find left interval
  lIdx <- findInterval(my, mx, rightmost.closed=FALSE, all.inside=TRUE)
  ## find right interval
  rIdx <- lIdx+1L

  ## calculate relative differences for left and right nearest point
  lDiff <- abs(mx[lIdx]-my)/mx[lIdx]
  rDiff <- abs(mx[rIdx]-my)/mx[rIdx]

  return(pmin(lDiff, rDiff) < tolerance)
}

#' number of common peaks
#' @param x spectrum1 (MSnbase::Spectrum2)
#' @param y spectrum2 (MSnbase::Spectrum2)
#' @param tolerance double, allowed deviation
#' @return double, number of common peaks
.nCommonPeaks <- function(x, y, tolerance=25e-6) {
  return(sum(.commonPeaks(x, y, tolerance=tolerance)))
}


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


