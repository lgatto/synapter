## maybe these functions should become part of MSnbase

#' common peaks
#' @param x spectrum1 (MSnbase::Spectrum2)
#' @param y spectrum2 (MSnbase::Spectrum2)
#' @param tolerance double, allowed deviation
#' @return logical, common peaks in y
.commonPeaks <- function(x, y, tolerance=25e-6) {
  mx <- mz(x)
  my <- mz(y)

  if (length(mx) == 0 || length(my) == 0) {
    return(logical(length(my)))
  } else if (length(mx) == 1) {
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

