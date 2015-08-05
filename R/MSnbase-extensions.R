## these functions are only wrapper around some MSnbase functions and
## should never be part of MSnbase

#' get spectra from MSnExp objects
#' same functionality like "get" but doesn't throw an error if the key is
#' not found
#' @param key character, key
#' @param msnexp MSnExp object
#' @return Spectrum2
#' @noRd
.getSpectrum <- function(key, msnexp) {
  if (!is.character(key)) {
    key <- as.character(key)
  }

  if (exists(key, envir=assayData(msnexp))) {
    return(msnexp[[key]])
  } else {
    return(.createEmptyMsnbaseSpectrum2(key=key))
  }
}

.sumAllPeakCounts <- function(msexp) {
  sum(unlist(peaksCount(msexp)))
}

#' handle plgs/synapter agreement see issue 73:
#' https://github.com/lgatto/synapter/issues/73
#' TODO: should be a method for MSnSets in future
#' @noRd
synapterPlgsAgreement <- function(msnset) {

  i <- grep("synapterPlgsAgreement", colnames(fData(msnset)))

  ## metric a: how often identified
  fData(msnset)$nIdentified <- rowSums(!is.na(exprs(msnset)))

  ## metric b: how often PLGS/synapter agree
  fData(msnset)$nAgree <- rowSums(fData(msnset)[, i] == "agree")

  ## metric c: how often PLGS/synapter disagree
  fData(msnset)$nDisagree <- rowSums(fData(msnset)[, i] == "disagree")

  ## calculate ratio for easier interpretation
  fData(msnset)$synapterPlgsAgreementRatio <-
    fData(msnset)$nAgree/(fData(msnset)$nAgree + fData(msnset)$nDisagree)

  msnset
}

