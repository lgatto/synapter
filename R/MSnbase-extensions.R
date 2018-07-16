#' Synapter/PLGS Agreement
#'
#' This method checks the agreement between synapter analysis and PLGS results.
#'
#' @details
#' Each synapter object has \code{synapterPlgsAgreement} column in its
#' \code{MatchedEMRTs} \code{data.frame} (see \code{\link{writeMatchedEMRTs}}).
#' After converting the synapter object into an \code{\linkS4class{MSnSet}}
#' instance via \code{as(synapterobject, "MSnSet"} this column could be find in
#' the feature data (\code{fData(msnset)$synapterPlgsAgreement}).\cr
#' In the \code{synapterPlgsAgreement} each peptide is classified as:
#' \itemize{
#'  \item{\code{"agree"}: EMRT identified in identification and quantitation run
#'  by PLGS and same EMRT matched in synapter's grid search.}
#'  \item{\code{"disagree"}: EMRT identified in identification and quantitation
#'  run by PLGS and a different EMRT was matched in synapter's grid search.}
#'  \item{\code{"no_plgs_id"}: EMRT was \emph{not} identified in the
#'  quantitation run by PLGS but matched in synapter's grid search.}
#'  \item{\code{"no_synapter_transfer"}: EMRT was identified in the
#'  identification and quantitation run by PLGS but not matched in synapter's
#'  grid search.}
#'  \item{\code{"no_id_or_transfer"}: EMRT was \emph{not} identified in the
#'  quantitation run by PLGS and not matched in synapter's grid search.}
#'  \item{\code{"multiple_ident_matches"}: a single quantitation EMRT was
#'  matched by synapter to multiple identification EMRTs found by PLGS (could
#'  happen if the grid search parameters are too relaxed).}
#' }
#' After combining multiple \code{\linkS4class{MSnSet}} the method
#' \code{synapterPlgsAgreement} adds additional columns to the feature data:
#' \itemize{
#'  \item{\code{nIdentified}: how often a peptide was identified across
#'  multiple runs?}
#'  \item{\code{nAgree}: how often a peptide was identified by PLGS and synapter
#'  across multiple runs (counts \code{"agree"} entries)?}
#'  \item{\code{nDisagree}: how often a peptide was differently identified by
#'  PLGS and synapter across multiple runs (counts \code{"disagree"} entries)?}
#'  \item{\code{synapterPlgsAgreementRatio}: \code{nAgree/(nAgree +
#'  nDisagree)}.}
#' }
#'
#' @usage
#' \S4method{synapterPlgsAgreement}{MSnSet}(object, \dots)
#'
#' @param object An \code{\linkS4class{MSnSet}} object.
#' @param \dots further arguments, not used yet.
#' @return \code{\linkS4class{MSnSet}} where the columns \code{nIdentified},
#' \code{nAgree}, \code{nDisagree} and \code{synapterPlgsAgreementRatio} were
#' added to the feature data.
#'
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @references
#' See discussion on github: \url{https://github.com/lgatto/synapter/issues/73}
#' @seealso MSnSet documentation: \code{\linkS4class{MSnSet}}
#' @aliases synapterPlgsAgreement
#' @rdname synapterPlgsAgreement
setMethod("synapterPlgsAgreement", signature(object="MSnSet"),
          function(object, ...) .synapterPlgsAgreement(object))

.synapterPlgsAgreement <- function(msnset) {

  i <- grep("synapterPlgsAgreement", fvarLabels(msnset))

  ## metric a: how often identified
  fData(msnset)$nIdentified <- rowSums(!is.na(exprs(msnset)))

  ## metric b: how often PLGS/synapter agree
  fData(msnset)$nAgree <-
    rowSums(fData(msnset)[, i, drop=FALSE] == "agree", na.rm=TRUE)

  ## metric c: how often PLGS/synapter disagree
  fData(msnset)$nDisagree <-
    rowSums(fData(msnset)[, i, drop=FALSE] == "disagree", na.rm=TRUE)

  ## calculate ratio for easier interpretation
  fData(msnset)$synapterPlgsAgreementRatio <-
    fData(msnset)$nAgree/(fData(msnset)$nAgree + fData(msnset)$nDisagree)

  msnset
}

#' Correct Intensity
#'
#' This method corrects the intensity values of an \code{\linkS4class{MSnSet}}
#' object by applying the intensity model built in the synapter workflow by
#' \code{\link{modelIntensity}}.
#'
#' @param object An \code{\linkS4class{MSnSet}} object.
#' @param method Correct (\code{*}) or undo correction (\code{/}).
#' @param \dots further arguments, not used yet.
#' @return A \code{\linkS4class{MSnSet}} with corrected intensity values
#' (\code{exprs}).
#'
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @seealso MSnSet documentation: \code{\linkS4class{MSnSet}}
#' @aliases correctIntensity
#' @noRd
.correctIntensity <- function(msnset, method=c("correct", "undo")) {

  method <- match.arg(method)

  i <- grep("intensityCorrectionFactor", fvarLabels(msnset))

  if (length(i)) {
    m <- as.matrix(fData(msnset)[, i, drop=FALSE])

    if (method == "undo") {
      exprs(msnset) <- exprs(msnset) / m
    } else {
      exprs(msnset) <- exprs(msnset) * m
    }
  }

  msnset
}

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
