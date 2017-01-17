#' Intensity requantification
#'
#' This method tries to remove saturation effects from the intensity counts.
#'
#' @details
#' Currently \code{requantify} supports 3 (5) different requantification
#' methods. \cr
#'
#' \code{"sum"} is the simplest requantification method. All ions of a
#' peptide below the saturation threshold are summed to get the new intensity.
#' This method accept an additional argument, namely \code{onlyCommonIsotopes}
#' If \code{onlyCommonIsotopes=TRUE} (default) all ions that are not seen in
#' all runs are removed and only the common seen ions are summed. In contrast
#' \code{onlyCommonIsotopes=FALSE} sums all ions regardless they are present
#' in all runs. \cr
#'
#' In \code{"reference"} the run that has the most unsaturated ions in common
#' with all the other runs. If there are more than one run, the most intense is
#' used as reference. The other runs are corrected as
#' follows:
#' \itemize{
#'  \item{Find common ions between current and the reference run.}
#'  \item{Divide the intensities of the common ions and calculate the mean of
#'  these quotients as a run specific scaling factor.}
#'  \item{Multiply the unsaturated ions of the current run by the scaling
#'  factor and replace the saturated ones by the product of the scaling
#'  factor and the intensities of their corresponding ions in the reference
#'  run.}
#'  \item{Sum the rescaled ion intensities.}
#' }
#'
#' The \code{"th.*"} methods are nearly identical. All of them calculate the
#' theoretical isotopic distribution for the given sequence of the peptide.
#' Subsequently the unsaturated ions are divided by their theoretical
#' proportion and the \code{mean}/\code{median}/\code{weighted.mean}
#' (proportions are used as weights) of these intensities are calculated per
#' charge state. The sum of the charge state values is used as requantified
#' intensity for this peptide. \cr
#' If \code{requantifyAll=FALSE} (default) just peptides with at least one
#' saturated ion are requantified (unsaturated peptides are unaffected). If
#' \code{requantify=TRUE} all peptides even these where all ions are below
#' \code{saturationThreshold} are requantified by their theoretical
#' distribution.
#'
#' @usage
#' \S4method{requantify}{MSnSet}(object, saturationThreshold,
#' method=c("sum", "reference", "th.mean", "th.median",
#' "th.weighted.mean"), \ldots)
#'
#' @param object An \code{\linkS4class{MSnSet}} object.
#' @param saturationThreshold \code{double}, intensity of an ion (isotope of a
#' given charge state) at which saturation is starting to occur.
#' @param method \code{character}, requantification method, please see details
#' section.
#' @param \ldots further arguments passed to internal functions. Currently
#' \code{onlyCommonIsotopes} for \code{method="sum"} and \code{requantifyAll}
#' for \code{method=c("th.mean", "th.median", "th.weighted.mean")} are supported.
#' @return \code{\linkS4class{MSnSet}} where the
#' \code{assayData} are requantified.
#'
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de} and Pavel Shliaha
#' @references
#' See discussion on github: \url{https://github.com/lgatto/synapter/issues/39}
#' @seealso MSnSet documentation: \code{\linkS4class{MSnSet}}
#' @aliases requantify
#' @rdname requantify
setMethod("requantify", signature(object="MSnSet"),
          function(object, saturationThreshold,
                   method=c("sum", "reference",
                            "th.mean", "th.median", "th.weighted.mean"),
                   ...) {
            .requantify(object, saturationThreshold=saturationThreshold,
                        method=method, ...)
})

.requantify <- function(msnset, saturationThreshold,
                        method=c("sum", "reference",
                                 "th.mean", "th.median", "th.weighted.mean"),
                        ...) {
  i <- grep("isotopicDistr", fvarLabels(msnset))

  if (!length(i)) {
    stop("Could not find any isotopic distribution information.")
  }

  msnset <- .correctIntensity(msnset, method="undo")

  e <- exprs(msnset)
  f <- as.matrix(fData(msnset)[, i])
  ## don't introduce new values for missing entries (could happen if there is
  ## isotopicDistribution but no Counts)
  ## see https://github.com/lgatto/synapter/issues/39#issuecomment-200355965
  isNA <- is.na(e)
  f[isNA] <- NA_character_

  method <- match.arg(method)

  if (substring(method, 1, 2) != "th") {
    fun <- switch(method,
           "sum" = {
             .requantifySum
           },
           "reference" = {
             .requantifyReferenceRun
           }
    )
    e <- t(apply(f, 1, function(x) {
      fun(.isotopicDistr2matrix(x),
          saturationThreshold=saturationThreshold, ...)
    }))
  } else {
    if (!requireNamespace("BRAIN")) {
      stop("Please install the BRAIN package via 'biocLite(\"BRAIN\")' to use this method.")
    }
    ## remove "th." prefix
    method <- substring(method, 4, nchar(method))
    rn <- rownames(e)

    for (i in 1:nrow(e)) {
      e[i, ] <- .requantifyTheoreticalDistribution(
        .isotopicDistr2matrix(f[i, ]), sequence=rn[i],
        saturationThreshold=saturationThreshold, method=method, ...)
    }
  }

  exprs(msnset)[!isNA] <- e[!isNA]
  .correctIntensity(msnset, method="correct")
}

.isUnsaturatedIsotope <- function(x, saturationThreshold=Inf) {
  x < saturationThreshold
}

.runsUnsaturated <- function(x, saturationThreshold=Inf) {
  colSums(!.isUnsaturatedIsotope(x, saturationThreshold=saturationThreshold),
          na.rm=TRUE) == 0L
}

.isCommonIsotope <- function(x) {
  colSums(is.na(x)) == 0L & nrow(x) > 0L
}

.filterMissingRuns <- function(x) {
  x[rowSums(!is.na(x)) != 0L,, drop=FALSE]
}

.names2chargesIsotopes <- function(x) {
  cn <- as.numeric(unlist(strsplit(x, "_", fixed=TRUE)))
  sel <- seq(1L, length(cn), by=2L)
  list(charges=cn[sel], isotopes=cn[sel + 1L])
}

.referenceRun <- function(x, unsat) {
  unsat[is.na(unsat)] <- FALSE
  ref <- unsat %*% t(unsat)
  diag(ref) <- NA_real_
  refrs <- rowSums(ref, na.rm=TRUE)
  which.max((refrs == max(refrs)) * rowSums(x, na.rm=TRUE))
}

.requantifySum <- function(x, saturationThreshold=Inf,
                           onlyCommonIsotopes=TRUE) {
  i <- .runsUnsaturated(x, saturationThreshold=saturationThreshold)

  if (onlyCommonIsotopes) {
    i <- i & .isCommonIsotope(.filterMissingRuns(x))
  }
  rs <- rowSums(x[, i, drop=FALSE], na.rm=TRUE)
  rs[rs == 0L] <- NA_real_
  rs
}

.requantifyReferenceRun <- function(x, saturationThreshold=Inf) {
  unsat <- .isUnsaturatedIsotope(x, saturationThreshold=saturationThreshold)

  runSums <- rowSums(x, na.rm=TRUE)
  runSat <- rowSums(!unsat, na.rm=TRUE) == 0L
  runSums[!runSat] <- NA_real_

  ref <- .referenceRun(x, unsat)

  for (i in seq(along=runSat)) {
    if (!runSat[i]) {
      co <- .isCommonIsotope(x[c(ref, i),, drop=FALSE]) & unsat[i,]
      if (sum(co)) {
        f <- mean(x[i, co]/x[ref, co], na.rm=TRUE)
        # use ref run as estimator for saturated intensities
        runSums[i] <- sum(f*x[ref, ]*!unsat[i, ], x[i, ]*unsat[i, ],
                          na.rm=TRUE)
      }
    }
  }
  runSums
}

.requantifyTheoreticalDistribution <- function(x, sequence,
                                               saturationThreshold=Inf,
                                               method=c("mean", "median",
                                                        "weighted.mean"),
                                               requantifyAll=FALSE) {
  unsat <- .isUnsaturatedIsotope(x, saturationThreshold=saturationThreshold)
  x <- x * unsat

  ci <- .names2chargesIsotopes(colnames(x))
  ci$isotopes <- ci$isotopes + 1L

  probs <- BRAIN::calculateIsotopicProbabilities(
    BRAIN::getAtomsFromSeq(sequence), nrPeaks=max(ci$isotopes))

  th <- t(t(x)/probs[ci$isotopes])
  th[th == 0L] <- NA_real_

  fun <- switch(match.arg(method),
                "mean" = function(x)colMeans(x, na.rm=TRUE),
                "median" = function(x)apply(x, 2, median, na.rm=TRUE),
                "weighted.mean" = function(x)apply(x, 2,
                                    function(xx)weighted.mean(xx, w=probs[seq_along(xx)], na.rm=TRUE)))

  r <- t(MSnbase:::utils.applyColumnwiseByGroup(t(th), groupBy = ci$charges, FUN = fun))
  r <- rowSums(r, na.rm=TRUE)

  if (!requantifyAll) {
    runUnsat <- which(rowSums(!unsat, na.rm=TRUE) == 0L)
    if (length(runUnsat)) {
      r[runUnsat] <- rowSums(x[runUnsat, , drop=FALSE], na.rm=TRUE)
    }
  }
  r
}

#' Rescale for TOP3
#'
#' This method rescales the intensity values of an \code{\linkS4class{MSnSet}}
#' object to be suiteable for TOP3 quantification.
#'
#' @details
#' If an \code{\linkS4class{MSnSet}} object was requantified using the
#' \code{method="sum"} requantification method
#' (see \code{\link{requantify,MSnSet-method}}) TOP3 is not valid anymore
#' because the most abundant proteins are penalised by removing high intensity
#' isotopes.
#'
#' To overcome this \code{rescaleForTop3} takes the proportion
#' \code{isotope/sum(isotopes} for each requantified peptide and calculates a
#' correction factor by comparing these proportions against the unsaturated
#' isotopes before requantification. The new rescale intensity values of the
#' isotopes are the mean correction factor multiplied with the corrected
#' intensity values (see
#' \url{https://github.com/lgatto/synapter/issues/39#issuecomment-207987278} for
#' the complete explanation/discussion).
#'
#' @usage
#' \S4method{rescaleForTop3}{MSnSet,MSnSet}(before, after, saturationThreshold,
#' onlyForSaturatedRuns=TRUE, \ldots)
#'
#' @param before An \code{\linkS4class{MSnSet}} object before requantification.
#' @param after The same \code{\linkS4class{MSnSet}} object as \code{before} but
#' after requantification.
#' @param saturationThreshold \code{double}, intensity of an ion (isotope of a
#' given charge state) at which saturation is starting to occur.
#' @param onlyForSaturatedRuns \code{logical}, rescale just runs where at least
#' one isotope is affected by saturation.
#' @param \ldots further arguments passed to internal functions. Currently
#' ignored.
#' @return \code{\linkS4class{MSnSet}} where the
#' \code{assayData} are requantified.
#'
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de} and Pavel Shliaha
#' @references
#' See discussion on github: \url{https://github.com/lgatto/synapter/issues/39}
#' @seealso MSnSet documentation: \code{\linkS4class{MSnSet}}
#' \url{https://github.com/lgatto/synapter/issues/39#issuecomment-207987278}
#' @aliases rescaleForTop3
#' @rdname rescaleForTop3

setMethod("rescaleForTop3", signature(before="MSnSet", after="MSnSet"),
          function(before, after, saturationThreshold, onlyForSaturatedRuns=TRUE, ...) {
            .rescaleForTop3(before=before, after=after,
                            saturationThreshold=saturationThreshold,
                            onlyForSaturatedRuns=onlyForSaturatedRuns, ...)
})

.rescaleForTop3 <- function(before, after, saturationThreshold, onlyForSaturatedRuns=TRUE) {
  i <- grep("isotopicDistr", fvarLabels(before))

  if (!length(i)) {
    stop("Could not find any isotopic distribution information.")
  }

  eBefore <- exprs(before)
  eAfter <- exprs(after)
  f <- as.matrix(fData(before)[, i])
  ## don't introduce new values for missing entries (could happen if there is
  ## isotopicDistribution but no Counts)
  ## see https://github.com/lgatto/synapter/issues/39#issuecomment-200355965
  isNA <- is.na(eBefore)
  f[isNA] <- NA_character_

  unsat <- t(apply(f, 1, function(x) {
    .runsUnsaturated(t(.isotopicDistr2matrix(x)),
                     saturationThreshold=saturationThreshold)
  }))

  eNew <- eBefore
  eNew[is.na(eAfter)] <- NA_real_
  prop <- eAfter/rowSums(eAfter, na.rm=TRUE)
  cf <- eBefore/prop
  cf[unsat] <- NA_real_

  if (!onlyForSaturatedRuns) {
    unsat[rowSums(!unsat & !is.na(cf), na.rm=TRUE) != 0L, ] <- FALSE
  }
  eNew[!unsat] <- (rowMeans(cf, na.rm=TRUE) * prop)[!unsat]

  exprs(after) <- eNew
  after
}
