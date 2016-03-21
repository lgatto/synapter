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
#' In \code{"reference"} the run where all ions are unsaturated and has the
#' highest intensity is used as reference. The other runs are corrected as
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
#' (proportions are used as weights) of these intensities are used as
#' requantified intensity for this peptide.
#'
#' @usage
#' \S4method{requantify}{MSnSet}(object, saturationThreshold,
#' method=c("sum", "reference", "th.mean", "th.median",
#' "th.weighted.mean"), \ldots)
#'
#' @param object An \code{\linkS4class{MSnSet}} object.
#' @param saturationThreshold \code{double}, all intensities above this
#' threshold are considered as saturated.
#' @param method \code{character}, requantification method, please see details
#' section.
#' @param \ldots further arguments passed to internal functions. Currently just
#' \code{onlyCommonIsotopes} for \code{method="sum"} is supported.
#' @return \code{\linkS4class{MSnSet}} where the
#' \code{assayData} are requantified.
#'
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de} and Pavel Shliaha
#' @references
#' See discussion on github: \url{https://github.com/lgatto/synapter/issues/39}
#' @seealso MSnSet documentation: \code{\linkS4class{MSnSet}}
#' @aliases requantify requantify-method,MSnSet
#' @rdname requantify
setMethod("requantify", signature(object="MSnSet"),
          function(object, saturationThreshold,
                   method=c("sum", "reference",
                            "th.mean", "th.median", "th.weighted.mean"),
                   ...) .requantify(object, ...))

.requantify <- function(msnset, saturationThreshold,
                        method=c("sum", "reference",
                                 "th.mean", "th.median", "th.weighted.mean"),
                        ...) {
  cn <- fvarLabels(msnset)
  i <- grep("isotopicDistr", cn)

  if (!length(i)) {
    stop("Could not find any isotopic distribution information.")
  }

  f <- fData(msnset)[, i]
  e <- exprs(msnset)

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
        saturationThreshold=saturationThreshold, method=method)
    }
  }

  dimnames(e) <- dimnames(exprs(msnset))
  exprs(msnset) <- e
  msnset
}

.isUnsaturatedIsotope <- function(x, saturationThreshold=Inf) {
  x < saturationThreshold
}

.runsUnsaturated <- function(x, saturationThreshold=Inf) {
  colSums(!.isUnsaturatedIsotope(x, saturationThreshold=saturationThreshold),
          na.rm=TRUE) == 0L
}

.isCommonIsotope <- function(x) {
  colSums(is.na(x)) == 0L
}

.sumIsotopes <- function(x) {
  isoNames <- factor(gsub("^[0-9]+_([0-9]+)$", "\\1", colnames(x)))
  m <- matrix(NA_real_, ncol=length(levels(isoNames)), nrow=nrow(x),
              dimnames=list(rownames(x), levels(isoNames)))
  for (i in seq(along=levels(isoNames))) {
    m[, i] <- rowSums(x[, isoNames == levels(isoNames)[i], drop=FALSE],
                      na.rm=TRUE)
  }
  m
}

.requantifySum <- function(x, saturationThreshold=Inf,
                           onlyCommonIsotopes=TRUE) {
  i <- .runsUnsaturated(x, saturationThreshold=saturationThreshold)

  if (onlyCommonIsotopes) {
    i <- i & .isCommonIsotope(x)
  }
  rowSums(x[, i, drop=FALSE], na.rm=TRUE)
}

.requantifyReferenceRun <- function(x, saturationThreshold=Inf) {
  unsat <- .isUnsaturatedIsotope(x, saturationThreshold=saturationThreshold)

  runSat <- rowSums(!unsat, na.rm=TRUE) == 0L
  rs <- rowSums(x, na.rm=TRUE)
  runSums <- rs * runSat

  ref <- which.max(runSums)

  for (i in seq(along=runSat)) {
    if (!runSat[i]) {
      co <- .isCommonIsotope(x[c(ref, i),, drop=FALSE]) & unsat[i,]
      if (sum(co)) {
        f <- mean(x[i, co]/x[ref, co], na.rm=TRUE)
        # use ref run as estimator for saturated intensities
        runSums[i] <- sum(f*x[ref, ]*!unsat[i, ], x[i, ]*unsat[i, ],
                          na.rm=TRUE)
      } else {
        return(rep.int(NA_real_, length(runSums)))
      }
    }
  }
  runSums
}

.requantifyTheoreticalDistribution <- function(x, sequence,
                                               saturationThreshold=Inf,
                                               method=c("mean", "median",
                                                        "weighted.mean")) {
  unsat <- .isUnsaturatedIsotope(x, saturationThreshold=saturationThreshold)
  x <- .sumIsotopes(x * unsat)
  nc <- ncol(x)

  props <- BRAIN::calculateIsotopicProbabilities(
    BRAIN::getAtomsFromSeq(sequence), nrPeaks=nc)

  th <- t(t(x)/props)
  th[th == 0L] <- NA_real_

  switch(match.arg(method),
         "mean" = rowMeans(th, na.rm=TRUE),
         "median" = apply(th, 1, median, na.rm=TRUE),
         "weighted.mean" = apply(th, 1, weighted.mean, w=props,
                                 na.rm=TRUE))
}
