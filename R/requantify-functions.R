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
#' @aliases requantify requantify-method,MSnSet
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
        saturationThreshold=saturationThreshold, method=method, ...)
    }
  }

  dimnames(e) <- dimnames(exprs(msnset))
  ## don't introduce new values for missing entries (could happen if there is
  ## isotopicDistribution but no Counts)
  ## see https://github.com/lgatto/synapter/issues/39#issuecomment-200355965
  isNA <- is.na(exprs(msnset))
  exprs(msnset)[!isNA] <- e[!isNA]
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

.applyByChargeState <- function(x, charges, fun, ...) {
  stopifnot(length(charges) == ncol(x))
  ucharges <- unique(charges)
  m <- matrix(NA_real_, ncol=length(ucharges), nrow=nrow(x),
              dimnames=list(rownames(x), ucharges))
  for (i in seq(along=ucharges)) {
    m[, i] <- fun(x[, charges == ucharges[i], drop=FALSE], ...)
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
                                                        "weighted.mean"),
                                               requantifyAll=FALSE) {
  unsat <- .isUnsaturatedIsotope(x, saturationThreshold=saturationThreshold)
  x <- x * unsat

  cn <- as.numeric(unlist(strsplit(colnames(x), "_", fixed=TRUE)))
  sel <- as.logical(seq_along(cn) %% 2L)
  charges <- cn[sel]
  iIsotopes <- cn[!sel] + 1L
  nIsotopes <- max(iIsotopes)

  probs <- BRAIN::calculateIsotopicProbabilities(
    BRAIN::getAtomsFromSeq(sequence), nrPeaks=nIsotopes)

  th <- t(t(x)/probs[iIsotopes])
  th[th == 0L] <- NA_real_

  r <- switch(match.arg(method),
    "mean" = .applyByChargeState(th, charges, rowMeans, na.rm=TRUE),
    "median" = .applyByChargeState(th, charges, function(x)apply(x, 1, median,
                                                                 na.rm=TRUE)),
    "weighted.mean" = .applyByChargeState(th, charges, function(x)apply(x, 1,
                          function(xx)weighted.mean(xx, w=probs[seq_along(xx)],
                                                    na.rm=TRUE))))
  r <- rowSums(r, na.rm=TRUE)

  if (!requantifyAll) {
    runUnsat <- which(rowSums(!unsat, na.rm=TRUE) == 0L)
    if (length(runUnsat)) {
      r[runUnsat] <- rowSums(x[runUnsat, , drop=FALSE], na.rm=TRUE)
    }
  }
  r
}
