requantify <- function(msnset, saturationThreshold,
                       method=c("sum", "reference",
                                "th.mean", "th.median", "th.weighted.mean")) {
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
          saturationThreshold=saturationThreshold)
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

  return(switch(match.arg(method),
                "mean" = rowMeans(th, na.rm=TRUE),
                "median" = apply(th, 1, median, na.rm=TRUE),
                "weighted.mean" = apply(th, 1, weighted.mean, w=props, na.rm=TRUE)))
}
