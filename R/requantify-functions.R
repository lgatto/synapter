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

.requantifySum <- function(x, saturationThreshold=Inf,
                           onlyCommonIsotopes=FALSE) {
  i <- .runsUnsaturated(x, saturationThreshold=saturationThreshold)

  if (onlyCommonIsotopes) {
    i <- i & .isCommonIsotope(x)
  }
  rowSums(x[, i], na.rm=TRUE)
}

.requantifyReferenceRun <- function(x, saturationThreshold=Inf) {
  unsat <- .isUnsaturatedIsotope(x, saturationThreshold=saturationThreshold)

  runSat <- rowSums(!unsat, na.rm=TRUE) == 0L
  rs <- rowSums(x, na.rm=TRUE)
  runSums <- rs * runSat

  ref <- which.max(runSums)

  for (i in seq(along=runSat)) {
    if (!runSat[i]) {
      co <- .isCommonIsotope(x[c(ref, y),, drop=FALSE]) & unsat[y,]
      f <- mean(x[y, co]/x[ref, co], na.rm=TRUE)
      # use ref run as estimator for saturated intensities
      runSums[i] <- sum(f*x[ref, ]*!unsat[y, ], x[y,]*unsat[y,], na.rm=TRUE)
    }
  }
  runSums
}

