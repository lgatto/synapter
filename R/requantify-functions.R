requantify <- function(msnset, saturationThreshold,
                       method=c("sum", "th.mean", "th.median", "th.weighted.mean")) {
  cn <- fvarLabels(msnset)
  i <- grep("isotopicDistr", cn)

  if (!length(i)) {
    stop("Could not find any isotopic distribution information.")
  }

  f <- fData(msnset)[, i]
  e <- exprs(msnset)

  method <- match.arg(method)

  if (method == "sum") {
    e <- t(apply(f, 1, .requantifySum,
                 saturationThreshold=saturationThreshold))
  } else {
    if (!requireNamespace("BRAIN")) {
      stop("Please install the BRAIN package via 'biocLite(\"BRAIN\")' to use this method.")
    }
    ## remove "th." prefix
    method <- substring(method, 4, nchar(method))
    rn <- rownames(e)

    for (i in 1:nrow(e)) {
      e[i, ] <- .requantifyTheoreticalDistribution(f[i, ], sequence=rn[i],
                                                   saturationThreshold=saturationThreshold,
                                                   method=method)
    }
  }

  dimnames(e) <- dimnames(exprs(msnset))
  exprs(msnset) <- e
  msnset
}

.requantifySum <- function(x, saturationThreshold=Inf) {
  x <- .splitIsotopicDistr(unlist(x))
  commonnm <- .commonIsotopes(x, saturationThreshold)

  if (length(commonnm)) {
    return(sapply(x, function(x)sum(x[commonnm])))
  }

  rep.int(NA, length(x))
}

.requantifyTheoreticalDistribution <- function(x, sequence,
                                               saturationThreshold=Inf,
                                               method=c("mean", "median", "weighted.mean")) {
  x <- .splitIsotopicDistr(unlist(x))
  commonnm <- .commonIsotopes(x, saturationThreshold)

  if (length(commonnm)) {
    distr <- lapply(x, function(y).sumIsotopes(y[commonnm]))
    n <- max(sapply(distr, length))
    props <- BRAIN::calculateIsotopicProbabilities(BRAIN::getAtomsFromSeq(sequence),
                                                   nrPeaks=n)
    distr <- lapply(distr, function(y)y/props)

    return(switch(match.arg(method),
                  "mean" = sapply(distr, mean, na.rm=TRUE),
                  "median" = sapply(distr, median, na.rm=TRUE),
                  "weighted.mean" = sapply(distr, weighted.mean, w=props, na.rm=TRUE)))
  }

  rep.int(NA, length(x))
}

.commonIsotopes <- function(x, saturationThreshold=Inf) {
  distr <- lapply(x, function(y)y[y<saturationThreshold])
  nm <- lapply(distr, names)
  Reduce(intersect, nm)
}

.sumIsotopes <- function(x) {
  nm <- factor(sapply(strsplit(names(x), "_", fixed=TRUE), "[", 2))
  tapply(x, nm, sum)
}
