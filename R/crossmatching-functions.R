#' @param spectra list, 4 MSnbase::Spectrum2 objects
#' @param sequences list, 4 character vectors containing the peptide sequences
#' for spectra and fragments
#' @param fragments list, 2 character vectors containing the fragment.str
#' @param tolerance double, allowed deviation
#' @param ... passed to MSnbase:::.plotSingleSpectrum
#' @noRd
.plotSpectraVsFragments <- function(spectra, sequences, fragments,
                                    tolerance=25e-6, ...) {
  spectra <- lapply(spectra, normalize, method="precursor")

  if (missing(sequences)) {
    sequences <- character(4)
  }

  mass <- unlist(lapply(spectra, mz))
  xlim <- c(min(mass, na.rm=TRUE), max(mass, na.rm=TRUE))

  inten <- unlist(lapply(spectra, intensity))
  maxInten <- max(c(0, inten), na.rm=TRUE)
  ylim <- c(-maxInten, maxInten)

  oldPar <- par(no.readonly=TRUE)
  on.exit(par(oldPar))
  par(mfrow=c(1, 2))

  par(mar=c(2, 2, 2, 0.5)) #c(bottom, left, top, right)
  .plotSpectrumVsSpectrum(
    spectra[1:2],
    sequences=unlist(sequences[1:2]),
    common=list(MSnbase:::commonPeaks(spectra[[1]],
                                      spectra[[4]],
                                      tolerance=tolerance),
                MSnbase:::commonPeaks(spectra[[2]],
                                      spectra[[3]],
                                      tolerance=tolerance)),
    fragments=list(data.frame(), data.frame()),
    main="spectra", xlim=xlim, ylim=ylim, ...)
  par(mar=c(2, 0.5, 2, 2)) #c(bottom, left, top, right)
  .plotSpectrumVsSpectrum(
    spectra[3:4],
    sequences=unlist(sequences[3:4]),
    common=list(MSnbase:::commonPeaks(spectra[[3]],
                                      spectra[[4]],
                                      tolerance=tolerance),
                MSnbase:::commonPeaks(spectra[[4]],
                                      spectra[[3]],
                                      tolerance=tolerance)),
    fragments=list(data.frame(mz=mz(spectra[[3]]), ion=fragments[[1]],
                              stringsAsFactors=FALSE),
                   data.frame(mz=mz(spectra[[4]]), ion=fragments[[2]],
                              stringsAsFactors=FALSE)),
    main="fragments", xlim=xlim, ylim=ylim, yaxt="n", ...)
  axis(4)
}

#' plot spectrum1 vs spectrum2
#' @param spectra list, 2 MSnbase::Spectrum2 objects
#' @param sequences character vector (length==2) containing the peptide
#' sequences for both spectra
#' @param common list (length==2), containing logical vector for common peaks
#' @param fragments list (length==2), containing fragments data.frames
#' @param gridSearchResult character
#' @param fragments.cex cex for fragments
#' @param legend.cex cex for legend
#' @param ... passed to MSnbase:::.plotSingleSpectrum
#' @noRd
.plotSpectrumVsSpectrum <- function(spectra,
                                    sequences,
                                    common,
                                    fragments=vector(mode="list", length=2),
                                    gridSearchResult,
                                    xlim, ylim,
                                    legend.cex=1, ...) {
  orientation <- c(1, -1)
  add <- c(FALSE, TRUE)
  legend.pos <- c("topleft", "bottomleft")
  legend.prefix <- c("ident", "quant")
  ## colors: ColorBrewer RdYlBu c(9, 11, 3, 1)
  cols <- c("#74ADD1", "#313695", "#F46D43", "#A50026")
  pch <- c(NA, 19)

  for (i in seq(along=spectra)) {
    if (!length(common[[i]])) {
      common[[i]] <- logical(peaksCount(spectra[[i]]))
    }
    MSnbase:::.plotSingleSpectrum(spectra[[i]],
                                  sequence=sequences[[i]],
                                  fragments=fragments[[i]],
                                  orientation=orientation[i],
                                  add=add[i], xlim=xlim, ylim=ylim,
                                  col=cols[(i-1)*2+common[[i]]+1],
                                  pch=pch[common[[i]]+1], ...)

    label <- paste0(legend.prefix[i], " prec leID: ",
                    precScanNum(spectra[[i]]))

    if (peaksCount(spectra[[i]])) {
      label <- paste0(label, ", prec mass: ", round(precursorMz(spectra[[i]]), 3),
                             ", prec z: ", precursorCharge(spectra[[i]]),
                             ", # common: ", sum(common[[i]]),
                             ", match: ", gridSearchResult)
      if (nchar(sequences[[i]])) {
        label <- paste0(label, ",\nseq: ", sequences[[i]])
      }
    }

    legend(legend.pos[i], legend=label, bty="n", cex=legend.cex)
  }
}

#' crossmatching, workhorse function
#' compares ident fragments vs quant product spectrum and
#' quant fragments vs ident product spectrum
#' @param flatEmrts flattened EMRTs (see flatMatchedEMRTs)
#' @param spectra list of 4 MSnExp containing the spectra and fragment spectra
#' @param tolerance double, allowed deviation to consider a m/z as equal
#' @param verbose verbose output?
#' @return data.frame, extend/flatted matchedEmrts df with additional columns:
#' gridSearchResult,
#' spectrum.identXfragments.ident, spectrum.quantXfragments.quant,
#' spectrum.identXfragments.quant, spectrum.quantXfragments.ident
#' sorry for the names
#' @noRd
crossmatching <- function(flatEmrts, spectra, tolerance=25e-6, verbose=TRUE) {
  prefixes <- paste(rep(c("spectrum", "fragments"), each=2),
                    rep(c("ident", "quant"), times=2), sep=".")
  # "spectrum.ident"   "spectrum.quant"   "fragments.ident" "fragments.quant"

  keys <- as.character(rep(unlist(flatEmrts[, c("precursor.leID.ident",
                                                "matched.quant.spectrumIDs")]),
                           times=2))
  keysm <- matrix(keys, ncol=4)

  #cmb <- list(c(1, 3), c(2, 4),
  #            c(1, 4), c(2, 3))
  cmb <- list(c(1, 4), c(2, 3), c(3, 4))

  cols <- sapply(cmb, function(x)paste0(prefixes[x], collapse="X"))
  # "spectrum.identXfragments.ident" ...

  if (verbose) {
    message("Look for common peaks")
    pb <- txtProgressBar(0, length(cmb)*nrow(flatEmrts), style=3)
  }

  flatEmrts[, cols] <- lapply(cmb, function(i) {
    apply(keysm[, i], 1, function(k) {
      if (verbose) {
        setTxtProgressBar(pb, pb$getVal()+1)
      }
      MSnbase:::numberOfCommonPeaks(
        .getSpectrum(k[1], spectra[[i[1]]]),
        .getSpectrum(k[2], spectra[[i[2]]]),
        tolerance=tolerance)
    })
  })

  if (verbose) {
    close(pb)
    message("Calculate Differences")
  }

  flatEmrts <- .crossMatchingDifferences(flatEmrts)

  return(flatEmrts)
}

#' calculate difference between first highest and second highest number of
#' common peaks in a non-unique match group
#' @param cx data.frame
#' @param mcol column name of the matching results (e.g.
#' fragments.identXfragments.quant)
#' @return modified cx data.frame
#' @noRd
.crossMatchingDifferences <- function(cx,
                                      mcol="spectrum.quantXfragments.ident") {
  # use only non-unique matches
  idx <- grep("^non-unique", cx$gridSearchResult)

  dcol <- paste0(mcol, ".diff")
  rcol <- paste0(mcol, ".rank")
  cx[, dcol] <- NA
  cx[, rcol] <- NA
  ## we use the second argument (precursor.leID.quant) is just needed to allow a
  ## return value with 2 columns
  cx[idx, c(dcol, rcol)] <- ave(cx[idx, c("precursor.leID.quant", mcol)],
                                cx$precursor.leID.ident[idx],
                                FUN=function(x) {
    sorted <- sort(x[, 2], decreasing=TRUE)
    r <- rank(-x[,2])
    d <- diff(sorted[2:1])
    cbind(d, r)
  })

  return(cx)
}

#' plot crossmatching
#' @param obj synapter object
#' @param key all rows with a matching "key" in "column" are plotted.
#' If missing everything is plotted.
#' @param column in which column we should look for "key"
#' @param verbose verbose output?
#' @noRd
.plotCrossMatching <- function(obj, key, column="peptide.seq", verbose=TRUE,
                              ...) {
  cx <- obj$CrossMatching

  if (!column %in% colnames(cx)) {
    stop(sQuote("column"), "=", column, " is not a valid column of the ",
         "cross matching data.frame!")
  }

  if (missing(key) && verbose) {
    message(sQuote("key"), " is missing. Plotting everything.")
  }

  if (!missing(key)) {
    cx <- cx[which(cx[, column] %in% key), ]

    if (!nrow(cx)) {
      stop("No row matches your criteria!")
    }
  }

  if (verbose) {
    pb <- txtProgressBar(0, nrow(cx), style=3)
  }

  for (i in 1:nrow(cx)) {
    .plotCxRow(cxrow=cx[i, ],
               spectra=list(obj$IdentSpectrumData,
                            obj$QuantSpectrumData,
                            obj$IdentFragmentData,
                            obj$QuantFragmentData),
               ...)

   if (verbose) {
      setTxtProgressBar(pb, i)
    }
  }

  if (verbose) {
    close(pb)
  }
}

#' @param cxrow a single row of the crossmatching df
#' @param spectra list of 4 MSnExp containing the spectra and fragment spectra
#' @param legend.cex cex for legend text
#' @param fragments.cex cex for fragment text
#' @noRd
.plotCxRow <- function(cxrow, spectra,
                       legend.cex=0.6, fragments.cex=0.5,
                       ...) {
  keys <- as.character(rep(c(cxrow$precursor.leID.ident,
                             cxrow$matched.quant.spectrumIDs), times=2))
  sequences <- mapply(function(x, k) {
    ## avoid partial matching of rows
    ## ( [.data.frame always use partial match for rows; exact=TRUE works only
    ## for column names)
    return(fData(x)[match(k, rownames(fData(x))), "peptide.seq"])
  }, x=spectra, k=keys, SIMPLIFY=FALSE, USE.NAMES=FALSE)

  fragments <- mapply(function(x, k) {
    ## avoid partial matching (see above)
    fragment.str <- fData(x)[match(k, rownames(fData(x))), "fragment.str"]
    return(na.omit(MSnbase:::utils.ssv2vec(fragment.str)))
  }, x=spectra[3:4], k=keys[3:4], SIMPLIFY=FALSE, USE.NAMES=FALSE)

  .plotSpectraVsFragments(spectra=.getSpectra(keys, spectralist=spectra),
                          sequences=sequences,
                          fragments=fragments,
                          gridSearchResult=cxrow$gridSearchResult,
                          legend.cex=legend.cex,
                          fragments.cex=fragments.cex,
                          ...)
}

#' plot cx performance
#' @param cx cx data.frame
#' @return invisible list with two elements (unique,nonunique) each containing a
#' matrix of 3 columns threshold, true, false
.plotCrossMatchingPerformance <- function(cx,
                                          mcol="spectrum.quantXfragments.ident") {
  l <- setNames(vector(mode="list", length=2), c("unique", "nonunique"))
  what <- c("unique", "non-unique")
  xlab <- c("# of common peaks", "delta common peaks")

  par(mfcol=c(1, 2))
  for (i in seq(along=l)) {
    l[[i]] <- .crossMatchingContingencyMatrix(cx, what=what[i], mcol=mcol)

    ylim <- range(l[[i]][, c("tp", "fp")])
    plot(l[[i]][, 1], l[[i]][, "tp"],
         type="b", lty=1, pch=ifelse(l[[i]][, "tp"], 19, 1),
         ylim=ylim, col=4, xlab=xlab[i], ylab="# of peptides",
         main=paste("performance", what[i]))
    lines(l[[i]][, 1], l[[i]][, "fp"], type="b",
          pch=ifelse(l[[i]][, "fp"], 19, 1), col=2)
    grid()
    legend("topright", legend=c("true matches", "false matches"),
           col=c(4, 2), lwd=1, pch=19, bty="n")
  }
  par(mfcol=c(1, 1))

  invisible(l)
}

#' @param cx cross matching df, result of cross matching
#' @param what character unique/non-unique
#' @param mcol column name of the matching results (e.g.
#' fragments.identXfragments.quant)
#' @return matrix with cols threshold,
#' @noRd
.crossMatchingConfusionMatrix <- function(cx, what=c("unique", "non-unique"),
                                          mcol) {
  what <- match.arg(what)

  cx <- cx[grep(paste0("^", what), cx$gridSearchResult), ]

  rcol <- paste0(mcol, ".rank")

  if (what == "unique") {
    cx[, rcol] <- 1
  } else {
    mcol <- paste0(mcol, ".diff")
  }

  train <- tapply(1:nrow(cx), cx$precursor.leID.ident, function(i) {
    ytrain <- any(grepl("true$", cx$gridSearchResult[i]) & cx[i, rcol] == 1)
    xtrain <- cx[i[which.min(cx[i, rcol])], mcol]
    cbind(xtrain, ytrain)
  }, simplify=FALSE)
  train <- do.call(rbind, train)
  xtrain <- train[, 1]
  ytrain <- train[, 2]

  thresholds <- 0:max(xtrain, na.rm=TRUE)

  confusion <- t(sapply(thresholds, function(th) {
    tp <- sum(xtrain >= th & ytrain)
    fp <- sum(xtrain >= th & !ytrain)
    tn <- sum(xtrain < th & !ytrain)
    fn <- sum(xtrain < th & ytrain)
    return(c(tp=tp, fp=fp, tn=tn, fn=fn))
  }))
  confusion <- cbind(thresholds, confusion)

  colnames(confusion)[1] <- ifelse(what == "unique", "ncommon", "deltacommon")
  rownames(confusion) <- NULL

  return(confusion)
}

#' @param cx cross matching df, result of crossmatching
#' @param what character unique/non-unique
#' @param mcol column name of the matching results (e.g.
#' fragments.identXfragments.quant)
#' @return matrix with cols tp, fp, tn, fn
#' @noRd
.crossMatchingContingencyMatrix <- function(cx, what, mcol) {
  confusion <- .crossMatchingConfusionMatrix(cx, what, mcol)
  return(cbind(confusion, fdr=diagnosticErrors(confusion)[, "fdr"]))
}

#' filter unique matches
#' @param obj synapter object
#' @param mincommon a correct match must have at least "mincommon" common peaks
#' @param mcol column name of the matching results (e.g.
#' fragments.identXfragments.quant)
#' @return filtered emrts data.frame
#' @noRd
.filterUniqueMatches <- function(obj, mincommon,
                                 mcol="spectrum.quantXfragments.ident") {
  emrts <- obj$MatchedEMRTs
  cx <- obj$CrossMatching

  if ("Function.1" %in% colnames(emrts)) {
    warning("This function should not used multiple times. It removes your ",
            "original ", sQuote("Function"), " column.", immediate.=TRUE)
  }

  cx <- cx[cx$precursor.leID.ident %in% emrts$precursor.leID.ident &
           cx$Function == 1, ]

  keep <- cx[, mcol] >= mincommon

  if (!any(keep)) {
    stop("No EMRT match your criteria! Try a lower threshold")
  }

  rows <- match(cx$precursor.leID.ident, emrts$precursor.leID.ident)

  exclude <- rows[!keep]

  ## backup function
  emrts$Function.1 <- emrts$Function
  emrts$Function[rows] <- 1

  if (length(exclude)) {
    # comment the next line and uncomment the over next line to *not* remove the
    # filtered emrts
    #emrts <- emrts[-exclude, ]
    emrts$Function[exclude] <-  -1
  }
  return(emrts)
}

#' filter nonunique matches
#' @param obj synapter object
#' @param mindelta a correct match must have at least "nmin" common peaks
#' @param mcol column name of the matching results (e.g.
#' fragments.identXfragments.quant)
#' @return filtered emrts data.frame
#' @noRd
.filterNonUniqueMatches <- function(obj, mindelta,
                                    mcol="spectrum.quantXfragments.ident") {
  emrts <- obj$MatchedEMRTs
  cx <- obj$CrossMatching

  if ("Function.2" %in% colnames(emrts)) {
    warning("This function should not used multiple times. It removes your ",
            "original ", sQuote("Function"), " column.", immediate.=TRUE)
  }

  rcol <- paste0(mcol, ".rank")
  mcol <- paste0(mcol, ".diff")

  cx <- cx[cx$precursor.leID.ident %in% emrts$precursor.leID.ident &
           cx$Function > 1, ]

  keep <- cx[, mcol] >= mindelta & cx[, rcol] == 1

  if (!any(keep)) {
    stop("No EMRT match your criteria! Try a lower threshold")
  }

  rows <- match(cx$precursor.leID.ident, emrts$precursor.leID.ident)

  exclude <- rows[!keep]

  ## backup function
  emrts$Function.2 <- emrts$Function

  if (length(exclude)) {
    cols <- c("matched.quant.spectrumIDs", "precursor.leID.quant", "spectrumID")
    emrts[rows, cols] <- cx[, cols]
    # comment the next line and uncomment the over next line to *not* remove the
    # filtered emrts
    #emrts <- emrts[-exclude, ]

    ## this order is important because include could contain row numbers that
    ## are present in keep as well (because cx has many duplicated
    ## precursor.leID.idents)
    emrts$Function[exclude] <- -2
    emrts$Function[rows[keep]] <- 1
  }
  return(emrts)
}

