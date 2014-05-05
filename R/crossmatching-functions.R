#' @param spectra list, 4 MSnbase::Spectrum2 objects
#' @param fragments list, 2 character vectors containing the fragment.str
#' @param tolerance double, allowed deviation
#' @param ... passed to MSnbase:::.plotSingleSpectrum
#' @noRd
.plotSpectraVsFragments <- function(spectra, fragments, tolerance=25e-6, ...) {
  spectra <- lapply(spectra, normalize, method="precursor")

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
#' @param matchType character
#' @param fragments.cex cex for fragments
#' @param legend.cex cex for legend
#' @param ... passed to MSnbase:::.plotSingleSpectrum
#' @noRd
.plotSpectrumVsSpectrum <- function(spectra,
                                    sequences,
                                    common,
                                    fragments=vector(mode="list", length=2),
                                    matchType,
                                    xlim, ylim,
                                    legend.cex=1, ...) {
  if (missing(sequences)) {
    sequences <- character(2)
  }

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
                             ", match: ", matchType)
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
#' matchType,
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
  }

  return(flatEmrts)
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
                          matchType=cxrow$matchType,
                          legend.cex=legend.cex,
                          fragments.cex=fragments.cex,
                          ...)
}

#' create an equal sized sample of true and false matches
#' @param cx cross matching df, result of crossmatching
#' @return list with trueIdx and falseIdx
#' @noRd
.groundTruthIndices <- function(cx) {
  ## select "ground-truth"
  trueIdx <- grep("true", cx$matchType)
  falseIdx <- grep("false", cx$matchType)

  n <- min(c(length(trueIdx), length(falseIdx)))

  if (length(trueIdx) > n) {
    trueIdx <- sample(trueIdx, n)
  }

  if (length(falseIdx) > n) {
    falseIdx <- sample(falseIdx, n)
  }
  return(list(trueIdx=trueIdx, falseIdx=falseIdx))
}

#' @param cx cross matching df, result of cross matching
#' @param mcol column name of the matching results (e.g.
#' fragments.identXfragments.quant)
#' @return matrix with cols tp, fp, tn, fn
#' @noRd
.crossMatchingConfusionMatrix <- function(cx, mcol) {
  trueIdx <- grep("true", cx$matchType)
  falseIdx <- grep("false", cx$matchType)

  ytrain <- rep(c(TRUE, FALSE), c(length(trueIdx), length(falseIdx)))
  xtrain <- c(cx[trueIdx, mcol], cx[falseIdx, mcol])
  thresholds <- 0:max(xtrain, na.rm=TRUE)

  confusion <- t(sapply(thresholds, function(th) {
    tp <- sum(xtrain >= th & ytrain)
    fp <- sum(xtrain >= th & !ytrain)
    tn <- sum(xtrain < th & !ytrain)
    fn <- sum(xtrain < th & ytrain)
    return(c(tp=tp, fp=fp, tn=tn, fn=fn))
  }))
  confusion <- cbind(ncommon=thresholds, confusion)
  rownames(confusion) <- NULL

  return(confusion)
}

#' @param cx cross matching df, result of crossmatching
#' @param mcol column name of the matching results (e.g.
#' fragments.identXfragments.quant)
#' @return matrix with cols tp, fp, tn, fn
#' @noRd
.crossMatchingContingencyMatrix <- function(cx, mcol) {
  confusion <- .crossMatchingConfusionMatrix(cx, mcol)
  return(cbind(confusion, diagnosticErrors(confusion)))
}

#' plot cross matching F1/FDR/Confusion/Boxplots
#' @param cx cross matching df, result of crossmatching
#' @param mcol column name of the matching results (e.g.
#' fragments.identXfragments.quant)
#' @return invisible matrix with rows tp, fp, tn, fn, accuracy, precision,
#' recall, fdr, f1
#' @noRd
.plotCrossMatchingSummary <- function(cx,
                                      mcol="spectrum.quantXfragments.ident") {
  sampled <- .groundTruthIndices(cx)
  contengency <- .crossMatchingContingencyMatrix(cx[unlist(sampled), ], mcol)

  par(mfcol=c(1, 3))
  x <- (1:nrow(contengency))-1
  matplot(x, contengency[, c("fdr", "f1"), drop=FALSE], type="b", lty=1,
          xlab="# of common peaks", ylab="performance",
          main="cross matching performance", pch=19)
  legend("topright", legend=c("FDR", "F1"),
         col=1:2, lwd=1, pch=19, bty="n")
  grid()
  matplot(x, contengency[, c("tp", "fp", "tn", "fn"), drop=FALSE],
          type="b", lty=1, pch=19,
          xlab="# of common peaks", ylab="# of peptides",
          main="cross matching confusion")
  grid()
  legend("right", legend=c("TP", "FP", "TN", "FN"),
         col=1:4, lwd=1, pch=19, bty="n")

  trueIdx <- grep("true", cx$matchType)
  falseIdx <- grep("false", cx$matchType)

  l <- list()

  l[["all true matches"]] <- cx[trueIdx, mcol]

  if (length(trueIdx) > length(sampled$trueIdx)) {
    l[["sampled true matches"]] <- cx[sampled$trueIdx, mcol]
  }

  l[["all false matches"]] <- cx[falseIdx, mcol]

  if (length(falseIdx) > length(sampled$falseIdx)) {
    l[["sampled false matches"]] <- cx[sampled$falseIdx, mcol]
  }

  jitteredBoxplot(l, jitter.factor=1,
                  main="cross matching distribution", ylab="# of common peaks")

  par(mfcol=c(1, 1))

  invisible(contengency)
}

#' filter non-unique matches
#' @param obj synapter object
#' @param nmin a correct match must have at least "nmin" common peaks
#' @param mcol column name of the matching results (e.g.
#' fragments.identXfragments.quant)
#' @return invisible data.frame, columns: diff, highest
#' @noRd
.filterMatchedEMRTsUsingCrossMatching <- function(obj, nmin,
                                                  mcol="spectrum.quantXfragments.ident") {
  cx <- obj$CrossMatching
  emrts <- obj$MatchedEMRTs
  pep3d <- obj$QuantPep3DData

  if ("Function.0" %in% colnames(emrts)) {
    warning("You already filtered using the CrossMatching results! ",
            "The results may be wrong!", immediate. = TRUE)
  }

  ## backup old Function column
  emrts$Function.0 <- emrts$Function

  cx <- cx[cx$Function > 1 & cx[, mcol] >= nmin, ]

  ## if we find duplicates (multiple matches even after applying the cross
  ## matching rule) don't transfer them into unique ones.
  nonDuplicated <- !(duplicated(cx$precursor.leID.ident) |
                     rev(duplicated(rev(cx$precursor.leID.ident))))
  cx <- cx[nonDuplicated, ]

  if (nrow(cx)) {
   rowsInEMRTs <- match(cx$precursor.leID.ident, emrts$precursor.leID.ident)
   rowsInPep3D <- match(cx$matched.quant.spectrumIDs, pep3d$spectrumID)

   columns <- intersect(colnames(cx), colnames(pep3d))

   cx[, c("precursor.leID.quant", columns)] <-
     pep3d[rowsInPep3D, c("spectrumID", columns)]
   cx$matchType <- "unique-true"
   cx$idSource <- "crossmatching"
   cx$Function <- 1
   emrts[rowsInEMRTs, columns] <- cx[, columns]
   emrts[rowsInEMRTs, c("idSource", "precursor.leID.quant", "spectrumID")] <-
     cx[, c("idSource", "precursor.leID.quant", "spectrumID")]
  }
  return(emrts)
}

