#' @param spectra list, 2 MSnbase::Spectrum2 objects
#' @param sequences list, 2 character vectors containing the peptide sequences
#' for spectra and fragments
#' @param fragments, characters vector containing the fragment.str
#' @param tolerance double, allowed deviation
#' @param ... passed to MSnbase:::.plotSingleSpectrum
#' @noRd
.plotFragmentsVsSpectrum <- function(spectra, sequences, fragments,
                                    tolerance=25e-6, ...) {
  spectra <- lapply(spectra, normalize, method="precursor")

  if (missing(sequences)) {
    sequences <- character(2)
  }

  mass <- unlist(lapply(spectra, mz))

  ## if both spectra are removed because of NeutralLoss=TRUE
  if (length(mass)) {
    xlim <- c(min(mass, na.rm=TRUE), max(mass, na.rm=TRUE))
  } else {
    xlim <- c(0, 0)
  }

  inten <- unlist(lapply(spectra, intensity))
  maxInten <- max(c(0, inten), na.rm=TRUE)
  ylim <- c(-maxInten, maxInten)

  oldPar <- par(no.readonly=TRUE)
  on.exit(par(oldPar))

  par(mar=c(2, 2, 2, 0.5)) #c(bottom, left, top, right)
  .plotSpectrumVsSpectrum(
    spectra,
    sequences=unlist(sequences),
    common=list(MSnbase:::commonPeaks(spectra[[1]],
                                      spectra[[2]],
                                      tolerance=tolerance),
                MSnbase:::commonPeaks(spectra[[2]],
                                      spectra[[1]],
                                      tolerance=tolerance)),
    fragments=list(data.frame(mz=mz(spectra[[1]]), ion=fragments,
                              stringsAsFactors=FALSE),
                   data.frame()),
    main="Identification Fragments vs Quantitation Spectrum",
    xlim=xlim, ylim=ylim, ...)
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
#' compares ident fragments vs quant product spectrum
#' @param flatEmrts flattened EMRTs (see flatMatchedEMRTs)
#' @param identFragments MSnExp of the identification fragments
#' @param quantSpectra MSnExp of the quantitation spectra
#' @param tolerance double, allowed deviation to consider a m/z as equal
#' @param verbose verbose output?
#' @return data.frame, extend/flatted matchedEmrts df with additional columns:
#' gridSearchResult, CrossMatching
#' @noRd
crossmatching <- function(flatEmrts, identFragments, quantSpectra,
                          tolerance=25e-6, verbose=TRUE) {
  if (verbose) {
    message("Look for common peaks")
    pb <- txtProgressBar(0, nrow(flatEmrts), style=3)
  }

  for (i in 1:nrow(flatEmrts)) {
    flatEmrts[i, "CrossMatching"] <-
      MSnbase:::numberOfCommonPeaks(
        .getSpectrum(flatEmrts$precursor.leID.ident[i], identFragments),
        .getSpectrum(flatEmrts$matched.quant.spectrumIDs[i], quantSpectra),
        tolerance=tolerance)

    if (verbose) {
      setTxtProgressBar(pb, i)
    }
  }

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
#' @return modified cx data.frame
#' @noRd
.crossMatchingDifferences <- function(cx) {
  # use only non-unique matches
  idx <- grep("^non-unique", cx$gridSearchResult)

  cx$CrossMatchingRank <- cx$CrossMatchingDiff <- NA

  cx$CrossMatchingRank[idx] <- ave(cx$CrossMatching[idx],
                                   cx$precursor.leID.ident[idx],
                                   FUN=function(x)rank(-x))

  cx$CrossMatchingDiff[idx] <- ave(cx$CrossMatching[idx],
                                   cx$precursor.leID.ident[idx],
                                   FUN=function(x)diff(sort(x, decreasing=TRUE)[2:1]))

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
               identFragments=obj$IdentFragmentData,
               quantSpectra=obj$QuantSpectrumData,
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
#' @param identFragments MSnExp of the identification fragments
#' @param quantSpectra MSnExp of the quantitation spectra
#' @param legend.cex cex for legend text
#' @param fragments.cex cex for fragment text
#' @noRd
.plotCxRow <- function(cxrow, identFragments, quantSpectra,
                       legend.cex=0.6, fragments.cex=0.5,
                       ...) {

  keys <- unlist(cxrow[,c("precursor.leID.ident", "matched.quant.spectrumIDs")])
  spectra <- list(ident=.getSpectrum(keys[1], identFragments),
                  quant=.getSpectrum(keys[2], quantSpectra))

  sequences <- mapply(function(x, k) {
    ## avoid partial matching of rows
    ## ( [.data.frame always use partial match for rows; exact=TRUE works only
    ## for column names)
    return(fData(x)[match(k, rownames(fData(x))), "peptide.seq"])
  }, x=list(identFragments, quantSpectra), k=keys,
  SIMPLIFY=FALSE, USE.NAMES=FALSE)

  fragments <- na.omit(MSnbase:::utils.ssv2vec(
    fData(identFragments)[match(keys[1], rownames(fData(identFragments))), "fragment.str"]))

  .plotFragmentsVsSpectrum(spectra=spectra,
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
.plotCrossMatchingPerformance <- function(cx) {
  l <- setNames(vector(mode="list", length=2), c("unique", "nonunique"))
  what <- c("unique", "non-unique")
  xlab <- c("# of common peaks", "delta common peaks")

  par(mfcol=c(1, 2))
  for (i in seq(along=l)) {
    l[[i]] <- .crossMatchingContingencyMatrix(cx, what=what[i])

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
#' @return matrix with cols threshold,
#' @noRd
.crossMatchingConfusionMatrix <- function(cx, what=c("unique", "non-unique")) {
  what <- match.arg(what)

  cx <- cx[grep(paste0("^", what), cx$gridSearchResult), ]

  rcol <- "CrossMatchingRank"

  if (what == "unique") {
    cx[, rcol] <- 1
    mcol <- "CrossMatching"
  } else {
    mcol <- "CrossMatchingDiff"
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
.crossMatchingContingencyMatrix <- function(cx, what) {
  confusion <- .crossMatchingConfusionMatrix(cx, what)
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

