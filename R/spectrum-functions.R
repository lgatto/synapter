#' parse synapter Spectrum.xml files
#'
#' Stupid parser for one of the ugliest pseudo-xml formats I have ever seen.
#'
#' @param file file path
#' @param ms1 should we read ms1 spectra?
#' @param encoding file encoding, seems to be windows specific
#' @param verbose verbose output?
#'
#' @return a list of length == 3; names: ms1, m2, assignments
#'  ms1: matrix of ms1 spectra information
#'  ms2: matrix of ms2 spectra information
#'  assignments: environment (key == leID), values == he_id
#' @noRd
.readSynapterSpectrumXml <- function(file, ms1=FALSE,
                                     encoding="Windows-1252",
                                     verbose=TRUE) {
  stopifnot(file.exists(file))

  ## helper functions
  .extractName <- function(x) {
    gsub("^.*NAME=\\\"([[:alpha:]_]+)\\\".*$", "\\1", x)
  }
  .extractLEID <- function(x) {
    gsub("^.*LE_ID=\\\"([[:digit:]]+)\\\".*$", "\\1", x)
  }
  .extractHEID <- function(x) {
    gsub("^.*HE_ID=\\\"([[:digit:],]+)\\\".*$", "\\1", x)
  }

  .createMsMatrix <- function(x, header) {
    ms <- do.call(rbind, strsplit(x, "[[:space:]]+"))
    ms <- ms[, -1]
    mode(ms) <- "double"
    colnames(ms) <- header
    return(ms)
  }

  if (verbose) {
    message("Reading ", file)
  }
  content <- readLines(file, encoding=encoding, warn=FALSE)

  if (verbose) {
    message("Search line numbers for Header information")
  }
  linesField <- grep("<FIELD", content, fixed=TRUE)
  splitPoint <- which(diff(linesField) > 1)
  linesFieldMs1 <- linesField[1:splitPoint]
  linesFieldMs2 <- linesField[(splitPoint+1):length(linesField)]

  if (verbose) {
    message("Read Header information")
  }

  header1 <- .extractName(content[linesFieldMs1])
  header2 <- .extractName(content[linesFieldMs2])

  if (ms1) {
    if (verbose) {
      message("Search line numbers for MS1 Data")
    }
    linesMs1 <- c(grep("<DATA", content, fixed=TRUE)+1,
                  grep("</DATA", content, fixed=TRUE)-1)
  }

  if (verbose) {
    message("Search line numbers for MS2 Data")
  }
  linesMs2 <- c(grep("<HE_DATA", content, fixed=TRUE)+1,
                grep("</HE_DATA", content, fixed=TRUE)-1)

  if (ms1) {
    if (verbose) {
      message("Read MS1 Data (", paste(linesMs1, collapse=":"), ")")
    }
    range <- seq(linesMs1[1], linesMs1[2])
    ms1 <- .createMsMatrix(content[range], header1)
  } else {
    ms1 <- matrix()
  }

  if (verbose) {
    message("Read MS2 Data (", paste(linesMs2, collapse=":"), ")")
  }
  range <- seq(linesMs2[1], linesMs2[2])
  ms2 <- .createMsMatrix(content[range], header2)

  if (verbose) {
    message("Search line numbers for MS1 to MS2 assignment")
  }
  linesAssignment <- c(grep("<PRECURSOR_PRODUCT_BIN", content, fixed=TRUE)+1,
                       grep("</PRECURSOR_PRODUCT_BIN", content, fixed=TRUE)-1)

  if (verbose) {
    message("Read Assignment Data (", paste(linesMs2, collapse=":"), ")")
  }
  range <- seq(linesAssignment[1], linesAssignment[2])

  le_ids <- .extractLEID(content[range])
  he_ids <- lapply(MSnbase:::utils.ssv2list(
                     .extractHEID(content[range]), sep=","), as.numeric)

  assignments <- new.env(hash=TRUE, parent=emptyenv(), size=length(le_ids))
  for (i in seq(along=le_ids)) {
    assign(le_ids[i], he_ids[[i]], envir=assignments)
  }

  return(list(ms1=ms1, ms2=ms2, assignments=assignments))
}

#' read spectrum.xml and turn data into MSnbase::Spectrum2 objects
#' @param df corresponding df from the synapter object ({Ident,Quant}PeptideData)
#' @param file filename
#' @param assaydata environment
#' @param storeAll should all spectra stored? or only the needed ones?
#' @param verbose verbose output
#' @return modified assaydata
#' @noRd
.spectrumXml2spectra <- function(df, file, storeAll=TRUE, verbose=TRUE) {
  xml <- .readSynapterSpectrumXml(file, ms1=TRUE, verbose=verbose)

  peptideinfo <- df[!duplicated(df$precursor.leID), ]

  if (storeAll) {
    uleID <- union(df$precursor.leID, ls(xml$assignments))
  } else {
    uleID <- intersect(df$precursor.leID, ls(xml$assignments))
  }

  if (verbose) {
    message("Convert spectra data.frames to MSnbase::Spectrum2 objects")
    pb <- txtProgressBar(0, length(uleID), style=3)
  }

  sequences <- rep(NA, length(uleID))

  assaydata <- new.env(hash=TRUE, parent=emptyenv(), size=length(uleID))

  for (i in seq(along=uleID)) {
    assign(as.character(uleID[i]),
           .createMs2SpectrumFromSpectrumXml(uleID[i], xml$ms1, xml$ms2,
                                             xml$assignments),
           envir=assaydata)
    if (i <= nrow(df)) {
      ## in the xml files are no sequence information that's why we can
      ## only use the information from the df
      sequences[i] <- peptideinfo$peptide.seq[i]
    }

    if (verbose) {
      setTxtProgressBar(pb, i)
    }
  }

  if (verbose) {
    close(pb)
    message("Create MSnExp object")
  }

  process <- new("MSnProcess",
                 processing=paste("Data loaded:", date()),
                 files=file)

  nm <- ls(assaydata)
  pdata <- new("NAnnotatedDataFrame", data=data.frame(sampleNames=1L,
                                                      type="spectrum"))
  fdata <- data.frame(spectrum=1:length(nm),
                      leID=uleID,
                      peptide.seq=sequences,
                      row.names=uleID,
                      stringsAsFactors=FALSE)
  ## reorder according to assaydata
  fdata <- fdata[nm, ]

  fdata <- new("AnnotatedDataFrame", data=fdata)

  msnexp <- new("MSnExp", assayData=assaydata,
                          phenoData=pdata,
                          featureData=fdata,
                          processingData=process)
  return(msnexp)
}

#' create MS2 spectrum from spectrum.xml data
#' @param leID precursor.leID
#' @param ms1 ms1 matrix generated by .readSynapterSpectrumXml
#' @param ms2 ms2 matrix generated by .readSynapterSpectrumXml
#' @param assignments env generated by .readSynapterSpectrumXml
#' @return Spectrum2 object
#' @noRd
.createMs2SpectrumFromSpectrumXml <- function(leID, ms1, ms2, assignments) {
  key <- as.character(leID)

  if (exists(key, envir=assignments)) {
    i <- get(key, envir=assignments)
    leID <- as.integer(leID)

    return(.createMsnbaseSpectrum2(leID=leID,
                                   precursor.mhp=ms1[leID, "Mass"],
                                   precursor.inten=ms1[leID, "Intensity"],
                                   precursor.z=ms1[leID, "Z"],
                                   precursor.retT=ms1[leID, "RT"],
                                   mass=ms2[i, "Mass"],
                                   intensity=ms2[i, "Intensity"]))
  } else {
    return(.createEmptyMsnbaseSpectrum2(leID=leID))
  }
}

#' create MS2 spectrum from final_fragment.csv data
#' @param leID precursor.leID
#' @param fragments data.frame generated by readFragments
#' @param assignments env
#' @return Spectrum2 object
#' @noRd
.createMs2SpectrumFromFragments <- function(leID, fragments, assignments) {
  key <- as.character(leID)

  if (exists(key, envir=assignments)) {
    i <- get(key, envir=assignments)
    j <- i[1]

    return(.createMsnbaseSpectrum2(leID,
                                   precursor.mhp=fragments$precursor.mhp[j],
                                   precursor.inten=fragments$precursor.inten[j],
                                   precursor.z=fragments$precursor.z[j],
                                   precursor.retT=fragments$precursor.retT[j],
                                   mass=fragments$product.mhp[i],
                                   intensity=fragments$product.inten[i]))
  } else {
    return(.createEmptyMsnbaseSpectrum2(leID=leID))
  }
}

#' create an instance of an MSnbase::Spectrum2 object
#' @param leID precursor.leID
#' @param precursor.mhp precursor mass
#' @param precursor.inten precursor intensity
#' @param precursor.z precursor charge
#' @param precursor.retT precursor retention time
#' @param mass mass/mz
#' @param intensity intensity
#' @return Spectrum2 object
#' @noRd
.createMsnbaseSpectrum2 <- function(leID, precursor.mhp, precursor.inten,
                                    precursor.z, precursor.retT,
                                    mass, intensity) {
  ## just to be sure, you can't trust any information in the csv/xml files
  o <- order(mass)
  new("Spectrum2",
      precScanNum=as.integer(leID),
      precursorMz=precursor.mhp,
      precursorIntensity=precursor.inten,
      precursorCharge=as.integer(precursor.z),
      rt=precursor.retT,
      centroided=TRUE,
      # tic=precursor.inten,
      peaksCount=length(mass),
      mz=mass[o], intensity=intensity[o],
      fromFile=1L)
}

#' create an empty instance of an MSnbase::Spectrum2 object
#' @param leID precursor.leID
#' @param key character (used in enviroments to access the correct spetrum)
#' @noRd
.createEmptyMsnbaseSpectrum2 <- function(leID, key="spectrum:-1") {
  if (missing(leID)) {
    leID <- as.integer(tail(strsplit(key, ":")[[1]], 1))
  }
  new("Spectrum2", precScanNum=as.integer(leID))
}

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
    cx <- cx[which(key == cx[, column]), ]

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

#' plot crossmatching FDR
#' @param cx crossmatching df, result of crossmatching
#' @param tolerance double, allowed deviation
#' @param verbose verbose output?
#' @return invisible matrix with rows TP, FP, TN, FN, FDR
#' @noRd
plotCrossmatchingFDR <- function(cx, tolerance=25e-6, verbose=TRUE) {

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


  ytrain <- rep(c(TRUE, FALSE), each=n)
  xtrain <- c(cx$fragments.identXfragments.quant[trueIdx],
              cx$fragments.identXfragments.quant[falseIdx])
  thresholds <- 0:max(xtrain, na.rm=TRUE)

  contingencies <- sapply(thresholds, function(th) {
    tp <- sum(xtrain > th & ytrain)
    fp <- sum(xtrain > th & !ytrain)
    tn <- sum(xtrain < th & !ytrain)
    fn <- sum(xtrain < th & ytrain)
    fdr <- fp/(tp+fp)
    return(c(tp=tp, fp=fp, tn=tn, fn=fn, fdr=fdr))
  })
  colnames(contingencies) <- thresholds

  par(mfcol=c(1, 2))
  plot(thresholds, contingencies[5, ], type="b",
       xlab="# of common peaks", ylab="FDR",
       main="crossmatching FDR", pch=19)
  matplot(t(contingencies[1:4, ]), type="l", lty=1,
          xlab="# of common peaks", ylab="# of peptides",
          main="crossmatching details")
  legend("bottomright", legend=c("TP", "FP", "TN", "FN"),
         col=1:4, lwd=1, bty="n", inset=0.05)
  par(mfcol=c(1, 1))

  invisible(contingencies)
}

#' calculate and plot the difference of the first and the second highest
#' number of peaks in a non-unique-match group
#' @param cx crossmatching df, result of crossmatching
#' @param tolerance double, allowed deviation
#' @param verbose verbose output?
#' @return invisible data.frame, columns: diff, true
#' @noRd
plotCrossmatchingNonuniqueDiff <- function(cx, breaks=50) {
  cx <- cx[grep("non-unique", cx$matchType), ]
  cx$true <- grepl("true", cx$matchType)

  d <- tapply(1:nrow(cx), cx$precursor.leID.ident, function(i) {
    s <- cx[i, ]
    com <- sort(s$fragments.identXfragments.quant,
                decreasing=TRUE, index.return=TRUE)
    return(c(diff=com$x[1]-com$x[2], true=isTRUE(s$true[com$ix[1]])))
  })
  d <- as.data.frame(do.call(rbind, d))
  d$true <- as.logical(d$true)
  d$precursor.leID.ident <- as.numeric(rownames(d))

  hist(d$diff[d$true], col=3,
       main="Common Peak Number Differences for Non-Unique Matches",
       xlab="# of common peaks", ylab="# peptides", breaks=breaks)
  hist(d$diff[!d$true], col=2, add=TRUE, breaks=breaks)
  legend("topright", legend=c("true match", "false match"),
         col=3:2, pch=15, bty="n")

  invisible(d)
}

