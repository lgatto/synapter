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
#'
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
#' @param prefix character, prefix for the keyvalue in assaydata
#' (e.g.,"ident.spectra")
#' @param assaydata environment
#' @param storeAll should all spectra stored? or only the needed ones?
#' @param verbose verbose output
#' @return modified assaydata
.spectrumXml2spectra <- function(df, file, prefix, storeAll=TRUE, 
                                 verbose=TRUE) {
  stopifnot(!missing(prefix))

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

  keys <- paste(prefix, uleID, sep=":")
  sequences <- rep(NA, length(uleID))

  assaydata <- new.env(hash=TRUE, parent=emptyenv(), size=length(uleID))

  for (i in seq(along=uleID)) {
    assign(keys[i],
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
                      row.names=keys,
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
.createEmptyMsnbaseSpectrum2 <- function(leID, key="spectrum:-1") {
  if (missing(leID)) {
    leID <- as.integer(tail(strsplit(key, ":")[[1]], 1))
  }
  new("Spectrum2", precScanNum=as.integer(leID))
}

#' @param spectra list, 4 MSnbase::Spectrum2 objects
#' @param norm normalise spectra?
#' @param fragments list, 2 character vectors containing the fragment.str
#' @param tolerance double, allowed deviation
#' @param ... passed to .plotIdentVsQuantSpectra
.plotSpectraVsFragments <- function(spectra, norm=TRUE, fragments,
                                    tolerance=25e-6, ...) {
  if (norm) {
    spectra <- lapply(spectra, normalize, method="precursor")
  }

  mass <- unlist(lapply(spectra, mz))
  xlim <- c(min(mass), max(mass))

  inten <- unlist(lapply(spectra, intensity))
  maxInten <- max(inten)
  ylim <- c(-maxInten, maxInten)

  oldPar <- par(no.readonly=TRUE)
  on.exit(par(oldPar))
  par(mfrow=c(1, 2))

  par(mar=c(2, 2, 2, 0.5)) #c(bottom, left, top, right)
  .plotIdentVsQuantSpectra(spectra[1:2], main="spectra",
                           common=list(.commonPeaks(spectra[[4]],
                                                    spectra[[1]],
                                                    tolerance),
                                       .commonPeaks(spectra[[3]],
                                                    spectra[[2]],
                                                    tolerance)),
                           xlim=xlim, ylim=ylim, ...)
  par(mar=c(2, 0.5, 2, 2)) #c(bottom, left, top, right)
  .plotIdentVsQuantSpectra(spectra[3:4], main="fragments", fragments=fragments,
                           common=list(.commonPeaks(spectra[[4]],
                                                    spectra[[3]],
                                                    tolerance),
                                       .commonPeaks(spectra[[3]],
                                                    spectra[[4]],
                                                    tolerance)),
                           xlim=xlim, ylim=ylim, yaxt="n", ...)
  axis(4)
}

#' plot ident vs quant spectra
#' @param spectra list, 2 MSnbase::Spectrum2 objects
#' @param sequences list, 2 character vectors containing the peptide.seq
#' @param fragments list, 2 character vectors containing the fragment.str
#' @param common list, 2 logical vectors
#' @param xlab label for x-axis
#' @param ylab label for y-axis
#' @param xlim limits for x-axis
#' @param ylim limits for y-axis
#' @param fragments.cex cex for fragments
#' @param legend.cex cex for legend
.plotIdentVsQuantSpectra <- function(spectra,
                                     sequences,
                                     fragments=list(character(), character()),
                                     common,
                                     main=character(),
                                     xlab="mhp", ylab="intensity",
                                     xlim, ylim,
                                     fragments.cex=0.5,
                                     legend.cex=0.5, ...) {
  if (missing(xlim)) {
    mass <- unlist(lapply(spectra, mz))
    xlim <- c(min(c(mass, 0)), max(c(mass, 0)))
  }

  if (missing(ylim)) {
    inten <- unlist(lapply(spectra, intensity))
    maxInten <- max(c(inten, 0))
    ylim <- c(-maxInten, maxInten)
  }

  if (missing(common)) {
    common <- lapply(spectra, function(x)logical(peaksCount(x)))
  }

  plot(NA, type="h", col=1,
       main=main,
       xlab=xlab, ylab=ylab,
       xlim=xlim, ylim=ylim,
       ...)
  abline(h=0, col="#808080")

  orientation <- c(1, -1)
  text.pos <- c(3, 1)
  legend.pos <- c("topleft", "bottomleft")
  legend.prefix <- c("ident", "quant")
  pal <- brewer.pal(11, "RdYlBu")
  cols <- c(pal[c(9, 11)], pal[c(3, 1)])
  pch <- c(NA, 19)

  for (i in seq(along=spectra)) {
    lines(mz(spectra[[i]]), orientation[i]*intensity(spectra[[i]]),
          type="h", col=cols[(i-1)*2+common[[i]]+1L], lwd=1.5)
    points(mz(spectra[[i]]), orientation[i]*intensity(spectra[[i]]),
           col=cols[(i-1)*2+common[[i]]+1L], pch=pch[common[[i]]+1L],
           cex=0.5)

    if (length(fragments[[i]]) && !anyNA(fragments[[i]])) {
      text(mz(spectra[[i]]), orientation[i]*intensity(spectra[[i]]),
           fragments[[i]], pos=text.pos[i], offset=0.25,
           cex=fragments.cex, col="#808080")
    }

    label <- paste0(legend.prefix[i], ".leID: ", precScanNum(spectra[[i]]))

    if (peaksCount(spectra[[i]])) {
      label <- paste0(label, ", mhp: ", precursorMz(spectra[[i]]),
                             ", z: ", precursorCharge(spectra[[i]]),
                             ", seq: ", sequences[[i]],
                             ", cx: ", sum(common[[i]]))
    }

    legend(legend.pos[i], legend=label, bty="n", cex=legend.cex)
  }
}

#' readSpectraAndFragments
#' @param obj synapter object
#' @param filenames named list of filenames list(identspectrum, quantspectrum,
#' @param removeNeutralLoss remove rows with neutral loss != "none"?
#' @param verbose verbose output?
#' @return  env with MSnbase::Spectrum2 objects
readSpectraAndFragments <- function(obj, filenames, removeNeutralLoss=TRUE,
                                    verbose=TRUE) {
  stopifnot(all(names(filenames) %in% c("identspectrum", "quantspectrum",
                                        "identfragments", "quantfragments")))
  filenames <- as.list(filenames)

  spectra <- vector(mode="list", length=length(filenames))
  names(spectra) <- c("spectra.ident", "spectra.quant",
                      "fragments.ident", "fragments.quant")

  spectra$spectra.ident <- 
    .spectrumXml2spectra(df=obj$IdentPeptideData,
                         file=filenames$identspectrum,
                         prefix="spectra.ident",
                         storeAll=FALSE,
                         verbose=verbose)
  spectra$spectra.quant <- 
    .spectrumXml2spectra(df=obj$QuantPeptideData,
                         file=filenames$quantspectrum,
                         prefix="spectra.quant",
                         storeAll=TRUE,
                         verbose=verbose)

  spectra$fragments.ident <- 
    .finalFragment2spectra(df=obj$IdentPeptideData,
                           file=filenames$identfragments,
                           prefix="fragments.ident",
                           storeAll=FALSE,
                           removeNeutralLoss=removeNeutralLoss,
                           verbose=verbose)
  spectra$fragments.quant <- 
    .finalFragment2spectra(df=obj$QuantPeptideData,
                           file=filenames$quantfragments,
                           prefix="fragments.quant",
                           storeAll=TRUE,
                           removeNeutralLoss=removeNeutralLoss,
                           verbose=verbose)
  return(spectra)
}

#' crossmatching, wrapper function
#' @param obj synapter object
#' @param spectra list of 4 MSnExp containing the spectra and fragment spectra
#' @param tolerance double, allowed deviation to consider a m/z as equal
#' @param verbose verbose output?
#' @return data.frame, extend/flatted matchedEmrts df with additional columns:
#' matchType,
#' spectra.identXfragments.ident, spectra.quantXfragments.quant,
#' spectra.identXfragments.quant, spectra.quantXfragments.ident
#' sorry for the names
crossmatching <- function(obj, spectra, tolerance=25e-6, verbose=TRUE) {
  if (verbose) {
    message("create flat EMRTs data.frame")
  }
  emrts <- flatMatchedEMRTs(obj$MatchedEMRTs)

  return(.cossmatching(flatEmrtsemrts=emrts, spectra=spectra,
                       tolerance=tolerance, verbose=verbose))
}

#' .crossmatching, workhorse function
#' compares ident fragments vs quant product spectrum and
#' quant fragments vs ident product spectrum
#' @param flatEmrts flattened EMRTs (see flatMatchedEMRTs)
#' @param spectra list of 4 MSnExp containing the spectra and fragment spectra
#' @param tolerance double, allowed deviation to consider a m/z as equal
#' @param verbose verbose output?
#' @return data.frame, extend/flatted matchedEmrts df with additional columns:
#' matchType,
#' spectra.identXfragments.ident, spectra.quantXfragments.quant,
#' spectra.identXfragments.quant, spectra.quantXfragments.ident
#' sorry for the names
.crossmatching <- function(flatEmrts, spectra, tolerance=25e-6, verbose=TRUE) {
  prefixes <- paste(rep(c("spectra", "fragments"), each=2),
                    rep(c("ident", "quant"), times=2), sep=".")
  # "spectra.ident"   "spectra.quant"   "fragments.ident" "fragments.quant"

  keys <- paste(rep(prefixes, each=nrow(flatEmrts)),
                rep(unlist(flatEmrts[, c("precursor.leID.ident",
                                         "precursor.leID.quant")]), times=2),
                sep=":")
  keysm <- matrix(keys, ncol=4)

  cmb <- list(c(1, 3), c(2, 4),
              c(1, 4), c(2, 3))

  cols <- sapply(cmb, function(x)paste0(prefixes[x], collapse="X"))
  # "spectra.identXfragments.ident" ...

  if (verbose) {
    message("Look for common peaks")
    pb <- txtProgressBar(0, 4*nrow(flatEmrts), style=3)
  }

  flatEmrts[, cols] <- lapply(cmb, function(i) {
    apply(keysm[, i], 1, function(k) {
      if (verbose) {
        setTxtProgressBar(pb, pb$getVal()+1)
      }
      .nCommonPeaks(.getSpectrum(k[1], spectra[[i[1]]]),
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
#' @param cx crossmatching df, result of crossmatching
#' @param spectra list of 4 MSnExp containing the spectra and fragment spectra
#' @param legend.cex cex for legend text
#' @param tolerance double, allowed deviation
#' @param verbose verbose output?
#' @return data.frame, extend/flatted matchedEmrts df with additional columns:
#' cxIdentSxQuantF, cxIdentFxQuantS, matchType
plotCrossmatching <- function(cx, spectra,
                              legend.cex=0.75, tolerance=25e-6,
                              verbose=TRUE) {
  prefixes <- paste(rep(c("spectra", "fragments"), each=2),
                    rep(c("ident", "quant"), times=2), sep=".")
  keys <- paste(rep(prefixes, each=nrow(cx)),
                rep(unlist(cx[, c("precursor.leID.ident",
                                  "matched.quant.spectrumIDs")]), times=2),
                sep=":")
  keysm <- matrix(keys, ncol=4)

  if (verbose) {
    pb <- txtProgressBar(0, nrow(cx), style=3)
  }

  for (i in 1:nrow(cx)) {

    sequences <- mapply(function(x, k) {
      ## avoid partial matching of rows 
      ## ( [.data.frame always use partial match for rows; exact=TRUE works only
      ## for column names)
      return(fData(x)[match(k, rownames(fData(x))), "peptide.seq"])
    }, x=spectra, k=keysm[i, ], SIMPLIFY=FALSE, USE.NAMES=FALSE)

    fragments <- mapply(function(x, k) {
      ## avoid partial matching (see above)
      fragment.str <- fData(x)[match(k, rownames(fData(x))), "fragment.str"]
      return(MSnbase:::utils.ssv2vec(fragment.str))
    }, x=spectra[3:4], k=keysm[i, 3:4], SIMPLIFY=FALSE, USE.NAMES=FALSE)

    .plotSpectraVsFragments(spectra=.getSpectra(keysm[i, ], spectralist=spectra),
                            sequences=sequences, fragments=fragments,
                            tolerance=tolerance, legend.cex=legend.cex)
    if (verbose) {
      setTxtProgressBar(pb, i)
    }
  }

  if (verbose) {
    close(pb)
  }
}

