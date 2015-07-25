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

  leRows <- new.env(hash=TRUE, parent=emptyenv(), size=length(le_ids))
  heRows <- new.env(hash=TRUE, parent=emptyenv(), size=length(le_ids))

  le_rows <- findInterval(le_ids, ms1[, "LE_ID"])
  he_rows <- relist(findInterval(unlist(he_ids), ms2[, "HE_ID"]), he_ids)

  for (i in seq(along=le_ids)) {
    assign(le_ids[i], le_rows[[i]], envir=leRows)
    assign(le_ids[i], he_rows[[i]], envir=heRows)
  }

  return(list(ms1=ms1, ms2=ms2, leRows=leRows, heRows=heRows))
}

#' read spectrum.xml and turn data into MSnbase::Spectrum2 objects
#' @param df corresponding df from the synapter object ({Ident,Quant}PeptideData)
#' @param file filename
#' @param assaydata environment
#' @param storeAll should all spectra stored? or only the needed ones?
#' @param removePrecursor remove precursor ion from fragments
#' @param tolerance tolerance to look for the precursor ion
#' @param verbose verbose output
#' @return modified assaydata
#' @noRd
.spectrumXml2spectra <- function(df, file, storeAll=TRUE, removePrecursor=TRUE,
                                 tolerance=25e-6, verbose=TRUE) {
  xml <- .readSynapterSpectrumXml(file, ms1=TRUE, verbose=verbose)

  peptideinfo <- df[!duplicated(df$precursor.leID), ]

  if (storeAll) {
    uleID <- union(df$precursor.leID, xml$ms1[, "LE_ID"])
  } else {
    uleID <- intersect(df$precursor.leID, xml$ms1[, "LE_ID"])
  }

  if (verbose) {
    message("Convert spectra data.frames to MSnbase::Spectrum2 objects")
    pb <- txtProgressBar(0, length(uleID), style=3)
  }

  sequences <- rep(NA_character_, length(uleID))

  assaydata <- new.env(hash=TRUE, parent=emptyenv(), size=length(uleID))

  for (i in seq(along=uleID)) {
    assign(as.character(uleID[i]),
           .createMs2SpectrumFromSpectrumXml(uleID[i], xml$ms1, xml$ms2,
                                             xml$leRows, xml$heRows,
                                             removePrecursor=removePrecursor,
                                             tolerance=tolerance),
           envir=assaydata)
    if (i <= nrow(df)) {
      ## in the xml files are no sequence information that's why we have to use
      ## information from the df
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
#' @param leRows env generated by .readSynapterSpectrumXml
#' @param heRows env generated by .readSynapterSpectrumXml
#' @param removePrecursor remove precursor ion from fragments
#' @param tolerance tolerance to look for the precursor ion
#' @return Spectrum2 object
#' @noRd
.createMs2SpectrumFromSpectrumXml <- function(leID, ms1, ms2, leRows, heRows,
                                              removePrecursor=TRUE,
                                              tolerance=25e-6) {
  key <- as.character(leID)

  if (exists(key, envir=leRows)) {
    lrow <- get(key, envir=leRows)
    hrows <- get(key, envir=heRows)

    return(.createMsnbaseSpectrum2(leID=leID,
                                   precursor.mhp=ms1[lrow, "Mass"],
                                   precursor.inten=ms1[lrow, "Intensity"],
                                   precursor.z=ms1[lrow, "Z"],
                                   precursor.retT=ms1[lrow, "RT"],
                                   mass=ms2[hrows, "Mass"],
                                   intensity=ms2[hrows, "Intensity"],
                                   removePrecursor=removePrecursor,
                                   tolerance=tolerance))
  } else {
    return(.createEmptyMsnbaseSpectrum2(leID=leID))
  }
}

#' create MS2 spectrum from final_fragment.csv data
#' @param leID precursor.leID
#' @param fragments data.frame generated by readFragments
#' @param assignments env
#' @param removePrecursor remove precursor ion from fragments
#' @param tolerance tolerance to look for the precursor ion
#' @return Spectrum2 object
#' @noRd
.createMs2SpectrumFromFragments <- function(leID, fragments, assignments,
                                            removePrecursor=TRUE,
                                            tolerance=25e-6) {
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
                                   intensity=fragments$product.inten[i],
                                   removePrecursor=removePrecursor,
                                   tolerance=25e-6))
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
#' @param removePrecursor remove precursor ion from fragments
#' @param tolerance tolerance to look for the precursor ion
#' @return Spectrum2 object
#' @noRd
.createMsnbaseSpectrum2 <- function(leID, precursor.mhp, precursor.inten,
                                    precursor.z, precursor.retT,
                                    mass, intensity,
                                    removePrecursor=TRUE, tolerance=25e-6) {
  ## just to be sure, you can't trust any information in the csv/xml files
  o <- order(mass)

  if (removePrecursor) {
    m <- !as.logical(MSnbase:::relaxedMatch(x=mass[o], table=precursor.mhp,
                                            nomatch=0,
                                            relative=TRUE,
                                            tolerance=tolerance))
    o <- o[m]

    if (!length(o)) {
      return(.createEmptyMsnbaseSpectrum2(leID=leID))
    }
  }

  new("Spectrum2",
      precScanNum=as.integer(leID),
      precursorMz=precursor.mhp,
      precursorIntensity=precursor.inten,
      precursorCharge=as.integer(precursor.z),
      rt=precursor.retT,
      centroided=TRUE,
      # tic=precursor.inten,
      peaksCount=length(o),
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

  new("Spectrum2", precScanNum=as.integer(leID), fromFile=1L)
}

#' @param msexp MSnExp object to filter
#' @param minIntensity minimal intensity
#' @param maxNumber maximal number of fragments to preserve
#' @param verbose verbose output?
#' @return modified MSnExp object
#' @noRd
.filterIntensity <- function(msexp, minIntensity=NULL, maxNumber=NULL,
                             verbose=TRUE) {
  if(is.null(minIntensity) && is.null(maxNumber)) {
    stop("At least one of arguments ", sQuote("minIntensity"), " and ",
         sQuote("maxNumber"), " must be given!")
  }

  minIntensity2 <- 0

  if (!is.null(maxNumber)) {
    cs <- .cumsumIntensities(msexp)
    minIntensity2 <- as.double(names(cs)[which(cs > (tail(cs, 1)-maxNumber))[1]])
  }

  minIntensity <- max(minIntensity, minIntensity2)

  if (verbose) {
    message("Set peaks with an intensity < ", minIntensity, " to zero")
    pcBefore <- .sumAllPeakCounts(msexp)
  }

  msexp <- removePeaks(msexp, t=minIntensity, verbose=verbose)

  ## remove corresponding fragment information if available
  stopifnot(fData(msexp)$leID == precScanNum(msexp))
  if ("fragment.str" %in% colnames(fData(msexp))) {
    fragments <- MSnbase:::utils.ssv2list(as.character(fData(msexp)$fragment.str))
    fData(msexp)$fragment.str <- MSnbase:::utils.list2ssv(mapply(function(f, i) {
      if (length(f) && length(i)) {
        as.character(f[i > 0])
      } else {
        character()
      }
    }, f=fragments, i=intensity(msexp), SIMPLIFY=FALSE, USE.NAMES=FALSE))

    stopifnot(validObject(msexp))
  }

  if (verbose) {
    message("Remove peaks with zero intensity")
    pcBefore <- .sumAllPeakCounts(msexp)
  }

  msexp <- clean(msexp, all=TRUE, verbose=verbose)

  if (verbose) {
    pcAfter <- .sumAllPeakCounts(msexp)
    message(pcBefore - pcAfter, " peaks removed; new total number of peaks: ",
            pcAfter)
  }
  msexp
}

#' @param msexp MSnExp object to filter
#' @noRd
.plotIntensityVsNumber <- function(msexp, what) {
  cs <- .cumsumIntensities(msexp)
  plot(as.double(names(cs)), cs, log="x", type="l",
       main=paste0("Cumulative Number of Fragments (", what, ")"),
       xlab="Intensity", ylab="# of Fragments")
  grid()
}

.cumsumIntensities <- function(msexp) {
  cumsum(table(unlist(intensity(msexp))))
}

