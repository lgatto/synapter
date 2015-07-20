#' filter Neutral Loss
#' @param df data.frame (needs a column Neutral.LossType)
#' @return a data.frame that contains only rows with Neutral.LossType == "None"
#' @noRd
.filterNeutralLoss <- function(df) {
  df[which(df$Neutral.LossType == "None"), ]
}

#' filter empty fragments
#' @param df data.frame (needs a column Neutral.LossType)
#' @return a data.frame that contains only rows with Neutral.LossType != ""
#' @noRd
.filterEmptyFragments <- function(df) {
  df[which(nchar(df$Neutral.LossType) > 0), ]
}

#' read final_fragments.csv
#' @param file file path
#' @param removeNeutralLoss remove rows with neutral loss != "none"?
#' @param verbose verbose output?
#' @return a data.frame
#' @noRd
.readFragments <- function(file, removeNeutralLoss=FALSE, verbose=TRUE) {
  stopifnot(file.exists(file))

  if (verbose) {
    message("Reading ", file)
  }

  df <- readcsv(file, keepCols = c("precursor.leID" = "d",
                                   "peptide.seq" = "c",
                                   "fragment.str" = "c",
                                   "precursor.mhp" = "d",
                                   "precursor.inten" = "d",
                                   "precursor.retT" = "d",
                                   "precursor.z" = "d",
                                   "product.mhp" = "d",
                                   "product.inten" = "d",
                                   "Neutral.LossType" = "c"))
  df <- .filterEmptyFragments(df)

  ## convert non-ascii characters to _
  df$fragment.str <- iconv(x=df$fragment.str, to="ASCII", sub="_")

  if (removeNeutralLoss) {
    return(.filterNeutralLoss(df))
  } else {
    return(df)
  }
}

#' read final_fragment.csv and turn data into MSnbase::Spectrum2 objects
#' @param df corresponding df from the synapter object ({Ident,Quant}PeptideData)
#' @param file filename
#' @param storeAll should all spectra stored? or only the needed ones?
#' @param removeNeutralLoss remove rows with neutral loss != "none"?
#' @param removePrecursor remove precursor ion from fragments
#' @param tolerance tolerance to look for the precursor ion
#' @param verbose verbose output
#' @return modified assaydata
#' @noRd
.finalFragment2spectra <- function(df, file, storeAll=TRUE,
                                   removeNeutralLoss=TRUE,
                                   removePrecursor=TRUE,
                                   tolerance=25e-6, verbose=TRUE) {
  fragments <- .readFragments(file, removeNeutralLoss=removeNeutralLoss,
                              verbose=verbose)
  .fragments2spectra(df=df, fragments=fragments, file=file, storeAll=storeAll,
                     removePrecursor=removePrecursor, tolerance=tolerance,
                     verbose=verbose)
}

#' create a MSnbase::Spectrum2 objects out of a final_fragments data.frame
#' @param df corresponding df from the synapter object ({Ident,Quant}PeptideData)
#' @param fragments final_fragments data.frame
#' @param file filename (just for the metadata)
#' @param storeAll should all spectra stored? or only the needed ones?
#' @param removePrecursor remove precursor ion from fragments
#' @param tolerance tolerance to look for the precursor ion
#' @param verbose verbose output
#' @return modified assaydata
#' @noRd
.fragments2spectra <- function(df, fragments, file, storeAll=TRUE,
                               removePrecursor=TRUE, tolerance=25e-6,
                               verbose=TRUE) {

  assignments <- new.env(hash=TRUE, parent=emptyenv())
  idx <- split(1:nrow(fragments), f=fragments$precursor.leID)
  idx_keys <- names(idx)

  for (i in seq(along=idx)) {
    assign(idx_keys[i], idx[[i]], envir=assignments)
  }

  if (storeAll) {
    uleID <- unique(fragments$precursor.leID)
  } else {
    uleID <- intersect(df$precursor.leID, fragments$precursor.leID)
  }

  if (verbose) {
    message("Convert fragment data.frame to MSnbase::Spectrum2 objects")
    pb <- txtProgressBar(0, length(uleID), style=3)
  }

  fragment.str <- rep(NA, length(uleID))
  sequences <- rep(NA, length(uleID))

  assaydata <- new.env(hash=TRUE, parent=emptyenv(), size=length(uleID))

  for (i in seq(along=uleID)) {
    assign(as.character(uleID[i]),
           .createMs2SpectrumFromFragments(uleID[i],
                                           fragments=fragments,
                                           assignments=assignments,
                                           removePrecursor=removePrecursor,
                                           tolerance=tolerance),
           envir=assaydata)
    fragment.str[i] <- .getFragmentStrFromFragments(uleID[i],
                                                    fragments=fragments,
                                                    assignments=assignments,
                                                    removePrecursor=removePrecursor,
                                                    tolerance=tolerance)
    sequences[i] <- .getPeptideSeqFromFragments(uleID[i],
                                                fragments=fragments,
                                                assignments=assignments)
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
                                                      type="fragment"))
  fdata <- data.frame(spectrum=1:length(nm),
                      leID=uleID,
                      fragment.str=fragment.str,
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

#' fetch fragment.str in correct order
#' @param leID precursor.leID
#' @param fragments data.frame generated by readFragments
#' @param assignments env
#' @param removePrecursor remove precursor ion from fragments
#' @param tolerance tolerance to look for the precursor ion
#' @return character vector
#' @noRd
.getFragmentStrFromFragments <- function(leID, fragments, assignments,
                                         removePrecursor=TRUE,
                                         tolerance=25e-6) {
  key <- as.character(leID)

  if (exists(key, envir=assignments)) {
    i <- get(key, envir=assignments)

    o <- order(fragments$product.mhp[i])

    if (removePrecursor) {
      m <- !as.logical(MSnbase:::relaxedMatch(x=fragments$product.mhp[i[o]],
                                              table=fragments$precursor.mhp[i[1]],
                                              nomatch=0,
                                              relative=TRUE,
                                              tolerance=tolerance))
      o <- o[m]

    }

    fragment.str <- MSnbase:::utils.vec2ssv(fragments$fragment.str[i][o])

    if (nzchar(fragment.str)) {
      return(fragment.str)
    }
  }
  return(NA)
}

#' fetch peptide.seq
#' @param leID precursor.leID
#' @param fragments data.frame generated by readFragments
#' @param assignments env
#' @return character vector
#' @noRd
.getPeptideSeqFromFragments <- function(leID, fragments, assignments) {
  key <- as.character(leID)

  if (exists(key, envir=assignments)) {
    i <- get(key, envir=assignments)

    peptide.seq <- fragments$peptide.seq[i[1]]
    if (!is.null(peptide.seq)) {
      return(peptide.seq)
    }
  }
  return(NA)
}

