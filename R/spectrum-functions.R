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
#'  assignments: data.frame that holds the connection between ms1 and ms2
#'
.readSynapterSpectrumXml <- function(file, ms1=FALSE,
                                     encoding="Windows-1252",
                                     verbose=TRUE) {
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

  assignments <- data.frame(le_id=.extractLEID(content[range]),
                            he_id=.extractHEID(content[range]),
                            stringsAsFactors=FALSE)

  return(list(ms1=ms1, ms2=ms2, assignments=assignments))
}

