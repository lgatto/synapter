#' fetch colnames from file
#' @param file
#' @param sep separator that splits the columns
#' @return character
#' @noRd
fetchColnames <- function(file, sep = ",") {
  header <- scan(file, what = character(), nlines = 1L, quiet = TRUE)
  header <- gsub("\"", "", header)
  trimws(unlist(strsplit(paste(header, collapse=" "), split=sep)))
}

#' create col_types argument for readr::read_csv
#' @param file
#' @param keepCols named character vector for columns of interest
#' @return character
#' @examples
#'  createReadrColTypes("final_peptide.csv", c("peptide.seq" = "c",
#'                                             "peptide.score" = "d",
#'                                             "peptide.matchType" = "c",
#'                                             "protein.dataBaseType" = "c"))
#' @noRd
createReadrColTypes <- function(file, keepCols) {
  cols <- fetchColnames(file)
  cols <- setNames(rep("_", length(cols)), cols)
  keepCols <- keepCols[names(keepCols) %in% names(cols)]
  cols[names(keepCols)] <- keepCols
  paste0(cols, collapse="")
}

#' wrapper around readr::read_csv to read just needed columns
#' @param file
#' @param keepCols named character vector for columns of interest
#' @param verbose
#' @return tbl_df
#' @examples
#'  readcsv("final_peptide.csv", c("peptide.seq" = "c",
#'                                 "peptide.score" = "d",
#'                                 "peptide.matchType" = "c",
#'                                 "protein.dataBaseType" = "c"))
#' @noRd
readcsv <- function(file, keepCols, verbose=interactive()) {
  as.data.frame(read_csv(file, col_types = createReadrColTypes(file, keepCols),
                progress=verbose))
}

#' @noRd
readFinalPeptides <- function(file, verbose=interactive()) {
  readcsv(file, keepCols=c("protein.Accession" = "c",
                           "protein.Description" = "c",
                           "protein.dataBaseType" = "c",
                           "protein.falsePositiveRate" = "d",
                           "peptide.matchType" = "c",
                           "peptide.mhp" = "d",
                           "peptide.seq" = "c",
                           "peptide.score" = "d",
                           "precursor.leID" = "d",
                           "precursor.mhp" = "d",
                           "precursor.retT" = "d",
                           "precursor.inten" = "d",
                           "precursor.z" = "d",
                           "precursor.mz" = "d",
                           "precursor.Mobility" = "d") , verbose=verbose)
}

#' @noRd
readPep3D <- function(file, verbose=interactive()) {
  pep3d <- readcsv(file, keepCols=c("Function" = "d",
                                    "spectrumID" = "d",
                                    "rt_min" = "d",
                                    "mwHPlus" = "d",
                                    "charge" = "d",
                                    "Intensity" = "d",
                                    "Counts" = "d",
                                    "clust_drift" = "d",
                                    "isFid" = "d",
                                    "ion_ID" = "d",
                                    "ion_z" = "d",
                                    "ion_iso" = "d",
                                    "ion_area" = "d",
                                    "ion_counts" = "d",
                                    "Model" = "d"), verbose=verbose)
  ## rename Function into matchedEMRTs; see issue #67
  colnames(pep3d)[1L] <- "matchedEMRTs"
  pep3d
}
