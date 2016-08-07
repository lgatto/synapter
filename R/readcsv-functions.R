#' wrapper around readr::read_csv to read just needed columns
#' @param file filename
#' @param keepCols named character vector for columns of interest
#' @param verbose verbose output?
#' @return data.frame
#' @examples
#'  readcsv("final_peptide.csv", c("peptide.seq" = "c",
#'                                 "peptide.score" = "d",
#'                                 "peptide.matchType" = "c",
#'                                 "protein.dataBaseType" = "c"))
#' @noRd
readcsv <- function(file, keepCols, verbose=interactive()) {
  as.data.frame(read_csv(file,
                         col_types = do.call(cols_only, as.list(keepCols)),
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
