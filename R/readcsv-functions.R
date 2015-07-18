#' fetch colnames from file
#' @param file
#' @param sep separator that splits the columns
#' @return character
#' @noRd
fetchColnames <- function(file, sep = ",") {
  header <- scan(file, what = character(), nlines = 1L, quiet = TRUE)
  unlist(strsplit(paste0(header, collapse=""), split=sep))
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
  cols[names(keepCols)] <- keepCols
  paste0(cols, collapse="")
}

#' wrapper around readr::read_csv to read just needed columns
#' @param file
#' @param keepCols named character vector for columns of interest
#' @return tbl_df
#' @examples
#'  createReadrColTypes("final_peptide.csv", c("peptide.seq" = "c",
#'                                             "peptide.score" = "d",
#'                                             "peptide.matchType" = "c",
#'                                             "protein.dataBaseType" = "c"))
#' @noRd
readcsv <- function(file, keepCols) {
  read_csv(file, col_types = createReadrColTypes(file, keepCols))
}
