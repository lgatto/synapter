##' This function takes a final peptide file and
##' returns information about the \emph{unique} peptide scores
##' and their number in each peptide matchType (PepFrag1 and
##' PepFrag2) by protein dataBaseType (Random and Regular)
##' category. The function is a lightweight verion of
##' \code{\link{getPepNumbers}} and \code{\link{plotPepScores}}
##' (to be used with \code{Synapter} instances) for
##' individual files.
##'
##' @title Inspect peptide scores.
##' @param filename The name of a final peptide file.
##' @return A table of peptide counts in each peptide matchType *
##' protein dataBasteType category. Also plots the distribution of
##' respective peptide scores.
##' @author Laurent Gatto
inspectPeptideScores <- function(filename) {
  protein.dataBaseType <- NULL ## no visible binding note
  x <- readcsv(filename, c("peptide.seq" = "c",
                           "peptide.score" = "d",
                           "peptide.matchType" = "c",
                           "protein.dataBaseType" = "c"))
  x <- x[x$peptide.matchType %in% c("PepFrag1", "PepFrag2"), ]
  x <- x[!duplicated(x$peptide.seq), ]
  p <- densityplot( ~ peptide.score | peptide.matchType,
                   data = x, groups = protein.dataBaseType,
                   plot.points = FALSE, ref = TRUE,
                   main = basename(filename),
                   auto.key = list(columns = 2))
  tab <- unlist(by(x, x$peptide.matchType,
                   function(.x) table(.x$protein.dataBaseType)))
  print(p)
  return(tab)
}

##' Opens the relevant vignette; a shortcut to using \code{vignette}.
##' \code{synapterGuide()} gives access to the main overview vignette.
##'
##' @title Opens the 'synapter' vignette
##' @usage
##' synapterGuide()
##' @author Laurent Gatto
synapterGuide <- function()
  vignette("synapter", package = "synapter")

synapterTinyData <- function() {
  f <- system.file("extdata", "04_test_database.fasta", package = "synapter")
  data(synapterTiny, envir = .GlobalEnv)
  .GlobalEnv$synapterTiny$DbFastaFile <- f
}
