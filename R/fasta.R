dbUniquePeptideSet <- function(file, missedCleavages = 0,
                               IisL = FALSE, verbose = TRUE) {
  if (tolower(file_ext(file)) == "rds") {
    peptides <- readRDS(file)

    noeffectMsg <- paste0("\nYour current setting has not any effect! ",
                          "Recreate your RDS file if you want to change ",
                          "these settings.")

    if (!all(attr(peptides, "missedCleavages") == missedCleavages)) {
      warning("The RDS file was created with missedCleavages=",
              attr(peptides, "missedCleavages"), noeffectMsg, immediate.=TRUE)
    }
    if (!all(attr(peptides, "IisL") == IisL)) {
      warning("The RDS file was created with IisL=",
              attr(peptides, "IisL"), noeffectMsg, immediate.=TRUE)
    }
  } else {
    peptides <- .dbUniquePeptideSet(file, missedCleavages=missedCleavages,
                                    IisL=IisL, verbose=verbose)
  }

  if (verbose) {
      topics <- c("Peptide database file: ",
                  "Allowed mis-cleavages: ",
                  "I/L treatment: ")
      values <- c(file,
                  attr(peptides, "missedCleavages"),
                  ifelse(attr(peptides, "IisL"), "I == L", "I != L"))

      topics <- format(topics, justify="left")

      message(paste0(topics, values, collapse="\n"))
  }

  return(peptides)
}

.dbUniquePeptideSet <- function(fastafile, missedCleavages = 0,
                                IisL = FALSE, verbose = TRUE) {

  missedCleavages <- max(missedCleavages)
  stopifnot(missedCleavages >= 0)

  proteins <- as.character(readAAStringSet(fastafile))

  ## using default cleaver cleavages rules for trypsin, defined:
  ## http://web.expasy.org/peptide_cutter/peptidecutter_enzymes.html#Tryps
  peptides <- unlist(cleave(proteins, enzym="trypsin",
                            missedCleavages=0:missedCleavages),
                     use.names=FALSE)

  ## remove non-unique peptides
  ## Note that this does however not take pre/post-fix peptides into account
  ## like "DLIELTESLFR" and "HNPEFTMMELYMAYADYHDLIELTESLFR",
  ##      "YYGYTGAFR" and "EGYYGYTGAFR"
  ## where the first ones are NOT unique!
  ## Such cases are handled by filtering duplicates in the peptide data
  if (IisL) {
    peptides <- gsub(pattern="[IL]", replacement="-", x=peptides)
    upeptides <- peptides[!.duplicated2(.peptides)]
  } else {
    upeptides <- peptides[!.duplicated2(peptides)]
  }

  attr(upeptides, "missedCleavages") <- missedCleavages
  attr(upeptides, "IisL") <- IisL

  if (verbose) {
    topics <- c("Cleavage results:      ", "")
    values <- c(paste0(length(proteins), " proteins"),
                paste0(length(upeptides), " out of ", length(peptides),
                " tryptic peptides are proteotypic."))

    topics <- format(topics, justify="left")
    message(paste0(topics, values, collapse="\n"))
  }

  return(upeptides)
}

#' This function creates an RDS file to store the \dQuote{Unique Peptides
#' Database} for future runs of \code{\link{synergise}} or
#' \code{\link{Synapter}}.
#'
#' @title Create an RDS file for the 'Unique Peptides Database'
#' @param fastaFile file path of the input fasta file
#' @param outputFile file path of the target RDS file; must have the file
#' extension ".rds"
#' @param missedCleavages Number of maximal allowed missed cleavages. Default
#' is 0.
#' @param IisL If \code{TRUE} Isoleucin and Leucin are treated as identical.
#' In this case sequences like "ABCI", "ABCL" are removed because they
#' are not unqiue. If \code{FALSE} (default) "ABCI" and "ABCL" are reported as
#' unique.
#' @param verbose If \code{TRUE} a verbose output is provied.
#' @author Sebastian Gibb <mail@@sebastiangibb.de>
#' @seealso \code{\link{Synapter}} for details about the cleavage procedure.
#' @examples
#' \dontrun{
#' createUniquePeptideDbRds("uniprot.fasta", "uniprot.fasta.rds")
#' }
createUniquePeptideDbRds <- function(fastaFile,
                                     outputFile = paste0(fastaFile, ".rds"),
                                     missedCleavages = 0,
                                     IisL = FALSE,
                                     verbose = TRUE) {
  if (!file.exists(fastaFile)) {
    stop("File ", sQuote(fastaFile), " does not exists!")
  }

  if (tolower(file_ext(outputFile)) != "rds") {
    stop("outputFile must have the file extention .rds!")
  }

  peptides <- dbUniquePeptideSet(fastaFile, missedCleavages=missedCleavages,
                                 IisL=IisL, verbose=verbose)
  if (verbose) {
    message("Save unique peptides to ", sQuote(outputFile))
  }
  saveRDS(peptides, file=outputFile)
}

