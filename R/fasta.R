dbUniquePeptideSet <- function(file, missedCleavages = 0, PLGS = TRUE,
                               IisL = TRUE, verbose = TRUE) {
  if (tolower(getExtension(file)) == "rds") {
    peptides <- readRDS(file)

    noeffectMsg <- paste0("\nYour current setting has not any effect! ",
                          "Recreate your RDS file if you want to change ",
                          "these settings.")

    if (attr(peptides, "PLGS") != PLGS) {
      warning("The RDS file was created with PLGS=", attr(peptides, "PLGS"),
              noeffectMsg, immediate.=TRUE)
    }
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
                                    PLGS=PLGS, verbose=verbose)
  }

  if (verbose) {
      topics <- c("Peptide database file: ",
                  "Allowed mis-cleavages: ",
                  "Used cleavage rule: ",
                  "I/L treatment: ")
      values <- c(file,
                  attr(peptides, "missedCleavages"),
                  ifelse(attr(peptides, "PLGS"), "PLGS", "cleaver"),
                  ifelse(attr(peptides, "IisL"), "I == L", "I != L"))

      topics <- format(topics, justify="left")

      message(paste0(topics, values, collapse="\n"))
  }

  return(peptides)
}

.dbUniquePeptideSet <- function(fastafile, missedCleavages = 0, PLGS = TRUE,
                                IisL = TRUE, verbose = TRUE) {

    missedCleavages <- max(missedCleavages)
    stopifnot(missedCleavages >= 0)

    proteins <- as.character(readAAStringSet(fastafile))

    if (PLGS) {
      ## we modify the standard "cleaver" rules to emulate PLGS' behaviour
      ## according to the discussion on:
      ## https://github.com/sgibb/cleaver/issues/4 and
      ## https://github.com/lgatto/synapter/issues/53
      ##
      ## 1. We cleave after every K or R that is *not* followed by P (default
      ##    trypsin rule), or D, E, K, R (PLGS rule)
      ##    [perfect cleavage, no mis-cleavages allowed].
      ## 2. Subsequently we count the consecutive K or R in each peptide and
      ##    consider these numbers as maximal possible mis-cleavages to get
      ##    every peptide of interest.
      ## 3. We compare the numbers from 2 with the user-defined threshold
      ##    "missedCleavages" and only allow 0 to
      ##    min(number from 2, user-defined threshold).
      ## 4. We cleave after every K or R that is *not* followed by P (default
      ##    trypsin rule) and allow maximal n mis-cleavages. n is determined
      ##    in 3.
      ##
      ## See also:
      ##
      ## Glatter, Timo, et al. Large-scale quantitative assessment of
      ## different in-solution protein digestion protocols reveals superior
      ## cleavage efficiency of tandem Lys-C/trypsin proteolysis over trypsin
      ## digestion.
      ## Journal of proteome research 11.11 (2012): 5145-5156.
      ## http://dx.doi.org/10.1021/pr300273g
      ##
      ## Rodriguez, Jesse, et al. Does trypsin cut before proline?.
      ## Journal of proteome research 7.01 (2007): 300-305.
      ## http://dx.doi.org/10.1021/pr0705035
      ##
      ## Brownridge, Philip, and Robert J. Beynon. The importance of the digest:
      ## proteolysis and absolute quantification in proteomics.
      ## Methods 54.4 (2011): 351-360.
      ## http://dx.doi.org/10.1016/j.ymeth.2011.05.005
      ##

      ## 1) "perfect split" without any mis-cleavage
      peptides <- unlist(cleave(proteins, custom="[K|R](?=[^PDEKR])",
                                missedCleavages=0), use.names=FALSE)

      ## 2) guess maximal number of maximal needed mis-cleavages
      rx <- gregexpr("[KR]+", peptides)
      maxMissedCleavages <- sapply(rx, function(x)max(attr(x, "match.length")))

      ## We reduce maxMissedCleavages by one because we assume that all peptides
      ## that are cleaved in 1) have K, R at their right end (and we don't want to
      ## cleave something like "AAAK" again).
      maxMissedCleavages <- maxMissedCleavages-1

      ## 3) determine maximal allowed mis-cleavages
      maxMissedCleavages <- pmin(maxMissedCleavages, missedCleavages)

      ## only cleave peptides that contain at least one K or R
      toCleave <- maxMissedCleavages > 0

      ## 4) final cleaving
      peptides2 <- unlist(mapply(function(x, m) {
        cleave(x, custom="[K|R](?=[^P])", missedCleavages=(0:m))
      }, x=peptides[toCleave], m=maxMissedCleavages[toCleave],
      SIMPLIFY=FALSE, USE.NAMES=FALSE), use.names=FALSE)

      peptides <- c(peptides, peptides2)
    } else {
      ## using default cleaver cleavages rules for trypsin, defined:
      ## http://web.expasy.org/peptide_cutter/peptidecutter_enzymes.html#Tryps
      peptides <- unlist(cleave(proteins, enzym="trypsin",
                                missedCleavages=0:missedCleavages),
                         use.names=FALSE)
    }

    ## remove non-unique peptides
    ## Note that this does however not take pre/post-fix peptides into account
    ## like "DLIELTESLFR" and "HNPEFTMMELYMAYADYHDLIELTESLFR",
    ##      "YYGYTGAFR" and "EGYYGYTGAFR"
    ## where the first ones are NOT unique!
    ## Such cases are handled by filtering duplicates in the peptide data
    if (IisL) {
      .peptides <- gsub(pattern="[IL]", replacement="-", x=peptides)
      uniqueIdx <- !(duplicated(.peptides) |
                            duplicated(.peptides, fromLast=TRUE))
    } else {
      uniqueIdx <- !(duplicated(peptides) |
                            duplicated(peptides, fromLast=TRUE))
    }

    upeptides <- peptides[uniqueIdx]

    attr(upeptides, "missedCleavages") <- missedCleavages
    attr(upeptides, "PLGS") <- PLGS
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
#' @param PLGS If \code{TRUE} (default) try to emulate PLGS' peptide cleavage
#' rules. Otherwise use the default rules from the \code{cleaver} package. See
#' \code{\link{Synapter}} for references.
#' @param IisL If \code{TRUE} (default) Isoleucin and Leucin are treated as
#' equal. In this case sequences like "ABCI", "ABCL" are removed because they
#' are not unqiue. If \code{FALSE} "ABCI" and "ABCL" are reported as unique.
#' @param verbose If \code{TRUE} a verbose output is provied.
#' @author Sebastian Gibb <mail@@sebastiangibb.de>
#' @seealso \code{\link{Synapter}} for details about the cleavage procedure.
#' @examples
#' \dontrun{
#' createUniquePeptideDbRds("uniprot.fasta", "uniprot.fasta.rds", PLGS=TRUE)
#' }
createUniquePeptideDbRds <- function(fastaFile,
                                     outputFile = paste0(fastaFile, ".rds"),
                                     missedCleavages = 0,
                                     PLGS = TRUE,
                                     IisL = TRUE,
                                     verbose = TRUE) {
  if (!file.exists(fastaFile)) {
    stop("File ", sQuote(fastaFile), " does not exists!")
  }

  if (tolower(getExtension(outputFile)) != "rds") {
    stop("outputFile must have the file extention .rds!")
  }

  peptides <- dbUniquePeptideSet(fastaFile, missedCleavages=missedCleavages,
                                 PLGS=PLGS, IisL=IisL, verbose=verbose)
  if (verbose) {
    message("Save unique peptides to ", sQuote(outputFile))
  }
  saveRDS(peptides, file=outputFile)
}

