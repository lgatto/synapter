loadIdentOnly <- function(pepfile,
                          fdr = 0.01,
                          method = "BH",
                          verbose = interactive()) {
  if (verbose)
    message("Processing ", basename(pepfile))
  x1 <- .Synapter$new()
  x1$IdentPeptideFile <- pepfile
  x1$IdentPeptideData <- readFinalPeptides(x1$IdentPeptideFile,
                                           verbose = verbose)
  ## x1$DbFastaFile <- fastafile
  x1$IdentPeptideData$errorppm <-
    error.ppm(obs = x1$IdentPeptideData$precursor.mhp,
              theo = x1$IdentPeptideData$peptide.mhp)
  x1$IdentPeptideData$protein.falsePositiveRate <-
    x1$IdentPeptideData$protein.falsePositiveRate / 100
  x1$filterMatchType()
  x1$addIdentIdStats()
  x1$setPepScoreFdr(fdr)
  x1$filterIdentPepScore(fdrMethod = method)
  ## x1$setProtFpr(fpr)
  ## x1$filterIdentProtFpr()
  return(x1)
}


setMethod("show", "MasterPeptides",
          function(object) {
            if (!length(object@masters)) {
              cat("Empty object of class \"", class(object), "\"\n", sep = "")
            } else {
              cat("Object of class \"", class(object), "\"\n", sep = "")
              f <- basename(object@pepfiles)[object@orders[[1]]]
              f <- paste0(" [", 1:length(f), "] ", f)
              o1 <- paste(1:length(object@pepfiles), collapse = " ")
              o2 <- paste(c(2:length(object@pepfiles), 1), collapse = " ")
              cat(" 1st Master [", o1, "] has", nrow(object@masters[[1]]), "peptides \n")
              cat(" 2nd Master [", o2, "] has", nrow(object@masters[[2]]), "peptides \n")
              cat(paste(f, collapse = "\n"), "\n")
            }

              if (.hasSlot(object, "fragmentlibrary")) {
                  if (!length(object@fragmentlibrary)) {
                      cat("\n No fragment library.\n")
                  } else {
                      cat("\n Fragment library contains fragments for ",
                          length(object@fragmentlibrary), " peptides.\n", sep="")
                      cat(paste0(" [", seq_along(object@fragmentfiles), "] ",
                                 basename(object@fragmentfiles), "\n"), sep="")
                  }
              }
            invisible(NULL)
          })

setMethod("writeMasterPeptides", c("MasterPeptides", "character"),
          function(x, file, ...) {
  filename <- sprintf("%02i-%s", seq_along(x@masters), basename(file))
  path <- dirname(file)
  for (i in seq(along=x@masters)) {
    write.csv(x@masters[[i]], file = file.path(path, filename[i]), row.names = FALSE, ...)
  }
})

setMethod("writeFragmentLibrary", c("MasterPeptides", "character"),
          function(x, file, ...) {
  write.csv(x@fragments, file = file, row.names = FALSE, ...)
})

##' This function takes all possible combination of \code{pepfiles}
##' of length greater or equal than 2 and computes the number of
##' estimated incorrect peptides, the number of unique peptides,
##' the number of unique protetypic peptides and the
##' false discovery rate after merging for each combination.
##' The best combination has an fdr lower than \code{masterFdr}
##' and the highest number of unique (proteotypic) peptides.
##'
##' The false discovery rate for the master (merged) file is calcualted
##' by summing the number of estimated false discoveries for each
##' individual final peptide file (number of unique peptides in that file
##' multiplied by \code{fdr}) divided by the total number of unique
##' peptides for that specific combination.
##'
##' The function returns an instance of the class
##' \code{"\linkS4class{MasterFdrResults}"}.
##'
##' @title Computes FDR for all possible final peptide combinations
##' @param pepfiles A \code{list} of \code{vector} of final peptide
##' filenames.
##' @param fastafile A \code{character} with the fasta filename.
##' @param masterFdr A \code{numeric} indicating the maximum merged
##' false discovery to be allowed.
##' @param fdr Peptide FDR level for individual peptide files filtering.
##' @param proteotypic Logical. Should number proteotypic peptides be
##' used to choose best combination and plot results or total number
##' of unique peptides.
##' @param missedCleavages Number of maximal missed cleavage sites. Default is 0.
##' @param IisL If \code{TRUE} Isoleucin and Leucin are treated as identical.
##' In this case sequences like "ABCI", "ABCL" are removed because they
##' are not unqiue. If \code{FALSE} (default) "ABCI" and "ABCL" are reported as
##' unique.
##' @param maxFileComb A \code{numeric} to limit the accepted file
##' combinations to reduce computation time. Default is
##' \code{length(pepfiles)} meaning no limit set.
##' @param verbose Should progress messages be printed?
##' @return An instance of class \code{"\linkS4class{MasterFdrResults}"}.
##' See details above.
##' @author Laurent Gatto
##' @references Bond N. J., Shliaha P.V., Lilley K.S. and Gatto L.,
##' (2013) J. Prot. Research.
##' @seealso The \code{\link{makeMaster}} function to combine
##' the peptide data as suggested by \code{estimateMasterFdr} into
##' one single \emph{master} peptide file.
##'
##' The vignette, accessible with \code{synapterGuide()} illustrates a
##' complete pipeline using \code{estimateMasterFdr} and
##' \code{makeMaster}.
estimateMasterFdr <- function(pepfiles,
                              fastafile,
                              masterFdr = 0.025,
                              fdr = 0.01,
                              proteotypic = TRUE,
                              missedCleavages = 0,
                              IisL = FALSE,
                              maxFileComb = length(pepfiles),
                              verbose = interactive()) {
  hdmseList <- lapply(pepfiles,
                      loadIdentOnly,
                      fdr = fdr,
                      verbose = verbose)
  if (verbose) {
    message("Generating unique proteotypic peptides...")
  }
  proteotyptic <- dbUniquePeptideSet(fastafile,
                                     missedCleavages = missedCleavages,
                                     IisL = IisL,
                                     verbose = FALSE)

  cmbs <- .calculatePeptideFileCombinations(length(pepfiles),
                                            maxFileComb = maxFileComb,
                                            verbose = verbose)

  uniquePeptidesList <- lapply(hdmseList, function(x)
                               unique(x$IdentPeptideData$peptide.seq))

  if (verbose) {
    message("Calculating...")
  }

  ans <- .masterFdrSummary(cmbs, proteotyptic, uniquePeptidesList, fdr)

  b <- which.max(ans[ans$fdr < masterFdr, ifelse(proteotypic,
                                                 "proteotypic",
                                                 "unique")])
  new("MasterFdrResults",
      all = ans,
      best = ans[which(ans$fdr < masterFdr)[b], ],
      files = pepfiles,
      masterFdr = masterFdr)
}

#' Calculate the possible combinations of pepfiles for master creation.
#' @param n number of files
#' @param maxFileComb maximal accepted file combinations
#' @param verbose verbose output?
#' @return list, with numeric vectors containing all combinations
#' @noRd
.calculatePeptideFileCombinations <- function(n, maxFileComb=n,
                                              verbose=interactive()) {
  cmbs <- lapply(2:maxFileComb, function(k)combn(n, k, simplify = FALSE))
  cmbs <- Reduce(c, cmbs)

  if (verbose) {
    msg <- paste(n, "peptide files available -", length(cmbs), "combinations")
    if (maxFileComb < n) {
      msg <- paste0(msg, " (limited by maxFileComb=", maxFileComb, ")")
    }
    message(msg)
  }
  cmbs
}

#' Create bookholder matrix for unique peptides in specific ident files
#' @param x list of unique peptides per file (one element per file)
#' @return matrix, with unique peptides in rows and files in columns
#' @noRd
.peptideMatrix <- function(x) {
  uniquePeptides <- unique(unlist(x))

  l <- lapply(x, "%in%", x=uniquePeptides)
  m <- do.call(cbind, l)
  rownames(m) <- uniquePeptides
  colnames(m) <- seq_along(x)
  storage.mode(m) <- "integer"
  m
}

#' Create summary for combinations needed for MasterFdrResults
#' @param cmbs list of combinations
#' @param proteotyptic character, proteotypic peptides
#' @param uniquePeptidesList list of unique peptides per file
#' (one element per file)
#' @param fdr fdr
#' @return data.frame
#' @noRd
.masterFdrSummary <- function(cmbs, proteotyptic, uniquePeptidesList, fdr) {
  contMat <- .peptideMatrix(uniquePeptidesList)

  nIncorrect <- round(colSums(contMat) * fdr, 0L)

  ans <- lapply(cmbs, function(cmb) {
    uniquePeptides <- unique(unlist(uniquePeptidesList[cmb]))
    c(incorrect = round(sum(contMat[, cmb]) * fdr, 0L),
      unique = sum(rowSums(contMat[, cmb]) != 0L),
      proteotypic = sum(uniquePeptides %in% proteotyptic))
  })

  ans <- as.data.frame(do.call(rbind, ans))
  ans$fdr <- ans$incorrect / ans$unique
  ans$combination <- cmbs
  ans$nbSample <- lengths(cmbs)
  ans
}

setMethod("show", "MasterFdrResults",
          function(object) {
            cat(length(object@files), "files -", nrow(object@all), "combinations\n")
            cat("Best combination:", object@best$combination[[1]], "\n")
            cat(" -", object@best$proteotypic, "proteotypic peptides\n")
            cat(" -", object@best$unique, "unique peptides\n")
            cat(" -", round(object@best$fdr,3), "FDR\n")
          })

setMethod("plot", c("MasterFdrResults", "missing"),
          function(x,y, proteotypic = TRUE) {
            .x <- x@all$fdr
            ifelse(proteotypic,
                   { .y <- x@all[, "proteotypic"]; .ylab <- "Number of unique proteotypic peptides"},
                   { .y <- x@all[, "unique"]; .ylab <- "Number of unique peptides"})
            plot(.y ~ .x, type = "n",
                 xlab = "Total FDR", ylab = .ylab)
            text(.x, .y,
                 as.character(x@all$nbSample),
                 cex = .6,
                 col = ifelse(x@all$fdr < x@masterFdr,
                   "#0000FF80",
                   "#00000080"))
            abline(v = x@masterFdr, lty = "dotted")
            ifelse(proteotypic,
                   .x <- x@best[, c("fdr", "proteotypic")],
                   .x <- x@best[, c("fdr", "unique")])
            points(.x, pch = 19, cex = 2, col = "#FF000020")
            points(.x, pch = 19, cex = 2, col = "#FF000020")
          })



setMethod("bestComb", "MasterFdrResults",
          function(object) object@best$combination[[1]])

setMethod("fileNames", "MasterFdrResults",
          function(object) object@files)

setMethod("masterFdr", "MasterFdrResults",
          function(object) object@masterFdr)

setMethod("allComb", "MasterFdrResults",
          function(object) object@all)





##' This function combines a list of peptide final peptide files into
##' one single \emph{master} file that is obtained by merging
##' the unique peptides from the filtered original peptide files. \cr
##' Additionally it can combine multiple final fragment files into a fragment
##' library.
##'
##' The merging process is as follows:
##' \enumerate{
##' \item Each individual peptide final peptide file is filtered to retain
##' (i) non-duplicated unique tryptic peptides, (ii) peptides with a
##' false discovery rate <= \code{fdr} and (iii) proteins with a false
##' positive rate <= \code{fpr}.
##'
##' \item The filtered peptide files are ordered (1) according to their total
##' number of peptides (for example [P1, P2, P3]) and (2) as before with the first
##' item is positioned last ([P2, P3, P1] in the previous example).
##' The peptide data are then combined in pairs in these respective orders.
##' The first one is called the \emph{master} file.
##'
##' \item For each (master, slave) pair, the slave peptide file
##' retention times are modelled according to the (original)
##' master's retention times and slave peptides, not yet present in the
##' master file are added to the master file.
##'
##' \item The final \emph{master} datasets, containing their own peptides and
##' the respective slave specific retention time adjusted peptides are returned
##' as a \code{MasterPeptides} instance.
##' }
##'
##' The resulting \code{MasterPeptides} instance can be further used
##' for a complete master vs. peptides/Pep3D analysis, as described in
##' \code{\link{Synapter}}, \code{\link{synergise}} or using the GUI
##' (\code{\link{synapterGUI}}). To do so, it must be serialised (using the
##' \code{saveRDS} function) with a \code{.rds} file
##' extension, to be recognised (and loaded) as a \code{R} object.
##'
##' When several quantitation (or identification) files are combined as a master set
##' to be mapped back against the inidividual final peptide files,
##' the second master [P2, P3, P1] is used when analysing the peptide data
##' that was first selected in the master generation (P1 above).
##' This is to avoid aligning two identical sets of peptides (those of P1)
##' and thus not being able to generate a valid retention time model.
##' This is detected automatically for the user.
##'
##' The two master peptides dataframes can be exported to disk as
##' two \code{csv} files with \code{writeMasterPeptides}. The
##' \code{MasterPeptides} object returned by \code{makeMaster} can be
##' saved to disk (with \code{save} or \code{saveRDS}) and later reloaded
##' (with \code{load} or \code{readRDS}) for further analysis.
##'
##' The fragment library generation works as follows:
##' \enumerate{
##' \item Each individual final fragment file is imported and only peptides
##' present in the \emph{master} dataset are used.
##'
##' \item The fragments are combined based on their precursor ions.
##'
##' \item The intensities of identical fragments (seen in different runs)
##' is summed and divided by the summed precursor intensity (of the same
##' peptide in different runs).
##'
##' \item Afterwards the intensities are normalized to the average precursor
##' intensity of the different runs.
##'
##' \item Finally a \code{\linkS4class{MSnExp}} object is created.
##' }
##'
##' The fragment library dataframe can be exported to disk as
##' \code{csv} file with \code{writeFragmentLibrary}.
##'
##' @title Merges final peptide files
##' @param pepfiles A \code{character} vector of final peptide file
##' names to be merged.
##' @param fragmentfiles A \code{character} vector of final fragment file
##' names to be combined into an fragment library. These files should be
##' from the same runs as the final peptide files used in \code{pepfiles}.
##' @param fdr A \code{numeric} indicating the peptide false
##' discovery  rate limit.
##' @param method A \code{character} indicating the p-value adjustment to
##' be used. One of \code{BH} (default), \code{Bonferroni} or \code{qval}.
##' @param span.rt A \code{numeric} with the loess span parameter value
##' to be used for retention time modelling.
##' @param span.int A \code{numeric} with the loess span parameter value
##' to be used for intensity modelling.
##' @param maxDeltaRt A \code{double} value that sets a maximum limit for
##' the retention time deviaton between master and slave run to be included
##  in the retention time modelling.
##' @param removeNeutralLoss A \code{logical}, if \code{TRUE} peptides with
##' neutral loss are removed from the fragment library.
##' @param removePrecursor A \code{logical}, if \code{TRUE} precursor ions are
##' removed from the fragment spectra.
##' @param tolerance A \code{double} value that determines the tolerance used
##' to look for the precursor ions.
##' @param verbose A \code{logical} indicating if information should be
##' printed out.
##' @return An instance of class \code{"\linkS4class{MasterPeptides}"}.
##' @author Laurent Gatto, Sebastian Gibb
##' @aliases writeMasterPeptides
##' writeMasterPeptides,MasterPeptides,character-method
##' writeFragmentLibrary writeFragmentLibrary,MasterPeptides,character-method
##' @references Shliaha P.V., Bond N. J., Lilley K.S. and Gatto L., in prep.
##' @seealso See the \code{\link{Synapter}} class manual page for
##' detailed information on filtering and modelling and the general
##' algorithm implemented in the \code{synapter} package.
##'
##' The \code{\link{estimateMasterFdr}} function allows to control
##' false dicovery rate when combining several peptide files while
##' maximising the number of identifications and suggest which
##' combination of peptide files to use.
##'
##' The vignette, accessible with \code{synapterGuide()}
##' illustrates a complete pipeline using \code{estimateMasterFdr} and
##' \code{makeMaster}.
makeMaster <- function(pepfiles,
                       fragmentfiles,
                       fdr = 0.01,
                       ## fpr
                       method = c("BH", "Bonferroni", "qval"),
                       span.rt = 0.05,
                       span.int = 0.05,
                       maxDeltaRt = Inf,
                       removeNeutralLoss = TRUE,
                       removePrecursor = TRUE,
                       tolerance = 25e-6,
                       verbose = interactive()) {
  method <- match.arg(method)
  hdmseList <- lapply(pepfiles, loadIdentOnly,
                      fdr = fdr,
                      method = method,
                      verbose = verbose)
  hdmseList <- lapply(hdmseList, .filterDuplicatedPeptideSequences)
  orders <- .orderForMasterModels(hdmseList)
  mergedList <- lapply(orders,
                       function(o).mergeMaster(hdmseList[o],
                                               span.rt = span.rt,
                                               span.int = span.int,
                                               maxDeltaRt = maxDeltaRt,
                                               verbose = verbose))

  mergedList <- .regeneratePrecursorLeId(mergedList, verbose = verbose)

  ## create fragment library
  if (!missing(fragmentfiles)) {
    if (verbose) {
      message("Creating fragment library.")
    }

    fragments <- .createFragmentLibrary(master = mergedList[[1L]],
                                        files = fragmentfiles,
                                        removeNeutralLoss = removeNeutralLoss,
                                        verbose = verbose)
    lib <- .fragments2spectra(df = NULL, # df is not needed for storeAll=TRUE, because we don't use
                                         # ID intersection
                              fragments = fragments,
                              file = fragmentfiles,
                              storeAll = TRUE,
                              removePrecursor = removePrecursor,
                              tolerance = tolerance,
                              verbose = verbose)
  } else {
    lib <- new("MSnExp")
    fragmentfiles <- character()
    fragments <- data.frame()
  }

  new("MasterPeptides",
      masters = mergedList,
      pepfiles = pepfiles,
      fdr = fdr,
      method = method,
      orders = orders,
      fragmentfiles = fragmentfiles,
      fragments = fragments,
      fragmentlibrary = lib)
}

#' Remove rows with douplicated peptide.seq entries.
#' @param x hdmse list generated by lapply(..., loadIdentOnly)
#' @noRd
.filterDuplicatedPeptideSequences <- function(x) {
  x$IdentPeptideData <-
    x$IdentPeptideData[!duplicated(x$IdentPeptideData$peptide.seq), ]
  x
}

#' Use the data set with the highest number of peptides for model creation
#' @param x hdmse list generated by lapply(..., loadIdentOnly)
#' @noRd
.orderForMasterModels <- function(x) {
  n <- vapply(x, function(xx)nrow(xx$IdentPeptideData), double(1L))
  o <- order(n, decreasing = TRUE)
  list(o, c(o[-1L], o[1L]))
}

#' Merge a peptide data (master and a single slave)
#' @param master master df
#' @param slave df that should be merged
#' @param maxDeltaRt maximal allowed deviation in retTime; see issue #107
#' @param verbose verbose output?
#' @noRd
.mergePeptideData <- function(master, slave, maxDeltaRt = Inf,
                              verbose = interactive()) {
  m <- merge(master, slave, by.x = "peptide.seq", by.y = "peptide.seq",
             suffixes = c(".master", ".slave"))
  m$deltaRt <- m$precursor.retT.master - m$precursor.retT.slave

  m$intenRatio <- log2(m$precursor.inten.master / m$precursor.inten.slave)

  keep <- abs(m$deltaRt) < maxDeltaRt

  if (verbose && sum(!keep)) {
    message(" |     Ignoring ", sum(!keep), " features because ",
            "deltaRt > ", maxDeltaRt,
            " (using ", sum(keep), " for rt model).")
  }

  m[keep, ]
}

#' Merge a master and its corresponding slaves.
#' @param x hdmse list generated by lapply(..., loadIdentOnly)
#' @param span.rt loess span for retention time model
#' @param span.int loess span for intensity model
#' @param maxDeltaRt maximal allowed deviation in retTime; see issue #107
#' @param verbose verbose output?
#' @noRd
.mergeMaster <- function(x, span.rt = 0.05, span.int = 0.05, maxDeltaRt = Inf,
                         verbose = interactive()) {
  if (length(x) < 2L) {
    stop("To create a master at least two identification data are needed.")
  }

  master <- x[[1L]]
  slaves <- x[-1L]
  merged <- master$IdentPeptideData

  if (verbose) {
    message("Master: ", basename(master$IdentPeptideFile),
            " (", nrow(master$IdentPeptideData), " peptides)")
  }

  for (i in seq(along=slaves)) {

    if (verbose) {
      message(" +- Merging master and ",
              basename(slaves[[i]]$IdentPeptideFile),
              " (", nrow(slaves[[i]]$IdentPeptideData), " peptides)")
    }

    mergedPeptideData <- .mergePeptideData(master$IdentPeptideData,
                                           slaves[[i]]$IdentPeptideData,
                                           maxDeltaRt = maxDeltaRt,
                                           verbose = verbose)
    rtModel <- loessModel(mergedPeptideData$precursor.retT.master,
                          mergedPeptideData$deltaRt, span=span.rt)

    slaves[[i]]$IdentPeptideData$precursor.retT <-
      slaves[[i]]$IdentPeptideData$precursor.retT +
      predict(rtModel, slaves[[i]]$IdentPeptideData$precursor.retT)

    intModel <- loessModel(mergedPeptideData$precursor.retT.slave,
                           mergedPeptideData$intenRatio, span=span.int)

    slaves[[i]]$IdentPeptideData$precursor.inten <-
      slaves[[i]]$IdentPeptideData$precursor.inten *
      2L^(predict(intModel, slaves[[i]]$IdentPeptideData$precursor.retT))

    isMissing <- !(slaves[[i]]$IdentPeptideData$peptide.seq %in%
                   merged$peptide.seq)
    isNotNa <- !is.na(slaves[[i]]$IdentPeptideData$precursor.retT)

    merged <- rbind(merged,
                    slaves[[i]]$IdentPeptideData[isMissing & isNotNa, ])

    if (verbose) {
      sep <- if (i == length(slaves)) { "\\" } else { "|" }
      message(" ", sep, "--- (", i, ") Merged: ", nrow(merged), " features.")
    }
  }
  rownames(merged) <- NULL
  merged
}

#' Create new values for the precursor.leID column.
#' @param x list generated by lapply(..., .mergeMaster)
#' @param verbose verbose output?
#' @noRd
.regeneratePrecursorLeId <- function(x, verbose = interactive()) {
  if (verbose) {
    message("Regenerate precursor.leIDs.")
  }

  peptide.seqs <- unique(x[[1L]]$peptide.seq)
  precursor.leID <- setNames(seq_along(peptide.seqs), peptide.seqs)

  lapply(x, function(xx) {
    xx$precursor.leID <- precursor.leID[xx$peptide.seq]
    xx
  })
}
