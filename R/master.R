loadIdentOnly <- function(pepfile, 
                          fdr = 0.01,
                          method = "BH",
                          verbose = TRUE) {
  if (verbose)
    message("Processing ", basename(pepfile))
  x1 <- .Synapter$new()
  x1$IdentPeptideFile <- pepfile
  x1$IdentPeptideData <- read.csv(x1$IdentPeptideFile, stringsAsFactors = FALSE)
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
            if (length(object@masters) == 0) {
              cat("Empty object of class \"", class(object), "\"\n",sep = "")
            } else  {              
              cat("Object of class \"", class(object), "\"\n",sep = "")              
              f <- basename(object@pepfiles)[object@orders[[1]]]
              f <- paste0(" [", 1:length(f), "] ", f)
              o1 <- paste(1:length(object@pepfiles), collapse = " ")
              o2 <- paste(c(2:length(object@pepfiles), 1), collapse = " ")
              cat(" 1st Master [", o1, "] has", nrow(object@masters[[1]]), "peptides \n")
              cat(" 2nd Master [", o2, "] has", nrow(object@masters[[2]]), "peptides \n")
              cat(paste(f, collapse = "\n"), "\n")
            }
            invisible(NULL)
          })



writeMasterPeptides <- function(x, file, ...) {
  filename <- basename(file)
  path <- dirname(file)
  file1 <- paste0("01-", filename)
  file2 <- paste0("02-", filename)
  write.csv(x@masters[[1]], file = file.path(path, file1), ...)
  write.csv(x@masters[[2]], file = file.path(path, file2), ...)
}


##' This function takes all possible combination of \code{pepfiles}
##' of length greater or equal than 2 and computes the number of
##' estimated incorrect petides, the number of unique peptides,
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
##' @param verbose Should progress messages be printed?
##' @param PLGS If \code{TRUE} (default) try to emulate PLGS' peptide cleavage
##' rules. Otherwise use the default rules from the \code{cleaver} package. See
##' \code{\link{Synapter}} for references.
##' @param IisL If \code{TRUE} Isoleucin and Leucin are treated as
##' equal. In this case sequences like "ABCI", "ABCL" are removed because they
##' are not unqiue. If \code{FALSE} "ABCI" and "ABCL" are reported as unique.
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
                              PLGS = TRUE,
                              IisL,
                              verbose = TRUE) {
    pepfile <- unlist(pepfiles)
    m <- length(pepfiles)
    cmbs <- lapply(2:m, function(.k) combn(m, .k, simplify = FALSE))
    cmbs <- Reduce(c, cmbs)
    k <- length(cmbs)
    if (verbose) {
        message(m, " Peptide files available - ", k, " combinations")
    }
    hdmseList <- lapply(pepfiles,
                        loadIdentOnly,
                        fdr = fdr,
                        verbose = verbose)
    if (verbose)
        message("Generating unique proteotypic peptides...")
    proteotyptic <- dbUniquePeptideSet(fastafile,
                                       missedCleavages = missedCleavages,
                                       PLGS = PLGS,
                                       IisL = IisL,
                                       verbose = FALSE)
    if (verbose)
        message("Calculating...")
    uniquePepList <-
        lapply(hdmseList, function(.x) .x$IdentPeptideData$peptide.seq)
    uniquePeps <- unique(unlist(uniquePepList))  
    n <- length(uniquePeps)  
    contMat <- matrix(0L, nrow = n, ncol = m)
    rownames(contMat) <- uniquePeps
    colnames(contMat) <- 1:m 
    for (i in 1:ncol(contMat)) {
        x <- uniquePeps[uniquePeps %in% uniquePepList[[i]]]
        contMat[x, i] <- 1L
    }
    nbIncorrect <- round(colSums(contMat) * fdr, 0)
    ans <- sapply(cmbs, function(.cmb) {
        x <- unique(unlist(uniquePepList[.cmb]))
        .nbProteotypic <- sum(x %in% proteotyptic)
        mat <- contMat[, .cmb]
        .nbIncorrect <- sum(nbIncorrect[.cmb])
        .nbUnique <- sum(rowSums(contMat[, .cmb]) != 0)
        .fdr <- .nbIncorrect/.nbUnique
        c(incorrect = .nbIncorrect,
          unique = .nbUnique,
          proteotypic = .nbProteotypic,
          fdr = .fdr)
    })  
    ans <- data.frame(t(ans))
    ans$combinationAsText <- sapply(cmbs, paste, collapse = ".")
    ans$combination <- cmbs
    ans$nbSample <- sapply(cmbs, length)
    ifelse(proteotypic,
           b <- which.max(ans[ans$fdr < masterFdr, "proteotypic"]),
           b <- which.max(ans[ans$fdr < masterFdr, "unique"]))
    best <- ans[which(ans$fdr < masterFdr)[b], ]
    ret <- new("MasterFdrResults",
               all = ans,
               best = best,
               files = pepfiles,
               masterFdr = masterFdr)
    return(ret)
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
##' the unique peptides from the filtered original peptide files.
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
##' extension, to be recognised (and loded) as a \code{R} object.
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
##' @title Merges final peptide files
##' @param pepfiles A \code{list} of peptide final peptide file names
##' to be merged.
##' @param fdr A \code{numeric} indicating the preptide false
##' discovery  rate limit.
##' @param method A \code{character} indicating the p-value adjustment to
##' be used. One of \code{BH} (default), \code{Bonferroni} or \code{qval}.
##' @param span A \code{numeric} with the loess span parameter value
##' to be used for retention time modelling.
##' @param verbose A \code{logical} indicating information should be
##' printed out.
##' @return An instance of class \code{"\linkS4class{MasterPeptides}"}.
##' @author Laurent Gatto
##' @aliases writeMasterPeptides
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
                       fdr = 0.01,
                       ## fpr
                       method = c("BH", "Bonferroni", "qval"),
                       span = 0.05,
                       verbose = TRUE) {
  method <- match.arg(method)
  n <- length(pepfiles)
  hdmseList <- lapply(pepfiles, loadIdentOnly,
                      fdr = fdr,
                      method = method,
                      verbose = verbose)
  npeps <- sapply(hdmseList, function(.x) nrow(.x$IdentPeptideData))
  o1 <- order(npeps, decreasing = TRUE)
  o2 <- c(o1[-1], o1[1])
  orders <- list(o1, o2)
  mergedList <- 
    lapply(orders, function(o) {
      ._hdmseList <- hdmseList[o]
      master <- ._hdmseList[[1]]
      merged <- master$IdentPeptideData
      if (verbose) {
        message("Master: ", basename(master$IdentPeptideFile),
                " (", nrow(master$IdentPeptideData), " peptides)")
        message(" |--- (1) Merged: ", nrow(merged), " features.")
      }
      for (i in 2:n) {
        slave <- ._hdmseList[[i]]
        if (verbose)                  
          message(" +- Merging master and ", basename(slave$IdentPeptideFile),
                  " (", nrow(slave$IdentPeptideData), " peptides)")
        mergedPeptideData <- merge(master$IdentPeptideData,
                                   slave$IdentPeptideData,
                                   by.x = "peptide.seq",
                                   by.y = "peptide.seq",
                                   suffixes = c(".master", ".slave"))
        RtModel <- loess(precursor.retT.master ~ precursor.retT.slave, data = mergedPeptideData, span = span)
        predictedslaveRtime <- predict(RtModel,
                                       newdata = data.frame(precursor.retT.slave = slave$IdentPeptideData$precursor.retT))
        slave$IdentPeptideData$precursor.retT <- predictedslaveRtime
        selPepsToAdd <- !(slave$IdentPeptideData$peptide.seq %in% merged$peptide.seq)
        merged <- rbind(merged, slave$IdentPeptideData[selPepsToAdd, ])
        if (verbose) {
          if (i == n) {
            message(" \\--- (", i, ") Merged: ", nrow(merged), " features.")
          } else {
            message(" |--- (", i, ") Merged: ", nrow(merged), " features.")
          }
        }
      }
      sel <- is.na(merged$precursor.retT)
      merged <- merged[!sel, ]    
      return(merged)
    })
  master <- new("MasterPeptides",
                masters = mergedList,
                pepfiles = pepfiles,
                fdr = fdr,
                method = method,
                orders = orders)
  return(master)
}

