## Constructor
Synapter <- function(filenames, master = FALSE) {
  xx <- .Synapter$new()
  xx$setMaster(master)
  if (missing(filenames)) {
    xx$loadFiles()
  } else {
    if (!all(names(filenames) %in% c("identpeptide", "quantpeptide", "quantpep3d", "fasta", "identfragments", "quantspectra")))
      stop("File names must be provided as a named list with names 'identpeptide','quantpeptide', 'quantpep3d' and 'fasta' [optional 'identfragments' and 'quantspectra'].")
    flength <- sapply(filenames, length)
    if (!all(flength == 1))
        stop("This interface only accepts single MSe/HDMSe files as input. See '?synapterGUI' and '?synergise' for alternative interfaces.")
    ftest <- sapply(filenames, file.exists)
    if (any(!ftest))
      stop(paste(filenames[!ftest], collapse = " "),
           " do(es) not exists.")
    xx$IdentPeptideFile <- filenames$identpeptide
    xx$QuantPeptideFile <- filenames$quantpeptide
    xx$QuantPep3DFile <- filenames$quantpep3d
    xx$DbFastaFile <- filenames$fasta

    ## add optional Fragment/Spectra files
    if (length(filenames$identfragments)) {
      xx$IdentFragmentFile <- filenames$identfragments
    } else {
      xx$IdentFragmentFile <- character()
    }
    if (length(filenames$quantspectra)) {
      xx$QuantSpectrumFile <- filenames$quantspectra
    } else {
      xx$QuantSpectrumFile <- character()
    }
  }

  if (xx$Master) {
    xx$loadMasterData()
  } else {
    xx$loadData()
  }
  return(xx)
}

setMethod(show, "Synapter",
          function(object) {
            'Textual display of the instance.'
            cat("Object of class", classLabel(class(object)), "\n")
            cat("Class version", object$ClassVersion, "\n")
            cat("Package version", object$Version, "\n")
            cat("Data files:\n")
            cat(" + Quantitation pep file:",
                basename(object$QuantPeptideFile), "\n")
            cat(" + Identification pep file:",
                basename(object$IdentPeptideFile), "\n")
            cat(" + Quantitation Pep3DAMRT file:",
                basename(object$QuantPep3DFile), "\n")
            cat(" + Fasta file:",
                basename(object$DbFastaFile), "\n")
            cat("Log:\n")
            l <- length(object$SynapterLog)
            if (l > 4) {
              msg <- getLog(object)
              print(msg[1:2])
              cat("[", l-4 ,"lines ]\n")
              cat("[", l-1, "] \"",  msg[l-1], "\"\n", sep = "")
              cat("[", l, "] \"",  msg[l], "\"\n", sep = "")
            } else {
              print(object$SynapterLog)
            }
            invisible(NULL)
          })

setMethod("dim", "Synapter",
          function(x) {
            dims <- list(IdentPeptideData = dim(x$IdentPeptideData),
                         QuantPeptideData = dim(x$QuantPeptideData),
                         MergedPeptides = dim(x$MergedFeatures),
                         MatchedEMRTs = dim(x$MatchedEMRTs))
            dims
          })

setMethod(inputFiles, "Synapter",
          function(object)
          c(identpeptide = object$IdentPeptideFile,
            quantpeptide = object$QuantPeptideFile,
            quantpep3d = object$QuantPep3DFile,
            fasta = object$DbFastaFile,
            identfragments = object$IdentFragmentFile,
            quantspectra = object$QuantSpectrumFile))

setMethod(getLog, "Synapter",
          function(object) object$SynapterLog)


setMethod(mergePeptides, "Synapter",
          function(object) {
            object$mergePeptides()
          })

setMethod(modelRt, "Synapter",
          function(object, span) {
            if (missing(span))
              span <- object$LowessSpan
            if (length(span) == 0) {
              object$setLowessSpan()
              span <- object$LowessSpan
            }
            object$modelRetentionTime(span)
          })

setMethod(findEMRTs, "Synapter",
          function(object, ppm, nsd, imdiff) {
            if (!missing(ppm))
              object$setPpmError(ppm)
            if (!missing(nsd))
              object$setRtNsd(nsd)
            if (!missing(imdiff))
              object$setImDiff(imdiff)
            object$findEMRTs()
          })

setMethod(rescueEMRTs, "Synapter",
          function(object, method = c("rescue", "copy")) {
            if (!nrow(object$MatchedEMRTs)) {
              stop("You have to run ", sQuote("findEMRTs"), " first!")
            }
            object$rescueEMRTs(method)
          })

setMethod(searchGrid, "Synapter",
          function(object,
                   ppms,
                   nsds,
                   imdiffs,
                   subset,
                   n,
                   verbose = TRUE) {
            if (missing(ppms))
              ppms <- seq(2, 20, 2)
            names(ppms) <- ppms
            if (missing(nsds))
              nsds <- seq(0.5, 5, 0.5)
            names(nsds) <- nsds
            if (missing(imdiffs))
              imdiffs <- seq(0.2, 2, 0.2)
            if (!missing(n) & !missing(subset))
              stop("Use either 'n' or 'subset', not both.")
            if (missing(n) & missing(subset))
              subset <- 1
            if (!missing(subset) && (subset > 1 | subset <= 0))
              subset <- 1
            object$searchGrid(ppms = ppms,
                              nsds = nsds,
                              imdiffs = imdiffs,
                              subset = subset,
                              n = n,
                              verbose = verbose)
          })

setMethod(getGrid, "Synapter",
          function(object, digits = 3) {
            lapply(object$Grid, round, digits = digits)
          })

setMethod(getGridDetails, "Synapter",
          function(object) object$GridDetails)

setMethod(getBestGridValue, "Synapter",
          function(object) object$getBestGridValue())

setMethod(getBestGridParams, "Synapter",
          function(object) object$getBestGridParams())

setMethod(setBestGridParams, "Synapter",
          function(object, what = c("auto", "model", "total", "details")) {
            what <- match.arg(what)
            object$setBestGridParams(what)
          })

setMethod(setPepScoreFdr, "Synapter",
          function(object, fdr = 0.01) object$setPepScoreFdr(fdr))

setMethod(getPepScoreFdr, "Synapter",
          function(object) object$PepScoreFdr)

setMethod(setProtFpr, "Synapter",
          function(object, fpr = 0.01) object$setProtFpr(fpr))

setMethod(getProtFpr, "Synapter",
          function(object) object$ProtFpr)

setMethod(setIdentPpmError, "Synapter",
          function(object, ppm = 10) object$setIdentPpmError(ppm))

setMethod(getIdentPpmError, "Synapter",
          function(object) object$IdentPpmError)

setMethod(setQuantPpmError, "Synapter",
          function(object, ppm = 10) object$setQuantPpmError(ppm))

setMethod(getQuantPpmError, "Synapter",
          function(object) object$QuantPpmError)

setMethod(setLowessSpan, "Synapter",
          function(object, span = 0.05) object$setLowessSpan(span))

setMethod(getLowessSpan, "Synapter",
          function(object) object$LowessSpan)

setMethod(setRtNsd, "Synapter",
          function(object, nsd = 2) object$setRtNsd(nsd))

setMethod(getRtNsd, "Synapter",
          function(object) object$RtNsd)

setMethod(setPpmError, "Synapter",
          function(object, ppm = 10) object$setPpmError(ppm))

setMethod(getPpmError, "Synapter",
          function(object) object$PpmError)

setMethod(setImDiff, "Synapter",
          function(object, imdiff = 0.5) object$setImDiff(imdiff))

setMethod(getImDiff, "Synapter",
          function(object) object$ImDiff)

setMethod(getPpmErrorQs, "Synapter",
          function(object,
                   qs = c(0.25, 0.5, 0.75, seq(0.9, 1, 0.01)),
                   digits = 3) {
            ## mass error quantile table.
            t <- rbind(round(getQs(object$IdentPeptideData$errorppm, qs)$y, digits),
                       round(getQs(object$QuantPeptideData$errorppm, qs)$y, digits))
            rownames(t) <- c("Ident", "Quant")
            return(t)
          })

setMethod(getRtQs, "Synapter",
          function(object,
                   qs = c(0.25, 0.5, 0.75, seq(0.9, 1, 0.01)),
                   digits = 3) {
            ## Retention time quantile table.
            diffs <- plotRetTimeDiffs(object, plot=FALSE)
            return(round(getQs(diffs, qs)$y, digits))
          })

setMethod(getPepNumbers, "Synapter",
          function(object) {
            if (object$Master) {
              quant <- unlist(by(object$.QuantPeptideScores,
                                 object$.QuantPeptideScores$peptide.matchType,
                                 function(x) table(x$protein.dataBaseType)))
              return(quant)
            } else {
              ident <- unlist(by(object$.IdentPeptideScores,
                                 object$.IdentPeptideScores$peptide.matchType,
                                 function(x) table(x$protein.dataBaseType)))
              quant <- unlist(by(object$.QuantPeptideScores,
                               object$.QuantPeptideScores$peptide.matchType,
                               function(x) table(x$protein.dataBaseType)))
            return(rbind(ident, quant))
            }
          })

setMethod(setFragmentMatchingPpmTolerance, "Synapter",
          function(object, ppm = 25)
            object$setFragmentMatchingPpmTolerance(ppm))

setMethod(getFragmentMatchingPpmTolerance, "Synapter",
          function(object) object$FragmentMatchingPpmTolerance)

setMethod(showFdrStats, "Synapter",
          function(object,
                   k = c(0.001, 0.01, 0.05, 0.1)) {
            names(k) <- as.character(k)
            ident <- list(pval = object$IdentPeptideData$pval,
                          BH = object$IdentPeptideData$BH,
                          Bonf = object$IdentPeptideData$Bonferroni,
                          qval = object$IdentPeptideData$qval)
                     quant <- list(pval = object$QuantPeptideData$pval,
                                   BH = object$QuantPeptideData$BH,
                                   Bonf = object$QuantPeptideData$Bonferroni,
                                   qval = object$QuantPeptideData$qval)
            .f <- function(x, k)
              sapply(k, function(.k) round(100 * sum(x<=.k)/length(x), 2))
            ans1 <- sapply(ident, .f, k)
            ans2 <- sapply(quant, .f, k)
            return(list(Identification = ans1,
                        Quantitation = ans2))
          })



## filter prior to merging

setMethod(filterPeptideLength, "Synapter",
          function(object, l = 7) {
              object$filterPeptideLength(l)
          })


setMethod(filterUniqueDbPeptides, "Synapter",
          function(object, missedCleavages = 0, IisL = FALSE, verbose = TRUE) {
            object$filterUniqueSeq()
            object$filterUniqueDbPeptides(object$DbFastaFile,
                                          what = c("ident", "quant"),
                                          missedCleavages = missedCleavages,
                                          IisL = IisL,
                                          verbose = verbose)
          })

setMethod(filterUniqueQuantDbPeptides, "Synapter",
          function(object, missedCleavages = 0, IisL = FALSE, verbose = TRUE) {
            object$filterUniqueQuantSeq()
            object$filterUniqueQuantDbPeptides(object$DbFastaFile,
                                               missedCleavages = missedCleavages,
                                               IisL = IisL,
                                               verbose = verbose)
          })

setMethod(filterUniqueIdentDbPeptides, "Synapter",
          function(object, missedCleavages = 0, IisL = FALSE, verbose = TRUE) {
            object$filterUniqueIdentSeq()
            object$filterUniqueIdentDbPeptides(object$DbFastaFile,
                                               missedCleavages = missedCleavages,
                                               IisL = IisL,
                                               verbose = verbose)
          })

setMethod(filterQuantPepScore, "Synapter",
          function(object, fdr,
                   method = c("BH", "Bonferroni", "qval")) {
            method <- match.arg(method)
            if (!missing(fdr))
              object$setPepScoreFdr(fdr)
            object$filterQuantPepScore(method)
          })

setMethod(filterIdentPepScore, "Synapter",
          function(object, fdr,
                   method = c("BH", "Bonferroni", "qval")) {
            method <- match.arg(method)
            if (!missing(fdr))
              object$setPepScoreFdr(fdr)
            object$filterIdentPepScore(method)
          })


setMethod(filterQuantProtFpr, "Synapter",
          function(object, fpr) {
            if (!missing(fpr))
              object$setProtFpr(fpr)
            object$filterQuantProtFpr()
          })

setMethod(filterIdentProtFpr, "Synapter",
          function(object, fpr) {
            if (!missing(fpr))
              object$setProtFpr(fpr)
            object$filterIdentProtFpr()
          })

setMethod(filterIdentPpmError, "Synapter",
          function(object, ppm) {
            if (!missing(ppm))
              object$setIdentPpmError(ppm)
            object$filterIdentPpmError(ppm)
          })

setMethod(filterQuantPpmError, "Synapter",
          function(object, ppm) {
            if (!missing(ppm))
              object$setQuantPpmError(ppm)
            object$filterQuantPpmError()
          })

# filter post merging
setMethod(filterFragments, "Synapter",
          function(object, what, minIntensity = NULL, maxNumber = NULL,
                   verbose = TRUE) {
            object$filterFragments(what = what,
                                   minIntensity = minIntensity,
                                   maxNumber = maxNumber,
                                   verbose = verbose)
          })

setMethod(filterUniqueMatches, "Synapter",
          function(object, minNumber) {
            object$filterUniqueMatches(minNumber)
          })

setMethod(filterNonUniqueMatches, "Synapter",
          function(object, minDelta) {
            object$filterNonUniqueMatches(minDelta)
          })

## Plotting
setMethod(plotPpmError, "Synapter",
          function(object,
                   what = c("Quant", "Ident", "both")) {
            what <- match.arg(what)
            switch(what,
                   Ident = qPlot(
                     object$IdentPeptideData$errorppm,
                     ylab = expression(Identification~Mass~Error~(ppm))),
                   Quant = qPlot(
                     object$QuantPeptideData$errorppm,
                     ylab = expression(Quantitation~Mass~Error~(ppm))),
                   both = {
                     par(mfrow=c(1,2))
                     qPlot(object$IdentPeptideData$errorppm,
                           ylab = expression(Identification~Mass~Error~(ppm)))
                     qPlot(object$IdentPeptideData$errorppm,
                           ylab = expression(Quantitation~Mass~Error~(ppm)))
                   })
          })

setMethod(plotRt, "Synapter",
          function(object,
                   what = c("data", "model"),
                   f = structure( ## for data
                     c(2/3, 1/2, 1/4, 1/10, 1/16, 1/25, 1/50),
                     names = c("2/3", "1/2", "1/4", "1/10", "1/16", "1/25", "1/50")),
                   nsd = c(1, 3, 5), ## for model
                   ... ) {
            what <- match.arg(what)
            if (length(nsd) > 3) {
              nsd <- nsd[1:3]
              warning("Using first 3 nsds.")
            }
            switch(what,
                   data = plotLowess(object$MergedFeatures, f = f),
                   model = plotLowess2(object$MergedFeatures,
                     object$RtModel,
                     nsd,
                     ...))
          })


setMethod(plotPepScores, "Synapter",
          function(object) {
            if (object$Master) {
              xx <- rbind(cbind(object$.QuantPeptideScores, data = "Quantitation"))
            } else {
              xx <- rbind(cbind(object$.IdentPeptideScores, data ="Identification"),
                          cbind(object$.QuantPeptideScores, data = "Quantitation"))
            }
            p <- (densityplot(~ peptide.score | peptide.matchType * data,
                              data = xx,
                              groups = protein.dataBaseType,
                              plot.points = FALSE, ref = TRUE))
            print(p)
            invisible(p)
          })


setMethod(plotFdr, "Synapter",
          function(object,
                   method = c("BH", "Bonferroni", "qval")) {
            ## Graphical display of qvalues (adapted from qvalue package).
            method <- match.arg(method)
            .qplot <- function(pepdata, rng = c(0, 0.1), ...) {
              pv1 <- pepdata[pepdata$peptide.matchType == "PepFrag1", "pval"]
              qv1 <- pepdata[pepdata$peptide.matchType == "PepFrag1", method]
              qv1 <- qv1[order(pv1)]
              pv2 <- pepdata[pepdata$peptide.matchType == "PepFrag2", "pval"]
              qv2 <- pepdata[pepdata$peptide.matchType == "PepFrag2", method]
              qv2 <- qv2[order(pv2)]
              if (min(c(qv1, qv2)) > rng[2])
                rng <- c(min(c(qv1, qv2)),
                         quantile(c(qv1, pv2), 0.1))
              plot(qv1[qv1 >= rng[1] & qv1 <= rng[2]],
                   (1 + sum(qv1 < rng[1])):sum(qv1 <= rng[2]),
                   type = "l", xlab = "FDR cut-off",
                   ylab = "significant peptides",
                   col = "red", ...)
              lines(qv2[qv2 >= rng[1] & qv2 <= rng[2]],
                    (1 + sum(qv2 < rng[1])):sum(qv2 <= rng[2]),
                    col = "steelblue")
              legend("bottomright", c("PepFrag1", "PepFrag2"),
                     col = c("red", "steelblue"), lty = 1,
                     bty = "n", cex = .6)
              plot((1 + sum(qv1 < rng[1])):sum(qv1 <= rng[2]),
                   qv1[qv1 >= rng[1] & qv1 <= rng[2]] * (1 + sum(qv1 < rng[1])):sum(qv1 <= rng[2]),
                   type = "l", xlab = "significant peptides",
                   ylab = "expected false positives",
                   col = "red", ...)
              lines((1 + sum(qv2 < rng[1])):sum(qv2 <= rng[2]),
                    qv2[qv2 >= rng[1] & qv2 <= rng[2]] * (1 + sum(qv2 < rng[1])):sum(qv2 <= rng[2]),
                    col = "steelblue")
              legend("topleft", c("PepFrag1", "PepFrag2"),
                     col = c("red", "steelblue"), lty = 1,
                     bty = "n", cex = 0.6)
            }
            if (object$Master) {
              par(mfrow=c(1,2))
              .qplot(object$QuantPeptideData,  main="Quantitation")
              par(mfrow=c(1,1))
            } else {
              par(mfrow=c(2,2))
              .qplot(object$IdentPeptideData, main="Identification")
              .qplot(object$QuantPeptideData,  main="Quantitation")
              par(mfrow=c(1,1))
            }
          })

setMethod(plotFeatures, "Synapter",
          function(object,
                   what = c("all", "some"),
                   xlim = c(40, 60),
                   ylim = c(1160, 1165),
                   ionmobility = FALSE) {
            what <- match.arg(what)
            switch(what,
                   all = plot.all.features(
                     object$MergedFeatures,
                     object$QuantPep3DData,
                     ionmobility=ionmobility),
                   some = {
                     if (length(object$PpmError) == 0) {
                       warning("Ppm error for EMRTs matching is not set. Using default value.")
                       object$setPpmError()
                     }
                     if (length(object$RtNsd) == 0) {
                       warning("Number of retention time stdevs for EMRTs matching is not set. Using default value.")
                       object$setRtNsd()
                     }
                     plot.some.features(object$MergedFeatures,
                                        object$IdentPeptideData,
                                        object$QuantPep3DData,
                                        object$RtModel,
                                        object$IdentPpmError,
                                        object$RtNsd,
                                        xlim = xlim,
                                        ylim = ylim)
                   })
          })

setMethod(plotFragmentMatching, "Synapter",
          function(object, key, column="peptide.seq", verbose=TRUE, ...) {
            if (!nrow(object$FragmentMatching)) {
              stop("You have to run ", sQuote("fragmentMatching"), " first!")
            }

            .plotFragmentMatching(object, key, column=column, verbose=verbose,
                               tolerance=object$FragmentMatchingPpmTolerance/1e6,
                               ...)
          })

setMethod(plotFragmentMatchingPerformance, "Synapter",
          function(object, ...) {
            if (!nrow(object$FragmentMatching)) {
              stop("You have to run ", sQuote("fragmentMatching"), " first!")
            }

            invisible(.plotFragmentMatchingPerformance(object$FragmentMatching))
          })

setMethod(plotCumulativeNumberOfFragments, "Synapter",
          function(object, what = c("fragments.ident",
                                    "spectra.quant")) {
            what <- match.arg(what)
            msexp <- switch(what,
                            "fragments.ident" = object$IdentFragmentData,
                            "spectra.quant" = object$QuantSpectrumData)
            .plotIntensityVsNumber(msexp, what = what)
          })

setMethod(getEMRTtable, "Synapter",
          function(object) table(object$MatchedEMRTs$matchedEMRTs))


setMethod(plotEMRTtable, "Synapter",
          function(object) {
            p <- barchart(table(object$MatchedEMRTs$matchedEMRTs), horizontal=FALSE)
            print(p)
            invisible(p)
          })

setMethod(performance, "Synapter",
          function(object, verbose = TRUE) {
            if (nrow(object$MergedFeatures) == 0)
              stop("Merging required before estimating performance.")
            if (nrow(object$MatchedEMRTs) == 0)
              stop("Matching required before estimating performance.")
            ## synapter results
            S <- object$MatchedEMRTs[object$MatchedEMRTs$matchedEMRTs == 1,
                                     "spectrumID"]
            nS <- length(S)
            uS <- unique(S)
            ## Ident peptides
            I <- object$IdentPeptideData$precursor.leID
            nI <- length(I)
            uI <- unique(I)
            ## Quant peptides
            Q <- object$QuantPeptideData$precursor.leID
            nQ <- length(Q)
            uQ <- unique(Q)

            e <- 100 * (nS - nQ) / nQ
            w <- c(length(setdiff(uQ, uS)),
                   length(setdiff(uS, uQ)),
                   length(intersect(uS, uQ)))
            names(w) <- c("Q", "S", "QS")

            ans <- list(nS, nI, nQ,
                        w, e)
            names(ans) <- c("Synapter", "Ident", "Quant",
                            "VennCounts", "Enrichment")
            if (verbose){
              cat("(S) Synapter: ", ans$Synapter, " EMRTs uniquely matched.\n", sep = "")
              cat("(I) Ident: ", ans$Ident, " peptides.\n", sep = "")
              cat("(Q) Quant: ", ans$Quant, " peptides.\n", sep = "")
              cat("Enrichment (S/Q): ", round(ans$Enrichment, 2), "%\n", sep = "")
              cat("Overlap:\n")
              print(ans$VennCounts)
            }
            invisible(ans)
          })

setMethod(performance2, "Synapter",
          function(object, verbose = TRUE) {
              id.source <- object$MatchedEMRTs$idSource
              counts <- object$MatchedEMRTs$Counts
              na.counts <- is.na(counts)
              ans <- table(id.source, na.counts)
              print(ans)
              invisible(list(id.source = id.source, counts = counts))
          })

setMethod(plotRtDiffs, "Synapter",
          function(object, ...) {
            diffs <- plotRetTimeDiffs(object, plot=TRUE, ...)
            invisible(diffs)
          })

setMethod(plotGrid, "Synapter",
          function(object, what = c("total", "model", "details"),
                   maindim = c("im", "rt", "mz")) {
            ## Plots the grid search results.
            if ( length(object$Grid) == 0 )
              stop("No grid search result to plot.")
            what <- match.arg(what)
            if (what == "total") {
              grd <- object$Grid[[1]]
              main <- "Percentage of total ident peptides uniquely matched."
            } else if (what == "model") {
              grd <- object$Grid[[2]]
              main <- "Percentage of modelled peptides matched."
            } else {  ## details
              grd <- object$Grid[[3]]
              main <- "Percentage of correct unique assignments."
            }

            maindim <- match.arg(maindim)
            if (maindim == "im") {
              xlab <- "nsd"
              ylab <- "ppm"
            } else if (maindim == "rt") {
              grd <- aperm(grd, c(2, 3, 1))
              xlab <- "ppm"
              ylab <- "imdiff"
            } else {
              grd <- aperm(grd, c(1, 3, 2))
              xlab <- "nsd"
              ylab <- "imdiff"
            }

            p <- levelplot(grd, xlab = xlab, ylab = ylab, main = main)
            print(p)
            invisible(p)
          })

setMethod(fragmentMatching, "Synapter",
          function(object, ppm, verbose=TRUE) {
            if (!missing(ppm)) {
              setFragmentMatchingPpmTolerance(object, ppm)
            }
            object$fragmentMatching(verbose=verbose)
          })

setMethod(getIdentificationFragments, "Synapter",
          function(object) {
            object$IdentFragmentData
          })

setMethod(getQuantitationSpectra, "Synapter",
          function(object) {
            object$QuantSpectrumData
          })

## Results to csv
setMethod(writeIdentPeptides, "Synapter",
          function(object,
                   file = "Res-IdentPeptides.csv",
                   ...) {
            write.csv(object$IdentPeptideData, file = file, ...)
          })

setMethod(writeQuantPeptides, "Synapter",
          function(object,
                   file = "Res-QuantPeptides.csv",
                   ...) {
            write.csv(object$QuantPeptideData, file = file, ...)
          })


setMethod(writeMergedPeptides, "Synapter",
          function(object,
                   file = "Res-MergedPeptides.csv",
                   what = c("light", "full"),
                   ...) {
            ## Writes merged peptides to a csv file.
            what <- match.arg(what)
            switch(what,
                   full = write.csv(object$MergedFeatures, file = file, ...),
                   light = write.csv(lightMergedFeatures(object$MergedFeatures),
                     file = file, ...))
          })


setMethod(writeMatchedEMRTs, "Synapter",
          function(object,
                   file = "Res-MatchedEMRTs.csv",
                   what = c("light", "full"),
                   ...) {
            ## Writes matched EMRTs to a csv file.
            what <- match.arg(what)
            switch(what,
                   full = write.csv(object$MatchedEMRTs, file = file, ...),
                   light = write.csv(lightMatchedEMRTs(object$MatchedEMRTs),
                     file = file, ...))
          })

setAs("Synapter", "MSnSet",
      function (from) {
        cols <- c("peptide.seq",
                  "protein.Accession",
                  "protein.Description",
                  "protein.falsePositiveRate",
                  "peptide.matchType",
                  "peptide.mhp",
                  "peptide.score",
                  "precursor.mhp",
                  "precursor.retT",
                  "precursor.inten",
                  "precursor.Mobility",
                  "spectrumID",
                  "Intensity",
                  "ion_ID",
                  "ion_area",
                  "ion_counts",
                  "pval",
                  "Bonferroni",
                  "BH",
                  "qval",
                  "isotopicDistr",
                  "synapterPlgsAgreement")
        eset <- matrix(from$MatchedEMRTs$Counts)
        colnames(eset) <- sub("_IA_final_peptide$", "",
                              basename(file_path_sans_ext(from$QuantPeptideFile,
                                                          compression=TRUE)))
        obj <- new("MSnSet",
                   exprs = eset,
                   processingData = new("MSnProcess",
                     processing = "Coerced from a 'Synapter' object."),
                   annotation = "No annotation",
                   featureData = new("AnnotatedDataFrame",
                     data = from$MatchedEMRTs[, cols]))
        fnames <- fData(obj)$peptide.seq
        if (any(duplicated(fnames)))
          fnames <- make.unique(fnames)
        featureNames(obj) <- fnames
        obj <- updateFvarLabels(obj, sampleNames(obj)[1L])
        if (validObject(obj))
          return(obj)
      })

as.MSnSet.Synapter <- function(x) as(x,"MSnSet")

## check class version/updates
setMethod(isCurrent, "Synapter",
          function(object) {
            .isCurrent(object)
})

.validSynapterObject <- function(object) {
    msg <- NULL

    if (!isCurrent(object)) {
      msg <- validMsg(msg,
                      paste0("Your Synapter object is out of date. ",
                             "Please run ",
                             sQuote("object <- updateObject(object)"), "."))
    }

    if (is.null(msg)) TRUE else msg
}
setValidity("Synapter", .validSynapterObject)

setMethod(updateObject, "Synapter",
          function(object, ..., verbose = TRUE) {
            .updateSynapterObject(object, ..., verbose=verbose)
})

