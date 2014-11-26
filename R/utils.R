filterFunction <- function(x) { 
  ## keep only function 1 in Pep3DAMRT
  return(x[x$Function == 1,])
}

filterPeptideMatchType <- function(x) {
  ## keep peptide.matchType PepFrag1 or pepFrag2
  return(x[x$peptide.matchType %in% c("PepFrag1", "PepFrag2"),])
}


## pepScore filtering
filterPepScore <- function(dfr,
                           fdr = 0.01,
                           method = c("BH", "qval", "Bonferroni")) {
  method <- match.arg(method)
  ## Assumes that peptide data has been
  ## 
  sel <- vector("logical", length=nrow(dfr))
  sel1 <- dfr$peptide.matchType == "PepFrag1"
  sel2 <- dfr$peptide.matchType == "PepFrag2"
  ## signif1 <- (dfr[sel1, "qval"] <= fdr)
  ## signif2 <- (dfr[sel2, "qval"] <= fdr)
  ## v 0.7.7 - switching 
  signif1 <- (dfr[sel1, method] <= fdr)
  signif2 <- (dfr[sel2, method] <= fdr)
  if (any(is.na(signif1)) | any(is.na(signif2))) {
    stop("Filtering NA qvalues out.")
    ## signif1[is.na(signif1)] <- FALSE
    ## signif2[is.na(signif2)] <- FALSE
  }    
  sel[sel1] <- signif1
  sel[sel2] <- signif2
  return(dfr[sel,])
}


filterProtFpr <- function(pepdata, fpr) {
  pepdata[pepdata$protein.falsePositiveRate < fpr, ]
}

filterCommonSeq <- function(quantpep, identpep) {
  ## keep common sequences
  cmnseq <- intersect(quantpep$peptide.seq, identpep$peptide.seq)
  quantpep <- quantpep[quantpep$peptide.seq %in% cmnseq, ]
  identpep <- identpep[identpep$peptide.seq %in% cmnseq, ]
  return(list(quantpep = quantpep,
              identpep = identpep))
}


filterKeepUniqueSeq <- function(quantpep, identpep) {
  ## DEBUG: some sequences are not unique
  ident2keep <- names(which(table(identpep$peptide.seq)==1))
  quant2keep   <- names(which(table(quantpep$peptide.seq)==1))
  ## DEBUD: save csv for diagnostic
  ## write.csv(quantpep[!quantpep$peptide.seq %in% quant2keep, ],
  ##          file="quant_pep_file_duplicated_entries.csv")
  ## write.csv(identpep[!identpep$peptide.seq %in% ident2keep, ],
  ##          file="ident_pep_file_duplicated_entries.csv")
  ## DEBUG: remove these from data
  identpep <- identpep[identpep$peptide.seq %in% ident2keep, ] 
  quantpep <- quantpep[quantpep$peptide.seq %in% quant2keep, ]
  return(list(quantpep = quantpep,
              identpep = identpep))
}


filterKeepUniqueProt <- function(quantpep, identpep) {
  ## DEBUG: there is still one protein that is
  ##        not common and unique in each data set
  ## if (debug) {
  ##   venn <- Venn(list(quant = quantpep$peptide.seq,
  ##                     ident = identpep$peptide.seq))
  ##   return(venn)
  ## }
  quantspecific <- setdiff(quantpep$peptide.seq, identpep$peptide.seq)
  identpep <- identpep[!identpep$peptide.seq %in% quantspecific, ]
  quantpep <-quantpep[!quantpep$peptide.seq %in% quantspecific, ]
  return(list(quantpep = quantpep,
              identpep = identpep))
}  


error.ppm <- function(obs, theo) {
  ## Estimating MS mass accuracy
  ## Error(ppm) = 1e6 * (observed m/z - exact m/z) / (exact m/z)
  ## with observed m/z = xx$precursor.mhp.[hd]mse
  ##      thoeretical m/z = xx$peptide.mhp.[hd]mse -- calculated from sequence
  ((obs - theo)/theo) * 1e6
}

filter.error.ppm <- function(x, colname, ppm = 2) {
  ## t is the +/- ppm error threshold
  ## this will also be a percentage of points
  ## we want to keep and will be estimated from
  ## calculated mean +/- sd
  ## using ... ?
  if (length(colname) != 1)
    stop("Require exactly 1 column name.")
  if (!colname %in% names(x))
    stop(paste(colname, "not found!"))
  sel <- abs(x[,colname]) < ppm
  return(x[sel,])
}

modelRetTime <- function(xx, span) {
  ## delta = hdmse.rt - mse.rt
  ## hdmse' = mse.rt + delta
  ##        = mse.rt + hdmse.rt - mse.rt
  ##        = hdmse.rt
  ## mse' = hdmse - delta
  ##      = hdmse - (hdmse.rt - mse.rt)
  ##      = hdmse - hdmse.rt + mse.rt
  ##      = mse
  lo <- loess(deltaRt ~ precursor.retT.ident, data = xx, span = span,
              degree = 1, family = "symmetric", iterations = 4, surface = "direct")
  o <- order(xx$precursor.retT.ident)
  pp <- predict(lo, data.frame(precursor.retT.ident = xx$precursor.retT.ident), se=TRUE)
  sd <- pp$se.fit * sqrt(lo$n) ## get sd from se
  stopifnot(all.equal(pp$fit, fitted(lo), check.attributes = FALSE))
  return(list(lo = lo,
              o = o,
              preds = pp,
              sd = sd))
}

doHDMSePredictions <- function(identpep, model, nsd) {
  ## changes in v 0.4.6 - if available, ans is retrieved
  ## from the input data.frame, else it is computed using
  ## the model. 
  if (all(c("predictedRt", "sdRt") %in% names(identpep))) {
    ## get from identpep dataframe
    .fitted <- NA
    .predicted <- identpep$predictedRt
    .sd <- identpep$sdRt
  } else {
    ## compute from model 
    .o <- order(identpep$precursor.retT)
    .allpreds <- predict(model$lo, data.frame(precursor.retT.ident = identpep$precursor.retT), se=TRUE)
    .sd <- .allpreds$se.fit[.o] * sqrt(model$lo$n) ## get sd from se
    .predicted <- identpep$precursor.retT - .allpreds$fit[.o]
    .fitted <- .allpreds$fit[.o]
  }
  if (!missing(nsd)) {
    .lower <- .predicted - (nsd * .sd)
    .upper <- .predicted + (nsd * .sd)
    stopifnot(all(.lower <= .upper))
  } else {
    .lower <- .upper <- NA
  }
  ans <- list(fitted = .fitted,
              predicted = .predicted,
              lower = .lower,
              upper = .upper,
              mass = identpep$peptide.mhp,
              sd = .sd)
  stopifnot(length(ans$predicted) == length(ans$sd))
  return(ans)
}

findRtIndices <- function(sortedPep3d, lowerHDMSeRt, upperHDMSeRt) {
  stopifnot(all(lowerHDMSeRt <= upperHDMSeRt))
  stopifnot(length(lowerHDMSeRt) == length(upperHDMSeRt))

  lowerIdx <- findInterval(lowerHDMSeRt, sortedPep3d$rt_min)+1
  upperIdx <- findInterval(upperHDMSeRt, sortedPep3d$rt_min)

  ## just to be sure
  lowerIdx <- pmin(lowerIdx, upperIdx)

  ## create matrix with boundaries (col 1 and 2) 
  return(cbind(lower=lowerIdx, upper=upperIdx))
}

findMSeEMRTMatches <- function(observedMass,
                               hdmseMass,
                               rtIndices,
                               ppmthreshold) {
  n <- length(hdmseMass)
  stopifnot(nrow(rtIndices) == n)

  res <- apply(cbind(rtIndices, hdmseMass), 1, function(x) {
    selPpm <- which((abs(error.ppm(obs = observedMass[x[1]:x[2]] ,
                                   theo = x[3])) < ppmthreshold))
    (x[1]:x[2])[selPpm]
  })
}

calculateGridPerformance <- function(identpep, sortedPep3d, mergedpep, matches) {
  ## Those that match *1* spectumIDs will be transferred
  ## BUT there is no guarantee that with *1* unique match,
  ##     we get the correct one, even for those that were
  ##     part of the matched high confidence ident+quant
  ##     identified subset! 

  ## #############################################################

  n <- length(matches)
  k <- sapply(matches, length)
  k1 <- which(k == 1)
  k2 <- which(k > 1)

  idx <- match(mergedpep$precursor.leID.ident, identpep$precursor.leID)
  precursor.leID.quant <- rep(NA, n)
  precursor.leID.quant[idx] <- mergedpep$precursor.leID.quant

  spectrumID <- rep(NA, n)
  spectrumID[k1] <- sortedPep3d$spectrumID[unlist(matches[k1])]
  
  multipleMatchedSpectrumIDs <- vector(mode="list", length=n)
  multipleMatchedSpectrumIDs[k2] <- matches[k2]

  ## grd1: number of unique matches divided by total number of matches
  ## => sum(k==1)/length(k) == mean(k==1)
  grd1 <- mean(k == 1)

  notNaIdx <- which(!is.na(precursor.leID.quant))

  ## grd2: number of correct matched 
  ## (MSe peptide's leID == MSe Pep3D's spectrumID)
  ## divided by number of features used in model
  ## non-NAs are those used for modelling
  grd2 <- sum(precursor.leID.quant == spectrumID, na.rm=TRUE)/length(notNaIdx)

  details <- integer(n)
  details[k1] <- ifelse(spectrumID[k1] == precursor.leID.quant[k1], 1, -1)
  details[k2] <- ifelse(unlist(lapply(k2, function(i) { 
    precursor.leID.quant[i] %in% multipleMatchedSpectrumIDs[[i]]})), 2, -2)

  ## exclude all values where precursor.leID.quant == NA
  details <- details[notNaIdx]

  ## tabulate needs positive integer values
  details <- setNames(tabulate(1-min(details)+details), sort(unique(details)))

  return(list(grd1=grd1, grd2=grd2, details=details))
}

findMSeEMRTs <- function(identpep, 
                         pep3d,
                         mergedpep,
                         nsd,
                         ppmthreshold, 
                         model,
                         mergedEMRTs) {
  sortedPep3d <- pep3d
  sortedPep3d <- sortedPep3d[order(pep3d$rt_min),]

  hdmseData <- doHDMSePredictions(identpep, model, nsd)
  rtIdx <- findRtIndices(sortedPep3d, hdmseData$lower, hdmseData$upper)
  res <- findMSeEMRTMatches(sortedPep3d$mwHPlus, hdmseData$mass, rtIdx,
                            ppmthreshold)
  k <- sapply(res, length)
 
  ## Those that match *1* spectumIDs will be transferred
  ## BUT there is no guarantee that with *1* unique match,
  ##     we get the correct one, even for those that were
  ##     part of the matched high confidence ident+quant
  ##     identified subset! 

  ## #############################################################

  n <- length(k)
  m <- ncol(pep3d)
  ## to initialise the new pep3d2 with with n rows 
  ## and same nb of columns than pep3d
  pep3d2 <- matrix(NA, nrow=n, ncol=m, dimnames=list(c(), colnames(pep3d)))
  pep3d2[, 1] <- k
  k1 <- which(k == 1)
  pep3d2[k1, ] <- unlist(sortedPep3d[unlist(res[k1]), ])

  ## #############################################################
  
  ans <- cbind(identpep, pep3d2)
  
  ans$matched.quant.spectrumIDs <- sapply(res, paste, collapse = ",")
  
  ans$precursor.leID.quant <- NA  
  idx <- match(mergedpep$precursor.leID.ident, ans$precursor.leID)
  
  ans$precursor.leID.quant[idx] <- mergedpep$precursor.leID.quant
  i <- grep("precursor.leID$", names(ans))
  names(ans)[i] <- "precursor.leID.ident" ## to avoid any confusion

  ans$idSource <- "transfer"
  
  if (mergedEMRTs == "rescue") {
    ## these are those that were in the merged data set but that
    ## did NOT get transferred because they did NOT uniquely matched
    ## a pep3D EMRT
    lost <- ans$Function != 1 & ans$precursor.leID.ident %in% mergedpep$precursor.leID.ident
    ## rescue these by adding their quant straight from QuantPeptideData
    lostids <- ans$precursor.leID.ident[lost]    
    ans[lost, "Counts"] <-
      mergedpep[match(lostids, mergedpep$precursor.leID.ident), "precursor.inten.quant"]
    ans[lost, "idSource"] <- "rescue"
  } else if (mergedEMRTs == "copy") {    
    allmerged <- ans$precursor.leID.ident %in% mergedpep$precursor.leID.ident
    allmergedids <- ans$precursor.leID.ident[allmerged]
    ans[allmerged, "Counts"] <-
      mergedpep[match(allmergedids, mergedpep$precursor.leID.ident), "precursor.inten.quant"]
    ans[allmerged, "idSource"] <- "copy"
  } ## else if (fromQuant == "transfer") keep as is
  
  ## since v 0.5.0 - removing multiply matched EMRTs
  ## dupIDs <- ans$spectrumID[ans$Function == 1 & duplicated(ans$spectrumID)]
  ## dupRows <- ans$spectrumID %in% dupIDs
  ## ans[dupRows, (ncol(identpep)+1) : ncol(ans)] <- -1 
  
  return(ans)
}


estimate.mass.range <- function(Mhdmse, ppmt) {
  mass.ranges <- sapply(Mhdmse, function(m)
                        c(( (ppmt * m) / 1e6) + m,
                          ((-ppmt * m) / 1e6) + m))
  return(t(mass.ranges))
}


lightMergedFeatures <- function(x) {
  cols <- c("peptide.seq",
            "protein.Accession.ident",
            "protein.Description.ident",
            "protein.dataBaseType.ident",
            "protein.falsePositiveRate.ident", 
            "peptide.matchType.ident", 
            "peptide.mhp.ident", 
            "peptide.score.ident", 
            "precursor.mhp.ident", 
            "precursor.retT.ident", 
            "precursor.inten.ident",
            "pval.ident",
            "Bonferroni.ident",
            "BH.ident",
            "qval.ident",
            ## "precursor.Mobility", ## not present in MSe master
            "peptide.matchType.quant", 
            "peptide.score.quant", 
            "precursor.retT.quant",
            "precursor.inten.quant", 
            "deltaRt", 
            "errorppm.ident", 
            "errorppm.quant",
            "pval.quant",
            "Bonferroni.quant",
            "BH.quant",
            "qval.quant")
  return(x[,cols])
}

lightMatchedEMRTs <- function(x) {
  cols <- c("peptide.seq",
            "protein.Accession",
            "protein.Description",
            "protein.dataBaseType",
            "protein.falsePositiveRate", 
            "peptide.matchType", 
            "peptide.mhp", 
            "peptide.score", 
            "precursor.mhp", 
            "precursor.retT", 
            "precursor.inten", 
            ## "precursor.Mobility", ## not present in MSe master
            "spectrumID",
            "Intensity",
            "Counts",
            "ion_ID",
            "ion_area",
            "ion_counts")
  return(x[,cols])
}


gridSearch2 <- function(model,
                        identpep,
                        pep3d,
                        mergedPeptides,
                        ppms,
                        nsds,
                        verbose = TRUE) {
  ## As initial gridSearch, but now returns a list
  ## with two grids; first one as before, percent of
  ## uniquely matched features; second is percent of
  ## correctly assigned features (based on marged feautres
  ## used to model rt).
  ##  
  n <- length(nsds)
  m <- length(ppms)  
  grd1 <- grd2 <- outer(nsds, ppms)
  N <- n*m
  details <- vector("list", length=N)

  ._k <- 0
  if (verbose) {
    pb <- txtProgressBar(min=0, max=N, style=3)
  }

  sortedPep3d <- pep3d[order(pep3d$rt_min),]

  for (i in 1:n) {
    hdmseData <- doHDMSePredictions(identpep, model, nsds[i])
    rtIdx <- findRtIndices(sortedPep3d, hdmseData$lower, hdmseData$upper)

    for (j in 1:m) {
      ._k <- ._k + 1
      matches <- findMSeEMRTMatches(sortedPep3d$mwHPlus,
                                    hdmseData$mass,
                                    rtIdx,
                                    ppms[j])
      
      gridDetails <- calculateGridPerformance(identpep, sortedPep3d, mergedPeptides, matches)
      grd1[i, j] <- gridDetails$grd1
      grd2[i, j] <- gridDetails$grd2
      details[[._k]] <- gridDetails$details

      if (verbose) {
        setTxtProgressBar(pb, ._k)
      }
    }
  }
  if (verbose) {
    close(pb)
  }
  nms <- paste(rep(nsds, each = m) ,
               rep(ppms, n),
               sep=":")
  names(details) <- nms
  return(list(prcntTotal = grd1,
              prcntModel = grd2,
              details = details))
}

makeFigurePath <- function(dirname, filename) {
  full <- paste(dirname, "/", filename, sep="")
  rel <- filename
  return(list(full=full, relative=rel))
}

getQs <- function(x, qtls) {
    xi <- length(x) * qtls
    yi <- quantile(sort(abs(x)), qtls)
    return(list(x=xi, y=yi))
}

keepUniqueSpectrumIds <- function(pep3d) {
  ## The EMRTs spectra are duplicated
  ## for different charge states, isotopes, ...
  ## This functions removes the duplicated
  ## spectrum ids and keeps the first instances
  ## only
  pep3d[!duplicated(pep3d$spectrumID),]
}

score2pval <- function(xx) {  
  srnd <- xx[xx$protein.dataBaseType == "Random", "peptide.score"]
  sreg <- xx[xx$protein.dataBaseType == "Regular", "peptide.score"]
  nrnd <- length(srnd)
  sapply(sreg, function(x) sum(srnd >= x)/nrnd)
}

getIdStats <- function(pepdata) {
  dbtypes <- sort(unique(pepdata$protein.dataBaseType))
  if (any(dbtypes != c("Random", "Regular"))) 
    stop("[Regular and Random] not found.")
  ## new in v 0.8.1 -- keeping only unique (as in unique()) 
  ## peptide entries to calculate statistics, to avoid to biasing
  ## towards repeated identical peptides (from different proteins
  ## for instance).
  ## fixed in v 0.8.3
  pepseqs <- as.character(pepdata$peptide.seq)
  res <- pepdata[, c("peptide.seq",
                     "peptide.score",
                     "peptide.matchType",
                     "protein.dataBaseType")]
  rownames(res) <- make.unique(pepseqs)    
  res$pval <- res$qval <- res$Bonferroni <- res$BH <- NA  
  ## name subsets
  pf1 <- rownames(res)[res$peptide.matchType == "PepFrag1"]
  pf2 <- rownames(res)[res$peptide.matchType == "PepFrag2"]
  pf1reg <- rownames(res)[res$peptide.matchType == "PepFrag1" &
                          res$protein.dataBaseType == "Regular"]
  pf2reg <- rownames(res)[res$peptide.matchType == "PepFrag2" &
                          res$protein.dataBaseType == "Regular"]
  ## PepFrag1  
  .pv1 <- score2pval(res[pf1, ])
  res[pf1reg, "pval"] <- .pv1  
  res[pf1reg, "qval"] <- qvalue(.pv1)$qvalues
  .adj1 <- mt.rawp2adjp(.pv1)
  res[pf1reg, "BH"] <- .adj1$adjp[order(.adj1$index), "BH"]
  res[pf1reg, "Bonferroni"] <- .adj1$adjp[order(.adj1$index), "Bonferroni"]    
  ## PepFrag2
  .pv2 <- score2pval(res[pf2, ])
  res[pf2reg, "pval"] <- .pv2
  res[pf2reg, "qval"] <- qvalue(.pv2)$qvalues
  .adj2 <- mt.rawp2adjp(.pv2)
  res[pf2reg, "BH"] <- .adj2$adjp[order(.adj2$index), "BH"]
  res[pf2reg, "Bonferroni"] <- .adj2$adjp[order(.adj2$index), "Bonferroni"]    
  ## checking
  testsel <- res$peptide.matchType %in% c("PepFrag1", "PepFrag2") & res$protein.dataBaseType == "Regular"
  if (any(is.na(res[testsel, "pval"])))
    stop("NA p-values generated")
  ## 
  ans <- res[, c("pval", "qval", "BH", "Bonferroni")]
  rownames(ans) <- rownames(pepdata)
  return(ans)                             
}

## closes #42
isCorrespondingPep3DataFile <- function(quant, pep3d) {
  idx <- match(quant$precursor.leID, pep3d$spectrumID)

  ## We disable the comparison of the intensity values for the file check
  ## temporarly because some files exist that have entries with different
  ## intensity values.
  ## See https://github.com/lgatto/synapter/issues/42 for details
  return(!anyNA(idx))
#  return(!anyNA(idx) && all(quant$precursor.inten == pep3d$Counts[idx]))
}

getExtension <- function (filename) {
  x <- strsplit(filename, "\\.")[[1]]
  ext <- x[length(x)]
  return(ext)
}

