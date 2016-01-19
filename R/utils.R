filterFunction <- function(x) {
  ## keep only function 1 in Pep3DAMRT
  ## was Function column before (see issue #67)
  return(x[x$matchedEMRTs == 1,])
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
  if (anyNA(signif1) | anyNA(signif2)) {
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

rawRetTimeModel <- function(retT, deltaRt, span) {
  loess(deltaRt ~ retT, span = span, degree = 1, family = "symmetric",
        iterations = 4, surface = "direct")
}

modelRetTime <- function(retT, deltaRt, span) {
  ## delta = hdmse.rt - mse.rt
  ## hdmse' = mse.rt + delta
  ##        = mse.rt + hdmse.rt - mse.rt
  ##        = hdmse.rt
  ## mse' = hdmse - delta
  ##      = hdmse - (hdmse.rt - mse.rt)
  ##      = hdmse - hdmse.rt + mse.rt
  ##      = mse
  lo <- rawRetTimeModel(retT, deltaRt, span)
  o <- order(retT)
  pp <- predict(lo, data.frame(retT = retT), se=TRUE)
  sd <- pp$se.fit * sqrt(lo$n) ## get sd from se
  stopifnot(all.equal(pp$fit, fitted(lo), check.attributes = FALSE))
  list(lo = lo,
       o = o,
       preds = pp,
       sd = sd)
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
    .allpreds <- predict(model$lo, data.frame(retT = identpep$precursor.retT),
                         se=TRUE)
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
  stopifnot(length(.predicted) == length(.sd))
  list(fitted = .fitted,
       predicted = .predicted,
       lower = .lower,
       upper = .upper,
       mass = identpep$peptide.mhp,
       sd = .sd)
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

findMzIndices <- function(pep3dMz, hdmseMz, rtIndices, ppmthreshold) {
  stopifnot(nrow(rtIndices) == length(hdmseMz))
  apply(cbind(rtIndices, hdmseMz), 1, function(x) {
    (x[1]:x[2])[which((abs(error.ppm(obs=pep3dMz[x[1]:x[2]] ,
                                     theo=x[3])) < ppmthreshold))]
  })
}

findImIndices <- function(pep3dIm, hdmseIm, mzIndices, imthreshold) {
  if (imthreshold == Inf) {
    return(mzIndices)
  }
  stopifnot(length(mzIndices) == length(hdmseIm))
  mapply(function(mzIdx, im) {
    mzIdx[which( abs(pep3dIm[mzIdx] - im) < imthreshold )]
  }, mzIdx=mzIndices, im=hdmseIm, SIMPLIFY=FALSE, USE.NAMES=FALSE)
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
  k1 <- which(k == 1L)
  k2 <- which(k > 1L)

  idx <- match(mergedpep$precursor.leID.ident, identpep$precursor.leID)
  precursor.leID.quant <- rep(NA_real_, n)
  precursor.leID.quant[idx] <- mergedpep$precursor.leID.quant

  spectrumID <- rep(NA_real_, n)
  spectrumID[k1] <- sortedPep3d$spectrumID[unlist(matches[k1])]

  multipleMatchedSpectrumIDs <- vector(mode="list", length=n)
  multipleMatchedSpectrumIDs[k2] <- matches[k2]

  ## grd1: number of unique matches divided by total number of matches
  ## => sum(k==1)/length(k) == mean(k==1)
  grd1 <- mean(k == 1L)

  notNaIdx <- which(!is.na(precursor.leID.quant))

  ## grd2: number of correct matched
  ## (MSe peptide's leID == MSe Pep3D's spectrumID)
  ## divided by number of features used in model
  ## non-NAs are those used for modelling
  grd2 <- sum(precursor.leID.quant == spectrumID, na.rm=TRUE)/length(notNaIdx)

  details <- integer(n)
  details[k1] <- -1L
  details[k1][spectrumID[k1] == precursor.leID.quant[k1]] <- 1L
  details[k2] <- -2L
  details[k2][unlist(lapply(k2, function(i) {
    precursor.leID.quant[i] %in% multipleMatchedSpectrumIDs[[i]]}))] <- 2L

  ## exclude all values where precursor.leID.quant == NA
  details <- details[notNaIdx]

  ## tabulate needs positive integer values
  details <- setNames(tabulate(1-min(details)+details), sort(unique(details)))

  ## make sure that we have alway 5 categories
  grddetails <- c("-2"=0, "-1"=0, "0"=0, "1"=0, "2"=0)
  grddetails[names(details)] <- details

  return(list(grd1=grd1, grd2=grd2, details=grddetails))
}

findMSeEMRTs <- function(identpep,
                         pep3d,
                         mergedpep,
                         nsd,
                         ppmthreshold,
                         imdiff,
                         model) {
  sortedPep3d <- pep3d
  sortedPep3d <- sortedPep3d[order(pep3d$rt_min),]

  hdmseData <- doHDMSePredictions(identpep, model, nsd)
  rtIdx <- findRtIndices(sortedPep3d, hdmseData$lower, hdmseData$upper)
  mzIdx <- findMzIndices(sortedPep3d$mwHPlus, hdmseData$mass, rtIdx,
                         ppmthreshold)
  res <- findImIndices(sortedPep3d$clust_drift, identpep$precursor.Mobility,
                       mzIdx, imdiff)

  k <- lengths(res)

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
  pep3d2 <- matrix(NA_real_, nrow=n, ncol=m, dimnames=list(c(), colnames(pep3d)))
  pep3d2[, 1] <- k
  k1 <- which(k == 1)
  ## convert to data.frame first to avoid conversion of matrix to character
  ## matrix (because of the new isotopicDistr column) that results in a
  ## data.frame full of factors
  pep3d2 <- as.data.frame(pep3d2, stringsAsFactors = FALSE)
  pep3d2[k1, ] <- sortedPep3d[unlist(res[k1]), ]

  ## #############################################################

  ans <- cbind(identpep, pep3d2)

  ans$matched.quant.spectrumIDs <- MSnbase:::utils.list2ssv(res, sep=";")

  ans$precursor.leID.quant <- NA_real_
  idx <- match(mergedpep$precursor.leID.ident, ans$precursor.leID)

  ans$precursor.leID.quant[idx] <- mergedpep$precursor.leID.quant
  i <- grep("precursor.leID$", names(ans))
  names(ans)[i] <- "precursor.leID.ident" ## to avoid any confusion

  ans$idSource <- "transfer"

  ## since v 0.5.0 - removing multiply matched EMRTs
  ## dupIDs <- ans$spectrumID[ans$Function == 1 & duplicated(ans$spectrumID)]
  ## dupRows <- ans$spectrumID %in% dupIDs
  ## ans[dupRows, (ncol(identpep)+1) : ncol(ans)] <- -1

  ans
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

gridSearch3 <- function(model,
                        identpep,
                        pep3d,
                        mergedPeptides,
                        ppms,
                        nsds,
                        imdiffs,
                        verbose = TRUE) {
  ## As initial gridSearch, but now returns a list
  ## with two grids; first one as before, percent of
  ## uniquely matched features; second is percent of
  ## correctly assigned features (based on merged features
  ## used to model rt).
  ##

  ## set imdiffs=Inf to disable 3D grid search

  n <- length(nsds)
  m <- length(ppms)
  o <- length(imdiffs)
  grd1 <- grd2 <- array(NA_integer_, dim=c(n, m, o),
                        dimnames=list(nsds, ppms, imdiffs))
  N <- n*m*o
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
      mzIdx <- findMzIndices(sortedPep3d$mwHPlus, hdmseData$mass, rtIdx,
                             ppms[j])
      for (k in 1:o) {
        ._k <- ._k + 1L

        matches <- findImIndices(sortedPep3d$clust_drift,
                                 identpep$precursor.Mobility, mzIdx,
                                 imdiffs[k])
        gridDetails <- calculateGridPerformance(identpep, sortedPep3d,
                                                mergedPeptides, matches)
        grd1[i, j, k] <- gridDetails$grd1
        grd2[i, j, k] <- gridDetails$grd2
        details[[._k]] <- gridDetails$details
      }

      if (verbose) {
        setTxtProgressBar(pb, ._k)
      }
    }
  }
  if (verbose) {
    close(pb)
  }
  nms <- paste(rep(nsds, each = m*o) ,
               replicate(n, rep(ppms, each = o)),
               rep(imdiffs, m*n),
               sep=":")
  names(details) <- nms
  list(prcntTotal = grd1,
       prcntModel = grd2,
       details = details)
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
  if (anyNA(res[testsel, "pval"]))
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

## duplicate rows if the have multiple matched.quant.spectrumIDs
## adds a columns "gridSearchResult" that could be "unique-true",
## "unique-false", "non-unique-true", "non-unique-false"
## please note that it would change the order of the data.frame
flatMatchedEMRTs <- function(emrts, pep3d, na.rm=TRUE, verbose=TRUE) {
  if (verbose) {
    message("create flat EMRT data.frame (one matched EMRT per row)")
  }

  if (na.rm) {
    emrts <- emrts[!is.na(emrts$precursor.leID.quant), ]
  }

  emrts$gridSearchResult <- "no_quant_id"

  ## unique matches
  k1 <- which(emrts$matchedEMRTs == 1)
  isCorrectMatch <- as.numeric(emrts$spectrumID[k1]) ==
                      as.numeric(emrts$precursor.leID.quant[k1])
  emrts$gridSearchResult[k1] <- ifelse(isCorrectMatch,
                                       "unique-true", "unique-false")
  emrts$matched.quant.spectrumIDs[k1] <- emrts$spectrumID[k1]

  ## non-unique matches
  k2 <- which(emrts$matchedEMRTs > 1)
  mIds <- matched.quant.spectrumIDs2numeric(emrts$matched.quant.spectrumIDs)

  flatEmrts <- lapply(k2, function(j) {
    ## duplicated current row
    curRow <- emrts[rep(j, length(mIds[[j]])), ]
    ## update matched.quant.spectrumIDs (contain only a single entry now)
    curRow$matched.quant.spectrumIDs <- mIds[[j]]
    ## update spectrumID
    curRow$spectrumID <- curRow$matched.quant.spectrumIDs

    curRow$gridSearchResult <- ifelse(mIds[[j]] %in%
                                      as.numeric(curRow$precursor.leID.quant),
                                      "non-unique-true", "non-unique-false")
    return(curRow)
  })

  ## build final df
  flatEmrts <- do.call(rbind, flatEmrts)
  flatEmrts <- rbind(emrts[k1, ], flatEmrts)
  flatEmrts$matched.quant.spectrumIDs <-
    as.numeric(flatEmrts$matched.quant.spectrumIDs)
  rownames(flatEmrts) <- NULL

  ## refill pep3d information
  rows <- flatEmrts$matched.quant.spectrumIDs

  ## find corresponding columns (but exclude spectrumID and Function
  cols <- intersect(colnames(flatEmrts), colnames(pep3d)[-c(1:2)])

  flatEmrts[, cols] <- pep3d[rows, cols]

  return(flatEmrts)
}

## see issue 73 for details
## https://github.com/lgatto/synapter/issues/73
.findSynapterPlgsAgreement <- function(emrts) {
  ## single match
  k1 <- emrts$matchedEMRTs == 1

  k1Idx <- which(k1)

  agreement <- rep.int("no_synapter_transfer", nrow(emrts))

  ## agree (same EMRT by database search and synapter)
  ## disagree (different EMRT by database search and synapter)
  mqsid <- lapply(matched.quant.spectrumIDs2numeric(emrts$matched.quant.spectrumIDs[k1Idx]), "[", 1L)
  agreement[k1Idx] <- ifelse(mqsid == as.numeric(emrts$precursor.leID.quant[k1Idx]), "agree", "disagree")

  ## no_plgs_id (transferred by synapter not ided by PLGS)
  agreement[is.na(emrts$precursor.leID.quant)] <- "no_plgs_id"
  ## no_id_or_transfer (not transferred by synapter not ided by PLGS)
  agreement[(is.na(emrts$precursor.leID.quant) & is.na(emrts$matched.quant.spectrumIDs)) |
            (!k1 & is.na(emrts$precursor.leID.quant))] <- "no_id_or_transfer"

  agreement
}

.appendFragmentMatchingColumn <- function(emrts, fm) {
  idx <- match(sort(unique(fm$precursor.leID.ident)),
               emrts$precursor.leID.ident)

  emrts[idx, "FragmentMatching"] <- MSnbase:::utils.list2ssv(
    split(fm[, "FragmentMatching"], fm$precursor.leID.ident))

  return(emrts)
}

matched.quant.spectrumIDs2numeric <- function(x) {
  lapply(MSnbase:::utils.ssv2list(x), as.numeric)
}

diagnosticErrors <- function(x) {
  tp <- x[, "tp"]
  fp <- x[, "fp"]
  tn <- x[, "tn"]
  fn <- x[, "fn"]

  acc <- (tp+tn)/(tp+fp+tn+fn)
  pre <- tp/(tp+fp)
  rec <- tp/(tp+fn)
  f1 <- 2*pre*rec/(pre+rec)
  return(cbind(accuracy=acc, precision=pre, recall=rec, fdr=1-pre, f1=f1))
}

# concatenate isotope distribution information
# @param df data.frame (e.g. pep3d data)
# @param idcol id column (e.g. "spectrumID")
# @param zcol charge column name (e.g. "ion_z")
# @param isocol isotope column name (e.g. "ion_iso")
# @param intcol intensity column name (e.g. "ion_counts")
# @return a character vector of length == nrow(df)
# @noRd
.catIsotopeDistr <- function(df,
                             idcol = "spectrumID",
                             zcol = "ion_z",
                             isocol = "ion_iso",
                             intcol = "ion_counts") {
  ## cat each row
  r <- paste0(df[, zcol], "_", df[, isocol], ":",  df[, intcol])
  ## combine rows with equal ids
  ave(r, df[, idcol], FUN=function(x)paste0(x, collapse=";"))
}

# convert isotope distribution information to matrix
# @param x vector of isotopicDistr strings, generated by .catIsotopeDistr
# @return a matrix of isotopic counts; rownames: runs;
#   colnames: isotopic names
# @noRd
.isotopicDistr2matrix <- function(x) {
  if (is.data.frame(x)) {
    x <- unlist(x)
  }
  nm <- names(x)
  n <- length(x)

  x <- strsplit(x, ":|;")
  ne <- lengths(x)

  ## if input contains NA we remove it manually
  ne1 <- which(ne == 1L)
  x[ne1] <- NULL
  ne[ne1] <- 0L
  nall <- sum(ne)
  x <- unlist(x, use.names=FALSE)

  cn <- x[seq(from=1L, to=nall, by=2L)]
  ucn <- unique(cn)

  m <- matrix(NA_real_, nrow=n, ncol=length(ucn),
              dimnames=list(nm, ucn))

  rows <- rep(1:n, ne/2L)
  cols <- match(cn, ucn)
  i <- rows + (cols - 1L) * n

  m[i] <- as.double(x[seq(from=2L, to=nall, by=2L)])
  m
}

.rescueEMRTs <- function(matchedEMRTs, mergedFeatures,
                         method=c("rescue", "copy")) {
  method <- match.arg(method)

  if (method == "rescue") {
    ## these are those that were in the merged data set but that
    ## did NOT get transferred because they did NOT uniquely matched
    ## a pep3D EMRT
    i <- matchedEMRTs$matchedEMRTs != 1 &
         matchedEMRTs$precursor.leID.ident %in%
           mergedFeatures$precursor.leID.ident
  } else {
    i <- matchedEMRTs$precursor.leID.ident %in%
          mergedFeatures$precursor.leID.ident
  }
  ids <- matchedEMRTs$precursor.leID.ident[i]
  matchedEMRTs$Counts[i] <- mergedFeatures$precursor.inten.quant[
                              match(ids, mergedFeatures$precursor.leID.ident)]
  matchedEMRTs$idSource[i] <- method
  matchedEMRTs
}
