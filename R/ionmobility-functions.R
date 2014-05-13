#' @param d ionmobility difference
#' @param label TRUE/FALSE matchType
#' @param thresholds
#' @return matrix with cols tp, fp, tn, fn
#' @noRd
.ionMobilityConfusionMatrix <- function(d, label,
                                        thresholds=seq(0, max(abs(d)), 0.5)) {
  xtrain <- abs(d)
  ytrain <- label

  confusion <- t(sapply(thresholds, function(th) {
    tp <- sum(xtrain <= th & ytrain)
    fp <- sum(xtrain <= th & !ytrain)
    tn <- sum(xtrain > th & !ytrain)
    fn <- sum(xtrain > th & ytrain)
    return(c(tp=tp, fp=fp, tn=tn, fn=fn))
  }))
  confusion <- cbind(ncommon=thresholds, confusion)
  rownames(confusion) <- NULL

  return(confusion)
}

#' @param d ionmobility difference
#' @param label TRUE/FALSE matchType
#' @param thresholds
#' @return matrix with cols tp, fp, tn, fn
#' @noRd
.ionMobilityContingencyMatrix <- function(d, label, thresholds) {
  confusion <- .ionMobilityConfusionMatrix(d, label, thresholds)
  return(cbind(confusion, diagnosticErrors(confusion)))
}

#' plot ionmobility F1/FDR/Confusion/Boxplots
#' @param obj synapter object
#' @param thresholds
#' @param equalCharge compare only equally charged ions
#' @return invisible matrix with cols tp, fp, tn, fn, accuracy, precision,
#' recall, fdr, f1
#' @noRd
.plotIonMobilitySummary <- function(obj, thresholds=seq(0, 100, 0.5),
                                    equalCharge=FALSE) {

  if (!"precursor.Mobility" %in% colnames(obj$IdentPeptideData) ||
      !"clust_drift" %in% colnames(obj$QuantPep3DData)) {
    stop("No ion mobility information available!")
  }

  emrts <- flatMatchedEMRTs(obj$MatchedEMRTs, obj$QuantPep3DData, verbose=TRUE)

  if (equalCharge) {
    isChargeEqual <- emrts$precursor.z == emrts$ion_z

    emrts <- emrts[isChargeEqual,]
  }

  sampled <- .groundTruthIndices(emrts)

  ## matched.quant.spectrumIDs contains only one ID per row
  ## (done in flatMatchedEMRTs)
  quantIds <- emrts$matched.quant.spectrumIDs

  emrts$ionmobdiff <- emrts$precursor.Mobility -
    obj$QuantPep3DData$clust_drift[quantIds]

  labels <- rep(c(TRUE, FALSE), times=c(length(sampled$trueIdx),
                                        length(sampled$falseIdx)))

  contengency <- .ionMobilityContingencyMatrix(emrts$ionmobdiff[unlist(sampled)],
                                               labels, thresholds)
  colnames(contengency)[1] <- "absdiff"

  par(mfcol=c(1, 3))
  x <- contengency[,1]
  matplot(x, contengency[, c("fdr", "f1"), drop=FALSE], type="b", lty=1,
          xlab="ion mobility difference", ylab="performance",
          main="ion mobility performance", pch=19)
  legend("topright", legend=c("FDR", "F1"),
         col=1:2, lwd=1, pch=19, bty="n")
  grid()
  matplot(x, contengency[, c("tp", "fp", "tn", "fn"), drop=FALSE],
          type="b", lty=1, pch=19,
          xlab="ion mobility difference", ylab="# of peptides",
          main="ion mobility confusion")
  grid()
  legend("right", legend=c("TP", "FP", "TN", "FN"),
         col=1:4, lwd=1, pch=19, bty="n")

  trueIdx <- grep("true", emrts$matchType)
  falseIdx <- grep("false", emrts$matchType)

  l <- list()

  l[["all true matches"]] <- emrts[trueIdx, "ionmobdiff"]

  if (length(trueIdx) > length(sampled$trueIdx)) {
    l[["sampled true matches"]] <- emrts[sampled$trueIdx, "ionmobdiff"]
  }

  l[["all false matches"]] <- emrts[falseIdx, "ionmobdiff"]

  if (length(falseIdx) > length(sampled$falseIdx)) {
    l[["sampled false matches"]] <- emrts[sampled$falseIdx, "ionmobdiff"]
  }

  jitteredBoxplot(l, jitter.factor=1,
                  main="ion mobility difference distribution",
                  ylab="ion mobility difference")

  par(mfcol=c(1, 1))

  invisible(contengency)
}

