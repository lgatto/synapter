## external script containing chunks for the synergise reports

## @knitr synergise.filtering.pepScores
plotPepScores(obj)

## @knitr synergise.filtering.pepNumbers
knitr::kable(getPepNumbers(obj))

## @knitr synergise.filtering.fdr
plotFdr(obj)

## @knitr synergise.filtering.data
if (uniquepep) {
    filterUniqueDbPeptides(obj, missedCleavages=missedCleavages,
                           IisL=IisL, verbose=verbose)
}
filterPeptideLength(obj, l=peplen)
setPepScoreFdr(obj, fdr=fdr)
filterQuantPepScore(obj, method=fdrMethod)
  if (!master)
    filterIdentPepScore(obj, method=fdrMethod)

## @knitr synergise.filtering.ppmErrorPre
par(mfcol=c(1, 2))
plotPpmError(obj, what="Ident")
plotPpmError(obj, what="Quant")
par(mfcol=c(1, 1))


## @knitr synergise.filtering.ppmErrorQsPre
knitr::kable(getPpmErrorQs(obj))

## @knitr synergise.filtering.ppm
filterQuantPpmError(obj, ppm=quantppm)
filterIdentPpmError(obj, ppm=identppm)

## @knitr synergise.filtering.ppmErrorPost
par(mfcol=c(1, 2))
plotPpmError(obj, what="Ident")
plotPpmError(obj, what="Quant")
par(mfcol=c(1, 1))

## @knitr synergise.filtering.ppmErrorQsPost
knitr::kable(getPpmErrorQs(obj))

## @knitr synergise.filtering.fpr
filterIdentProtFpr(obj, fpr=fpr)
filterQuantProtFpr(obj, fpr=fpr)

## @knitr synergise.merge
mergePeptides(obj)


## @knitr synergise.rtmodel.plotRt
plotRt(obj, what="data")

## @knitr synergise.rtmodel.modelling
setLowessSpan(obj, span.rt)
modelRt(obj)

## @knitr synergise.rtmodel.rtdiffs
par(mfcol=c(1, 2))
plotRtDiffs(obj)
plotRt(obj, what="model", nsd=1) ## better focus on model
par(mfcol=c(1, 1))

## @knitr synergise.rtmodel.rtdiffsTable
rtQs <- t(getRtQs(obj))
rownames(rtQs) <- c("Retention time difference:")
knitr::kable(rtQs)

## @knitr synergise.rtmodel.plotSomeFeatures
setRtNsd(obj, 2)     ## RtNsd and PpmError are used for detailed plot
setPpmError(obj, 10) ## if not set manually, default values are set automatically
plotFeatures(obj, what="some", xlim=c(30,50), ylim=c(1160, 1165))

## @knitr synergise.grid.bestParamaters
setBestGridParams(obj, what=grid.param.sel)

## @knitr synergise.emrt.table
tab <- as.data.frame(getEMRTtable(obj))
colnames(tab) <- c("Number of assigned EMRTs", "Freq")
knitr::kable(tab)

## @knitr synergise.performance
perf <- performance(obj, verbose=FALSE)
group <- c("(S) Synapter (uniquely matched EMRTs)",
           "(I) Identification", "(Q) Quantitation", "(S/Q) Enrichment")
values <- c(perf$Synapter, perf$Ident, perf$Quant, paste0(round(perf$Enrichment, 2), "%"))
knitr::kable(data.frame(group, values, stringsAsFactors=FALSE), align="lr")

overlap <- as.data.frame(perf$VennCounts)
colnames(overlap) <- "Number of peptides"
knitr::kable(overlap)

## @knitr synergise.export
writeMergedPeptides(obj, file=file.path(outputdir, "MergedPeptides.csv"))
writeMatchedEMRTs(obj, file=file.path(outputdir, "MatchedPeptides.csv"))
writeIdentPeptides(obj, file=file.path(outputdir, "IdentPeptides.csv"))
writeQuantPeptides(obj, file=file.path(outputdir, "QuantPeptides.csv"))
write.table(gridDetails, file=file.path(outputdir, "GridDetails.txt"))
saveRDS(obj, file=file.path(outputdir, "SynapterObject.rds"))

## @knitr synergise.log
log <- getLog(obj)
for (i in seq(along=log)){
  cat("  - ", log[i], "\n")
}

## @knitr synergise.sessionInfo
sessionInfo()
