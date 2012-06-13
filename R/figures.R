
plotErrorPpm <- function(xx) {
  layout(c(1,1,1,2))
  par(mar=c(0,4,3,2))
  hist(xx$errorppm.quant, breaks=100, freq=FALSE, xaxt="n", ylab="",
       main="Error(ppm)")
  rug(xx$errorppm.quant)
  legend("topright", c("Ident", "Quant", "normal"), lty=1,
         col=c("red", "blue", "green"), bty="n")
  lines(density(xx$errorppm.ident),col="red")
  lines(density(xx$errorppm.quant),col="blue")
  lines(density(rnorm(1e6, mean=mean(xx$errorppm.quant), sd=sd(xx$errorppm.quant))),
        col="green")
  abline(v=quantile(xx$errorppm.quant,
           c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)),
         lty="dotted")
  par(mar=c(5,4,1,2))
  boxplot(xx$errorppm.quant, horizontal=TRUE, xlab="Quantitation error [ppm]")
}

plotLowess <- function(xx,
                       f = structure(
                         c(2/3, 1/2, 1/4, 1/10, 1/16, 1/25, 1/50),
                         names = c("2/3", "1/2", "1/4", "1/10", "1/16", "1/25", "1/50"))
                       ) {
  par(mfrow=c(1,2))
  plot(xx$precursor.retT.ident, xx$deltaRt,
       col = "#8FBDDA80", ## col2hcl("steelblue", alpha=.5),
       pch = 19,
       ylab = expression(Identification - Quantitation),
       xlab = expression(Identification~retention~time))
  abline(h=0)
  grid()
  ## quick lowess to get range
  lws <- lowess(xx$precursor.retT.ident, xx$deltaRt, f = 1/2)
  plot(xx$precursor.retT.ident, xx$deltaRt,
       col = "#8FBDDA80", ## col2hcl("steelblue", alpha=.5),
       pch = 19,
       ylab = expression(Identification - Quantitation),
       ylim = range(lws$y) + c(-.2,.2),
       xlab = expression(Identification~retention~time))
  abline(h=0)
  grid()
  cls <- brewer.pal(length(f), "Set1")  
  sapply(1:length(f), function(i) {
    span <- fi <- f[i]
    ## using lowess
    lws <- lowess(xx$precursor.retT.ident, xx$deltaRt, f = fi)
    lines(lws$x, lws$y, col=cls[i])
    ## using loess -- much slower, but has predict, for later predictions
    ## lo <- loess(deltaRt ~ precursor.retT.ident, data = xx, span=span,
    ##             degree=1, family="symmetric", iterations=4, surface="direct")
    ## o <- order(xx$precursor.retT.ident)
    ## lines(xx$precursor.retT.ident[o], fitted(lo)[o], col=cls[i])
  })
  legend("bottomright",
         paste("span=", names(f), sep = ""),
         col = cls, lty = 1, bty = "n", cex=.7)
  par(mfrow=c(1,1))
}


plotLowess2 <- function(xx, model, nsd, ...) {
  o <- model$o
  lo <- model$lo
  pp <- model$preds
  sd <- model$sd
  
  sdlim1 <- sapply(nsd, function(i) pp$fit[o] - i * sd[o])
  sdlim2 <- sapply(nsd, function(i) pp$fit[o] + i * sd[o])
  ylim <- range(c(sdlim1, sdlim2))
  
  plot(xx$precursor.retT.ident, xx$deltaRt,
       col = "#8FBDDA80", ## col2hcl("steelblue", alpha=.5),
       pch = 19,
       ylim = ylim,
       ylab = expression(Identification - Quantitation),
       xlab = expression(Identification~retention~time),
       type = "n",
       ...)
  grid()
  polcls <- c("#8FBDDAFF", ## col2hcl("steelblue", alpha=1)
              "#8FBDDA80", ## col2hcl("steelblue", alpha=.5)
              "#8FBDDA33") ## col2hcl("steelblue", alpha=.2)
  for (i in ncol(sdlim1):1) {
    polygon(c(xx$precursor.retT.ident[o], rev(xx$precursor.retT.ident[o])),
            c(sdlim1[,i], rev(sdlim2[,i])),
            col = polcls[i],
            lty = 0)
  }
  abline(h = 0, lty = "dotted")
  points(xx$precursor.retT.ident, xx$deltaRt,
         col = "#0000004D", ## col2hcl("black", alpha=.3),
         pch = 19)
  lines(xx$precursor.retT.ident[o], fitted(lo)[o], col="red", lwd=2)
}

plotRetTimeDiffs <- function(object, plot = TRUE,
                             freq = FALSE,
                             xlab = "observed - fitted rt",
                             ...) {
  o <- object$RtModel$o
  lo <- object$RtModel$lo
  diffs <- object$MergedFeatures$deltaRt[o] - fitted(lo)[o]
  ## diffs <-  xx$deltaRt[o] - fitted(lo)[o]
  if (plot) {
    hist(diffs, freq = freq, xlab = xlab, ...)
    rug(diffs)
  }
  invisible(diffs)
}

plot.some.features <- function(xx,
                               mse,
                               model,
                               ppmthreshold,
                               nsd,                               
                               xlim = c(40, 60), 
                               ylim = c(1160, 1165)) {
  ## for EMRT mass, use [precursor|pretide].mhp.[hd]mse
  ## on the y axis, instead of calculating the mass
  ## with precursor.mz.[hd]mse * precursor.z.[hd]mse.
  ## The former takes into account relative proportions
  ## of different ions (1+, 2+, 3+, ...) whereas the
  ## latter uses the mz of the most abondant ion.
  ## Note: peptide.mhp.hdmse and peptide.mhp.mse are
  ## identical (for identical peptides), as these are
  ## the masses calculated from the peptide sequence.
  ## hdmse final peptide data
  plot(xx$precursor.retT.ident, xx$peptide.mhp.ident,
       col = "#D4D4D480", ##col2hcl("darkgrey", alpha=.5),
       cex = .4, pch = 19,
       xlim = xlim, ylim = ylim,
       xlab = "retention time", ylab = "precursor mass")
  text(xx$precursor.retT.ident, xx$peptide.mhp.ident,
       xx$peptide.seq, cex=.1, adj=c(1,0))
  grid()
  ## mse final peptide data
  ## observed
  points(xx$precursor.retT.quant, xx$peptide.mhp.quant,
         col = "#8FBDDA80", ## col2hcl("steelblue", alpha=.5),
         cex = 0.4, pch = 19)
  text(xx$precursor.retT.quant, xx$peptide.mhp.quant,
       xx$peptide.seq, cex=.1, adj=c(0,1))
  ## mse pep3d data
  points(mse$rt_min, mse$mwHPlus, 
         col = "#FF000080", ## col2hcl("red", alpha=.5),
         cex = 0.2, pch = 3)
  ## predicted
  lo <- model$lo
  ## nsd
  mass.ranges <- estimate.mass.range(xx$peptide.mhp.ident, ppmthreshold)
  ## from doHDMSePredictions 
  .o <- order(xx$precursor.retT.ident)
  .allpreds <- predict(model$lo, data.frame(precursor.retT.ident=xx$precursor.retT.ident), se=TRUE)
  .sd <- .allpreds$se.fit[.o] * sqrt(model$lo$n)
  .predicted <- xx$precursor.retT.ident - .allpreds$fit[.o]
  .lower <- .predicted - (nsd * .sd)
  .upper <- .predicted + (nsd * .sd)
  stopifnot(all(.lower <= .upper))
  rt.ranges <- cbind(.lower,
                     .upper)
  centers <- cbind(.predicted, xx$peptide.mhp.ident)
  ## points(centers, pch=".")
  rect(rt.ranges[,1], mass.ranges[,1],
       rt.ranges[,2], mass.ranges[,2],
       lwd=0.2)
  legend("topleft",
         c("Ident", "Quant", "EMRT", "Ident-EMRT matching"),
         col = c("#D4D4D480", "#8FBDDA80", "#FF000080", "#000000FF"),
         pch = c(19, 19, 3, 22),
         cex = .5,
         ncol = 4,
         bty = "n")
}

plot.all.features <- function(xx, mse) {
  par(mfrow=c(1,3))
  ylim <- range(range(xx$peptide.mhp.ident),
                range(xx$peptide.mhp.quant),
                range(mse$mwHPlus))
  xlim <- range(range(xx$precursor.retT.ident),
                range(xx$precursor.retT.quant),
                range(mse$rt_min))
  plot(xx$precursor.retT.ident, xx$peptide.mhp.ident,
       col = "#00000020", 
       pch = 19, xlim = xlim, ylim = ylim,
       xlab = "retention time", ylab = "precursor mass",
       main = "Identification final peptide")
  grid()
  plot(xx$precursor.retT.quant, xx$peptide.mhp.quant,
       col = "#00000020", 
       pch = 19, xlim = xlim, ylim = ylim,
       xlab = "retention time", ylab = "precursor mass",
       main = "Quantitation final peptide")
  grid()
  plot(mse$rt_min, mse$mwHPlus, 
       col = "#00000005", 
       pch = 19, xlim = xlim, ylim = ylim,
       xlab = "retention time", ylab = "precursor mass",
       main = "Quantitation Pep3D")
  grid()
  par(mfrow=c(1,1))
}

qPlot <- function(x, qtls=c(0.25, 0.5, 0.75, .90, 0.95, 0.99, 1), ...) {
  plot(sort(abs(x)),
       xaxt = "n", type = "l",
       xlab = "Percentage of data",
       ...)
  axis(1, at = seq(0,length(x), as.integer(length(x)/10)),
       labels = seq(0,100,10), cex.axis = .8)
  grid(nx = 0, ny = NULL)
  abline(v = seq(0,length(x), as.integer(length(x)/10)),
         col = "lightgray", lty = "dotted")
  qs <- getQs(x, qtls)
  segments(qs$x, 0, qs$x, qs$y, col = "red", lty = "dotted")
  segments(0, qs$y, qs$x, qs$y, col = "red", lty = "dotted")
}
