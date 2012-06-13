##' Performs a complete default analysis on the files defined
##' in \code{filenames}, creates a complete html report and
##' saves/exports all results as \code{csv} and \code{rda} files.
##' See details for a description of the pipeline and
##' \code{\link{Synapter}} for manual execution of individual steps.
##'
##' Data can be input as a \code{\link{Synapter}} object if available or
##' as a list of files (see \code{filenames}) that will be used to read the
##' data in. If none of \code{object} and \code{filenames} are provided, file
##' section menus are open to select input files. The html report and result
##' files will be created in the \code{outputdir} folder. If not provided,
##' the destination can be selected through a selection menu. All other
##' input parameters have default values.
##'
##' The data processing and analysis pipeline is as follows:
##' \enumerate{
##'  \item If \code{uniquepep} is set to TRUE (default), only unique
##'        proteotypic identification and quantitation peptides are retained.
##'  \item Peptides are filtered for a FDR <= \code{fdr} (default is 0.01)
##'        using the "BH" method (see \code{fdr} and \code{fdrMethod}
##'        parameters for details).
##'  \item Peptide with a mass tolerance > 20 ppms (see \code{quantppm} and
##'        \code{identppm}) are filtered out.
##'  \item Peptides with a protein false positive rate (as reported by the
##'        PLGS software) > \code{fpr} are filtered out.
##'  \item Common identification and quantitation peptides are merged and a retention time
##'        model is created using the Local Polynomial Regression Fitting
##'        (\code{\link{loess}} function for the \code{stats} package) using
##'        a default \code{span} value of 0.05.
##'  \item A grid search to optimise the width in retention time and mass
##'        tolerance for EMRTs matching is performed. The default grid
##'        search space is from 0.5 to 5 by 0.5 retention time model standard
##'        deviations (see \code{grid.nsd.from}, \code{grid.nsd.to} and
##'        \code{grid.nsd.by} parameters) and from 5 to 20 by 2 parts per
##'        million (ppm) for mass tolerance (see \code{grid.ppm.from},
##'        \code{grid.ppm.to} and \code{grid.ppm.by} parameters).
##'        The data can be subset using using an absolute number of features
##'        (see \code{grid.n}) or a fixed percentage (see \code{grid.subset}).
##'        The pair of optimal \code{nsd} and \code{ppm} is chosen
##'        (see \code{grid.param.sel} parameter).
##'  \item The quantitation EMRTs are matched using the optimised parameters.
##' }
##' If a master identification file is used (\code{master} is set to \code{TRUE}, default
##' is \code{FALSE}), the relevant actions that have already been executed when
##' the file was created with \code{\link{makeMaster}} are not repeated here.
##' 
##' @title Synergise identification and quantitation results
##' @param filenames A named \code{list} of file names to be load. The
##' names must be \code{identpeptide}, \code{quantpeptide},
##' \code{quantpep3d} and \code{fasta}.  If missing, a dialog box opens
##' to select files interactively. \code{identpeptide} can be a \code{csv}
##' final peptide file (from PLGS) or a saved \code{"\linkS4class{MasterPeptides}"}
##' data object as created by \code{\link{makeMaster}} if working with
##' \emph{master} peptide data. To serialise the \code{"\linkS4class{MasterPeptides}"}
##' instance, use the \code{saveRDS} function, and file extenstion \code{rds}.
##' @param master A \code{logical} indicating if the identification final
##' peptide files are master (see \code{\link{makeMaster}}) or
##' regular files. Default is \code{FALSE}.
##' @param object An instance of class \code{Synapter} that will be
##' copied,  processed and returned. If \code{filenames} are also
##' provided, the latter and \code{object}'s \code{inputFiles} will be
##' checked for equality. 
##' @param outputdir A \code{character} with the full path to an
##' existing directory.
##' @param fdr Peptide false discovery rate. Default is 0.01.
##' @param fdrMethod P-value adjustment method. One of \code{"BH"}
##' (default) for Benjamini and HochBerg (1995), \code{"Bonferroni"}
##' for Bonferroni's single-step adjusted p-values for strong control
##' of the FWER and \code{"qval"} from the \code{qvalue} package. See
##' \code{\link{Synapter}} for references.
##' @param fpr Protein false positive rate. Default is 0.01.
##' @param identppm Identification mass tolerance (in ppm). Default is 20.
##' @param quantppm Quantitation mass tolerance (in ppm). Default is 20.
##' @param uniquepep A \code{logical} is length 1 indicating if only
##' unique peptides in the identification and quantitation peptides as well as unique
##' tryptic peptides as defined in the fasta file. Default is
##' \code{TRUE}.
##' @param span The loess span parameter. Default is 0.05.
##' @param grid.ppm.from Mass tolerance (ppm) grid starting
##' value. Default is 5.
##' @param grid.ppm.to Mass tolerance (ppm) grid ending value. Default
##' is 20.
##' @param grid.ppm.by Mass tolerance (ppm) grid step value. Default
##' is 2.
##' @param grid.nsd.from Number of retention time stdev grid starting
##' value. Default is 0.5.
##' @param grid.nsd.to Number of retention time stdev grid ending
##' value. Default is 5.
##' @param grid.nsd.by Number of retention time stdev grid step
##' value. Default is 0.5.
##' @param grid.subset Percentage of features to be used for the grid
##' search. Default is 1.
##' @param grid.n Absolute number of features to be used for the grid
##' search. Default is 0, i.e ignored.
##' @param grid.param.sel Grid parameter selection method. One of
##' \code{auto} (default), \code{details}, \code{model} or
##' \code{total}. See \code{\link{Synapter}} for details on these
##' selection methods.
##' @param css An optional path to a custom css file. If \code{NULL}
##' (default), uses \code{synapter.css}.
##' @param verbose A \code{logical} indicating if progress output
##' should be printed to the console. Default is \code{TRUE}.
##' @return Invisibly returns an object of class \code{Synapter}.
##'         Used for its side effect of creating an html report of
##'         the run in \code{outputdir}.
##' @aliases synergize
##' @author Laurent Gatto
##' @references Shliaha P.V., Bond N. J., Lilley K.S. and Gatto L., in prep.
##' @examples
##' output <- tempdir() ## a temporary directory
##' synapterTinyData()
##' synergise(object = synapterTiny, outputdir = output, grid.subset = 0.2)
##' htmlReport <- paste0("file:///", file.path(output, "index.html")) ## the result report
##' if (interactive())
##'   browseURL(htmlReport) ## open the report with default browser
synergise <- function(filenames,
                      master = FALSE,
                      object,
                      outputdir,
                      ## filtering
                      fdr = 0.01,
                      fdrMethod = c("BH", "Bonferroni", "qval" ),
                      fpr = 0.01,
                      identppm = 20,
                      quantppm = 20,
                      uniquepep = TRUE,
                      ## modelling
                      span = 0.05,
                      ## grid
                      grid.ppm.from = 5,
                      grid.ppm.to = 20,
                      grid.ppm.by = 2,
                      grid.nsd.from = 0.5,
                      grid.nsd.to = 5,
                      grid.nsd.by = 0.5,
                      grid.subset = 1,
                      grid.n = 0,
                      grid.param.sel = c("auto", "model", "total", "details"),
                      css = NULL,
                      verbose = TRUE) {
  fdrMethod <- match.arg(fdrMethod)
  grid.param.sel <- match.arg(grid.param.sel)
  if (missing(outputdir)) {
    outputdir <- tk_choose.dir("", caption = "Select the output folder")
  }
  if (!file.exists(outputdir))
    stop("'outputdir' not found.")
  if (!missing(object)) {
    if (class(object) != "Synapter")
      stop("Input 'object' must be of class 'Synapter'; found ", class(object), ".")    
    if (!missing(filenames)) {
      if (any(unlist(filenames) != inputFiles(object)))
        stop("File names do not match for 'filenames' and 'object' inputs.")
    }
    obj <- object$copy()
    rm(object) ## will not be deleted in .GlobalEnv
    gc()
  } else {
    obj <- Synapter(filenames, master)
  }  
  htmlfile <- paste0(outputdir, "/index.html")
  p <- openPage(htmlfile,
                link.css = ifelse(is.null(css),
                  paste0(system.file(package="synapter", dir = "extdata"),
                         "/synapter.css"),
                  css))  
  hwrite("Synapter report", p, heading = 1)

  hwrite("Input", p, heading = 2)
  txt <- rbind(c("Identification peptide:", basename(obj$IdentPeptideFile), nrow(obj$IdentPeptideData)),
               c("Quantitation peptide:",   basename(obj$QuantPeptideFile), nrow(obj$QuantPeptideData)),
               c("Quantitation Pep3D:",     basename(obj$QuantPep3DFile),   nrow(obj$QuantPep3DData)),
               c("Fasta file:",    basename(obj$DbFastaFile), ""))
  colnames(txt) <- c("", "File name", "Number of features")
  hwrite(txt, p, br = TRUE, row.bgcolor = '#ffdc98')

  prms <- rbind(c("Master", master),
                c("Peptide FDR", fdr),
                c("Protein FPR", fpr),
                c("Identification mass tolerance (ppm)", identppm),
                c("Quantitation mass tolerance (ppm)", quantppm),
                c("Loess span (RT modelling)", span),
                c("Filtering unique peptides", as.character(uniquepep)))
  
  colnames(prms) <- c("Parameter name", "Value")
  hwrite(prms, p, br = TRUE, row.bgcolor = '#ffdc98')
  
  grd <- rbind(c("Mass tolerance (ppm) start", grid.ppm.from),
               c("Mass tolerance (ppm) end", grid.ppm.to),
               c("Mass tolerance (ppm) by", grid.ppm.by),
               c("Retention time stdev start", grid.nsd.from),
               c("Retention time stdev end", grid.nsd.to),
               c("Retention time stdev by", grid.nsd.by),
               c("Feature proportion ", grid.subset),
               c("Number of features", grid.n),
               c("parameters selection", grid.param.sel))
  
  colnames(grd) <- c("Grid parameter", "Value")
  hwrite(grd, p, br = TRUE, row.bgcolor = '#ffdc98')
  
  writeIdentPeptides(obj, file = paste0(outputdir, "/IdentPeptides_nonFiltered.csv"))
  writeQuantPeptides(obj, file = paste0(outputdir, "/QuantPeptides_nonFiltered.csv"))
  
  hwrite("Filtering", p, heading = 2)
  
  hwrite("Peptide false discovery rate", p, heading = 3)
  fig.jpg <- makeFigurePath(outputdir, "Fig-pepScores.jpg")
  jpeg(fig.jpg$full)
  plotPepScores(obj)
  dev.off()
  hwriteImage(fig.jpg$relative, p, link = fig.jpg$relative, br = TRUE)
  
  hwrite(getPepNumbers(obj), p, br = TRUE, row.bgcolor = '#ffdc98', col.bgcolor = '#ffdc98')
  
  fig.jpg <- makeFigurePath(outputdir, "Fig-qValues.jpg")
  jpeg(fig.jpg$full)
  plotFdr(obj)
  dev.off()
  
  hwriteImage(fig.jpg$relative, p, link = fig.jpg$relative, br = TRUE)
  
  ## filter data
  if (uniquepep) {
    if (verbose)
      message("Keeping unique peptides...")
    filterUniqueDbPeptides(obj, verbose) 
  }
  setPepScoreFdr(obj, fdr = fdr)
  filterQuantPepScore(obj, method = fdrMethod) 
  if (!master)
    filterIdentPepScore(obj, method = fdrMethod)
  
  hwrite("Mass tolerance", p, heading = 3)
  
  fig1.jpg <- makeFigurePath(outputdir, "Fig-Ident-ppmError1.jpg")
  jpeg(fig1.jpg$full)
  plotPpmError(obj, what="Ident")
  dev.off()
  
  fig2.jpg <- makeFigurePath(outputdir, "Fig-Quant-ppmError1.jpg")    
  jpeg(fig2.jpg$full)
  plotPpmError(obj, what="Quant")
  dev.off()
  
  hwriteImage(c(fig1.jpg$relative, fig2.jpg$relative),
              p, border = 0, 
              link = c(fig1.jpg$relative, fig2.jpg$relative),
              br = TRUE)
  
  hwrite(getPpmErrorQs(obj), p, br = TRUE, row.bgcolor = '#ffdc98', col.bgcolor = '#ffdc98')
  
  filterQuantPpmError(obj, ppm = quantppm) 
  filterIdentPpmError(obj, ppm = identppm) 
  
  fig1.jpg <- makeFigurePath(outputdir, "Fig-Ident-ppmError2.jpg")
  jpeg(fig1.jpg$full)
  plotPpmError(obj, what="Ident")
  dev.off()
  
  fig2.jpg <- makeFigurePath(outputdir, "Fig-Quant-ppmError2.jpg")    
  jpeg(fig2.jpg$full)
  plotPpmError(obj, what="Quant")
  dev.off()
  
  hwriteImage(c(fig1.jpg$relative, fig2.jpg$relative),
              p, border = 0, 
              link = c(fig1.jpg$relative, fig2.jpg$relative),
              br = TRUE)
  
  hwrite(getPpmErrorQs(obj), p, br = TRUE, row.bgcolor = '#ffdc98', col.bgcolor = '#ffdc98')
  
  filterIdentProtFpr(obj, fpr = fpr)
  filterQuantProtFpr(obj, fpr = fpr)
  
  ## (3) Merge peptide sequences
  mergePeptides(obj)
  
  hwrite("Retention time modelling", p, heading = 2)
  
  ## (4) Retention time modelling
  hwrite("Retention time model", p, heading = 3)
  fig.jpg <- makeFigurePath(outputdir, "Fig-rtData.jpg")
  jpeg(fig.jpg$full, width = 960)
  plotRt(obj, what="data")
  dev.off()
  
  hwriteImage(fig.jpg$relative, p, link = fig.jpg$relative, br = TRUE)
  
  setLowessSpan(obj, span)
  modelRt(obj) ## the actual modelling
    
  fig1.jpg <- makeFigurePath(outputdir, "Fig-rtDiff.jpg")
  jpeg(fig1.jpg$full)
  plotRtDiffs(obj)
  dev.off()
  
  fig2.jpg <- makeFigurePath(outputdir, "Fig-rtModel.jpg")
  jpeg(fig2.jpg$full)
  plotRt(obj, what="model", nsd=1) ## better focus on model
  dev.off()
  
  hwriteImage(c(fig1.jpg$relative, fig2.jpg$relative),
              p, border = 0, 
              link = c(fig1.jpg$relative, fig2.jpg$relative),
              br = TRUE)

  hwrite(c("Retention time difference: ", getRtQs(obj)), p, br = TRUE,
         row.bgcolor = '#ffdc98', col.bgcolor = '#ffdc98')
  
  hwrite("Feature space", p, heading = 3)
  
  fig.jpg <- makeFigurePath(outputdir, "Fig-all-features.jpg")
  jpeg(fig.jpg$full, width = 1200)
  plotFeatures(obj, what="all")
  dev.off()
  
  hwriteImage(fig.jpg$relative, p, link = fig.jpg$relative, br = TRUE)
  
  setRtNsd(obj, 2)     ## RtNsd and PpmError are used for detailed plot
  setPpmError(obj, 10) ## if not set manually, default values are set automatically
  
  fig.svg <- makeFigurePath(outputdir, "Fig-some-features.svg")
  svg(fig.svg$full)
  plotFeatures(obj, what="some", xlim=c(30,50), ylim=c(1160, 1165))
  dev.off()
  
  hwriteImage(fig.svg$relative, p, link = fig.svg$relative, br = TRUE)
  
  hwrite("Using nsd = 2 and ppm error = 10 as an example. ", p, br = FALSE)
  hwrite("Optimal values are selected below.", p, br = TRUE)
  
  hwrite("Grid optimisation", p, heading = 2)
  if (verbose)
    message("Running grid search...")
  
  ## (5) Grid search to optimise EMRT matching parameters
  if (grid.n == 0) {
    hwrite(paste0("Using ", grid.subset * 100, "% of the data."), p, br = TRUE)
    searchGrid(obj,
               ppms = seq(grid.ppm.from, grid.ppm.to, grid.ppm.by),
               nsds = seq(grid.nsd.from, grid.nsd.to, grid.nsd.by),
               subset = grid.subset,
               verbose = verbose)
  } else {
    hwrite(paste0("Using ", grid.n , "features.", p, br = TRUE))
    searchGrid(obj,
               ppms = seq(grid.ppm.from, grid.ppm.to, grid.ppm.by),
               nsds = seq(grid.nsd.from, grid.nsd.to, grid.nsd.by),
               n = grid.n,
               verbose = verbose)
  }
  
  
  fig1.jpg <- makeFigurePath(outputdir, "Fig-grid-total.jpg")
  jpeg(fig1.jpg$full)
  plotGrid(obj, what = "total") ## plot the grid
  dev.off()
  
  fig2.jpg <- makeFigurePath(outputdir, "Fig-grid-model.jpg")
  jpeg(fig2.jpg$full)
  plotGrid(obj, what = "model")
  dev.off()
  
  fig3.jpg <- makeFigurePath(outputdir, "Fig-grid-details.jpg")
  jpeg(fig3.jpg$full)
  plotGrid(obj, what = "details")  ## TODO
  dev.off()

  detfile <- paste0(outputdir, "/grid_details.txt")
  sink(file = detfile)
  for (k in names(getGridDetails(obj))) {
    cat(k, "\n");
    print(getGridDetails(obj)[[k]]);
    cat("------------\n")
  }
  sink()  

  img <- hwriteImage(c(fig1.jpg$relative, fig2.jpg$relative, fig3.jpg$relative),
                     link = c(fig1.jpg$full, fig2.jpg$full, fig3.jpg$full),
                     width = 300, 
                     table = FALSE)

  tab <- hwrite(c(hwrite(getGrid(obj)[[1]], row.bgcolor = '#ffdc98', col.bgcolor = '#ffdc98'),
                  hwrite(getGrid(obj)[[2]], row.bgcolor = '#ffdc98', col.bgcolor = '#ffdc98'),
                  hwrite(getGrid(obj)[[3]], row.bgcolor = '#ffdc98', col.bgcolor = '#ffdc98')),
                table = FALSE)

  prm <-   hwrite(c(hwrite(getBestGridParams(obj)[[1]], row.bgcolor = '#ffdc98'),
                    hwrite(getBestGridParams(obj)[[2]], row.bgcolor = '#ffdc98'),
                    hwrite(getBestGridParams(obj)[[3]], row.bgcolor = '#ffdc98')),
                  table = FALSE)

  val <-   hwrite(c(hwrite(round(getBestGridValue(obj)[1], 3), row.bgcolor = '#ffdc98'),
                    hwrite(round(getBestGridValue(obj)[2], 3), row.bgcolor = '#ffdc98'),
                    hwrite(round(getBestGridValue(obj)[3], 3), row.bgcolor = '#ffdc98')),
                  table = FALSE)

  mat <- rbind(img, tab, val, prm)
  rownames(mat) <- c("Grid image", "Grid table", "Best value", "Best parameter(s)")

  hwrite(mat, p, br = TRUE, border = 0)
  
  hwrite("Grid details", p, link = basename(detfile), br = TRUE)
  hwrite("", p, br = TRUE)
  
  setBestGridParams(obj, what = grid.param.sel)
  hwrite(paste0("Setting best grid parameters using '", grid.param.sel, "':\n"), p, br = TRUE)
  hwrite(paste0(" - Mass tolerance (ppm): ", obj$PpmError), p, br = TRUE)
  hwrite(paste0(" - Number of retention time stdev: ", obj$RtNsd), p, br = TRUE)
  
  hwrite("EMRT matching", p, heading = 2)
  if (verbose)
    message("Matching EMRTs...")
  
  ## (6) Matching ident peptides and quant EMRTs
  findEMRTs(obj)
  
  fig.jpg <- makeFigurePath(outputdir, "Fig-emrt-matching.jpg")
  jpeg(fig.jpg$full)    
  plotEMRTtable(obj)
  dev.off()
  
  hwriteImage(fig.jpg$relative, p, link = fig.jpg$relative, br = TRUE)
  
  hwrite(as.matrix(getEMRTtable(obj)), p, br = TRUE, col.bgcolor = '#ffdc98')
  
  hwrite("Performance", p, heading = 2)
  
  perf <- performance(obj, verbose = verbose)
  hwrite(paste0("(S) Synapter: ", perf$Synapter, " EMRTs uniquely matched."), p, br = TRUE)
  hwrite(paste0("(I) Identification: ", perf$Ident, " peptides."), p, br = TRUE)
  hwrite(paste0("(Q) Quantitation: ", perf$Quant, " peptides."), p, br = TRUE)
  hwrite(paste0("Enrichment (S/Q): ", round(perf$Enrichment, 2), "%"), p, br = TRUE)
  hwrite("Overlap:", p, br = TRUE)
  hwrite(as.matrix(perf$VennCounts), p, br = TRUE, col.bgcolor = '#ffdc98')
    
  hwrite("Exported result files", p, heading = 2)
  
  ## (7) Exporting data to csv spreadsheets
  writeMergedPeptides(obj, what = "light", file = paste0(outputdir, "/MergedPeptidesLight.csv"))
  writeMergedPeptides(obj, what = "full", file = paste0(outputdir, "/MergedPeptidesFull.csv"))
  writeMatchedEMRTs(obj, what = "light", file = paste0(outputdir, "/MatchedPeptidesLight.csv")) 
  writeMatchedEMRTs(obj, what = "full", file = paste0(outputdir, "/MatchedPeptidesFull.csv"))                       
  writeIdentPeptides(obj, file = paste0(outputdir, "/IdentPeptides_filtered.csv"))
  writeQuantPeptides(obj, file = paste0(outputdir, "/QuantPeptides_filtered.csv"))
  
  hwrite("Identification peptide file (non filtered)", link = "IdentPeptides_nonFiltered.csv", p, br = TRUE)
  hwrite("Quantitation peptide file (non filtered)", link = "QuantPeptides_nonFiltered.csv", p, br = TRUE)
  hwrite("Identification peptide file (filtered)", link = "IdentPeptides_filtered.csv", p, br = TRUE)
  hwrite("Quantitation peptide file (filtered)", link = "QuantPeptides_filtered.csv", p, br = TRUE)
  hwrite("Merged peptides (light)", link = "MergedPeptidesLight.csv", p, br = TRUE)
  hwrite("Merged peptides (full)", link = "MergedPeptidesFull.csv", p, br = TRUE)
  hwrite("Matched peptides (light)", link = "MatchedPeptidesLight.csv", p, br = TRUE)
  hwrite("Matched peptides (full)", link = "MatchedPeptidesFull.csv", p, br = TRUE)    
  
  save(obj, file = paste0(outputdir, "/SynapterObject.rda"))
  
  hwrite("Saved R object (binary)", link = "SynapterObject.rda", p, br = TRUE)    
  
  ## Write full log
  hwrite("Full analysis log", p, heading = 2)
  hwrite(cbind(1:length(getLog(obj)), getLog(obj)), p, br = TRUE, col.bgcolor = '#ffdc98')
  
  closePage(p)

  if (verbose)
    message("Report written to '", basename(outputdir), "'.")
  
  invisible(obj)
}


