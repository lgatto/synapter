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
##'        \code{grid.nsd.by} parameters) and from 2 to 20 by 2 parts per
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
##' \code{quantpep3d} and \code{fasta} (could be an RDS file created by
##' \code{link{createUniquePeptideDbRds}}). If missing, a dialog box opens
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
##' @param outputfile A \code{character} with the file name for the report.
##' @param fdrMethod P-value adjustment method. One of \code{"BH"}
##' (default) for Benjamini and HochBerg (1995), \code{"Bonferroni"}
##' for Bonferroni's single-step adjusted p-values for strong control
##' of the FWER and \code{"qval"} from the \code{qvalue} package. See
##' \code{\link{Synapter}} for references.
##' @param fpr Protein false positive rate. Default is 0.01.
##' @param peplen Minimum peptide length. Default is 7.
##' @param missedCleavages Number of allowed missed cleavages. Default
##' is 0.
##' @param IisL If \code{TRUE} Isoleucin and Leucin are treated as
##' equal. In this case sequences like "ABCI", "ABCL" are removed because they
##' are not unqiue. If \code{FALSE} (default) "ABCI" and "ABCL" are reported as
##' unique.
##' @param identppm Identification mass tolerance (in ppm). Default is 20.
##' @param quantppm Quantitation mass tolerance (in ppm). Default is 20.
##' @param uniquepep A \code{logical} is length 1 indicating if only
##' unique peptides in the identification and quantitation peptides as
##' well as unique tryptic peptides as defined in the fasta
##' file. Default is \code{TRUE}.
##' @param fdr Peptide false discovery rate. Default is 0.01.
##' @param span The loess span parameter. Default is 0.05.
##' @param grid.ppm.from Mass tolerance (ppm) grid starting
##' value. Default is 2.
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
##' @param mergedEMRTs One of \code{"rescue"} (default), \code{"copy"}
##' or \code{"transfer"}. See the documentation for the
##' \code{findEMRTs} function in \code{\link{Synapter}} for details.
##' @param template A \code{character} full path to Rmd template.
##' @param verbose A \code{logical} indicating if progress output
##' should be printed to the console. Default is \code{TRUE}.
##' @return Invisibly returns an object of class \code{Synapter}.
##'         Used for its side effect of creating an html report of
##'         the run in \code{outputdir}.
##' @aliases synergize synergise1 synergize1
##' @author Laurent Gatto, Sebastian Gibb
##' @references Bond N. J., Shliaha P.V., Lilley K.S. and Gatto
##' L. (2013) J. Prot. Research.
##' @examples
##' output <- tempdir() ## a temporary directory
##' synapterTinyData()
##' synergise(object = synapterTiny, outputdir = output, grid.subset = 0.2)
##' htmlReport <- paste0("file:///", file.path(output, "index.html")) ## the result report
##' \dontrun{
##' browseURL(htmlReport) ## open the report with default browser
##' }
synergise <-
synergise1 <-
synergize <-
synergize1 <- function(filenames,
                       master = FALSE,
                       object,
                       outputdir,
                       outputfile=paste0("synapter_report_",
                                         strftime(Sys.time(), "%Y%m%d-%H%M%S"),
                                         ".html"),
                       ## filtering
                       fdr = 0.01,
                       fdrMethod = c("BH", "Bonferroni", "qval" ),
                       fpr = 0.01,
                       peplen = 7,
                       missedCleavages = 0,
                       IisL = FALSE,
                       identppm = 20,
                       quantppm = 20,
                       uniquepep = TRUE,
                       ## modelling
                       span = 0.05,
                       ## grid
                       grid.ppm.from = 2,
                       grid.ppm.to = 20,
                       grid.ppm.by = 2,
                       grid.nsd.from = 0.5,
                       grid.nsd.to = 5,
                       grid.nsd.by = 0.5,
                       grid.subset = 1,
                       grid.n = 0,
                       grid.param.sel = c("auto", "model", "total", "details"),
                       mergedEMRTs = c("rescue", "copy", "transfer"),
                       template = system.file("reports", "synergise1.Rmd",
                                              package="synapter"),
                       verbose = interactive()) {

  fdrMethod <- match.arg(fdrMethod)
  grid.param.sel <- match.arg(grid.param.sel)
  mergedEMRTs <- match.arg(mergedEMRTs)
  gridDetails <- NULL

  if (missing(outputdir)) {
    stop("No 'outputdir' given.")
  }

  if (!file.exists(outputdir)) {
    stop("'outputdir' not found.")
  }

  if (!file.exists(template)) {
    stop("'template' not found.")
  }

  if (!missing(object)) {
    if (class(object) != "Synapter") {
      stop("Input 'object' must be of class 'Synapter'; found ", class(object), ".")
    }
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

  render(template, output_file=outputfile, output_dir=outputdir, quiet=!verbose)

  writeMergedPeptides(obj, file=file.path(outputdir, "MergedPeptides.csv"))
  writeMatchedEMRTs(obj, file=file.path(outputdir, "MatchedPeptides.csv"))
  writeIdentPeptides(obj, file=file.path(outputdir, "IdentPeptides.csv"))
  writeQuantPeptides(obj, file=file.path(outputdir, "QuantPeptides.csv"))
  write.table(gridDetails, file=file.path(outputdir, "GridDetails.txt"))
  save(obj, file=file.path(outputdir, "SynapterObject.rda"))

  invisible(obj)
}
