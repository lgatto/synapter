##' Performs a complete default analysis on the files defined
##' in \code{filenames}, creates a complete html report and
##' saves/exports all results as \code{csv} and \code{rds} files.
##' See details for a description of the pipeline and
##' \code{\link{Synapter}} for manual execution of individual steps.
##'
##' Data can be input as a \code{\link{Synapter}} object if available or
##' as a list of files (see \code{filenames}) that will be used to read the
##' data in. The html report and result files will be created in the
##' \code{outputdir} folder. All other input parameters have default values.
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
##' \code{link{createUniquePeptideDbRds}}). \code{identpeptide} can be a \code{csv}
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
##' synergise(object = synapterTiny, outputdir = output, outputfile = "synapter.html", grid.subset = 0.2)
##' htmlReport <- paste0("file:///", file.path(output, "synapter.html")) ## the result report
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
  .synergise(filenames=filenames, master=master, object=object,
             outputdir=outputdir, outputfile=outputfile,
             fdr=fdr, fdrMethod=fdrMethod, fpr=fpr,
             peplen=peplen, missedCleavages=missedCleavages, IisL=IisL,
             identppm=identppm, quantppm=quantppm, uniquepep=uniquepep, span.rt=span,
             grid.ppm.from=grid.ppm.from, grid.ppm.to=grid.ppm.to, grid.ppm.by=grid.ppm.by,
             grid.nsd.from=grid.nsd.from, grid.nsd.to=grid.nsd.to, grid.nsd.by=grid.nsd.by,
             grid.subset=grid.subset, grid.n=grid.n, grid.param.sel=grid.param.sel,
             mergedEMRTs=mergedEMRTs, template=template, verbose=verbose)
}


##' Performs a complete default analysis on the files defined
##' in \code{filenames}, creates a complete html report and
##' saves/exports all results as \code{csv} and \code{rds} files.
##' See details for a description of the pipeline and
##' \code{\link{Synapter}} for manual execution of individual steps.
##'
##' In contrast to \code{\link{synergise1}} \code{synergise2} extends the
##' default analysis and offers the follwing unique features:
##' \itemize{
##'   \item Performing 3D grid search (M/Z, Retention Time, Ion Mobility) for
##'         HDMSE data.
##'   \item Applying intensity correction.
##'   \item Filtering results by fragment matching.
##' }
##'
##' Data can be input as a \code{\link{Synapter}} object if available or
##' as a list of files (see \code{filenames}) that will be used to read the
##' data in. The html report and result files will be created in the
##' \code{outputdir} folder. All other input parameters have default values.
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
##'        a default \code{span.rt} value of 0.05.
##'  \item A grid search to optimise the width in retention time and mass
##'        tolerance (and ion mobility for HDMSE) for EMRTs matching is performed.
##'        The default grid search space is from 0.5 to 5 by 0.5 retention time
##'        model standard deviations (see \code{grid.nsd.from}, \code{grid.nsd.to} and
##'        \code{grid.nsd.by} parameters) and from 2 to 20 by 2 parts per
##'        million (ppm) for mass tolerance (see \code{grid.ppm.from},
##'        \code{grid.ppm.to} and \code{grid.ppm.by} parameters). If HDMSE data
##'        are used the search space is extended from ion mobility difference
##'        0.6 to 1.6 by 0.2 (see \code{grid.imdiffs.from}, \code{grid.imdiffs.to}
##'        and \code{grid.imdiffs.by}).
##'        The data can be subset using using an absolute number of features
##'        (see \code{grid.n}) or a fixed percentage (see \code{grid.subset}).
##'        The pair of optimal \code{nsd} and \code{ppm} is chosen
##'        (see \code{grid.param.sel} parameter).
##'  \item Fragment Matching is used to filter false-positive matches from the
##'        grid search using a default of 1 common peak for unique matches and
##'        at least a difference of 2 common peaks to choose the correct match
##'        out of non-unique matches (see \code{fm.minCommon} and
##'        \code{fm.minDelta}).
##'  \item The intensity is corrected by a Local Polynomial Regression Fitting
##'        (\code{\link{loess}} function for the \code{stats} package) using a
##'        default \code{span.int} value of 0.05.
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
##' \code{link{createUniquePeptideDbRds}}). \code{identpeptide} can be a \code{csv}
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
##' @param span.rt The loess span parameter for retention time correction.
##' Default is 0.05.
##' @param span.int The loess span parameter for intensity correction.
##' Default is 0.05.
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
##' @param grid.imdiffs.from Ion mobility difference grid starting value.
##' value. Default is 0.6.
##' @param grid.imdiffs.to Ion mobility difference grid ending value.
##' Default is 1.6.
##' @param grid.imdiffs.by Ion mobility difference grid step value.
##' Default is 0.2.
##' @param grid.subset Percentage of features to be used for the grid
##' search. Default is 1.
##' @param grid.n Absolute number of features to be used for the grid
##' search. Default is 0, i.e ignored.
##' @param grid.param.sel Grid parameter selection method. One of
##' \code{auto} (default), \code{details}, \code{model} or
##' \code{total}. See \code{\link{Synapter}} for details on these
##' selection methods.
##' @param fm.ppm Fragment Matching tolerance: Peaks in a range of \code{fm.ppm}
##' are considered as identical. Default is 25.
##' @param fm.minCommon Minimal number of peaks that unique matches need to have
##' in common. Default 1.
##' @param fm.minDelta Minimal difference in number of peaks that non-unique matches need
##' to have to be considered as true match. Default 2.
##' @param fm.fdr.unique Minimal FDR to select \code{fm.minCommon} automatically
##' (if both values are given the smaller is used). Default \code{Inf} means FDR
##' is not used.
##' @param fm.fdr.nonunique Minimal FDR to select \code{fm.minDelta} automatically
##' (if both values are given the smaller is used). Default \code{Inf} means FDR
##' is not used.
##' @param mergedEMRTs One of \code{"rescue"} (default), \code{"copy"}
##' or \code{"transfer"}. See the documentation for the
##' \code{findEMRTs} function in \code{\link{Synapter}} for details.
##' @param template A \code{character} full path to Rmd template.
##' @param verbose A \code{logical} indicating if progress output
##' should be printed to the console. Default is \code{TRUE}.
##' @return Invisibly returns an object of class \code{Synapter}.
##'         Used for its side effect of creating an html report of
##'         the run in \code{outputdir}.
##' @aliases synergise2 synergize2
##' @author Laurent Gatto, Sebastian Gibb
##' @references Bond N. J., Shliaha P.V., Lilley K.S. and Gatto
##' L. (2013) J. Prot. Research.
##' @examples
##' \dontrun{
##' library(synapterdata)
##' data(synobj2)
##' output <- tempdir() ## a temporary directory
##' synergise2(object = synobj2, outputdir = output, outputfile = "synapter.html")
##' htmlReport <- paste0("file:///", file.path(output, "synapter.html")) ## the result report
##' browseURL(htmlReport) ## open the report with default browser
##' }
synergise2 <-
synergize2 <- function(filenames,
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
                       span.rt = 0.05,
                       span.int = 0.05,
                       ## grid
                       grid.ppm.from = 2,
                       grid.ppm.to = 20,
                       grid.ppm.by = 2,
                       grid.nsd.from = 0.5,
                       grid.nsd.to = 5,
                       grid.nsd.by = 0.5,
                       grid.imdiffs.from = 0.6,
                       grid.imdiffs.to = 1.6,
                       grid.imdiffs.by = 0.2,
                       grid.subset = 1,
                       grid.n = 0,
                       grid.param.sel = c("auto", "model", "total", "details"),
                       ## fragment matching
                       fm.ppm = 25,
                       fm.minCommon = 1,
                       fm.minDelta = 2,
                       fm.fdr.unique = Inf,
                       fm.fdr.nonunique = Inf,
                       mergedEMRTs = c("rescue", "copy", "transfer"),
                       template = system.file("reports", "synergise2.Rmd",
                                              package="synapter"),
                       verbose = interactive()) {
  .synergise(filenames=filenames, master=master, object=object,
             outputdir=outputdir, outputfile=outputfile,
             fdr=fdr, fdrMethod=fdrMethod, fpr=fpr,
             peplen=peplen, missedCleavages=missedCleavages, IisL=IisL,
             identppm=identppm, quantppm=quantppm, uniquepep=uniquepep,
             span.rt=span.rt, span.int=span.int,
             grid.ppm.from=grid.ppm.from, grid.ppm.to=grid.ppm.to, grid.ppm.by=grid.ppm.by,
             grid.nsd.from=grid.nsd.from, grid.nsd.to=grid.nsd.to, grid.nsd.by=grid.nsd.by,
             grid.subset=grid.subset, grid.n=grid.n, grid.param.sel=grid.param.sel,
             fm.ppm=fm.ppm, fm.minCommon=fm.minCommon, fm.minDelta=fm.minDelta,
             fm.fdr.unique=fm.fdr.unique, fm.fdr.nonunique=fm.fdr.nonunique,
             mergedEMRTs=mergedEMRTs, template=template, verbose=verbose)
}



.synergise <- function(filenames,
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
                       span.rt = 0.05,
                       span.int = 0.05,
                       ## grid
                       grid.ppm.from = 2,
                       grid.ppm.to = 20,
                       grid.ppm.by = 2,
                       grid.nsd.from = 0.5,
                       grid.nsd.to = 5,
                       grid.nsd.by = 0.5,
                       grid.imdiffs.from = 0.6,
                       grid.imdiffs.to = 1.6,
                       grid.imdiffs.by = 0.2,
                       grid.subset = 1,
                       grid.n = 0,
                       grid.param.sel = c("auto", "model", "total", "details"),
                       ## fragment matching
                       fm.ppm = 25,
                       fm.minCommon = 1,
                       fm.minDelta = 2,
                       fm.fdr.unique = Inf,
                       fm.fdr.nonunique = Inf,
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
  } else {
    outputdir <- normalizePath(outputdir)
  }

  if (!file.exists(template)) {
    stop("'template' not found.")
  } else {
    template <- normalizePath(template)
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

  invisible(obj)
}
