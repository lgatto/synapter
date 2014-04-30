\name{Synapter}

% constructor
\alias{Synapter}
\alias{class:Synapter}

% descriptors
\alias{show,Synapter-method}
\alias{dim,Synapter-method}
\alias{inputFiles,Synapter-method}
\alias{inputFiles}
\alias{getLog,Synapter-method}
\alias{getLog}

% analysers
\alias{mergePeptides,Synapter-method}
\alias{mergePeptides}
\alias{modelRt,Synapter-method}
\alias{modelRt}
\alias{findEMRTs,Synapter-method}
\alias{findEMRTs}
\alias{searchGrid,Synapter-method}
\alias{searchGrid}
\alias{crossMatching,Synapter-method}
\alias{crossMatching}

% setters and getters
\alias{getGrid,Synapter-method}
\alias{getGrid}
\alias{getGridDetails,Synapter-method}
\alias{getGridDetails}
\alias{getBestGridValue,Synapter-method}
\alias{getBestGridValue}
\alias{getBestGridParams,Synapter-method}
\alias{getBestGridParams}
\alias{setBestGridParams,Synapter-method}
\alias{setBestGridParams}
\alias{setPepScoreFdr,Synapter-method}
\alias{setPepScoreFdr}
\alias{getPepScoreFdr,Synapter-method}
\alias{getPepScoreFdr}
\alias{setIdentPpmError,Synapter-method}
\alias{setIdentPpmError}
\alias{getIdentPpmError,Synapter-method}
\alias{getIdentPpmError}
\alias{setQuantPpmError,Synapter-method}
\alias{setQuantPpmError}
\alias{getQuantPpmError,Synapter-method}
\alias{getQuantPpmError}
\alias{setPpmError,Synapter-method}
\alias{setPpmError}
\alias{setLowessSpan,Synapter-method}
\alias{setLowessSpan}
\alias{getLowessSpan,Synapter-method}
\alias{getLowessSpan}
\alias{setRtNsd,Synapter-method}
\alias{setRtNsd}
\alias{getRtNsd,Synapter-method}
\alias{getRtNsd}
\alias{getPpmErrorQs,Synapter-method}
\alias{getPpmErrorQs}
\alias{getRtQs,Synapter-method}
\alias{getRtQs}
\alias{getPepNumbers,Synapter-method}
\alias{getPepNumbers}
\alias{showFdrStats,Synapter-method}
\alias{showFdrStats}
\alias{setProtFpr,Synapter-method}
\alias{setProtFpr}
\alias{getProtFpr,Synapter-method}
\alias{getProtFpr}
\alias{getEMRTtable,Synapter-method}
\alias{getEMRTtable}
\alias{setCrossMatchingPpmTolerance,Synapter-method}
\alias{setCrossMatchingPpmTolerance}
\alias{getCrossMatchingPpmTolerance,Synapter-method}
\alias{getCrossMatchingPpmTolerance}
\alias{setCrossMatchingMinimalNumberOfCommonPeaks,Synapter-method}
\alias{setCrossMatchingMinimalNumberOfCommonPeaks}
\alias{getCrossMatchingMinimalNumberOfCommonPeaks,Synapter-method}
\alias{getCrossMatchingMinimalNumberOfCommonPeaks}
\alias{performance,Synapter-method}
\alias{performance2,Synapter-method}
\alias{performance}
\alias{performance2}

% filters
\alias{filterPeptideLength,Synapter-method}
\alias{filterPeptideLength}
\alias{filterUniqueDbPeptides,Synapter-method}
\alias{filterUniqueDbPeptides}
\alias{filterUniqueQuantDbPeptides,Synapter-method}
\alias{filterUniqueQuantDbPeptides}
\alias{filterUniqueIdentDbPeptides,Synapter-method}
\alias{filterUniqueIdentDbPeptides}
\alias{filterQuantPepScore,Synapter-method}
\alias{filterQuantPepScore}
\alias{filterIdentPepScore,Synapter-method}
\alias{filterIdentPepScore}
\alias{filterIdentPpmError,Synapter-method}
\alias{filterIdentPpmError}
\alias{filterQuantPpmError,Synapter-method}
\alias{filterQuantPpmError}
\alias{filterIdentProtFpr,Synapter-method}
\alias{filterIdentProtFpr}
\alias{filterQuantProtFpr,Synapter-method}
\alias{filterQuantProtFpr}
\alias{filterMatchedEMRTsByCommonPeaks,Synapter-method}
\alias{filterMatchedEMRTsByCommonPeaks}

% plotting
\alias{plotEMRTtable,Synapter-method}
\alias{plotEMRTtable}
\alias{plotFeatures,Synapter-method}
\alias{plotFeatures}
\alias{plotGrid,Synapter-method}
\alias{plotGrid}
\alias{plotPepScores,Synapter-method}
\alias{plotPepScores}
\alias{plotPpmError,Synapter-method}
\alias{plotPpmError}
\alias{plotFdr,Synapter-method}
\alias{plotFdr}
\alias{plotRt,Synapter-method}
\alias{plotRt}
\alias{plotRtDiffs,Synapter-method}
\alias{plotRtDiffs}
\alias{plotCrossMatchingSummary,Synapter-method}
\alias{plotCrossMatchingSummary}

% importers
\alias{loadSpectrumXmlFiles,Synapter-method}
\alias{loadSpectrumXmlFiles}
\alias{loadFragmentCsvFiles,Synapter-method}
\alias{loadFragmentCsvFiles}

% exporters
\alias{writeMatchedEMRTs,Synapter-method}
\alias{writeMatchedEMRTs}
\alias{writeMergedPeptides,Synapter-method}
\alias{writeMergedPeptides}
\alias{writeIdentPeptides,Synapter-method}
\alias{writeIdentPeptides}
\alias{writeQuantPeptides,Synapter-method}
\alias{writeQuantPeptides}


% other
\alias{as.MSnSet.Synapter}

\docType{class}
\encoding{UTF-8}

\title{Class "Synapter"}


\description{
  A reference class to store, manage and process Synapt G2 data
  to combine identification and quantitation results.

  The data, intermediate and final results are stored together
  in such a ad-how container called a class. In the frame of the
  analysis of a set of 3 data files, namely as identification peptide,
  a quantitation peptide and a quantitation Pep3D, such a container is created
  and populated, updated according to the user's instructions
  and used to display and export results.

  The functionality of the \code{synapter} package implemented in the
  \code{Synapter} class in described in the \emph{Details} section
  below. Documentation for the individual methods is provided in the
  \emph{Methods} section. Finally, a complete example of an analysis is
  provided in the \emph{Examples} section, at the end of this document.

  See also papers by Shliaha et al. for details about ion mobility
  separation and the manuscript describing the \code{synapter}
  methodology.
}

\usage{
Synapter(filenames, master) ## creates an instance of class 'Synapter'
}

\arguments{
  \item{filenames}{A named \code{list} of file names to be load. The
    names must be 'identpeptide', 'quantpeptide', 'quantpep3d' and 'fasta'
    (could be an RDS file created by \code{link{createUniquePeptideDbRds}}).
    If missing, dialog boxes pop up to select the four files
    manually.
    \code{identpeptide} can be a \code{csv} final peptide file (from
    PLGS) or a saved \code{"\linkS4class{MasterPeptides}"} data object as
    created by \code{\link{makeMaster}} if working with \emph{master}
    peptide data. To serialise the \code{"\linkS4class{MasterPeptides}"}
    instance, use the \code{saveRDS} function, and file extenstion
    \code{rds}.
  }
  \item{master}{A \code{logical} that defines if the identification file
    is a \emph{master} file. See \code{\link{makeMaster}} for details
    about this strategy.
  }
}

\details{
  A \code{Synapter} object logs every operation that is applied to
  it. When displayed with \code{show} or when the name of the instance
  is typed at the R console, the original input file names, all
  operations and resulting the size of the respective data are
  displayed. This allows the user to trace the effect of respective
  operations.

  \subsection{Loading the data}{
    The construction of the data and analysis container, technically
    defined as an instance or object of class \code{Synapter}, is created
    with the \code{Synapter} constructor. This function opens four dialog
    boxes for the user to point to the input files,
    namely (and in that order), the identification final peptide csv file, the
    quantitation final peptide  csv file and the quantitation Pep3D csv file (as exported
    from the PLGS software) and the fasta file use for peptide
    identification. The latter could be an RDS file generated by
    \code{link{createUniquePeptideDbRds}}, too.
    The files are read and the data is stored in the newly
    created \code{Synapter} instance. The file names can also be
    specified as a named list with names 'identpeptide', 'quantpeptide'
    and 'quantpep3d' respectively.

    The final peptide files are filtered
    to retain peptides with \code{matchType} corresponding to
    \code{PepFrag1} and \code{PepFrag2}, corresponding to unmodified
    round 1 and 2 peptide identification. Other types, like
    \code{NeutralLoss_NH3}, \code{NeutralLoss_H20}, \code{InSource},
    \code{MissedCleavage} or \code{VarMod} are not considered in the rest
    of the analysis. The quantitation Pep3D data is filtered to retain
    \code{Function} equal to \code{1} and unique quantitation spectrum ids,
    i.e. unique entries for multiple charge states or isotopes of an EMRT
    (exact mass-retention time features).

    Then, p-values for \code{Regular} peptides are computed based on
    the \code{Regular} and \code{Random} database types score
    distributions, as described in  \cite{\enc{Käll}{Kall} et al.,
      2008a}. Only unique peptide sequences are taken into account:
    in case of duplicated peptides, only one entry is kept.
    Empirical p-values are adjusted using \cite{Bonferroni}
    and \cite{Benjamini and Hochberg, 1995} (\code{multtest} package)
    and q-values are computed using the \code{qvalue} package
    (\cite{Storey JD and Tibshirani R., 2003 and \enc{Käll}{Kall} et
      al., 2008b}). Only \code{Regular} entries are stored in the
    resulting data for subsequent analysis.

    The data tables can be exported as \code{csv} spreadsheets with the
    \code{writeIdentPeptides} and \code{writeQuantPeptides} methods.
  }

  \subsection{Filtering identification and quantitation peptide}{
    The first step of the analysis aims to match reliable peptide.
    The final peptide datasets are
    filtered based on the FDR (BH is default) using the
    \code{filterQuantPepScore} and  \code{filterIdentPepScore}
    methods. Several plots are provided to illustrate peptide score
    densities (from which p-values are estimated, \code{plotPepScores};
    use \code{getPepNumbers} to see how many peptides were available) and
    q-values (\code{plotFdr}).

    Peptides matching to multiple proteins in the fasta file (non-unique
    tryptic identification and quantitation peptides) can be
    discarded with the \code{filterUniqueDbPeptides} method. One can
    also filter on the peptide length using \code{filterPeptideLength}.

    Another filtering criterion is mass accuracy. Error tolerance
    quantiles (in ppm, parts per million) can be visualised with the
    \code{plotPpmError} method. The values can be retrieved with
    \code{getPpmErrorQs}. Filtering is then done separately for
    identification and quantitation peptide data using
    \code{filterIdentPpmError} and \code{filterQuantPpmError}
    respectively. The previous plotting functions can be used again to
    visualise the resulting distribution.

    Filtering can also be performed at the level of protein false
    positive rate, as computed by the PLGS application
    (\code{protein.falsePositiveRate} column), which counts the
    percentage of decoy proteins that have been identified prior to the
    regular protein of interest. This can be done with the
    \code{filterIdentProtFpr} and \code{filterQuantProtFpr} methods.
    Note that this field is erroneously called a false positive rate in
    the PLGS software and the associated manuscript; it is a false
    discovery rate.
  }

  \subsection{Merging identification and quantitation peptides}{
    Common and reliable identification and quantitation peptides are
    then matched based on their sequences and merged using the
    \code{mergePeptides} method.
  }

  \subsection{Retention time modelling}{
    Systematic differences between identification features and
    quantitation features retention times are modelled by
    fitting a local regression (see the \code{\link{loess}} function for
    details), using the \code{modelRt} method. The smoothing parameter,
    or number of neighbour data points used the for local fit, is
    controlled by the \code{span} parameter that can be set in the above
    method.

    The effect of this parameter can be observed with the \code{plotRt}
    method, specifying \code{what = "data"} as parameters. The resulting
    model can then be visualised with the above method specifying
    \code{what = "model"}, specifying up to 3 number of standard
    deviations to plot. A histogram of retention time differences can
    be produced with the \code{plotRtDiffs} method.

    Mention \code{plotFeatures} here.
  }

  \subsection{Grid search to optimise matching tolerances}{
    Matching of identification peptides and quantitation EMRTs is done
    within a mass tolerance in parts per million (ppm) and the modelled
    retention time +/- a certain number of standard deviations.
    To help in the choice of these two parameters, a grid search over a
    set of possible values is performed and performance metrics are
    recorded, to guide in the selection of a 'best' pair of parameters.

    The following metrics are computed:
    (1) the percentage of identification
    peptides that matched a single quantitation EMRT (called \code{prcntTotal}),
    (2) the percentage of identification peptides used in the retention time
    model that matched the quantitation EMRT corresponding to the
    correct quantitation peptide in ident/quant pair of the model
    (called \code{prcntModel})
    and
    (3) the detailed about the matching of the features used for
    modelling (accessible with \code{getGridDetails}) and the
    corresponding \code{details} grid that reports the percentage of
    correct unique assignments.
    The detailed grid results specify the number of non
    matched identification peptides (0), the number of correctly (1) or
    wrongly (-1) uniquely matched identification peptides, the number of
    identification peptides that matched 2 or more peptides including
    (2+) or excluding (2-) the correct quantitation equivalent are also
    available.

    See the next section for additional details about how matching.
    The search is performed with the \code{searchGrid} method, possibly
    on a subset of the data (see Methods and Examples sections for
    further details).

    The parameters used for matching can be set manually with
    \code{setPpmError} and \code{setRtNsd} respectively, or using
    \code{setBestGridParams} to apply best parameters as defined using
    the grid search. See example and method documentation for details.
  }

  \subsection{Identification transfer: matching identification peptides and quantitation EMRTs}{
    The identification peptide - quantitation EMRT matching, termed
    identification transfer, is performed using the best parameters, as
    defined above with a grid search, or using user-defined parameters.

    Matching is considered successful when one and only one EMRT is
    found in the mass tolerance/retention time window defined by the
    error ppm and number of retention time standard deviations
    parameters. The values of uniquely matched EMRTs are reported in the
    final \emph{matched} dataframe that can be exported (see below). If
    however, none or more than one EMRTs are matched, 0 or the number of
    matches are reported.

    As identification peptides are serially individually matched to 'close'
    EMRTs, it is possible for peptides to be matched the same EMRT
    independently. Such cases are reported as -1 in the results
    dataframes.

    The results can be assess using the \code{plotEMRTtable} (or
    \code{getEMRTtable} to retrieve the values) and \code{performace}
    methods. The former shows the number of identification peptides assigned to
    none (0), exactly 1 (1) or more (> 2) EMRTs.
    The latter method reports matched identification peptides, the number of
    (q-value and protein FPR filtered) identification and quantitation peptides.
    Matched EMRT and quantitation peptide numbers are then compared
    calculating the synapter enrichment (100 * ( synapter - quant ) / quant)
    and Venn counts.
  }

  \subsection{Importing Spectra and Fragments data}{
    \dQuote{Spectrum.xml} and \dQuote{final_fragment.csv} files could be
    imported via \code{loadSpectrumXmlFiles} respective \code{loadFragmentCsvFiles}.
  }

  \subsection{Cross Matching}{
    After importing spectra and fragment data peaks that matching between
    identification fragments and quantitation fragments, between
    identification spectrum and quantitation fragments, and between
    quantitation spectrum and identification fragments using the
    \code{crossMatching} method. Use \code{setCrossMatchingPpmTolerance} to set
    the maximal allowed tolerance for considering a peak as identical.
    There are two different methods to visualise the results of the cross
    matching procedure. \code{plotCrossMatching} plots the spectra and fragments
    for each considered cross matching pair. \code{plotCrossMatchingSummary}
    creates three plots showing the performance of the cross matching procedure.
    The left panel shows the F1-score and the FDR calculated for an equal sized
    subset of true and false matches. The panel in the middle shows the
    development of the true-positives, false-positives, true-negatives and
    false-negatives using an increasing number of common peaks as threshold. On
    the right panel boxplots for all matches and the sampled matches are shown.
    Use \code{setCrossMatchingMinimalNumberOfCommonPeaks} to set the number of
    common peaks between identification and quantitation fragments
    that have to be at least observed to be considered as unique match. The
    \code{filterMatchedEMRTsByCommonPeaks} method could be applied to turn non
    unique matches into unique ones in respect to the results of the cross
    matching procedure. If there are still multiple matches that have more
    common peaks than set by \code{setCrossMatchingMinimalNumberOfCommonPeaks}
    they are kept as non unique matches.
  }

  \subsection{Exporting and saving data}{
    The merged identification and quantitation peptides can be exported
    to csv using the  \code{writeMergedPeptides} method. Similarly, the
    matched identification peptides and quantitation EMRTs are exported
    with \code{writeMatchedEMRTs}.

    Complete \code{Synapter} instances can be serialised with
    \code{save}, as any R object, and reloaded with \code{load} for
    further analysis.
  }
}


\section{Methods}{

  \subsection{Analysis methods}{

    \describe{
      \item{mergePeptides}{\code{signature(object = "Synapter")}: Merges
	quantitation and identification final peptide data, used to
	perform retention time modelling (see \code{modelRt} below).
      }

      \item{modelRt}{\code{signature(object = "Synapter",
	  span = "numeric")}: Performs local polynomial regression
	fitting (see \code{\link{loess}}) retention time alignment using
	\code{span} parameter value to control the degree of smoothing.
      }

      \item{findEMRTs}{\code{signature(object = "Synapter", ppm =
	  "numeric", nsd = "numeric", mergedEMRTs = c("rescue", "copy",
	  "transfer"))}:
	Finds EMRTs matching identification peptides using \code{ppm}
	mass tolerance and \code{nsd} number of retention time standard
	deviations. The last two parameters are optional if previously
	set with \code{setPpmError} and \code{setRtNsd} or, better,
	\code{setBestGridParams} (see below). The \code{mergedEMRTs}
	parameter defined the behaviour for those high confidence
	features that where identified in both identification and
	quantitation acquisitions and used for the retention time model
	(see \code{mergePeptides}). Prior to version 1.1.1, these
	features were transferred from the quantitation pep3d file if
	unique matches were found, as any feature (\code{"transfer"}).
	As a result, those matching 0 or > 1 EMRTs were
	quantified as \code{NA}. The default is now to \code{"rescue"}
	the quantitation values of these by directly retrieving the data
	from the quantification peptide data. Alternatively, the
	quantitation values for these features can be directly taken
	from the quantitation peptide data using \code{"copy"}, thus
	effectively bypassing identification transfer.
      }

      \item{searchGrid}{\code{signature(object="Synapter",
	  ppms="character", nsds="character", subset="numeric", n =
	  "numeric", verbose="logical")}: Performs a grid search. The
	grid is defined by the \code{ppm} and \code{nsd} numerical
	vectors, representing the sequence of values to be
	tested. Default are \code{seq(5, 20, 2)} and \code{seq(0.5, 5,
	  0.5)} respectively. \code{subset} and \code{n} allow to use a
	randomly  chosen subset of the data for the grid search to
	reduce search time. \code{subset} is a numeric, between 0 and 1,
	describing the percentage of data to be used; \code{n} specifies
	the absolute number of feature to use. The default is to use all
	data. \code{verbose} controls whether textual output should be
	printed to the console. (Note, the \code{mergedEMRTs} value used
	in internal calls to \code{findEMRTs} is \code{"transfer"} - see
	\code{findEMRTs} for details).
      }

      \item{crossMatching}{\code{signature(object="Synapter",
      ppm = "numeric", verbose = "logical"}:
      Performs a cross matching between spectra and fragment data.
      The \code{ppm} argument controls the tolerance that is used to consider
      two peaks (MZ values) as identical.  If \code{verbose} is \code{TRUE}
      (default) a progress bar is shown.  }
  }

  \subsection{Methods to display, access and set data}{

    \describe{

      \item{show}{\code{signature(object = "Synapter")}: Display
	\code{object} by printing a summary to the console.
      }

      \item{dim}{\code{signature(x="Synapter")}: Returns a \code{list}
	of dimensions for the identification peptide, quantitation
	peptide, merged peptides and matched features data sets.
      }

      \item{inputFiles}{\code{signature(object="Synapter")}: Returns a
	\code{character} of length 4 with the names of the input files
	used as \code{identpeptide}, \code{quantpeptide},
	\code{quantpep3d} and \code{fasta}.
      }

      \item{getLog}{\code{signature(object="Synapter")}: Returns a
	\code{character} of variable length with a summary of processing
	undergone by \code{object}.
      }

      \item{getGrid}{\code{signature(object="Synapter", digits =
	  "numeric")}: Returns a named \code{list} of length 3 with the
	precent of total (\code{prcntTotal}), percent of model
	(\code{prcntModel}) and detailed (\code{details}) grid search
	results. The \code{details} grid search reports the proportion
	of correctly assigned features (+1) to all unique assignments
	(+1 and -1). Values are rounded to 3 \code{digits} by default.
      }

      \item{getGridDetails}{\code{signature(object="Synapter")}: Returns
	a \code{list} of number of ..., -2, -1, 0, +1, +2, ... results
	found for each of the \code{ppm}/\code{nsd} pairs tested during
	the grid search.
      }

      \item{getBestGridValue}{\code{signature(object="Synapter")}:
	Returns a named \code{numeric} of length 3 with best grid values
	for the 3 searches. Names are \code{prcntTotal},
	\code{prcntModel} and \code{details}.
      }

      \item{getBestGridParams}{\code{signature(object="Synapter")}:
	Returns a named \code{list} of matrices (\code{prcntTotal},
	\code{prcntModel} and \code{details}). Each matrix gives the
	\code{ppm} and \code{nsd} pairs that yielded the best grid
	values (see \code{getBestGridValue} above).
      }

      \item{setBestGridParams}{\code{signature(object="Synapter",
	  what="character")}: This methods set the best parameter pair,
	as determined by \code{what}. Possible values are \code{auto}
	(default), \code{model}, \code{total} and \code{details}. The 3
	last ones use the (first) best parameter values as reported by
	\code{getBestGridParams}. \code{auto} uses the best \code{model}
	parameters and, if several best pairs exists, the one that
	maximises \code{details} is selected.
      }

      \item{setPepScoreFdr}{\code{signature(object="Synapter", fdr =
	  "numeric")}: Sets the peptide score false discovery rate
	(default is 0.01) threshold used by \code{filterQuantPepScore}
	and \code{filterIdentPepScore}.
      }

      \item{getPepScoreFdr}{\code{signature(object="Synapter")}: Returns
	the peptide false discrovery rate threshold.
      }

      \item{setIdentPpmError}{\code{signature(object="Synapter", ppm =
	  "numeric")}: Set the identification mass tolerance to
	\code{ppm} (default 10).
      }

      \item{getIdentPpmError}{\code{signature(object="Synapter")}:
	Returns the identification mass tolerance.
      }

      \item{setQuantPpmError}{\code{signature(object="Synapter", ppm =
	  "numeric")}: Set the quantitation mass tolerance to \code{ppm}
	(default 10).
      }

      \item{getQuantPpmError}{\code{signature(object="Synapter")}:
	Returns the quantitation mass tolerance.
      }

      \item{setPpmError}{\code{signature(object="Synapter", ppm =
	  "numeric")}: Sets the identification and quantitation mass
	tolerance \code{ppm} (default is 10).
      }

      \item{setLowessSpan}{\code{signature(object="Synapter", span =
	  "numeric")}: Sets the \code{loess} \code{span} parameter;
	default is 0.05.
      }

      \item{getLowessSpan}{\code{signature(object="Synapter")}: Returns
	the \code{span} parameter value.
      }

      \item{setRtNsd}{\code{signature(object="Synapter", nsd =
	  "numeric")}: Sets the retention time tolerance \code{nsd},
	default is 2.
      }

      \item{getRtNsd}{\code{signature(object="Synapter")}: Returns the
	value of the retention time tolerance \code{nsd}.
      }

      \item{getPpmErrorQs}{\code{signature(object="Synapter", qs =
	  "numeric", digits = "numeric")}: Returns the mass tolerance
	\code{qs} quantiles (default is \code{c(0.25, 0.5, 0.75, seq(0.9,
	  1, 0.01)}) for the identification and quantitation
	peptides. Default is 3 \code{digits}.
      }

      \item{getRtQs}{\code{signature(object="Synapter", qs =
	  "numeric", digits = "numeric")}: Returns the retention time
	tolerance \code{qs} quantiles (default is \code{c(0.25, 0.5,
	  0.75, seq(0.9, 1, 0.01)}) for the identification and
	quantitation peptides. Default is 3 \code{digits}.
      }

      \item{getPepNumbers}{\code{signature(object="Synapter")}: Returns
	the number of regular and random quantitation and identification
	peptide considered for p-value calculation and used to plot the
	score densities (see \code{plotPepScores}). Especially the
	difference between random and regular entries are informative in
	respect with the confidence of the random scores distribution.
      }

      \item{setCrossMatchingPpmTolerance}{\code{signature(object="Synapter",
        ppm = "numeric")}: Sets the cross matching mass tolerance \code{ppm}
        (default is 25).
      }

      \item{getCrossMatchingPpmTolerance}{\code{signature(object="Synapter")}:
        Returns the cross matching mass tolerance in ppm.
      }

      \item{setCrossMatchingMinimalNumberOfCommonPeaks}{
        \code{signature(object="Synapter", n = "integer")}:
        Sets the minimal number of common peaks for cross matching.
      }

      \item{getCrossMatchingPpmTolerance}{\code{signature(object="Synapter")}:
        Returns the minimal number of common peaks for cross matching.
      }

      \item{showFdrStats}{\code{signature(object="Synapter", k =
	  "numeric")}: Returns a named \code{list} of length
	2 with the proportion of identification and quantitation
	peptides that are considered significant with a threshold of
	\code{k} (default is \code{c(0.001, 0.01, 0.5, 0.1)}) using raw
	and adjusted p-values/q-values.
      }

      \item{getEMRTtable}{\code{signature(object="Synapter")}: Returns a
	\code{table} with the number of 0, 1, 2, ... assigned EMRTs.
      }

      \item{performance}{\code{signatute(object="Synapter", verbose =
	  TRUE)}: Returns (and displays, if \code{verbose}) the
	performance of the synapter analysis.
      }

      \item{performance2}{\code{signatute(object="Synapter", verbose =
	  TRUE)}: Returns (and displays, if \code{verbose}) information
	about number of missing values and identification source of
	transfered EMRTs.
      }

    }
  }

  \subsection{Filters}{

    \describe{

      \item{filterUniqueDbPeptides}{\code{signature(object="Synapter",
          missedCleavages = 0, PLGS = TRUE, IisL = TRUE, verbose = TRUE)}:
	This method first digests the fasta database file and keeps
	unique tryptic peptides. (NOTE: since version 1.5.3, the tryptic
	digestion uses the \code{cleaver} package, replacing the more
	simplistic inbuild function. The effect of this change is
	documented in https://github.com/lgatto/synapter/pull/47).
  Since version 1.7.2 this method tries to emulate PLGS' peptide cleavage
  method. In fact PLGS performs a perfect cleavage after K or R if not followed
  by P, D, E, K, R first. In a second run it cleaves after K or R if not
  followed by P with mis-cleavages (NOTE: more details are documented in
  https://github.com/lgatto/synapter/pull/57)
  Set \code{PLGS = FALSE} to use the \code{cleaver} rules instead.
  See the references for more details about digest using trypsin.
  The number of maximal missed cleavages can be set as \code{missedCleavages}
	(default is 0).
  If \code{IisL = TRUE} (default) Isoleucin and Leucin are treated as
  the same aminoacid. In this case sequences like "ABCI", "ABCL" are removed
  because they are not unqiue anymore. If \code{IisL = FALSE} "ABCI" and
  "ABCL" are reported as unique.
  The peptide sequences are then used as a
	filter against the identification and quantitation peptides,
	where only unique proteotyptic instances (no miscleavage allowed
	by default) are eventually kept in the \code{object} instance.
	This method also removes any additional duplicated peptides,
	that would not match any peptides identified in the fasta
	database.  }

      \item{filterUniqueQuantDbPeptides}{\code{signature(object="Synapter",
	  missedCleavages = 0, PLGS = TRUE, IisL = TRUE, verbose = TRUE)}: As
	\code{filterUniqueDbPeptides} for quantitation peptides only.  }

      \item{filterUniqueIdentDbPeptides}{\code{signature(object="Synapter",
	  missedCleavages = 0, PLGS = TRUE, IisL = TRUE, verbose = TRUE)}: As
	\code{filterUniqueDbPeptides} for identification peptides
	only. }

      \item{filterQuantPepScore}{\code{signature(object="Synapter", fdr
	  = "numeric", method = "character")}: Filters the quantitation
	peptides using \code{fdr} false discovery rate. \code{fdr} is
	missing by default and is retrieved with \code{getPepScoreFdr}
	automatically. If not set, default value of 0.01 is
	used. \code{method} defines how to performe p-value adjustment;
	one of \code{BH}, \code{Bonferrone} or \code{qval}. See details
	section for more information.  }

      \item{filterIdentPepScore}{\code{signature(object="Synapter", fdr
	  = "numeric", method = "charactet")}: As
	\code{filterQuantPepScore}, but for identification peptides.  }

      \item{filterQuantProtFpr}{\code{signature(object="Synapter", fpr
	  = "numeric")}: Filters quantitation peptides using the protein
	false positive rate (erroneously defined as a FPR, should be
	FDR), as reported by PLGS, using threshold set by \code{fpr}
	(missing by default) or retrieved by \code{getProtFpr}.  }

      \item{filterIdentProtFpr}{\code{signature(object="Synapter", fpr =
	  "numeric")}: as \code{filterQuantProtFpr}, but for
	identification peptides.  }

      \item{filterQuantPpmError}{\code{signature(object="Synapter", ppm
	  = "numeric")}: Filters the quantitation peptides based on the
	mass tolerance \code{ppm} (default missing) provided or
	retrieved automatically using \code{getPpmError}.  }

      \item{filterIdentPpmError}{\code{signature(object="Synapter")}: as
	\code{filterQuantPpmError}, but for identification peptides.  }

    \item{filterMatchedEMRTsByCommonPeaks}{\code{signature(object="Synapter",
        matchColumn = "spectrum.quantXfragments.ident")}:
      Turns non unique matches into unique matches based on the cross matching
      results and a minimal number of common peaks between the identification
      and quantiation fragments. Use
      \code{\link{setCrossMatchingMinimalNumberOfCommonPeaks}} to set the
      minimal number of common peaks. If there are multiple non unique matches
      above this number they are kept as non unique matches.
      Use \code{matchColumn} to select which cross matching result should be
      used for filtering: \code{"spectrum.quantXfragments.ident"} for quantitation
      spectra vs identification fragments,
      \code{"spectrum.identXfragments.quant"} for identification spectra vs
      quantitation fragments, and \code{"fragments.identXfragments.quant"} for
      identification fragments vs quantitation fragments.
    }

    }
  }

  \subsection{Plotting}{

    \describe{

      \item{plotPpmError}{\code{signature(object="Synapter", what =
	  "character")}: Plots the proportion of data against the mass
	error tolerance in ppms. Depending on \code{what}, the data
	for identification (\code{what = "Ident"}), quantitation
	(\code{what = "Quant"}) or \code{"both"} is plotted.  }

      \item{plotRtDiffs}{\code{signature(object="Synapter", ...)}: Plots
	a histogram of retention time differences after
	alignments. \code{...} is passed to \code{hist}.  }

      \item{plotRt}{\code{signature(object="Synapter", what =
	  "character", f = "numeric", nsd = "numeric")}: Plots the
	Identification - Quantitation retention time difference as a
	function of the Identification retention time. If \code{what} is
	\code{"data"}, two plots are generated: one ranging the full range
	of retention time differences and one focusing on the highest data
	point density and showing models with various \code{span}
	parameter values, as defined by \code{f} (default is 2/3, 1/2,
	1/4, 1/10, 1/16, 1/25, 1/50, passed as a numed numeric). If
	\code{what} is \code{"model"}, a focused plot with the applied
	span parameter is plotted and areas of \code{nsd} (default is
	\code{x(1, 3, 5)} number of standard deviations are shaded around
      the model.
    }

    \item{plotPepScores}{\code{signature(object="Synapter")}: Plots the
      distribution of random and regular peptide scores for
      identification and quantitation features. This reflects how
      peptide p-values are computed. See also \code{getPepNumbers}.  }

    \item{plotFdr}{\code{signature(object="Synapter", method =
	"character")}: Displays 2 plots per identification and
      quantitation peptides, showing the number of significant
      peptides as a function of the FDR cut-off and the expected false
      number of false positive as a number of significant
      tests. PepFrag 1 and 2 peptides are illustrated on the same
      figures. These figures are adapted from \link[qvalue]{qplot}.
      \code{method}, one of \code{"BH"}, \code{"Bonferroni"} or
      \code{"qval"}, defines what identification statistics to use.
    }

    \item{plotEMRTtable}{\code{signature(object="Synapter")}: Plots
      the barchart of number or 0, 1, 2, ... assigned EMRTs (see
      \code{getEMRTtable}) .
    }

    \item{plotGrid}{\code{signature(object="Synapter", what =
	"character")}: Plots a heatmap of the respective grid search
      results. This grid to be plotted is controlled by \code{what}:
      \code{"total"}, \code{"model"} or \code{"details"} are
      available.
    }

    \item{plotFeatures}{\code{signature(object="Synapter", what =
	"character", xlim = "numeric", ylim = "numeric")}: Plots the
      retention time against precursor mass space. If \code{what} is
      \code{"all"}, three such plots are created side by side: for the
      identification peptides, the quantitation peptides and the
      quantitation Pep3D data. If \code{what} is \code{"some"}, a
      subset of the rt/mass space can be defined with \code{xlim}
      (default is \code{c(40, 60)}) and \code{ylim} (default is
      \code{c(1160, 1165)}) and identification peptide, quantitation
      peptides and EMRTs are presented on the same graph as grey
      dots, blue dots and red crosses respectively. In addition,
      rectangles based on the ppm and nsd defined tolerances (see
      \code{setPpmError} and \code{setNsdError}) are drawn and
      centered at the expected modelled retention time. This last
      figure allows to visualise the EMRT matching.
    }

    \item{plotCrossMatching}{\code{signature(object = "Synapter",
        key = "character", column = "character",
        verbose = "logical", \dots)}:
      Plots two spectra and fragments against each other. Please see
      \code{\link{plotCrossMatching}} for details.
    }

    \item{plotCrossMatchingSummary}{\code{signature(object = "Synapter",
        matchColumn = "spectrum.quantXfragments.ident")}:
      Creates three plots showing the performance of the cross matching
      procedure.  The left panel shows the F1-score and the FDR calculated for
      an equal sized subset of true and false matches.
      The panel in the middle shows the development of the true-positives,
      false-positives, true-negatives and false-negatives using an
      increasing number of common peaks as threshold. On the right panel
      boxplots for all matches and the sampled matches are shown.
      Use \code{matchColumn} to select which cross matching result should be
      visualized: \code{"spectrum.quantXfragments.ident"} for quantitation
      spectra vs identification fragments,
      \code{"spectrum.identXfragments.quant"} for identification spectra vs
      quantitation fragments, and \code{"fragments.identXfragments.quant"} for
      identification fragments vs quantitation fragments.
      This function returns a contingency matrix invisibly.
    }
  }
}

\subsection{Importers}{

    \describe{
      \item{loadSpectrumXmlFiles}{\code{signature(object="Synapter", filenames = "list",
          verbose = "logcial")}:
        Imports \sQuote{Spectrum.xml} files. The file path has to be given as
        named list (\code{filenames}). The names have to be
        \code{identspectrum} and \code{quantspectrum} for the identification and
        the quantitation spectrum file. If \code{verbose} is \code{TRUE}
        (default) a progress bar is shown.
      }
      \item{loadFragmentCsvFiles}{\code{signature(object="Synapter",
          filenames = "list", removeNeutralLoss = "logical",
          verbose = "logcial")}:
        Imports \sQuote{final_fragment.csv} files. The file path has to be given
        as named list (\code{filenames}). The names have to be
        \code{identfragments} and \code{quantfragments} for the identification
        and the quantitation fragment file. In the default setting rows
        with a neutrol loss value not equal to \sQuote{none} are removed.
        To keep these rows you have to set \code{removeNeutralLoss = FALSE}.
        If \code{verbose} is \code{TRUE} (default) a progress bar is shown.
      }
    }
  }

\subsection{Exporters}{

    \describe{
      \item{writeMergedPeptides}{\code{signature(object="Synapter", file
	  = "character", what = "character", ...)}: Exports the merged
	peptide data to a comma-separated \code{file} (default name is
	\code{"Res-MergedPeptides.csv"}). \code{what} can be
	\code{"light"} (default) or \code{"full"} and specifies if the
	full data or only selected columns are exported. \code{...} are
	passed to \code{\link{write.csv}}.
      }
      \item{writeMatchedEMRTs}{\code{signature(object="Synapter", file =
	  "character", what = "character", ...)}: As above, saving the
	matched EMRT table.
      }
      \item{writeIdentPeptides}{\code{signature(object="Synapter", file
	  = "character", ...)}: As above, exporting the identification
	peptide data.
      }
      \item{writeQuantPeptides}{\code{signature(object="Synapter", file
	  = "character", ...)}: A above, exporting the quantitation
	peptide data.
      }
    }
  }

  \subsection{Other}{

    \describe{
      \item{as(, "MSnSet")}{\code{signature(x = "Synapter")}: Coerce
	object from \code{Synapter} to \code{MSnSet} class. }
    }
  }
}
}


\references{
  \enc{Käll}{Kall} L, Storey JD, MacCoss MJ, Noble WS
  Posterior error probabilities and false discovery rates: two sides of
  the same coin.
  J Proteome Res. 2008a Jan; 7:(1)40-4

  Bonferroni single-step adjusted p-values for strong control of the
  FWER.

  Benjamini Y. and Hochberg Y.
  Controlling the false discovery rate: a practical and powerful
  approach to multiple testing.
  J. R. Statist. Soc. B., 1995, Vol. 57: 289-300.

  Storey JD and Tibshirani R.
  Statistical significance for genome-wide experiments.
  Proceedings of the National Academy of Sciences, 2003, 100: 9440-9445.

  \enc{Käll}{Kall}, Storey JD, MacCoss MJ, Noble WS
  Assigning significance to peptides identified by tandem mass
  spectrometry using decoy databases.
  J Proteome Res. 2008b Jan; 7:(1)29-34

  Improving qualitative and quantitative performance for MSE-based label
  free proteomics, N.J. Bond, P.V. Shliaha, K.S. Lilley and L. Gatto,
  Journal of Proteome Research, 2013, in press.

  The Effects of Travelling Wave Ion Mobility Separation on Data
  Independent Acquisition in Proteomics Studies, P.V. Shliaha, N.J. Bond,
  L. Gatto and K.S. Lilley, Journal of Proteome Research, 2013, in press.

  Trypsin cleavage:

  Glatter, Timo, et al.
  Large-scale quantitative assessment of different in-solution protein digestion
  protocols reveals superior cleavage efficiency of tandem Lys-C/trypsin
  proteolysis over trypsin digestion.
  Journal of proteome research 11.11 (2012): 5145-5156.
  \url{http://dx.doi.org/10.1021/pr300273g}

  Rodriguez, Jesse, et al.
  Does trypsin cut before proline?.
  Journal of proteome research 7.01 (2007): 300-305.
  \url{http://dx.doi.org/10.1021/pr0705035}

  Brownridge, Philip, and Robert J. Beynon.
  The importance of the digest: proteolysis and absolute quantification in
  proteomics.
  Methods 54.4 (2011): 351-360.
  \url{http://dx.doi.org/10.1016/j.ymeth.2011.05.005}

  cleaver's rules are taken from:
  \url{http://web.expasy.org/peptide_cutter/peptidecutter_enzymes.html#Tryps}
}


\examples{
library(synapter) ## always needed

\dontrun{
## (1) Construction - to create your own data objects
synapterTiny <- Synapter()
}

## let's use synapterTiny, shipped with the package
synapterTinyData() ## loads/prepares the data
synapterTiny ## show object

## (2) Filtering
## (2.1) Peptide scores and FDR

## visualise/explore peptide id scores
plotPepScores(synapterTiny)
getPepNumbers(synapterTiny)

## filter data
filterUniqueDbPeptides(synapterTiny) ## keeps unique proteotypic peptides
filterPeptideLength(synapterTiny, l = 7) ## default length is 7

## visualise before FDR filtering
plotFdr(synapterTiny)

setPepScoreFdr(synapterTiny, fdr = 0.01) ## optional
filterQuantPepScore(synapterTiny, fdr = 0.01) ## specifying FDR
filterIdentPepScore(synapterTiny) ## FDR not specified, using previously set value

## (2.2) Mass tolerance
getPpmErrorQs(synapterTiny)
plotPpmError(synapterTiny, what="Ident")
plotPpmError(synapterTiny, what="Quant")

setIdentPpmError(synapterTiny, ppm = 20) ## optional
filterQuantPpmError(synapterTiny, ppm = 20)
## setQuantPpmError(synapterTiny, ppm = 20) ## set quant ppm threshold below
filterIdentPpmError(synapterTiny, ppm=20)

filterIdentProtFpr(synapterTiny, fpr = 0.01)
filterQuantProtFpr(synapterTiny, fpr = 0.01)

getPpmErrorQs(synapterTiny) ## to be compared with previous output

## (3) Merge peptide sequences
mergePeptides(synapterTiny)

## (4) Retention time modelling
plotRt(synapterTiny, what="data")
setLowessSpan(synapterTiny, 0.05)
modelRt(synapterTiny) ## the actual modelling
getRtQs(synapterTiny)
plotRtDiffs(synapterTiny)
## plotRtDiffs(synapterTiny, xlim=c(-1, 1), breaks=500) ## pass parameters to hist()
plotRt(synapterTiny, what="model") ## using default nsd 1, 3, 5
plotRt(synapterTiny, what="model", nsd=0.5) ## better focus on model

plotFeatures(synapterTiny, what="all")
setRtNsd(synapterTiny, 3)     ## RtNsd and PpmError are used for detailed plot
setPpmError(synapterTiny, 10) ## if not set manually, default values are set automatically
plotFeatures(synapterTiny, what="some", xlim=c(36,44), ylim=c(1161.4, 1161.7))
## best plotting to svg for zooming

set.seed(1) ## only for reproducibility of this example

## (5) Grid search to optimise EMRT matching parameters
searchGrid(synapterTiny,
           ppms = 7:10,  ## default values are 5, 7, ..., 20
           nsds = 1:3,   ## default values are 0.5, 1,  ..., 5
           subset = 0.2) ## default is 1
## alternatively, use 'n = 1000' to use exactly
## 1000 randomly selected features for the grid search
getGrid(synapterTiny)  ## print the grid
getGridDetails(synapterTiny)  ## grid details
plotGrid(synapterTiny, what = "total")   ## plot the grid for total matching
plotGrid(synapterTiny, what = "model")   ## plot the grid for matched modelled feature
plotGrid(synapterTiny, what = "details") ## plot the detail grid
getBestGridValue(synapterTiny)  ## return best grid values
getBestGridParams(synapterTiny) ## return parameters corresponding to best values
setBestGridParams(synapterTiny, what = "auto") ## sets RtNsd and PpmError according the grid results
## 'what' could also be "model", "total" or "details"
## setPpmError(synapterTiny, 12) ## to manually set values
## setRtNsd(synapterTiny, 2.5)

## (6) Matching ident peptides and quant EMRTs
findEMRTs(synapterTiny)
plotEMRTtable(synapterTiny)
getEMRTtable(synapterTiny)
performance(synapterTiny)
performance2(synapterTiny)

## (7) Exporting data to csv spreadsheets
writeMergedPeptides(synapterTiny, what = "light") ## or what="full"
writeMergedPeptides(synapterTiny, file = "myresults.csv", what="light")
writeMatchedEMRTs(synapterTiny, what = "light")   ## or what="full"
writeMatchedEMRTs(synapterTiny, file = "myresults2.csv", what="light")
## These will export the filter peptide data
writeIdentPeptides(synapterTiny, file = "myIdentPeptides.csv")
writeQuantPeptides(synapterTiny, file = "myQuantPeptides.csv")
## If used right after loading, the non-filted data will be exported
}

\author{Laurent Gatto \email{lg390@cam.ac.uk}}

\keyword{classes}

