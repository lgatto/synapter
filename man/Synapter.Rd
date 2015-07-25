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
\alias{rescueEMRTs,Synapter-method}
\alias{rescueEMRTs}
\alias{searchGrid,Synapter-method}
\alias{searchGrid}
\alias{fragmentMatching,Synapter-method}
\alias{fragmentMatching}

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
\alias{setImDiff,Synapter-method}
\alias{setImDiff}
\alias{getImDiff,Synapter-method}
\alias{getImDiff}
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
\alias{setFragmentMatchingPpmTolerance,Synapter-method}
\alias{setFragmentMatchingPpmTolerance}
\alias{getFragmentMatchingPpmTolerance,Synapter-method}
\alias{getFragmentMatchingPpmTolerance}
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
\alias{filterFragments,Synapter-method}
\alias{filterFragments}
\alias{filterUniqueMatches,Synapter-method}
\alias{filterUniqueMatches}
\alias{filterNonUniqueMatches,Synapter-method}
\alias{filterNonUniqueMatches}

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
\alias{plotFragmentMatchingPerformance,Synapter-method}
\alias{plotFragmentMatchingPerformance}
\alias{plotCumulativeNumberOfFragments,Synapter-method}
\alias{plotCumulativeNumberOfFragments}

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
\alias{isCurrent,Synapter-method}
\alias{validObject,Synapter-method}
\alias{updateObject,Synapter-method}

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
    (could be an RDS file created by \code{link{createUniquePeptideDbRds}})
    and (optional 'identfragments' and 'quantspectra').
    If missing, dialog boxes pop up to select the six files
    manually.
    \code{identpeptide} can be a \code{csv} final peptide file (from
    PLGS) or a saved \code{"\linkS4class{MasterPeptides}"} data object as
    created by \code{\link{makeMaster}} if working with \emph{master}
    peptide data. To serialise the \code{"\linkS4class{MasterPeptides}"}
    instance, use the \code{saveRDS} function, and file extenstion
    \code{rds}. This master file could contain a fragment library as well. In
    this case the 'identfragments' argument would be ignored.
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
    with the \code{Synapter} constructor. This function opens six dialog
    boxes for the user to point to the input files,
    namely (and in that order), the identification final peptide csv file, the
    quantitation final peptide  csv file, the quantitation Pep3D csv file
    (as exported from the PLGS software), the fasta file use for peptide
    identification, the identification fragments csv file and the quantitation
    spectra xml file. You can choose 'cancel' to ignore the last two files and
    to skip the fragment matching later.
    The fasta file ('fasta') could be an RDS file generated by
    \code{link{createUniquePeptideDbRds}}, too.
    The files are read and the data is stored in the newly
    created \code{Synapter} instance. The file names can also be
    specified as a named list with names 'identpeptide', 'quantpeptide',
    'quantpep3d', 'identfragments' and 'quantspectra' respectively.

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
    \code{setPpmError}, \code{setRtNsd}, \code{setImDiff} respectively,
    or using \code{setBestGridParams} to apply best parameters as defined using
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

  \subsection{Remove Less Intense Peaks}{
    As an additional step it is possible to remove less intense peaks from the
    spectra and fragment data. Use \code{plotCumulativeNumberOfFragments} to
    plot the number of fragments vs the intensity and to find a good threshold.
    The \code{filterFragments} method could remove peaks if the intensity is
    below a specified threshold via the \code{minIntensity} argument. Set the
    \code{maxNumber} argument to keep only the \code{maxNumber} highest
    peaks/fragments. The \code{what} argument controls the data on which the
    filter is applied. Use \code{what = "fragments.ident"} for the
    identification fragments and \code{what = "spectra.quant"} for the
    quantiation spectra data.
  }

  \subsection{Fragment Matching}{
    After importing fragment and spectra data it is possible to
    match peaks between the identification fragments and the quantitation
    spectra using the \code{fragmentMatching} method.
    Use \code{setFragmentMatchingPpmTolerance} to set
    the maximal allowed tolerance for considering a peak as identical.
    There are two different methods to visualise the results of the fragment
    matching procedure. \code{plotFragmentMatching} plots the fragments and
    spectra for each considered pair.
    \code{plotFragmentMatchingPerformance} draws two plots. On the left panel
    you could see the performance of different thresholds for the number of
    common peaks for unique matches.  The right panel visualizes the performance
    of different differences (delta) of common peaks between the best match
    (highest number of common peaks) and the second best match in each non
    unique match group.
    Use \code{filterUniqueMatches} and \code{filterNonUniqueMatches} to remove
    unique or non unique matches below the threshold of common peaks
    respective the difference in common peaks from the MatchedEMRTs data.frame.
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
	  "numeric", nsd = "numeric", imdiff = "numeric")}:
	Finds EMRTs matching identification peptides using \code{ppm}
	mass tolerance, \code{nsd} number of retention time standard
	deviations and \code{imdiff} difference in ion mobility.
  The last three parameters are optional if previously
	set with \code{setPpmError}, \code{setRtNsd}, \code{setImDiff},
  or, better, \code{setBestGridParams} (see below).
      }

    \item{rescueEMRTs}{\code{signature(object = "Synapter",
    method = c("rescue", "copy"))}:
  The \code{method} parameter defined the behaviour for those
  high confidence features that where identified in both identification and
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
	  ppms="numeric", nsds="numeric", imdiffs = "numeric",
    subset="numeric", n = "numeric", verbose="logical")}:
  Performs a grid search. The
	grid is defined by the \code{ppm}, \code{nsd} and \code{imdiffs}
  numerical vectors, representing the sequence of values to be
	tested. Default are \code{seq(5, 20, 2)}, \code{seq(0.5, 5,
	  0.5)}, \code{seq(0.2, 2, 0.2)} respectively. To ignore ion mobility set
  \code{imdiffs = Inf}.
  \code{subset} and \code{n} allow to use a
	randomly  chosen subset of the data for the grid search to
	reduce search time. \code{subset} is a numeric, between 0 and 1,
	describing the percentage of data to be used; \code{n} specifies
	the absolute number of feature to use. The default is to use all
	data. \code{verbose} controls whether textual output should be
	printed to the console. (Note, the \code{mergedEMRTs} value used
	in internal calls to \code{findEMRTs} is \code{"transfer"} - see
	\code{findEMRTs} for details).
      }

      \item{fragmentMatching}{\code{signature(object="Synapter",
      ppm = "numeric", verbose = "logical"}:
      Performs a fragment matching between spectra and fragment data.
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
	\code{character} of length 6 with the names of the input files
	used as \code{identpeptide}, \code{quantpeptide},
	\code{quantpep3d}, \code{fasta}, \code{identfragments} and
  \code{quantspectra}.
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

      \item{setImDiff}{\code{signature(object="Synapter", imdiff =
	  "numeric")}: Sets the ion mobility tolerance \code{imdiff},
	default is 0.5.
      }

      \item{getImDiff}{\code{signature(object="Synapter")}: Returns the
	value of the ion mobility tolerance \code{imdiff}.
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

      \item{setFragmentMatchingPpmTolerance}{\code{signature(object="Synapter",
        ppm = "numeric")}: Sets the fragment matching mass tolerance \code{ppm}
        (default is 25).
      }

      \item{getFragmentMatchingPpmTolerance}{\code{signature(object="Synapter")}:
        Returns the fragment matching mass tolerance in ppm.
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
          missedCleavages = 0, IisL = TRUE, verbose = TRUE)}:
	This method first digests the fasta database file and keeps
	unique tryptic peptides. (NOTE: since version 1.5.3, the tryptic
	digestion uses the \code{cleaver} package, replacing the more
	simplistic inbuild function. The effect of this change is
	documented in https://github.com/lgatto/synapter/pull/47).
  The number of maximal missed cleavages can be set as \code{missedCleavages}
	(default is 0).
  If \code{IisL = TRUE} Isoleucin and Leucin are treated as the same aminoacid.
  In this case sequences like "ABCI", "ABCL" are removed
  because they are not unqiue anymore. If \code{IisL = FALSE} (default)
  "ABCI" and "ABCL" are reported as unique.
  The peptide sequences are then used as a
	filter against the identification and quantitation peptides,
	where only unique proteotyptic instances (no miscleavage allowed
	by default) are eventually kept in the \code{object} instance.
	This method also removes any additional duplicated peptides,
	that would not match any peptides identified in the fasta
	database.  }

      \item{filterUniqueQuantDbPeptides}{\code{signature(object="Synapter",
	  missedCleavages = 0, IisL = TRUE, verbose = TRUE)}: As
	\code{filterUniqueDbPeptides} for quantitation peptides only.  }

      \item{filterUniqueIdentDbPeptides}{\code{signature(object="Synapter",
	  missedCleavages = 0, IisL = TRUE, verbose = TRUE)}: As
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

    \item{filterFragments}{\code{signature(object = "Synapter",
        what = c("fragments.ident", "spectra.quant"),
        minIntensity = "numeric", maxNumber = "numeric", verbose = "logical")}:
      Filters the spectra/fragment data using a minimal intensity threshold
      (\code{minIntensity}) or a maximal number of peaks/fragments threshold
      (\code{maxNumber}). Please note that the maximal number is transfered to
      an intensity threshold and the result could contain less peaks than
      specified by \code{maxNumber}.
      If both arguments are given, the more aggressive one is chosen.
      Use the \code{what} argument to specify the data that should be filtered.
      Set \code{what = "fragments.ident"} for the identification fragment data
      or \code{what = "spectra.quant"} for the quantiation spectra.
      If \code{verbose} is \code{TRUE} (default) a progress bar is shown.
    }

    \item{filterUniqueMatches}{\code{signature(object="Synapter", minNumber =
        "numeric")}:
      Removes all unique matches that have less than \code{minNumber} of
      peaks/fragments in common. Use \code{plotFragmentMatchingPerformance}
      (left panel) to find an ideal threshold.
    }

    \item{filterNonUniqueMatches}{\code{signature(object="Synapter", minDelta =
        "numeric")}:
      Removes all non unique matches that have a difference between the best
      match (highest number of common peaks/fragments, treated as true match)
      and the second best match (second highest number of common
      peaks/fragments) less than \code{minDelta}. For the matches above the
      threshold only the one with the highest number of common peaks/fragments
      in each match group is kept.
      Use \code{plotFragmentMatchingPerformance} (right panel) to find an ideal
      threshold.
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
	"character"), maindim = "character"}:
      Plots a heatmap of the respective grid search
      results. This grid to be plotted is controlled by \code{what}:
      \code{"total"}, \code{"model"} or \code{"details"} are
      available. If ion mobility was used in the grid search you can use
      \code{maindim} to decided which dimensions should be shown. \code{maindim}
      could be one of \code{"im"} (default), \code{"rt"} and \code{"mz"}. If
      \code{maindim = "im"} a heatmap for each ion mobility threshold is drawn.
      For \code{maindim = "rt"} and \code{maindim} you get a heatmap for each
      retention time respective mass threshold.
    }

    \item{plotFeatures}{\code{signature(object="Synapter", what =
	"character", xlim = "numeric", ylim = "numeric", ionmobiltiy = "logical")}:
      Plots the retention time against precursor mass space.
      If \code{what} is \code{"all"}, three (six if ion mobility is available
      and \code{ionmobility = TRUE} (default is \code{FALSE});
      three additional plots with precursor mass against ion mobility)
      such plots are created side by side: for the
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

    \item{plotFragmentMatching}{\code{signature(object = "Synapter",
        key = "character", column = "character",
        verbose = "logical", \dots)}:
      Plots two spectra and fragments against each other. Please see
      \code{\link{plotFragmentMatching}} for details.
    }

    \item{plotFragmentMatchingPerformance}{\code{signature(object = "Synapter",
        matchColumn = "spectrum.quantXfragments.ident")}:
    Creates two plots. The left panel shows the performance of filtering the
    unique matches using a different number of common peaks. The right panel
    shows the performance of filtering the non unique matches using different
    differences (delta) in common peaks/fragments. These differences (delta) are
    calculated between the match with the highest number of common
    peaks/fragments and the second highest one.
    Use \code{filterUniqueMatches} and \code{filterNonUniqueMatches} to filter
    the \code{MatchedEMRT} \code{data.frame} using one of these thresholds.
    This function returns a \code{list} with two named elements (\code{unqiue}
    and \code{nonunqiue} invisibly. Each \code{list} element contains a
    \code{matrix} with six columns. The first column
    \code{ncommon}/\code{deltacommon} contains the
    thresholds. Column 2 to 6 are the true positives \code{tp}, false positives
    \code{fp}, true negatives \code{tn}, false negatives \code{fn} and FDR
    \code{fdr} for the current threshold (in that row).
    }
    \item{plotFragmentMatchingDiff}{\code{signature(object = "Synapter",
        matchColumn = "spectrum.quantXfragments.ident")}:
    Shows a boxplot for the differences of the first highest number of common
    peaks to the second one for true and false matches in each match group.
    }

    \item{plotCumulativeNumberOfFragments}{\code{signature(object = "Synapter",
        what = c("fragments.ident", "spectra.quant"))}:
    Plots the cumulative number of the fragments/peaks vs their intensity (log10
    scaled). Use the \code{what} argument to create this plot for the
    identification fragments (\code{what = "fragments.quant"}) or the
    the quantitation spectra (\code{what = "spectra.quant"}).
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
    \describe{
      \item{validObject}{\code{signature(object = "Synapter")}:
        Test whether a given \code{Synapter} object is valid.}
    }
    \describe{
      \item{updateObject}{\code{signature(object = "Synapter")}: Updates an old
       \code{Synapter} object.}
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

