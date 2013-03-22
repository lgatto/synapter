\name{synapter-package}
\alias{synapter-package}
\alias{synapter}
\docType{package}

\title{
  Combine label-free data for optimal identification and quantitation 
}

\description{
  The synapter package provides functionality to reanalyse label-free
  proteomics data acquired on a Synapt G2 mass spectrometer. One or several 
  runs, possibly processed with additional ion mobility separation to
  increase identification accuracy can be combined to other quantitation
  files to maximise identification and quantitation accuracy.  
}


\details{
  Three pipelines of variying flexibility are proposed the preform data
  analysis: (1) a graphical user interface is provided through the
  \code{synatperGUI} function, (2) the \code{synergise} function is a
  single entry function for a complete analysis and (3) low level,
  step-by-step data processing can be achieved as described in 
  \code{?Synapter}. 
  
  A high-level overview of the package and how to operate it can be
  found in the vignette, accessinble with
  \code{synapterGuide()}. Detailed information about the data processing 
  can be found in the respective function and class manual pages
  appropriately referenced in the vignette.

  For questions, use the Biocondcutor mailing list or contact the
  author. The vignette has a section with details on where/how to get
  help. 
}

\author{
Laurent Gatto, Pavel V. Shliaha and Nick J. Bond

Maintainer: Laurent Gatto <lg390@cam.ac.uk>
}

\references{
  \emph{Improving qualitative and quantitative performance for MSE-based label
    free proteomics}, N.J. Bond, P.V. Shliaha, K.S. Lilley and L. Gatto,
  Journal of Proteome Research, 2013, in press.
  
  \emph{The Effects of Travelling Wave Ion Mobility Separation on Data
    Independent Acquisition in Proteomics Studies}, P.V. Shliaha, N.J. Bond,
  L. Gatto and K.S. Lilley, Journal of Proteome Research, 2013, in press. 
}


\keyword{ package }

